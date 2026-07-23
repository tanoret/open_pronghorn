"""Fit and audit the effective MSR corrosion/plating correlation.

Terminology matters here.  This module estimates the 61-parameter correlation
from the workbook's active rows and reports in-sample fit diagnostics.  It does
not create an independent validation split.  OpenPronghorn's separate 0D replay
tests software consistency between this Python correlation and the C++ model.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .ingest import build_model_tables
from .model import (
    MoltenSaltBVModel,
    PARAMETER_SPECS,
    initial_parameter_vector,
    lower_bounds,
    priors_as_residuals,
    upper_bounds,
    vector_to_params,
)


# These roles enter the numerical objective.  Input-only and auxiliary workbook
# rows are retained in generated audit tables but must not influence the fit.
ACTIVE_FIT_ROLES = {"direct", "range", "upper", "lower", "ranking"}
DIRECT_FIT_ROLES = {"direct", "range"}


def _is_positive_finite(x: Any) -> bool:
    try:
        return bool(np.isfinite(float(x)) and float(x) > 0)
    except Exception:
        return False


def _sqrt_weight(row: pd.Series) -> float:
    w = row.get("quality_weight", 1.0)
    try:
        w = float(w)
    except Exception:
        w = 1.0
    if not np.isfinite(w) or w <= 0:
        w = 0.5
    return math.sqrt(w)


def residuals_for_targets(
    model: MoltenSaltBVModel, targets: pd.DataFrame
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Evaluate relation-aware, source-quality-weighted data residuals.

    A direct or range row contributes a signed log residual.  Upper and lower
    rows contribute only when their inequality is violated.  The single MSRE
    deposition-ranking row contributes two ordering residuals, so 38 active
    workbook rows produce 39 scalar data residuals.
    """
    residuals: list[float] = []
    details: list[dict[str, Any]] = []

    for _, row in targets.iterrows():
        role = str(row.get("fit_role", ""))
        kind = str(row.get("response_kind", ""))
        if role not in ACTIVE_FIT_ROLES:
            continue
        sigma = float(row.get("default_sigma_ln", 0.75) or 0.75)
        sigma = max(sigma, 0.15)
        weight = _sqrt_weight(row)

        if role == "ranking" and kind == "noble_metal_deposition_ranking":
            ranking = model.deposition_ranking(row)
            # Two constraints from MSRE: metal > graphite; turbulent metal > laminar metal.
            margin = math.log(1.10)
            r1 = max(0.0, margin - math.log(max(ranking["metal"], 1e-30) / max(ranking["graphite"], 1e-30))) / sigma * weight
            r2 = max(0.0, margin - math.log(max(ranking["turbulent_metal"], 1e-30) / max(ranking["laminar_metal"], 1e-30))) / sigma * weight
            residuals.extend([r1, r2])
            details.append({"measurement_id": row.get("measurement_id"), "response_kind": kind, "fit_role": role, "prediction": ranking["turbulent_metal"] / max(ranking["graphite"], 1e-30), "target_mid": np.nan, "residual_ln": max(r1, r2), "constraint_pass": r1 == 0 and r2 == 0})
            continue

        pred = model.predict_response(row)
        if not _is_positive_finite(pred):
            continue
        pred = float(pred)
        target_mid = row.get("target_mid")
        target_low = row.get("target_low")
        target_high = row.get("target_high")
        res = None
        constraint_pass = True

        if role in DIRECT_FIT_ROLES:
            if not _is_positive_finite(target_mid):
                continue
            target = float(target_mid)
            res = math.log(pred / target) / sigma * weight
            if role == "range" and _is_positive_finite(target_low) and _is_positive_finite(target_high):
                constraint_pass = float(target_low) <= pred <= float(target_high)
            else:
                constraint_pass = abs(math.log(pred / target)) <= math.log(2.0)
        elif role == "upper":
            upper = target_high if _is_positive_finite(target_high) else target_mid
            if not _is_positive_finite(upper):
                continue
            upper = float(upper)
            violation = math.log(pred / upper)
            res = max(0.0, violation) / sigma * weight
            constraint_pass = violation <= 0
        elif role == "lower":
            lower = target_low if _is_positive_finite(target_low) else target_mid
            if not _is_positive_finite(lower):
                continue
            lower = float(lower)
            violation = math.log(lower / pred)
            res = max(0.0, violation) / sigma * weight
            constraint_pass = violation <= 0

        if res is None:
            continue
        residuals.append(float(res))
        details.append(
            {
                "measurement_id": row.get("measurement_id"),
                "case_id": row.get("case_id"),
                "response_kind": kind,
                "fit_role": role,
                "prediction": pred,
                "target_low": target_low,
                "target_mid": target_mid,
                "target_high": target_high,
                "sigma_ln": sigma,
                "weight": row.get("quality_weight", 1.0),
                "residual_ln_weighted": float(res),
                "constraint_pass": bool(constraint_pass),
            }
        )

    return np.asarray(residuals, dtype=float), details


def objective_vector(vector: np.ndarray, targets: pd.DataFrame, prior_strength: float = 0.35) -> np.ndarray:
    """Return data residuals followed by regularizing parameter-prior residuals.

    The priors are scientifically material: there are more fitted coefficients
    than active experimental rows.  Removing or changing ``prior_strength`` is
    a change to the calibration definition, not merely an optimizer setting.
    """
    model = MoltenSaltBVModel.from_vector(vector)
    data_resid, _ = residuals_for_targets(model, targets)
    prior_resid = priors_as_residuals(model.params, scale=prior_strength)
    return np.concatenate([data_resid, prior_resid])


def fit_model(
    workbook_path: str | Path,
    output_dir: str | Path | None = None,
    prior_strength: float = 0.35,
    max_nfev: int = 4000,
) -> dict[str, Any]:
    """Fit the effective BV model to a local source workbook.

    Parameters
    ----------
    workbook_path
        Workbook containing ``Validation_Cases``, ``Detailed_Measurements``,
        and ``Source_Index`` sheets.
    output_dir
        Optional destination for generated CSV/JSON artifacts.  ``None`` keeps
        the calculation in memory, which is useful for tests.
    prior_strength
        Multiplier on all parameter-prior residuals.
    max_nfev
        Maximum number of nonlinear least-squares function evaluations.
    """
    workbook_path = Path(workbook_path)
    if not workbook_path.is_file():
        raise FileNotFoundError(f"Calibration workbook does not exist: {workbook_path}")

    tables = build_model_tables(workbook_path)
    targets = tables["targets"].copy()

    x0 = initial_parameter_vector()
    result = least_squares(
        objective_vector,
        x0,
        args=(targets, prior_strength),
        bounds=(lower_bounds(), upper_bounds()),
        loss="soft_l1",
        f_scale=1.0,
        max_nfev=max_nfev,
        x_scale="jac",
        verbose=0,
    )
    model = MoltenSaltBVModel.from_vector(result.x)
    target_predictions = model.predict_targets(targets)
    fit_resid, details = residuals_for_targets(model, targets)
    residual_details = pd.DataFrame(details)
    case_predictions = model.predict_cases(tables["case_features"])
    metrics = compute_metrics(target_predictions, residual_details)
    metrics["optimizer_success"] = bool(result.success)
    metrics["optimizer_message"] = str(result.message)
    metrics["optimizer_cost"] = float(result.cost)
    metrics["optimizer_nfev"] = int(result.nfev)
    metrics["prior_strength"] = float(prior_strength)
    metrics["n_parameters"] = len(PARAMETER_SPECS)

    param_table = model.parameter_table()
    # Add an approximate equivalent overpotential table for common redox classes.
    redox_rows = []
    for redox in [
        "purified_baseline",
        "msre_or_fuel_baseline",
        "oxidizing_fef2",
        "reducing_be",
        "impure_moisture",
        "chloride_unspecified",
        "tellurium",
        "stressed",
        "multi_alloy",
        "fission_product",
    ]:
        redox_rows.append(
            {
                "redox_class": redox,
                "log_rate_offset": model.redox_offset(redox),
                "equivalent_eta_mV_at_650C": 1000.0 * model.bv_overpotential_equivalent_V(redox, 650.0 + 273.15),
                "corrosion_multiplier_vs_baseline": math.exp(model.redox_offset(redox) - model.redox_offset("purified_baseline")),
            }
        )
    redox_table = pd.DataFrame(redox_rows)

    outputs = {
        "model": model,
        "tables": tables,
        "target_predictions": target_predictions,
        "case_predictions": case_predictions,
        "residual_details": residual_details,
        "parameter_table": param_table,
        "redox_table": redox_table,
        "metrics": metrics,
        "optimizer_result": result,
    }

    if output_dir is not None:
        save_outputs(outputs, output_dir)
    return outputs


def compute_metrics(predictions: pd.DataFrame, residual_details: pd.DataFrame) -> dict[str, Any]:
    """Compute in-sample diagnostics without relabeling them as test accuracy."""
    metrics: dict[str, Any] = {}
    active = predictions[predictions["fit_role"].isin(ACTIVE_FIT_ROLES)].copy()
    direct = active[active["fit_role"].isin(DIRECT_FIT_ROLES)].copy()
    direct = direct[_positive_mask(direct["prediction"]) & _positive_mask(direct["target_mid"])]
    metrics["n_measurements_total"] = int(len(predictions))
    metrics["n_active_constraints"] = int(len(active))
    metrics["n_direct_or_range_targets"] = int(len(direct))
    metrics["n_upper_lower_ranking_constraints"] = int(len(active) - len(direct))

    if len(direct):
        ln_err = np.log(direct["prediction"].astype(float).to_numpy() / direct["target_mid"].astype(float).to_numpy())
        metrics["log_rmse_direct"] = float(np.sqrt(np.mean(ln_err**2)))
        metrics["median_factor_error_direct"] = float(np.exp(np.median(np.abs(ln_err))))
        metrics["within_factor_2_direct"] = float(np.mean(np.abs(ln_err) <= math.log(2.0)))
        metrics["within_factor_5_direct"] = float(np.mean(np.abs(ln_err) <= math.log(5.0)))
    else:
        metrics["log_rmse_direct"] = None
        metrics["median_factor_error_direct"] = None
        metrics["within_factor_2_direct"] = None
        metrics["within_factor_5_direct"] = None

    if len(residual_details):
        metrics["constraint_pass_fraction"] = float(residual_details["constraint_pass"].mean())
        metrics["max_abs_weighted_residual"] = float(np.nanmax(np.abs(residual_details["residual_ln_weighted"])))
    else:
        metrics["constraint_pass_fraction"] = None
        metrics["max_abs_weighted_residual"] = None

    role_counts = predictions.groupby("fit_role", dropna=False).size().to_dict()
    metrics["fit_role_counts"] = {str(k): int(v) for k, v in role_counts.items()}
    kind_counts = predictions.groupby("response_kind", dropna=False).size().to_dict()
    metrics["response_kind_counts"] = {str(k): int(v) for k, v in kind_counts.items()}
    return metrics


def _positive_mask(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").gt(0) & np.isfinite(pd.to_numeric(series, errors="coerce"))


def save_outputs(outputs: dict[str, Any], output_dir: str | Path) -> None:
    """Write generated tables and fitted results below ``output_dir``.

    This function does not copy results into OpenPronghorn's production
    database.  Promoting a new calibration is intentionally a separate review
    step.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = output_dir / "data"
    results_dir = output_dir / "results"
    data_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    tables = outputs["tables"]
    for key in ["targets", "case_features", "usage_summary"]:
        tables[key].to_csv(data_dir / f"{key}.csv", index=False)
    outputs["target_predictions"].to_csv(results_dir / "validation_predictions.csv", index=False)
    outputs["case_predictions"].to_csv(results_dir / "case_predictions_all_76_cases.csv", index=False)
    outputs["residual_details"].to_csv(results_dir / "fit_residual_details.csv", index=False)
    outputs["parameter_table"].to_csv(results_dir / "calibrated_parameters.csv", index=False)
    outputs["redox_table"].to_csv(results_dir / "effective_redox_overpotentials.csv", index=False)
    with open(results_dir / "metrics.json", "w", encoding="utf-8") as f:
        json.dump(outputs["metrics"], f, indent=2)
    with open(results_dir / "parameters.json", "w", encoding="utf-8") as f:
        json.dump(outputs["model"].params, f, indent=2)


def load_model_from_json(path: str | Path) -> MoltenSaltBVModel:
    with open(path, "r", encoding="utf-8") as f:
        params = json.load(f)
    return MoltenSaltBVModel(params)
