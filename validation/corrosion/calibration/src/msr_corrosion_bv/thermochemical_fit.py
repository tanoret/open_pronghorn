"""Calibration, validation, comparison, and persistence utilities."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .mechanistic import (
    MECHANISTIC_ACTIVE_ROLES,
    MECHANISTIC_DIRECT_ROLES,
    _constraint_pass,
    _is_positive_finite,
    mechanistic_is_supported,
    mechanistic_scope_reason,
)
from .model import MoltenSaltBVModel
from .mstdb import MSTDBPair
from .thermochemical_data import (
    THERMOCHEMICAL_PARAMETER_SPECS,
    thermochemical_initial_vector,
    thermochemical_lower_bounds,
    thermochemical_upper_bounds,
)
from .thermochemical_model import MSTDBThermochemicalCorrosionModel


def _residual_weighting(row: pd.Series | Mapping[str, Any]) -> tuple[float, float]:
    """Return the positive sigma and square-root quality weight for one target."""
    values: dict[str, float] = {}
    for name, default in (("default_sigma_ln", 0.75), ("quality_weight", 1.0)):
        raw = row.get(name, None)
        if raw is None or bool(pd.isna(raw)):
            value = default
        else:
            value = float(raw)
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive when supplied")
        values[name] = value
    return max(values["default_sigma_ln"], 0.15), math.sqrt(values["quality_weight"])


def relation_aware_log_error(
    row: pd.Series | Mapping[str, Any], prediction: float
) -> float:
    """Signed log error to the admissible experimental set.

    Direct targets use their central value. Interval targets have no residual
    inside their reported interval and use the nearest bound outside it.
    One-sided targets contribute only when their bound is violated.
    """
    if not _is_positive_finite(prediction):
        return float("nan")
    role = str(row.get("fit_role", ""))
    low, mid, high = row.get("target_low"), row.get("target_mid"), row.get("target_high")
    pred = float(prediction)
    if role == "range" and _is_positive_finite(low) and _is_positive_finite(high):
        lower, upper = float(low), float(high)
        if lower <= pred <= upper:
            return 0.0
        return math.log(pred / (lower if pred < lower else upper))
    if role == "upper":
        bound = high if _is_positive_finite(high) else mid
        if _is_positive_finite(bound):
            return max(0.0, math.log(pred / float(bound)))
        return float("nan")
    if role == "lower":
        bound = low if _is_positive_finite(low) else mid
        if _is_positive_finite(bound):
            return min(0.0, math.log(pred / float(bound)))
        return float("nan")
    if _is_positive_finite(mid):
        return math.log(pred / float(mid))
    return float("nan")

def thermochemical_prior_residuals(params: Mapping[str, float], strength: float = 0.25) -> np.ndarray:
    return np.asarray(
        [
            strength * (float(params[spec.name]) - spec.prior) / spec.prior_sigma
            for spec in THERMOCHEMICAL_PARAMETER_SPECS
            if spec.prior_sigma > 0.0
        ],
        dtype=float,
    )


def thermochemical_residuals_for_targets(
    model: MSTDBThermochemicalCorrosionModel,
    targets: pd.DataFrame,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    residuals: list[float] = []
    details: list[dict[str, Any]] = []
    for _, row in targets.iterrows():
        if not mechanistic_is_supported(row):
            continue
        raw_prediction = model.predict_response(row)
        pred = float(raw_prediction) if _is_positive_finite(raw_prediction) else 1.0e-30
        role = str(row.get("fit_role", ""))
        sigma, weight = _residual_weighting(row)
        low, mid, high = row.get("target_low"), row.get("target_mid"), row.get("target_high")
        residual: float | None = None
        if role in MECHANISTIC_ACTIVE_ROLES:
            error = relation_aware_log_error(row, pred)
            if np.isfinite(error):
                residual = error / sigma * weight
        if residual is None:
            continue
        residuals.append(float(residual))
        details.append(
            {
                "measurement_id": row.get("measurement_id"),
                "case_id": row.get("case_id"),
                "response_kind": row.get("response_kind"),
                "fit_role": role,
                "prediction": pred,
                "target_low": low,
                "target_mid": mid,
                "target_high": high,
                "relation_aware_log_error": relation_aware_log_error(row, pred),
                "residual_ln_weighted": float(residual),
            }
        )
    return np.asarray(residuals, dtype=float), details



def thermochemical_physics_residuals(
    model: MSTDBThermochemicalCorrosionModel,
    targets: pd.DataFrame,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Qualitative mechanistic constraints documented by the source experiments."""
    residuals: list[float] = []
    details: list[dict[str, Any]] = []

    # The 316H FLiNaK loop reported an Fe-rich cold-leg deposit.  The total
    # mass-gain target alone cannot identify deposit composition, so preserve
    # the observed ordering as a one-sided qualitative constraint.
    m030 = targets[targets["measurement_id"].astype(str) == "M-030"]
    if len(m030):
        row = m030.iloc[0]
        fractions = model.deposition_species_fractions(row)
        margin = math.log(1.01)
        residual = max(0.0, margin - math.log(max(fractions["Fe"], 1.0e-30) / max(fractions["Cr"], 1.0e-30))) / 0.05
        residuals.append(residual)
        details.append(
            {
                "constraint_id": "P-001",
                "measurement_id": "M-030",
                "description": "FLiNaK cold-leg deposit is Fe-rich: Fe deposit fraction > Cr deposit fraction",
                "source_doi": "10.1016/j.jnucmat.2022.153551",
                "value": fractions["Fe"] / max(fractions["Cr"], 1.0e-30),
                "required_minimum": 1.0,
                "fit_margin": 1.01,
                "residual": residual,
                "constraint_pass": fractions["Fe"] > fractions["Cr"],
            }
        )

    # MSRE/Hastelloy-N evidence is selective Cr removal.  Keep Cr as the
    # dominant source species without prescribing an exact fraction.
    m003 = targets[targets["measurement_id"].astype(str) == "M-003"]
    if len(m003):
        row = m003.iloc[0]
        fractions = model.species_flux_fractions(row)
        for competitor in ("Fe", "Ni"):
            margin = math.log(1.05)
            ratio = fractions["Cr"] / max(fractions[competitor], 1.0e-30)
            residual = max(0.0, margin - math.log(ratio)) / 0.75
            residuals.append(residual)
            details.append(
                {
                    "constraint_id": f"P-002-{competitor}",
                    "measurement_id": "M-003",
                    "description": f"Selective Hastelloy-N dissolution: Cr flux fraction > {competitor} flux fraction",
                    "source_doi": None,
                    "value": ratio,
                    "required_minimum": 1.05,
                    "fit_margin": 1.05,
                    "residual": residual,
                    "constraint_pass": residual == 0.0,
                }
            )
    return np.asarray(residuals, dtype=float), details

def _objective(
    vector: np.ndarray,
    base_model: MoltenSaltBVModel,
    mstdb: MSTDBPair,
    targets: pd.DataFrame,
    prior_strength: float,
) -> np.ndarray:
    model = MSTDBThermochemicalCorrosionModel.from_vector(base_model, mstdb, vector)
    data, _ = thermochemical_residuals_for_targets(model, targets)
    physics, _ = thermochemical_physics_residuals(model, targets)
    priors = thermochemical_prior_residuals(model.params, prior_strength)
    return np.concatenate([data, physics, priors])


def fit_thermochemical_model(
    targets: pd.DataFrame,
    base_model: MoltenSaltBVModel,
    mstdb: MSTDBPair,
    *,
    prior_strength: float = 0.25,
    max_nfev: int = 5000,
) -> dict[str, Any]:
    result = least_squares(
        _objective,
        thermochemical_initial_vector(),
        args=(base_model, mstdb, targets, prior_strength),
        bounds=(thermochemical_lower_bounds(), thermochemical_upper_bounds()),
        loss="soft_l1",
        f_scale=1.0,
        x_scale="jac",
        max_nfev=max_nfev,
        verbose=0,
    )
    model = MSTDBThermochemicalCorrosionModel.from_vector(base_model, mstdb, result.x)
    residuals, details = thermochemical_residuals_for_targets(model, targets)
    physics_residuals, physics_details = thermochemical_physics_residuals(model, targets)
    return {
        "model": model,
        "optimizer_result": result,
        "residuals": residuals,
        "residual_details": pd.DataFrame(details),
        "physics_residuals": physics_residuals,
        "physics_residual_details": pd.DataFrame(physics_details),
        "parameter_table": model.parameter_table(),
        "prior_strength": float(prior_strength),
    }



def validate_held_out_measurements(
    targets: pd.DataFrame,
    base_model: MoltenSaltBVModel,
    mstdb: MSTDBPair,
    held_out_measurement_ids: Iterable[str],
    *,
    prior_strength: float = 0.25,
    max_nfev: int = 1000,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Refit the thermochemical layer without selected measurements.

    The shared effective model is supplied pre-calibrated and is *not* refit here,
    so this is a conditional sensitivity test rather than an independent holdout.
    It must not be cited as out-of-sample validation of the full model stack.
    """
    held_ids = {str(item) for item in held_out_measurement_ids}
    training = targets[~targets["measurement_id"].astype(str).isin(held_ids)].copy()
    held_out = targets[targets["measurement_id"].astype(str).isin(held_ids)].copy()
    held_out = held_out[held_out.apply(mechanistic_is_supported, axis=1)].copy()
    fit = fit_thermochemical_model(
        training,
        base_model,
        mstdb,
        prior_strength=prior_strength,
        max_nfev=max_nfev,
    )
    model = fit["model"]
    rows: list[dict[str, Any]] = []
    for _, row in held_out.iterrows():
        prediction = float(model.predict_response(row))
        item = row.to_dict()
        item["thermochemical_prediction"] = prediction
        item["constraint_pass"] = _constraint_pass(row, prediction)
        target_mid = row.get("target_mid")
        if _is_positive_finite(target_mid) and _is_positive_finite(prediction):
            target = float(target_mid)
            item["ln_error"] = math.log(prediction / target)
            item["factor_error"] = max(prediction / target, target / prediction)
        else:
            item["ln_error"] = np.nan
            item["factor_error"] = np.nan
        rows.append(item)
    predictions = pd.DataFrame(rows)
    direct = predictions[predictions["fit_role"].isin(MECHANISTIC_DIRECT_ROLES)].copy() if len(predictions) else predictions
    if len(direct):
        direct = direct[pd.to_numeric(direct["target_mid"], errors="coerce").gt(0)]
    errors = pd.to_numeric(direct.get("ln_error", pd.Series(dtype=float)), errors="coerce").dropna().to_numpy()
    optimizer = fit["optimizer_result"]
    metrics: dict[str, Any] = {
        "held_out_measurement_ids": sorted(held_ids),
        "n_training_rows": int(len(training)),
        "n_held_out_supported_constraints": int(len(predictions)),
        "n_held_out_direct_or_range_targets": int(len(direct)),
        "optimizer_success": bool(optimizer.success),
        "optimizer_message": str(optimizer.message),
        "optimizer_nfev": int(optimizer.nfev),
        "shared_effective_model_refit_on_training_only": False,
        "independent_holdout": False,
        "interpretation": (
            "Conditional thermochemical-layer sensitivity only: the shared effective "
            "model was calibrated using the full dataset, including these rows."
        ),
        "constraint_pass_fraction": float(predictions["constraint_pass"].mean()) if len(predictions) else None,
    }
    if len(errors):
        metrics.update(
            {
                "log_rmse_direct": float(np.sqrt(np.mean(errors**2))),
                "median_factor_error_direct": float(np.exp(np.median(np.abs(errors)))),
                "within_factor_2_direct": float(np.mean(np.abs(errors) <= math.log(2.0))),
                "within_factor_5_direct": float(np.mean(np.abs(errors) <= math.log(5.0))),
            }
        )
    return predictions, metrics

def compare_thermochemical_models(
    targets: pd.DataFrame,
    thermochemical_model: MSTDBThermochemicalCorrosionModel,
    effective_model: MoltenSaltBVModel,
    reduced_model: Any | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for _, row in targets.iterrows():
        if mechanistic_scope_reason(row) != "supported":
            continue
        thermo = float(thermochemical_model.predict_response(row))
        effective = float(effective_model.predict_response(row))
        reduced = float(reduced_model.predict_response(row)) if reduced_model is not None else np.nan
        item = row.to_dict()
        item.update(
            {
                "effective_prediction": effective,
                "reduced_mechanistic_prediction": reduced,
                "thermochemical_prediction": thermo,
                "effective_constraint_pass": _constraint_pass(row, effective),
                "reduced_mechanistic_constraint_pass": _constraint_pass(row, reduced) if _is_positive_finite(reduced) else False,
                "thermochemical_constraint_pass": _constraint_pass(row, thermo),
            }
        )
        mid = row.get("target_mid")
        if _is_positive_finite(mid):
            target = float(mid)
            for label, prediction in (("effective", effective), ("reduced_mechanistic", reduced), ("thermochemical", thermo)):
                item[f"{label}_factor_error"] = (
                    max(prediction / target, target / prediction) if _is_positive_finite(prediction) else np.nan
                )
                item[f"{label}_relation_aware_log_error"] = relation_aware_log_error(row, prediction)
        rows.append(item)
    comparison = pd.DataFrame(rows)
    direct = comparison[comparison["fit_role"].isin(MECHANISTIC_DIRECT_ROLES)].copy()
    direct = direct[pd.to_numeric(direct["target_mid"], errors="coerce").gt(0)]
    metrics: dict[str, Any] = {
        "n_thermochemical_constraints": int(len(comparison)),
        "n_direct_or_range_targets": int(len(direct)),
        "n_thermochemical_parameters": len(THERMOCHEMICAL_PARAMETER_SPECS),
        "mstdb": thermochemical_model.mstdb.metadata(),
    }
    columns = {
        "effective": "effective_prediction",
        "reduced_mechanistic": "reduced_mechanistic_prediction",
        "thermochemical": "thermochemical_prediction",
    }
    for label, column in columns.items():
        valid = direct[pd.to_numeric(direct[column], errors="coerce").gt(0)]
        if len(valid):
            error = np.log(valid[column].astype(float).to_numpy() / valid["target_mid"].astype(float).to_numpy())
            metrics[f"{label}_log_rmse_direct"] = float(np.sqrt(np.mean(error**2)))
            metrics[f"{label}_median_factor_error_direct"] = float(np.exp(np.median(np.abs(error))))
            metrics[f"{label}_within_factor_2_direct"] = float(np.mean(np.abs(error) <= math.log(2.0)))
            metrics[f"{label}_within_factor_5_direct"] = float(np.mean(np.abs(error) <= math.log(5.0)))
            relation_aware = np.asarray(
                [relation_aware_log_error(row, float(row[column])) for _, row in valid.iterrows()],
                dtype=float,
            )
            metrics[f"{label}_relation_aware_log_rmse_direct"] = float(
                np.sqrt(np.mean(relation_aware**2))
            )
            metrics[f"{label}_relation_aware_median_factor_error_direct"] = float(
                np.exp(np.median(np.abs(relation_aware)))
            )
            metrics[f"{label}_constraint_pass_fraction"] = float(comparison[f"{label}_constraint_pass"].mean())
    if len(direct):
        metrics["thermochemical_better_than_effective_count"] = int(
            (direct["thermochemical_factor_error"] < direct["effective_factor_error"]).sum()
        )
        if "reduced_mechanistic_factor_error" in direct:
            metrics["thermochemical_better_than_reduced_count"] = int(
                (direct["thermochemical_factor_error"] < direct["reduced_mechanistic_factor_error"]).sum()
            )
    return comparison, metrics


def _json_safe(value: Any) -> Any:
    """Convert NumPy/pandas scalars and non-finite values to strict JSON values."""
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def save_thermochemical_outputs(outputs: Mapping[str, Any], output_dir: str | Path) -> None:
    output_dir = Path(output_dir)
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    outputs["comparison"].to_csv(results_dir / "thermochemical_validation_comparison.csv", index=False)
    outputs["fit"]["parameter_table"].to_csv(results_dir / "thermochemical_parameters.csv", index=False)
    outputs["fit"]["physics_residual_details"].to_csv(
        results_dir / "thermochemical_physics_constraints.csv", index=False
    )
    outputs["diagnostics"].to_csv(results_dir / "thermochemical_species_inventory.csv", index=False)
    if "case_predictions" in outputs and outputs["case_predictions"] is not None:
        outputs["case_predictions"].to_csv(
            results_dir / "thermochemical_case_predictions_all_cases.csv", index=False
        )
    with open(results_dir / "thermochemical_parameters.json", "w", encoding="utf-8") as handle:
        json.dump(_json_safe(outputs["fit"]["model"].params), handle, indent=2, allow_nan=False)
    with open(results_dir / "thermochemical_comparison_metrics.json", "w", encoding="utf-8") as handle:
        json.dump(_json_safe(outputs["metrics"]), handle, indent=2, allow_nan=False)
