"""Calibration, comparison, and validation utilities for the DRIDN model.

The dynamic model is calibrated in physically separated stages rather than with a
single high-dimensional least-squares solve.  This both reduces optimizer cost and
makes parameter identifiability easier to audit: dissolution/redox, special
interfaces, inventory/deposition, mass recession, and grain-boundary penetration
are constrained by the response classes that actually inform them.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .dynamic_network import (
    DYNAMIC_FIT_PARAMETER_NAMES,
    DYNAMIC_PARAMETER_SPECS,
    DynamicCaseContext,
    DynamicRedoxInventoryDepletionModel,
)
from .mechanistic import (
    MECHANISTIC_DIRECT_ROLES,
    _constraint_pass,
    _is_positive_finite,
    mechanistic_is_supported,
    mechanistic_scope_reason,
)
from .thermochemical_fit import _json_safe, _residual_weighting, fit_thermochemical_model
from .thermochemical_model import MSTDBThermochemicalCorrosionModel

_DYNAMIC_SPEC_BY_NAME = {spec.name: spec for spec in DYNAMIC_PARAMETER_SPECS}


@dataclass(frozen=True)
class StagedOptimizerResult:
    """Small optimizer-compatible summary returned by the staged calibration."""

    success: bool
    message: str
    cost: float
    nfev: int
    stage_results: tuple[dict[str, Any], ...]


def dynamic_prior_residuals(
    params: Mapping[str, float],
    strength: float = 0.35,
    names: Iterable[str] | None = None,
) -> np.ndarray:
    selected = set(names) if names is not None else None
    values: list[float] = []
    for spec in DYNAMIC_PARAMETER_SPECS:
        if selected is not None and spec.name not in selected:
            continue
        if spec.prior_sigma > 0.0:
            values.append(strength * (float(params[spec.name]) - spec.prior) / spec.prior_sigma)
    return np.asarray(values, dtype=float)


def build_dynamic_contexts(
    targets: pd.DataFrame,
    model: DynamicRedoxInventoryDepletionModel,
) -> list[tuple[pd.Series, DynamicCaseContext]]:
    return [
        (row, model.build_context(row))
        for _, row in targets.iterrows()
        if mechanistic_is_supported(row)
    ]


def _prediction_from_context(
    model: DynamicRedoxInventoryDepletionModel,
    row: pd.Series,
    context: DynamicCaseContext,
) -> float:
    kind = str(row.get("response_kind", ""))
    if kind == "cr_diffusion_cm2_s":
        return model.thermochemical_model.cr_diffusion_cm2_s(row)
    if kind in {"redox_acceleration_ratio", "redox_acceleration_qualitative"}:
        oxidized = model.simulate_context(model.build_context(row, redox_override="oxidizing_fef2"))
        baseline = model.simulate_context(model.build_context(row, redox_override="purified_baseline"))
        return oxidized.corrosion_rate_um_y / max(baseline.corrosion_rate_um_y, 1.0e-30)
    result = model.simulate_context(context)
    if kind == "corrosion_rate_um_y":
        return result.corrosion_rate_um_y
    if kind == "corrosion_depth_um":
        return result.front_depth_um
    if kind == "mass_loss_mg_cm2":
        return result.mass_loss_mg_cm2
    if kind == "mass_gain_mg_cm2":
        return result.mass_gain_mg_cm2
    if kind == "salt_cr_ppm":
        return result.dissolved_ppm["Cr"]
    if kind == "salt_fe_decrease_ppm":
        return model.predict_response(row)
    if kind == "igc_depth_um":
        return result.igc_depth_um
    return np.nan


def range_aware_log_error(row: pd.Series | Mapping[str, Any], prediction: float) -> float:
    """Signed log error to the admissible experimental interval.

    Exact/direct targets are referenced to ``target_mid``.  Range targets have
    zero error inside their reported interval and are referenced to the nearest
    bound outside it.  This avoids treating an arbitrary range midpoint as more
    certain than the source data.
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
        reference = lower if pred < lower else upper
        return math.log(pred / reference)
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


def midpoint_log_error(row: pd.Series | Mapping[str, Any], prediction: float) -> float:
    if not _is_positive_finite(prediction) or not _is_positive_finite(row.get("target_mid")):
        return float("nan")
    return math.log(float(prediction) / float(row.get("target_mid")))


def dynamic_residuals_for_contexts(
    model: DynamicRedoxInventoryDepletionModel,
    contexts: list[tuple[pd.Series, DynamicCaseContext]],
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    residuals: list[float] = []
    details: list[dict[str, Any]] = []
    for row, context in contexts:
        raw_prediction = _prediction_from_context(model, row, context)
        pred = float(raw_prediction) if _is_positive_finite(raw_prediction) else 1.0e-30
        role = str(row.get("fit_role", ""))
        sigma, weight = _residual_weighting(row)
        low, mid, high = row.get("target_low"), row.get("target_mid"), row.get("target_high")
        residual: float | None = None
        if role == "direct" and _is_positive_finite(mid):
            residual = math.log(pred / float(mid)) / sigma * weight
        elif role == "range" and _is_positive_finite(low) and _is_positive_finite(high):
            error = range_aware_log_error(row, pred)
            residual = error / sigma * weight
        elif role == "upper":
            bound = high if _is_positive_finite(high) else mid
            if _is_positive_finite(bound):
                residual = max(0.0, math.log(pred / float(bound))) / sigma * weight
        elif role == "lower":
            bound = low if _is_positive_finite(low) else mid
            if _is_positive_finite(bound):
                residual = max(0.0, math.log(float(bound) / pred)) / sigma * weight
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
                "midpoint_log_error": midpoint_log_error(row, pred),
                "range_aware_log_error": range_aware_log_error(row, pred),
                "residual_ln_weighted": float(residual),
            }
        )
    return np.asarray(residuals, dtype=float), details


def dynamic_physics_residuals(
    model: DynamicRedoxInventoryDepletionModel,
    targets: pd.DataFrame,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    residuals: list[float] = []
    details: list[dict[str, Any]] = []

    m030 = targets[targets["measurement_id"].astype(str) == "M-030"]
    if len(m030):
        result = model.simulate(m030.iloc[0])
        fractions = result.deposit_fractions
        ratio = fractions["Fe"] / max(fractions["Cr"], 1.0e-30)
        margin = math.log(1.01)
        residual = max(0.0, margin - math.log(max(ratio, 1.0e-30))) / 0.08
        residuals.append(residual)
        details.append(
            {
                "constraint_id": "D-P001",
                "measurement_id": "M-030",
                "description": "FLiNaK cold-leg deposit remains Fe-rich relative to Cr",
                "value": ratio,
                "required_minimum": 1.0,
                "residual": residual,
                "constraint_pass": fractions["Fe"] > fractions["Cr"],
            }
        )

    m003 = targets[targets["measurement_id"].astype(str) == "M-003"]
    if len(m003):
        result = model.simulate(m003.iloc[0])
        fractions = result.cumulative_source_fractions
        for competitor in ("Fe", "Ni"):
            ratio = fractions["Cr"] / max(fractions[competitor], 1.0e-30)
            margin = math.log(1.05)
            residual = max(0.0, margin - math.log(max(ratio, 1.0e-30))) / 0.75
            residuals.append(residual)
            details.append(
                {
                    "constraint_id": f"D-P002-{competitor}",
                    "measurement_id": "M-003",
                    "description": f"Selective dissolution retains Cr source fraction > {competitor}",
                    "value": ratio,
                    "required_minimum": 1.05,
                    "residual": residual,
                    "constraint_pass": ratio > 1.05,
                }
            )
    return np.asarray(residuals, dtype=float), details


def _fit_stage(
    label: str,
    params: dict[str, float],
    parameter_names: tuple[str, ...],
    contexts: list[tuple[pd.Series, DynamicCaseContext]],
    thermochemical_model: MSTDBThermochemicalCorrosionModel,
    *,
    prior_strength: float,
    max_nfev: int,
    integration_steps: int,
) -> tuple[dict[str, float], dict[str, Any]]:
    if not contexts or not parameter_names:
        return params, {
            "stage": label,
            "success": True,
            "message": "No active contexts or parameters",
            "nfev": 0,
            "cost": 0.0,
            "parameters": list(parameter_names),
            "n_constraints": len(contexts),
        }
    x0 = np.asarray([params[name] for name in parameter_names], dtype=float)
    lower = np.asarray([_DYNAMIC_SPEC_BY_NAME[name].lower for name in parameter_names], dtype=float)
    upper = np.asarray([_DYNAMIC_SPEC_BY_NAME[name].upper for name in parameter_names], dtype=float)

    def objective(vector: np.ndarray) -> np.ndarray:
        candidate = dict(params)
        candidate.update({name: float(value) for name, value in zip(parameter_names, vector)})
        model = DynamicRedoxInventoryDepletionModel(
            thermochemical_model,
            candidate,
            integration_steps=integration_steps,
        )
        data_residuals, _ = dynamic_residuals_for_contexts(model, contexts)
        prior_residuals = dynamic_prior_residuals(
            candidate,
            strength=prior_strength,
            names=parameter_names,
        )
        return np.concatenate([data_residuals, prior_residuals])

    result = least_squares(
        objective,
        x0,
        bounds=(lower, upper),
        loss="soft_l1",
        f_scale=1.0,
        x_scale="jac",
        diff_step=5.0e-4,
        max_nfev=max_nfev,
        verbose=0,
    )
    updated = dict(params)
    updated.update({name: float(value) for name, value in zip(parameter_names, result.x)})
    return updated, {
        "stage": label,
        "success": bool(result.success),
        "message": str(result.message),
        "nfev": int(result.nfev),
        "cost": float(result.cost),
        "parameters": list(parameter_names),
        "n_constraints": len(contexts),
    }


def _default_dynamic_parameters() -> dict[str, float]:
    return {spec.name: float(spec.initial) for spec in DYNAMIC_PARAMETER_SPECS}


def fit_dynamic_network_model(
    targets: pd.DataFrame,
    thermochemical_model: MSTDBThermochemicalCorrosionModel,
    *,
    prior_strength: float = 0.35,
    max_nfev: int = 35,
    integration_steps: int = 24,
    initial_params: Mapping[str, float] | None = None,
    workers: int | None = None,
) -> dict[str, Any]:
    """Calibrate DRIDN using a response-informed staged optimization.

    ``workers`` is accepted for backward compatibility; the ODE solves are kept
    serial to avoid oversubscribing BLAS/LSODA in small validation datasets.
    """
    del workers
    params = _default_dynamic_parameters()
    if initial_params is not None:
        params.update(
            {
                name: float(value)
                for name, value in initial_params.items()
                if name in _DYNAMIC_SPEC_BY_NAME
            }
        )
    template = DynamicRedoxInventoryDepletionModel(
        thermochemical_model,
        params,
        integration_steps=integration_steps,
    )
    all_contexts = build_dynamic_contexts(targets, template)
    per_stage_nfev = max(6, min(int(max_nfev), 35))

    def select(predicate: Callable[[pd.Series], bool]) -> list[tuple[pd.Series, DynamicCaseContext]]:
        return [(row, context) for row, context in all_contexts if predicate(row)]

    baseline_kinds = {
        "corrosion_rate_um_y",
        "corrosion_depth_um",
        "redox_acceleration_ratio",
        "redox_acceleration_qualitative",
    }
    stages: tuple[
        tuple[str, tuple[str, ...], list[tuple[pd.Series, DynamicCaseContext]], float, int], ...
    ] = (
        (
            "baseline dissolution and redox",
            (
                "log_rate_scale",
                "log_surface_reservoir_um",
                "log_surface_replenishment_y_inv",
                "surface_availability_exponent",
                "log_redox_relaxation_y_inv",
                "redox_buffer_retention",
            ),
            select(
                lambda row: str(row.get("response_kind")) in baseline_kinds
                and str(row.get("redox_class")) not in {"stressed", "impure_moisture"}
            ),
            prior_strength,
            per_stage_nfev,
        ),
        (
            "stress and impurity interfacial activation",
            (
                "log_stress_interfacial_factor",
                "log_fluoride_impurity_interfacial_factor",
            ),
            select(
                lambda row: str(row.get("response_kind")) == "corrosion_rate_um_y"
                and str(row.get("redox_class")) in {"stressed", "impure_moisture"}
                and str(row.get("salt_class")) != "chloride"
            ),
            max(prior_strength, 0.35),
            min(per_stage_nfev, 20),
        ),
        (
            "inventory and cold-leg deposition",
            (
                "log_deposition_rate_y_inv_fuel",
                "log_deposition_rate_y_inv_flinak",
                "log_inventory_scale_msre",
                "log_inventory_scale_loop",
                "log_deposit_area_scale_fuel",
                "log_deposit_area_scale_flinak",
            ),
            select(
                lambda row: str(row.get("response_kind"))
                in {"salt_cr_ppm", "salt_fe_decrease_ppm", "mass_gain_mg_cm2"}
            ),
            max(prior_strength, 0.25),
            min(per_stage_nfev, 22),
        ),
        (
            "non-congruent mass recession",
            ("log_mass_loss_scale",),
            select(
                lambda row: str(row.get("response_kind"))
                in {"mass_loss_mg_cm2", "corrosion_rate_um_y"}
            ),
            max(prior_strength, 0.35),
            min(per_stage_nfev, 35),
        ),
        (
            "grain-boundary penetration",
            ("gb_dynamic_scale",),
            select(lambda row: str(row.get("response_kind")) == "igc_depth_um"),
            max(prior_strength, 0.30),
            min(per_stage_nfev, 18),
        ),
        (
            "joint endpoint polish",
            (
                "log_rate_scale",
                "log_mass_loss_scale",
                "gb_dynamic_scale",
                "log_stress_interfacial_factor",
                "log_fluoride_impurity_interfacial_factor",
            ),
            all_contexts,
            max(prior_strength, 0.40),
            min(per_stage_nfev, 35),
        ),
    )

    summaries: list[dict[str, Any]] = []
    for label, names, contexts, stage_prior, stage_nfev in stages:
        params, summary = _fit_stage(
            label,
            params,
            names,
            contexts,
            thermochemical_model,
            prior_strength=stage_prior,
            max_nfev=stage_nfev,
            integration_steps=integration_steps,
        )
        summaries.append(summary)

    model = DynamicRedoxInventoryDepletionModel(
        thermochemical_model,
        params,
        integration_steps=max(integration_steps, 36),
    )
    final_contexts = build_dynamic_contexts(targets, model)
    residuals, details = dynamic_residuals_for_contexts(model, final_contexts)
    physics_residuals, physics_details = dynamic_physics_residuals(model, targets)
    prior_residuals = dynamic_prior_residuals(
        model.params,
        strength=prior_strength,
        names=DYNAMIC_FIT_PARAMETER_NAMES,
    )
    final_vector = np.concatenate([residuals, physics_residuals, prior_residuals])
    optimizer = StagedOptimizerResult(
        success=all(bool(item["success"]) for item in summaries),
        message="; ".join(f"{item['stage']}: {item['message']}" for item in summaries),
        cost=float(0.5 * np.dot(final_vector, final_vector)),
        nfev=int(sum(int(item["nfev"]) for item in summaries)),
        stage_results=tuple(summaries),
    )
    return {
        "model": model,
        "optimizer_result": optimizer,
        "residuals": residuals,
        "residual_details": pd.DataFrame(details),
        "physics_residuals": physics_residuals,
        "physics_residual_details": pd.DataFrame(physics_details),
        "parameter_table": model.parameter_table(),
        "contexts": final_contexts,
        "prior_strength": float(prior_strength),
        "stage_results": pd.DataFrame(summaries),
    }


def _metric_bundle(comparison: pd.DataFrame, label: str) -> dict[str, float]:
    direct = comparison[comparison["fit_role"].isin(MECHANISTIC_DIRECT_ROLES)].copy()
    direct = direct[pd.to_numeric(direct["target_mid"], errors="coerce").gt(0)]
    midpoint_errors: list[float] = []
    interval_errors: list[float] = []
    for _, row in direct.iterrows():
        prediction = float(row[f"{label}_prediction"])
        midpoint_errors.append(midpoint_log_error(row, prediction))
        interval_errors.append(range_aware_log_error(row, prediction))
    midpoint = np.asarray(midpoint_errors, dtype=float)
    interval = np.asarray(interval_errors, dtype=float)
    return {
        "log_rmse_direct": float(np.sqrt(np.mean(midpoint**2))),
        "median_factor_error_direct": float(np.exp(np.median(np.abs(midpoint)))),
        "within_factor_2_direct": float(np.mean(np.abs(midpoint) <= math.log(2.0))),
        "within_factor_5_direct": float(np.mean(np.abs(midpoint) <= math.log(5.0))),
        "range_aware_log_rmse_direct": float(np.sqrt(np.mean(interval**2))),
        "range_aware_median_factor_error_direct": float(np.exp(np.median(np.abs(interval)))),
        "range_aware_within_factor_2_direct": float(np.mean(np.abs(interval) <= math.log(2.0))),
        "constraint_pass_fraction": float(comparison[f"{label}_constraint_pass"].mean()),
    }


def compare_dynamic_models(
    targets: pd.DataFrame,
    dynamic_model: DynamicRedoxInventoryDepletionModel,
    thermochemical_model: MSTDBThermochemicalCorrosionModel,
    reduced_model: Any,
    effective_model: Any,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    labels = ("effective", "reduced_mechanistic", "thermochemical", "dynamic_network")
    for _, row in targets.iterrows():
        if mechanistic_scope_reason(row) != "supported":
            continue
        predictions = {
            "effective": float(effective_model.predict_response(row)),
            "reduced_mechanistic": float(reduced_model.predict_response(row)),
            "thermochemical": float(thermochemical_model.predict_response(row)),
            "dynamic_network": float(dynamic_model.predict_response(row)),
        }
        item = row.to_dict()
        for label, prediction in predictions.items():
            item[f"{label}_prediction"] = prediction
            item[f"{label}_constraint_pass"] = _constraint_pass(row, prediction)
            item[f"{label}_midpoint_log_error"] = midpoint_log_error(row, prediction)
            item[f"{label}_range_aware_log_error"] = range_aware_log_error(row, prediction)
            if _is_positive_finite(row.get("target_mid")) and _is_positive_finite(prediction):
                target = float(row.get("target_mid"))
                item[f"{label}_factor_error"] = max(prediction / target, target / prediction)
            else:
                item[f"{label}_factor_error"] = np.nan
        rows.append(item)
    comparison = pd.DataFrame(rows)
    direct = comparison[comparison["fit_role"].isin(MECHANISTIC_DIRECT_ROLES)].copy()
    direct = direct[pd.to_numeric(direct["target_mid"], errors="coerce").gt(0)]
    metrics: dict[str, Any] = {
        "n_supported_constraints": int(len(comparison)),
        "n_direct_or_range_targets": int(len(direct)),
        "n_dynamic_parameters": len(DYNAMIC_PARAMETER_SPECS),
        "n_fitted_dynamic_parameters": len(DYNAMIC_FIT_PARAMETER_NAMES),
        "integration_steps": dynamic_model.integration_steps,
        "midpoint_metric_note": "Range targets are referenced to their reported midpoint for legacy comparability.",
        "range_aware_metric_note": "Range targets contribute zero error inside the reported interval and error to the nearest bound outside it.",
    }
    for label in labels:
        for key, value in _metric_bundle(comparison, label).items():
            metrics[f"{label}_{key}"] = value
    if len(direct):
        metrics["dynamic_better_than_effective_count"] = int(
            (direct["dynamic_network_factor_error"] < direct["effective_factor_error"]).sum()
        )
        metrics["dynamic_better_than_reduced_count"] = int(
            (direct["dynamic_network_factor_error"] < direct["reduced_mechanistic_factor_error"]).sum()
        )
        metrics["dynamic_better_than_thermochemical_count"] = int(
            (direct["dynamic_network_factor_error"] < direct["thermochemical_factor_error"]).sum()
        )
    return comparison, metrics


def validate_dynamic_held_out_measurements(
    targets: pd.DataFrame,
    thermochemical_model: MSTDBThermochemicalCorrosionModel,
    held_out_measurement_ids: Iterable[str],
    *,
    prior_strength: float = 0.35,
    max_nfev: int = 18,
    integration_steps: int = 18,
    refit_thermochemical: bool = True,
    initial_params: Mapping[str, float] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    held_ids = {str(item) for item in held_out_measurement_ids}
    training = targets[~targets["measurement_id"].astype(str).isin(held_ids)].copy()
    held_out = targets[targets["measurement_id"].astype(str).isin(held_ids)].copy()
    held_out = held_out[held_out.apply(mechanistic_is_supported, axis=1)].copy()
    training_thermochemical = thermochemical_model
    thermochemical_refit_result = None
    if refit_thermochemical:
        thermochemical_refit_result = fit_thermochemical_model(
            training,
            thermochemical_model.base_model,
            thermochemical_model.mstdb,
            prior_strength=0.25,
            max_nfev=max(250, max_nfev),
        )
        training_thermochemical = thermochemical_refit_result["model"]
    fit = fit_dynamic_network_model(
        training,
        training_thermochemical,
        prior_strength=prior_strength,
        max_nfev=max_nfev,
        integration_steps=integration_steps,
        initial_params=initial_params,
    )
    model = fit["model"]
    rows: list[dict[str, Any]] = []
    for _, row in held_out.iterrows():
        prediction = float(model.predict_response(row))
        item = row.to_dict()
        item["dynamic_network_prediction"] = prediction
        item["constraint_pass"] = _constraint_pass(row, prediction)
        item["midpoint_log_error"] = midpoint_log_error(row, prediction)
        item["range_aware_log_error"] = range_aware_log_error(row, prediction)
        item["factor_error"] = (
            max(prediction / float(row.get("target_mid")), float(row.get("target_mid")) / prediction)
            if _is_positive_finite(row.get("target_mid")) and _is_positive_finite(prediction)
            else np.nan
        )
        rows.append(item)
    predictions = pd.DataFrame(rows)
    direct = predictions[predictions["fit_role"].isin(MECHANISTIC_DIRECT_ROLES)].copy() if len(predictions) else pd.DataFrame()
    direct = direct[pd.to_numeric(direct.get("target_mid"), errors="coerce").gt(0)] if len(direct) else direct
    optimizer = fit["optimizer_result"]
    metrics: dict[str, Any] = {
        "held_out_measurement_ids": sorted(held_ids),
        "n_training_rows": int(len(training)),
        "n_held_out_supported_constraints": int(len(predictions)),
        "n_held_out_direct_or_range_targets": int(len(direct)),
        "optimizer_success": bool(optimizer.success),
        "optimizer_message": str(optimizer.message),
        "optimizer_nfev": int(optimizer.nfev),
        "thermochemical_refit_for_holdout": bool(refit_thermochemical),
        "thermochemical_refit_optimizer_success": (
            bool(thermochemical_refit_result["optimizer_result"].success)
            if thermochemical_refit_result is not None
            else None
        ),
        "thermochemical_refit_optimizer_nfev": (
            int(thermochemical_refit_result["optimizer_result"].nfev)
            if thermochemical_refit_result is not None
            else None
        ),
        "constraint_pass_fraction": float(predictions["constraint_pass"].mean()) if len(predictions) else None,
    }
    if len(direct):
        midpoint = direct["midpoint_log_error"].astype(float).to_numpy()
        interval = direct["range_aware_log_error"].astype(float).to_numpy()
        metrics.update(
            {
                "log_rmse_direct": float(np.sqrt(np.mean(midpoint**2))),
                "median_factor_error_direct": float(np.exp(np.median(np.abs(midpoint)))),
                "within_factor_2_direct": float(np.mean(np.abs(midpoint) <= math.log(2.0))),
                "within_factor_5_direct": float(np.mean(np.abs(midpoint) <= math.log(5.0))),
                "range_aware_log_rmse_direct": float(np.sqrt(np.mean(interval**2))),
                "range_aware_median_factor_error_direct": float(np.exp(np.median(np.abs(interval)))),
            }
        )
    return predictions, metrics


def save_dynamic_outputs(outputs: Mapping[str, Any], output_dir: str | Path) -> None:
    output_dir = Path(output_dir)
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    outputs["comparison"].to_csv(results_dir / "three_model_validation_comparison.csv", index=False)
    outputs["comparison"].to_csv(
        results_dir / "three_mechanistic_models_validation_comparison.csv",
        index=False,
    )
    outputs["fit"]["parameter_table"].to_csv(results_dir / "dynamic_network_parameters.csv", index=False)
    outputs["fit"]["residual_details"].to_csv(results_dir / "dynamic_network_residual_details.csv", index=False)
    outputs["fit"]["physics_residual_details"].to_csv(results_dir / "dynamic_network_physics_constraints.csv", index=False)
    outputs["fit"]["stage_results"].to_csv(results_dir / "dynamic_network_calibration_stages.csv", index=False)
    endpoint_audit = outputs["fit"].get("endpoint_state_audit")
    if isinstance(endpoint_audit, pd.DataFrame):
        endpoint_audit.to_csv(results_dir / "dynamic_network_endpoint_state_audit.csv", index=False)
        endpoint_audit.to_csv(results_dir / "dynamic_network_species_inventory.csv", index=False)
    with open(results_dir / "dynamic_network_parameters.json", "w", encoding="utf-8") as handle:
        json.dump(_json_safe(outputs["fit"]["model"].params), handle, indent=2, allow_nan=False)
    safe_metrics = _json_safe(outputs["metrics"])
    with open(results_dir / "three_model_comparison_metrics.json", "w", encoding="utf-8") as handle:
        json.dump(safe_metrics, handle, indent=2, allow_nan=False)
    with open(
        results_dir / "three_mechanistic_models_comparison_metrics.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(safe_metrics, handle, indent=2, allow_nan=False)
