#!/usr/bin/env python3
"""Validate and optionally refit the DRIDN state-evolving corrosion model."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from msr_corrosion_bv.calibrate import load_model_from_json
from msr_corrosion_bv.advanced_context import merge_advanced_case_context
from msr_corrosion_bv.dynamic_network import (
    DYNAMIC_FIXED_PARAMETER_NAMES,
    DynamicRedoxInventoryDepletionModel,
    ELEMENTS,
)
from msr_corrosion_bv.dynamic_network_fit import (
    StagedOptimizerResult,
    build_dynamic_contexts,
    compare_dynamic_models,
    dynamic_physics_residuals,
    dynamic_residuals_for_contexts,
    fit_dynamic_network_model,
    save_dynamic_outputs,
    validate_dynamic_held_out_measurements,
)
from msr_corrosion_bv.mechanistic import MechanisticCorrosionModel, mechanistic_is_supported
from msr_corrosion_bv.mstdb import MSTDBPair, resolve_mstdb_directory
from msr_corrosion_bv.publication_plots import make_all_publication_figures
from msr_corrosion_bv.thermochemical_model import MSTDBThermochemicalCorrosionModel


CPP_PARITY_MEASUREMENT_IDS = (
    "M-003",
    "M-005",
    "M-014",
    "M-018",
    "M-027",
    "M-029",
    "M-030",
    "M-038",
    "M-041",
)


def cpp_parity_cases(
    targets: pd.DataFrame,
    thermochemical: MSTDBThermochemicalCorrosionModel,
    model: DynamicRedoxInventoryDepletionModel,
) -> pd.DataFrame:
    """Freeze explicit C++ context inputs and Python endpoint references."""
    indexed = targets.set_index("measurement_id", drop=False)
    records: list[dict[str, Any]] = []
    for measurement_id in CPP_PARITY_MEASUREMENT_IDS:
        row = indexed.loc[measurement_id]
        context = model.build_context(row)
        endpoint = model.simulate_context(context)
        record: dict[str, Any] = {
            "measurement_id": measurement_id,
            "response_kind": row.get("response_kind"),
            "material": context.material,
            "salt_class": context.salt_class,
            "redox_class": context.redox_class,
            "position_class": context.position_class,
            "temperature_K": context.temperature_K,
            "cold_temperature_K": context.cold_temperature_K,
            "delta_T_C": context.delta_T_C,
            "flow_factor": context.flow_factor,
            "time_years": context.time_years,
            "density_g_cm3": context.density_g_cm3,
            "cr_fraction_ratio": context.cr_fraction_ratio,
            "inventory_scale": context.inventory_scale,
            "explicit_inventory_scale": context.explicit_inventory_scale,
            "deposition_closure": context.deposition_closure,
            "area_to_salt_mass_cm2_g": context.area_to_salt_mass_cm2_g,
            "inventory_coupling_factor": context.inventory_coupling_factor,
            "deposit_area_factor": context.deposit_area_factor,
            "selectivity_scale": context.selectivity_scale,
            "product_floor_ppm": context.product_floor_ppm,
            "initial_dissolved_Cr_ppm": context.initial_dissolved_ppm[0],
            "initial_dissolved_Fe_ppm": context.initial_dissolved_ppm[1],
            "initial_dissolved_Ni_ppm": context.initial_dissolved_ppm[2],
            "redox_shift_initial": context.redox_shift_initial,
            "log_charge_base_no_redox": context.log_charge_base_no_redox,
            "mass_transfer_rate_um_y": context.mass_transfer_rate_um_y,
            "inventory_capacity_ppm": context.inventory_capacity_ppm,
            "mass_loss_fraction": context.mass_loss_fraction,
            "cr_diffusion_cm2_s": context.cr_diffusion_cm2_s,
            "front_damage_multiplier": context.front_damage_multiplier,
            "gb_length_multiplier": context.gb_length_multiplier,
            "transient_redox": context.transient_redox,
            "stress_interfacial_activation": context.stress_interfacial_activation,
            "fluoride_impurity_interfacial_activation": (
                context.fluoride_impurity_interfacial_activation
            ),
            "chloride_salt": context.chloride_salt,
            "thermochemical_response_prediction": float(thermochemical.predict_response(row)),
            "dridn_response_prediction": float(model.predict_response(row)),
            "front_depth_um": endpoint.front_depth_um,
            "mass_recession_um": endpoint.mass_recession_um,
            "mass_loss_mg_cm2": endpoint.mass_loss_mg_cm2,
            "mass_gain_mg_cm2": endpoint.mass_gain_mg_cm2,
            "igc_depth_um": endpoint.igc_depth_um,
            "corrosion_rate_um_y": endpoint.corrosion_rate_um_y,
            "redox_shift_endpoint": endpoint.redox_shift,
            "mass_balance_relative_error": endpoint.mass_balance_relative_error,
        }
        for index, element in enumerate(ELEMENTS):
            record[f"mass_fraction_{element}"] = context.mass_fractions[index]
            record[f"log_exchange_offset_{element}"] = context.log_exchange_offsets[index]
            record[f"affinity_baseline_{element}"] = context.affinity_baseline[index]
            record[f"cold_capture_fraction_{element}"] = context.cold_capture_fraction[index]
            record[f"dissolved_{element}_ppm"] = endpoint.dissolved_ppm[element]
            record[f"cumulative_source_{element}_ppm"] = endpoint.cumulative_source_ppm[element]
            record[f"coupon_deposit_{element}_mg_cm2"] = endpoint.deposit_mg_cm2[element]
            record[f"bulk_captured_{element}_ppm"] = endpoint.bulk_captured_ppm[element]
            record[f"surface_{element}_availability"] = endpoint.surface_availability[element]
        records.append(record)
    base_context = model.build_context(indexed.loc["M-014"])
    boundary_context = replace(
        base_context,
        measurement_id="boundary_zero_elements",
        mass_fractions=np.asarray([1.0, 0.0, 0.0], dtype=float),
        initial_dissolved_ppm=np.asarray(
            [base_context.initial_dissolved_ppm[0], 0.0, 0.0], dtype=float
        ),
        cr_fraction_ratio=1.0 / 0.07,
    )
    boundary_endpoint = model.simulate_context(boundary_context)
    boundary = dict(records[CPP_PARITY_MEASUREMENT_IDS.index("M-014")])
    boundary.update(
        {
            "measurement_id": "boundary_zero_elements",
            "response_kind": "dridn_context_boundary",
            "material": "synthetic_cr_only",
            "cr_fraction_ratio": boundary_context.cr_fraction_ratio,
            "initial_dissolved_Fe_ppm": 0.0,
            "initial_dissolved_Ni_ppm": 0.0,
            "thermochemical_response_prediction": "",
            "dridn_response_prediction": "",
            "front_depth_um": boundary_endpoint.front_depth_um,
            "mass_recession_um": boundary_endpoint.mass_recession_um,
            "mass_loss_mg_cm2": boundary_endpoint.mass_loss_mg_cm2,
            "mass_gain_mg_cm2": boundary_endpoint.mass_gain_mg_cm2,
            "igc_depth_um": boundary_endpoint.igc_depth_um,
            "corrosion_rate_um_y": boundary_endpoint.corrosion_rate_um_y,
            "redox_shift_endpoint": boundary_endpoint.redox_shift,
            "mass_balance_relative_error": boundary_endpoint.mass_balance_relative_error,
        }
    )
    for index, element in enumerate(ELEMENTS):
        boundary[f"mass_fraction_{element}"] = boundary_context.mass_fractions[index]
        boundary[f"dissolved_{element}_ppm"] = boundary_endpoint.dissolved_ppm[element]
        boundary[f"cumulative_source_{element}_ppm"] = boundary_endpoint.cumulative_source_ppm[element]
        boundary[f"coupon_deposit_{element}_mg_cm2"] = boundary_endpoint.deposit_mg_cm2[element]
        boundary[f"bulk_captured_{element}_ppm"] = boundary_endpoint.bulk_captured_ppm[element]
        boundary[f"surface_{element}_availability"] = boundary_endpoint.surface_availability[element]
    records.append(boundary)
    return pd.DataFrame(records)


def _validation_only_fit(
    targets: pd.DataFrame,
    thermochemical: MSTDBThermochemicalCorrosionModel,
    params: dict[str, float],
    integration_steps: int,
) -> dict[str, Any]:
    model = DynamicRedoxInventoryDepletionModel(
        thermochemical,
        params,
        integration_steps=integration_steps,
    )
    contexts = build_dynamic_contexts(targets, model)
    residuals, details = dynamic_residuals_for_contexts(model, contexts)
    physics_residuals, physics_details = dynamic_physics_residuals(model, targets)
    optimizer = StagedOptimizerResult(
        success=True,
        message="Validated checked-in calibrated parameters; no refit requested",
        cost=float(0.5 * np.dot(residuals, residuals)),
        nfev=0,
        stage_results=(),
    )
    return {
        "model": model,
        "optimizer_result": optimizer,
        "residuals": residuals,
        "residual_details": pd.DataFrame(details),
        "physics_residuals": physics_residuals,
        "physics_residual_details": pd.DataFrame(physics_details),
        "parameter_table": model.parameter_table(),
        "contexts": contexts,
        "prior_strength": None,
        "stage_results": pd.DataFrame(
            [
                {
                    "stage": "checked-in calibrated parameter validation",
                    "success": True,
                    "message": optimizer.message,
                    "nfev": 0,
                    "cost": optimizer.cost,
                    "parameters": "all",
                    "n_constraints": len(contexts),
                }
            ]
        ),
    }


def endpoint_state_audit(
    targets: pd.DataFrame,
    model: DynamicRedoxInventoryDepletionModel,
) -> pd.DataFrame:
    """Audit endpoint states and elemental conservation for every supported case."""
    records: list[dict[str, Any]] = []
    supported = targets[targets.apply(mechanistic_is_supported, axis=1)].copy()
    for _, row in supported.iterrows():
        result = model.simulate(row)
        record: dict[str, Any] = {
            "measurement_id": row.get("measurement_id"),
            "case_id": row.get("case_id"),
            "response_kind": row.get("response_kind"),
            "front_depth_um": result.front_depth_um,
            "mass_recession_um": result.mass_recession_um,
            "mass_loss_mg_cm2": result.mass_loss_mg_cm2,
            "mass_gain_mg_cm2": result.mass_gain_mg_cm2,
            "igc_depth_um": result.igc_depth_um,
            "corrosion_rate_um_y": result.corrosion_rate_um_y,
            "redox_log_shift": result.redox_shift,
            "mass_balance_relative_error": result.mass_balance_relative_error,
        }
        for element in ELEMENTS:
            record[f"dissolved_{element}_ppm"] = result.dissolved_ppm[element]
            record[f"cumulative_source_{element}_ppm"] = result.cumulative_source_ppm[element]
            record[f"coupon_deposit_{element}_mg_cm2"] = result.deposit_mg_cm2[element]
            record[f"bulk_captured_{element}_ppm"] = result.bulk_captured_ppm[element]
            record[f"surface_{element}_availability"] = result.surface_availability[element]
        records.append(record)
    return pd.DataFrame(records)


def parameter_bound_audit(parameter_table: pd.DataFrame) -> dict[str, Any]:
    """Identify calibrated and fixed quantities that lie on admissible bounds."""
    lower_names: list[str] = []
    upper_names: list[str] = []
    for _, row in parameter_table.iterrows():
        value = float(row["value"])
        lower = float(row["lower"])
        upper = float(row["upper"])
        tolerance = 1.0e-8 * max(1.0, abs(lower), abs(upper), abs(upper - lower))
        if abs(value - lower) <= tolerance:
            lower_names.append(str(row["parameter"]))
        if abs(value - upper) <= tolerance:
            upper_names.append(str(row["parameter"]))
    fixed = set(DYNAMIC_FIXED_PARAMETER_NAMES)
    fitted_lower = [name for name in lower_names if name not in fixed]
    fitted_upper = [name for name in upper_names if name not in fixed]
    fixed_lower = [name for name in lower_names if name in fixed]
    fixed_upper = [name for name in upper_names if name in fixed]
    return {
        "fitted_parameters_at_lower_bound": fitted_lower,
        "fitted_parameters_at_upper_bound": fitted_upper,
        "fixed_parameters_at_lower_bound": fixed_lower,
        "fixed_parameters_at_upper_bound": fixed_upper,
        "n_fitted_parameters_at_any_bound": len(set(fitted_lower + fitted_upper)),
        "n_fixed_parameters_at_any_bound": len(set(fixed_lower + fixed_upper)),
    }


def numerical_convergence_check(
    targets: pd.DataFrame,
    thermochemical: MSTDBThermochemicalCorrosionModel,
    params: dict[str, float],
    steps: tuple[int, ...] = (30, 60, 120, 240),
) -> dict[str, Any]:
    supported = targets[targets.apply(mechanistic_is_supported, axis=1)].copy()
    predictions: dict[int, np.ndarray] = {}
    for count in steps:
        model = DynamicRedoxInventoryDepletionModel(
            thermochemical,
            params,
            integration_steps=count,
        )
        predictions[count] = np.asarray(
            [float(model.predict_response(row)) for _, row in supported.iterrows()],
            dtype=float,
        )
    reference = predictions[steps[-1]]
    comparisons: dict[str, Any] = {}
    for count in steps[:-1]:
        denominator = np.maximum(np.abs(reference), 1.0e-14)
        relative = np.abs(predictions[count] - reference) / denominator
        comparisons[f"{count}_vs_{steps[-1]}"] = {
            "maximum_relative_difference": float(np.max(relative)),
            "median_relative_difference": float(np.median(relative)),
            "p95_relative_difference": float(np.quantile(relative, 0.95)),
        }
    return {
        "integration_steps_evaluated": list(steps),
        "reference_integration_steps": steps[-1],
        "n_supported_predictions": int(len(supported)),
        "comparisons": comparisons,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate the DRIDN dynamic corrosion model.")
    parser.add_argument("--mstdb-dir", type=Path, default=None)
    parser.add_argument("--targets", type=Path, default=ROOT / "data" / "targets.csv")
    parser.add_argument("--case-context", type=Path, default=ROOT / "results" / "advanced" / "advanced_case_context.csv")
    parser.add_argument("--parameters", type=Path, default=ROOT / "data" / "parameters.json")
    parser.add_argument("--reduced-parameters", type=Path, default=ROOT / "results" / "advanced" / "mechanistic_parameters.json")
    parser.add_argument("--thermochemical-parameters", type=Path, default=ROOT / "results" / "advanced" / "thermochemical_parameters.json")
    parser.add_argument("--dynamic-parameters", type=Path, default=ROOT / "results" / "advanced" / "dynamic_network_parameters.json")
    parser.add_argument("--outdir", type=Path, default=ROOT / "generated" / "advanced")
    parser.add_argument("--prior-strength", type=float, default=0.35)
    parser.add_argument("--max-nfev", type=int, default=6, help="Per-stage evaluation cap used only with --refit")
    parser.add_argument("--integration-steps", type=int, default=120)
    parser.add_argument("--fit-integration-steps", type=int, default=16)
    parser.add_argument("--refit", action="store_true", help="Run the staged DRIDN calibration before validation")
    parser.add_argument("--run-holdout", action="store_true", help="Run the nested modern-family holdout stress test")
    parser.add_argument("--skip-plots", action="store_true")
    args = parser.parse_args()

    targets = merge_advanced_case_context(pd.read_csv(args.targets), args.case_context)
    effective = load_model_from_json(args.parameters)
    reduced = MechanisticCorrosionModel(
        effective,
        json.loads(args.reduced_parameters.read_text(encoding="utf-8")),
    )
    mstdb = MSTDBPair.from_calibration_directory(resolve_mstdb_directory(args.mstdb_dir))
    thermochemical = MSTDBThermochemicalCorrosionModel(
        effective,
        mstdb,
        json.loads(args.thermochemical_parameters.read_text(encoding="utf-8")),
    )
    initial_dynamic = json.loads(args.dynamic_parameters.read_text(encoding="utf-8"))

    if args.refit:
        fit = fit_dynamic_network_model(
            targets,
            thermochemical,
            prior_strength=args.prior_strength,
            max_nfev=args.max_nfev,
            integration_steps=args.fit_integration_steps,
            initial_params=initial_dynamic,
        )
        # Re-evaluate the fitted parameters with the requested production resolution.
        fit["model"] = DynamicRedoxInventoryDepletionModel(
            thermochemical,
            fit["model"].params,
            integration_steps=args.integration_steps,
        )
        fit["contexts"] = build_dynamic_contexts(targets, fit["model"])
        fit["residuals"], details = dynamic_residuals_for_contexts(fit["model"], fit["contexts"])
        fit["residual_details"] = pd.DataFrame(details)
        fit["physics_residuals"], physics_details = dynamic_physics_residuals(fit["model"], targets)
        fit["physics_residual_details"] = pd.DataFrame(physics_details)
        fit["parameter_table"] = fit["model"].parameter_table()
        calibration_mode = "staged refit"
    else:
        fit = _validation_only_fit(
            targets,
            thermochemical,
            initial_dynamic,
            args.integration_steps,
        )
        calibration_mode = "checked-in calibrated parameters"

    comparison, metrics = compare_dynamic_models(
        targets,
        fit["model"],
        thermochemical,
        reduced,
        effective,
    )
    optimizer = fit["optimizer_result"]
    physics = fit["physics_residual_details"]
    convergence = numerical_convergence_check(
        targets,
        thermochemical,
        fit["model"].params,
    )
    endpoint_audit = endpoint_state_audit(targets, fit["model"])
    fit["endpoint_state_audit"] = endpoint_audit
    balance_values = endpoint_audit["mass_balance_relative_error"].astype(float).to_numpy()
    bound_audit = parameter_bound_audit(fit["parameter_table"])
    metrics.update(
        {
            "optimizer_success": bool(optimizer.success),
            "optimizer_message": str(optimizer.message),
            "optimizer_cost": float(optimizer.cost),
            "optimizer_nfev": int(optimizer.nfev),
            "prior_strength": float(args.prior_strength) if args.refit else None,
            "calibration_mode": calibration_mode,
            "dynamic_model_name": "Dynamic Redox-Inventory-Depletion Network (DRIDN)",
            "ode_integrator": "SciPy solve_ivp/LSODA with nonnegative state projection in rate evaluations",
            "validation_interpretation": (
                "Calibrated like-for-like endpoint benchmark with explicit transient state evolution. "
                "The interval-aware score respects reported experimental ranges; numerical convergence, "
                "species ordering, and elemental inventory closure are audited independently."
            ),
            "dynamic_state_variables": [
                "surface Cr/Fe/Ni availability",
                "dissolved Cr/Fe/Ni inventory",
                "cumulative Cr/Fe/Ni dissolution source",
                "redox buffer coordinate",
                "reaction-front depth",
                "equivalent mass recession",
                "species-resolved measured-coupon deposit",
                "species-resolved uninstrumented/bulk capture",
                "grain-boundary penetration squared",
            ],
            "fixed_dynamic_parameters": list(DYNAMIC_FIXED_PARAMETER_NAMES),
            "endpoint_mass_balance": {
                "n_cases": int(len(balance_values)),
                "maximum_relative_error": float(np.max(balance_values)),
                "median_relative_error": float(np.median(balance_values)),
                "p95_relative_error": float(np.quantile(balance_values, 0.95)),
            },
            "parameter_bound_audit": bound_audit,
            "physics_constraints": {
                "n_constraints": int(len(physics)),
                "pass_fraction": float(physics["constraint_pass"].mean()) if len(physics) else None,
                "constraints": physics.to_dict(orient="records"),
            },
            "numerical_convergence": convergence,
        }
    )

    if args.run_holdout:
        modern_ids = {f"M-{index:03d}" for index in range(29, 39)}
        held_out_predictions, held_out_metrics = validate_dynamic_held_out_measurements(
            targets,
            thermochemical,
            modern_ids,
            prior_strength=args.prior_strength,
            max_nfev=args.max_nfev,
            integration_steps=args.fit_integration_steps,
            initial_params=None,
        )
        metrics["external_modern_family_holdout"] = held_out_metrics
        results_dir = args.outdir / "results"
        results_dir.mkdir(parents=True, exist_ok=True)
        held_out_predictions.to_csv(results_dir / "dynamic_network_external_modern_validation.csv", index=False)
        (results_dir / "dynamic_network_external_modern_validation_metrics.json").write_text(
            json.dumps(held_out_metrics, indent=2),
            encoding="utf-8",
        )

    results_dir = args.outdir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    convergence_rows: list[dict[str, Any]] = []
    for comparison_name, values in convergence["comparisons"].items():
        coarse_steps = int(comparison_name.split("_vs_")[0])
        convergence_rows.append(
            {
                "coarse_integration_steps": coarse_steps,
                "reference_integration_steps": convergence["reference_integration_steps"],
                **values,
            }
        )
    pd.DataFrame(convergence_rows).to_csv(
        results_dir / "dynamic_network_numerical_convergence.csv",
        index=False,
    )

    outputs = {"fit": fit, "comparison": comparison, "metrics": metrics}
    save_dynamic_outputs(outputs, args.outdir)
    cpp_parity_cases(targets, thermochemical, fit["model"]).to_csv(
        results_dir / "cpp_parity_cases.csv", index=False
    )

    if not args.skip_plots:
        manifest = make_all_publication_figures(
            comparison,
            metrics,
            fit["model"],
            targets,
            args.outdir / "figures",
            args.outdir / "results",
        )
        metrics["publication_figures"] = manifest
    else:
        # Preserve the auditable figure manifest when running a validation-only
        # or holdout job after publication figures have already been generated.
        manifest_path = args.outdir / "results" / "publication_figure_manifest.json"
        if manifest_path.exists():
            metrics["publication_figures"] = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )
    metrics_text = json.dumps(metrics, indent=2)
    (results_dir / "three_model_comparison_metrics.json").write_text(
        metrics_text,
        encoding="utf-8",
    )
    (results_dir / "three_mechanistic_models_comparison_metrics.json").write_text(
        metrics_text,
        encoding="utf-8",
    )

    print("DRIDN dynamic validation complete")
    print(f"Calibration mode: {calibration_mode}")
    print(f"Optimizer/validation status: {optimizer.success} ({optimizer.message})")
    for label in ("reduced_mechanistic", "thermochemical", "dynamic_network"):
        print(
            f"{label}: midpoint log RMSE={metrics[f'{label}_log_rmse_direct']:.6f}; "
            f"interval-aware log RMSE={metrics[f'{label}_range_aware_log_rmse_direct']:.6f}; "
            f"median factor={metrics[f'{label}_median_factor_error_direct']:.6f}; "
            f"within factor 2={metrics[f'{label}_within_factor_2_direct']:.1%}; "
            f"constraint pass={metrics[f'{label}_constraint_pass_fraction']:.1%}"
        )
    print(
        "Numerical convergence, 120 vs 240 steps: "
        f"max relative difference={convergence['comparisons']['120_vs_240']['maximum_relative_difference']:.3e}"
    )


if __name__ == "__main__":
    main()
