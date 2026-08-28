#!/usr/bin/env python3
"""Fit and validate the MSTDB-TC-grounded MSR corrosion model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["svg.fonttype"] = "none"
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from msr_corrosion_bv.calibrate import load_model_from_json
from msr_corrosion_bv.advanced_context import merge_advanced_case_context
from msr_corrosion_bv.mechanistic import MechanisticCorrosionModel, mechanistic_is_supported
from msr_corrosion_bv.mstdb import MSTDBPair, ThermochimicaRunner, resolve_mstdb_directory
from msr_corrosion_bv.thermochemical import (
    THERMOCHEMICAL_PARAMETER_SPECS,
    compare_thermochemical_models,
    fit_thermochemical_model,
    save_thermochemical_outputs,
    validate_held_out_measurements,
)

THERMOCHIMICA_REFERENCE_COMMIT = "0c35c8d7d1cf2084b4e2ca5d6608f7dcdf60adad"


def make_comparison_plot(comparison: pd.DataFrame, figure_dir: Path) -> None:
    direct = comparison[comparison["fit_role"].isin(["direct", "range"])].copy()
    direct = direct[pd.to_numeric(direct["target_mid"], errors="coerce").gt(0)]
    if direct.empty:
        return
    figure_dir.mkdir(parents=True, exist_ok=True)
    target = direct["target_mid"].astype(float).to_numpy()
    x = np.arange(len(direct), dtype=float)
    series = [
        ("effective BV", direct["effective_prediction"].astype(float).to_numpy() / target, "o", -0.20),
        ("reduced mechanistic", direct["reduced_mechanistic_prediction"].astype(float).to_numpy() / target, "x", 0.0),
        ("MSTDB thermochemical", direct["thermochemical_prediction"].astype(float).to_numpy() / target, "^", 0.20),
    ]
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    ax.axhline(1.0, linewidth=1.4, label="target")
    ax.axhline(2.0, linestyle="--", linewidth=0.9, label="factor 2")
    ax.axhline(0.5, linestyle="--", linewidth=0.9)
    for label, values, marker, shift in series:
        ax.scatter(x + shift, values, marker=marker, label=label)
    ax.set_yscale("log")
    ax.set_ylabel("Prediction / experimental target midpoint")
    ax.set_xlabel("Validation measurement")
    ax.set_title("MSTDB-TC thermochemical model compared with existing corrosion models")
    ax.set_xticks(x)
    ax.set_xticklabels(direct["measurement_id"].astype(str), rotation=60, ha="right")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(figure_dir / "fig09_thermochemical_model_comparison.png", dpi=220)
    fig.savefig(figure_dir / "fig09_thermochemical_model_comparison.svg")
    plt.close(fig)


def native_thermochimica_smoke_test(
    mstdb: MSTDBPair,
    thermochimica_executable: str | None = None,
) -> dict[str, object]:
    """Run a minimal FLiBe + trace-Cr equilibrium when a native executable exists.

    The calculation verifies executable/database compatibility and the JSON-output
    pathway.  It does not yet replace fitted activity closures in the corrosion
    model because most validation rows lack complete initial compositions.
    """
    candidate = thermochimica_executable or "InputScriptMode"
    resolved = (
        str(Path(candidate).expanduser().resolve())
        if Path(candidate).expanduser().is_file()
        else shutil.which(candidate)
    )
    if resolved is None:
        return {
            "available": False,
            "executed": False,
            "executable": None,
            "reason": "InputScriptMode was not found",
        }
    composition = {3: 2.0, 4: 1.0, 9: 4.0, 24: 1.0e-8}  # Li2BeF4 + trace Cr
    output = ThermochimicaRunner(resolved).run_elements(
        mstdb.fluoride.path,
        923.15,
        composition,
    )
    return {
        "available": True,
        "executed": True,
        "executable": resolved,
        "temperature_K": 923.15,
        "element_moles_by_atomic_number": {str(key): value for key, value in composition.items()},
        "output_top_level_keys": [str(key) for key in list(output)[:20]],
        "output_case_count": len(output),
    }


def database_checks(mstdb: MSTDBPair, thermochimica_executable: str | None = None) -> dict[str, object]:
    fluoride = mstdb.fluoride
    chloride = mstdb.chloride
    T = 923.15
    cr_u_reaction = {
        "CrF2_L1(liq)": 1.0,
        "UF3_L1(liq)": 2.0,
        "Cr_S1(s)": -1.0,
        "UF4_L1(liq)": -2.0,
    }
    fe_u_reaction = {
        "FeF2_L1(liq)": 1.0,
        "UF3_L1(liq)": 2.0,
        "Fe_bcc(s)": -1.0,
        "UF4_L1(liq)": -2.0,
    }
    ni_u_reaction = {
        "NiF2_L1(liq)": 1.0,
        "UF3_L1(liq)": 2.0,
        "Ni_fcc(s)": -1.0,
        "UF4_L1(liq)": -2.0,
    }
    cr_fef2_reaction = {
        "CrF2_L1(liq)": 1.0,
        "Fe_bcc(s)": 1.0,
        "Cr_S1(s)": -1.0,
        "FeF2_L1(liq)": -1.0,
    }
    cr_bef2_reaction = {
        "CrF2_L1(liq)": 1.0,
        "Be_S1(s)": 1.0,
        "Cr_S1(s)": -1.0,
        "BeF2_L1(liq)": -1.0,
    }
    lnK_cr_u = fluoride.equilibrium_log_constant(cr_u_reaction, T)
    lnK_fe_u = fluoride.equilibrium_log_constant(fe_u_reaction, T)
    lnK_ni_u = fluoride.equilibrium_log_constant(ni_u_reaction, T)
    lnK_cr_fef2 = fluoride.equilibrium_log_constant(cr_fef2_reaction, T)
    lnK_cr_bef2 = fluoride.equilibrium_log_constant(cr_bef2_reaction, T)
    crf2_hot = np.exp(
        (
            fluoride.standard_gibbs_J_mol("CrF2_P21/c_No.14(s)", 977.15)
            - fluoride.standard_gibbs_J_mol("CrF2_L1(liq)", 977.15)
        )
        / (8.31446261815324 * 977.15)
    )
    crf2_cold = np.exp(
        (
            fluoride.standard_gibbs_J_mol("CrF2_P21/c_No.14(s)", 811.15)
            - fluoride.standard_gibbs_J_mol("CrF2_L1(liq)", 811.15)
        )
        / (8.31446261815324 * 811.15)
    )
    native_smoke = native_thermochimica_smoke_test(mstdb, thermochimica_executable)
    checks = {
        "thermochimica_reference_commit": THERMOCHIMICA_REFERENCE_COMMIT,
        "native_thermochimica_available": bool(native_smoke["available"]),
        "native_thermochimica_smoke_test": native_smoke,
        "fluoride_species_records": len(fluoride.records),
        "chloride_species_records": len(chloride.records),
        "G_CrF2_liquid_J_mol_at_923_15K": fluoride.standard_gibbs_J_mol("CrF2_L1(liq)", T),
        "G_Cr_solid_J_mol_at_923_15K": fluoride.standard_gibbs_J_mol("Cr_S1(s)", T),
        "lnK_Cr_plus_2UF4_to_CrF2_plus_2UF3": lnK_cr_u,
        "lnK_Fe_plus_2UF4_to_FeF2_plus_2UF3": lnK_fe_u,
        "lnK_Ni_plus_2UF4_to_NiF2_plus_2UF3": lnK_ni_u,
        "lnK_Cr_plus_FeF2_to_CrF2_plus_Fe": lnK_cr_fef2,
        "lnK_Cr_plus_BeF2_to_CrF2_plus_Be": lnK_cr_bef2,
        "CrF2_ideal_saturation_activity_977_15K": float(crf2_hot),
        "CrF2_ideal_saturation_activity_811_15K": float(crf2_cold),
        "check_selective_Cr_over_Fe_over_Ni_in_U_buffer": bool(lnK_cr_u > lnK_fe_u > lnK_ni_u),
        "check_FeF2_oxidizes_Cr": bool(lnK_cr_fef2 > 0.0),
        "check_Be_buffer_suppresses_Cr_dissolution": bool(lnK_cr_bef2 < 0.0),
        "check_cold_leg_has_lower_CrF2_saturation": bool(crf2_cold < crf2_hot),
        "chloride_spot_check_G_CrCl2_liquid_J_mol_at_923_15K": chloride.standard_gibbs_J_mol("CrCl2_L1(liq)", T),
    }
    checks["all_thermochemical_sanity_checks_pass"] = bool(
        checks["check_selective_Cr_over_Fe_over_Ni_in_U_buffer"]
        and checks["check_FeF2_oxidizes_Cr"]
        and checks["check_Be_buffer_suppresses_Cr_dissolution"]
        and checks["check_cold_leg_has_lower_CrF2_saturation"]
    )
    return checks


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit the MSTDB-TC-grounded corrosion/species-inventory model.")
    parser.add_argument("--mstdb-dir", type=Path, default=None, help="Directory containing MSTDB-TC fluoride/chloride *_No_Func.dat files")
    parser.add_argument("--targets", type=Path, default=ROOT / "data" / "targets.csv")
    parser.add_argument("--case-context", type=Path, default=ROOT / "results" / "advanced" / "advanced_case_context.csv")
    parser.add_argument(
        "--cases",
        type=Path,
        default=None,
        help=(
            "Optional case table with the same strict explicit advanced context fields as --targets; "
            "the legacy case_features.csv is not evaluated implicitly"
        ),
    )
    parser.add_argument("--parameters", type=Path, default=ROOT / "data" / "parameters.json")
    parser.add_argument("--reduced-parameters", type=Path, default=ROOT / "results" / "advanced" / "mechanistic_parameters.json")
    parser.add_argument("--outdir", type=Path, default=ROOT / "generated" / "advanced")
    parser.add_argument(
        "--thermochimica-executable",
        default=None,
        help="Optional path to Thermochimica bin/InputScriptMode; also checked via PATH",
    )
    parser.add_argument("--prior-strength", type=float, default=0.25)
    parser.add_argument("--max-nfev", type=int, default=5000)
    args = parser.parse_args()

    mstdb_dir = resolve_mstdb_directory(args.mstdb_dir)
    mstdb = MSTDBPair.from_calibration_directory(mstdb_dir)
    targets = merge_advanced_case_context(pd.read_csv(args.targets), args.case_context)
    effective_model = load_model_from_json(args.parameters)
    reduced_model = None
    if args.reduced_parameters.is_file():
        reduced_model = MechanisticCorrosionModel(
            effective_model,
            json.loads(args.reduced_parameters.read_text(encoding="utf-8")),
        )

    fit = fit_thermochemical_model(
        targets,
        effective_model,
        mstdb,
        prior_strength=args.prior_strength,
        max_nfev=args.max_nfev,
    )
    comparison, metrics = compare_thermochemical_models(targets, fit["model"], effective_model, reduced_model)
    optimizer = fit["optimizer_result"]
    checks = database_checks(mstdb, args.thermochimica_executable)
    lower_hits: list[str] = []
    upper_hits: list[str] = []
    for spec in THERMOCHEMICAL_PARAMETER_SPECS:
        value = float(fit["model"].params[spec.name])
        tolerance = 1.0e-7 * max(1.0, abs(spec.lower), abs(spec.upper))
        if abs(value - spec.lower) <= tolerance:
            lower_hits.append(spec.name)
        if abs(value - spec.upper) <= tolerance:
            upper_hits.append(spec.name)

    modern_ids = {f"M-{index:03d}" for index in range(29, 39)}
    external_predictions, external_metrics = validate_held_out_measurements(
        targets,
        effective_model,
        mstdb,
        modern_ids,
        prior_strength=args.prior_strength,
        max_nfev=min(args.max_nfev, 1000),
    )
    physics_details = fit["physics_residual_details"]
    physics_metrics = {
        "n_constraints": int(len(physics_details)),
        "pass_fraction": float(physics_details["constraint_pass"].mean()) if len(physics_details) else None,
        "constraints": physics_details.to_dict(orient="records"),
    }
    metrics.update(
        {
            "optimizer_success": bool(optimizer.success),
            "optimizer_message": str(optimizer.message),
            "optimizer_cost": float(optimizer.cost),
            "optimizer_nfev": int(optimizer.nfev),
            "prior_strength": float(args.prior_strength),
            "parameters_at_lower_bound": lower_hits,
            "parameters_at_upper_bound": upper_hits,
            "thermochemical_backend": "mstdb_tc_standard_state_nernst",
            "native_thermochimica_equilibrium_executed": bool(
                checks["native_thermochimica_smoke_test"]["executed"]
            ),
            "validation_interpretation": (
                "Calibrated standard-state/Nernst benchmark. It is not a native SUBQ equilibrium result "
                "and is not, by itself, evidence of out-of-family predictive validity."
            ),
            "thermochemical_checks": checks,
            "species_mechanism_constraints": physics_metrics,
            "modern_family_conditional_sensitivity": external_metrics,
        }
    )
    supported = targets[targets.apply(mechanistic_is_supported, axis=1)]
    diagnostics = pd.DataFrame([fit["model"].diagnostic_row(row) for _, row in supported.iterrows()])
    case_predictions = None
    if args.cases is not None and args.cases.is_file():
        case_predictions = fit["model"].predict_cases(pd.read_csv(args.cases))
    outputs = {
        "fit": fit,
        "comparison": comparison,
        "metrics": metrics,
        "diagnostics": diagnostics,
        "case_predictions": case_predictions,
    }
    save_thermochemical_outputs(outputs, args.outdir)
    results_dir = args.outdir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "mstdb_thermochemistry_checks.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")
    external_predictions.to_csv(results_dir / "thermochemical_external_modern_validation.csv", index=False)
    (results_dir / "thermochemical_external_modern_validation_metrics.json").write_text(
        json.dumps(external_metrics, indent=2), encoding="utf-8"
    )
    make_comparison_plot(comparison, args.outdir / "figures")

    print("MSTDB-TC thermochemical validation complete")
    print(f"MSTDB version: {metrics['mstdb']['fluoride']['version']}")
    print(f"Optimizer success: {metrics['optimizer_success']} ({metrics['optimizer_message']})")
    print(f"Supported constraints: {metrics['n_thermochemical_constraints']}")
    print(f"Direct/range targets: {metrics['n_direct_or_range_targets']}")
    for label in ("effective", "reduced_mechanistic", "thermochemical"):
        if f"{label}_log_rmse_direct" in metrics:
            print(
                f"{label}: log RMSE={metrics[f'{label}_log_rmse_direct']:.6f}, "
                f"median factor={metrics[f'{label}_median_factor_error_direct']:.6f}, "
                f"within factor 2={metrics[f'{label}_within_factor_2_direct']:.1%}, "
                f"constraint pass={metrics[f'{label}_constraint_pass_fraction']:.1%}"
            )
    print(f"Thermochemical sanity checks: {checks['all_thermochemical_sanity_checks_pass']}")
    print(
        "Modern-family conditional sensitivity (shared base not held out): "
        f"median factor={external_metrics.get('median_factor_error_direct', float('nan')):.6f}, "
        f"within factor 2={external_metrics.get('within_factor_2_direct', float('nan')):.1%}, "
        f"within factor 5={external_metrics.get('within_factor_5_direct', float('nan')):.1%}"
    )
    print(f"Outputs written to: {args.outdir}")


if __name__ == "__main__":
    main()
