#!/usr/bin/env python3
"""Fit reduced-mechanistic MSR corrosion models and compare to effective BV."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from msr_corrosion_bv.calibrate import load_model_from_json
from msr_corrosion_bv.advanced_context import merge_advanced_case_context
from msr_corrosion_bv.mechanistic import (
    compare_mechanistic_to_effective,
    fit_mechanistic_model,
    save_mechanistic_outputs,
)


def make_comparison_plot(comparison: pd.DataFrame, figure_dir: Path) -> None:
    direct = comparison[comparison["fit_role"].isin(["direct", "range"])].copy()
    direct = direct[
        pd.to_numeric(direct["target_mid"], errors="coerce").gt(0)
        & pd.to_numeric(direct["effective_prediction"], errors="coerce").gt(0)
        & pd.to_numeric(direct["mechanistic_prediction"], errors="coerce").gt(0)
    ]
    if direct.empty:
        return
    figure_dir.mkdir(parents=True, exist_ok=True)
    target = direct["target_mid"].astype(float).to_numpy()
    effective_ratio = direct["effective_prediction"].astype(float).to_numpy() / target
    mechanistic_ratio = direct["mechanistic_prediction"].astype(float).to_numpy() / target
    x = np.arange(len(direct), dtype=float)

    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    ax.axhline(1.0, linewidth=1.4, label="target")
    ax.axhline(2.0, linestyle="--", linewidth=0.9, label="factor 2")
    ax.axhline(0.5, linestyle="--", linewidth=0.9)
    ax.scatter(x - 0.12, effective_ratio, marker="o", label="effective BV")
    ax.scatter(x + 0.12, mechanistic_ratio, marker="x", label="reduced mechanistic")
    ax.set_yscale("log")
    ax.set_ylabel("Prediction / experimental target midpoint")
    ax.set_xlabel("Validation measurement")
    ax.set_title("Mechanistic vs effective-model validation")
    ax.set_xticks(x)
    ax.set_xticklabels(direct["measurement_id"].astype(str), rotation=60, ha="right")
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(figure_dir / "fig08_mechanistic_vs_effective_parity.png", dpi=220)
    fig.savefig(figure_dir / "fig08_mechanistic_vs_effective_parity.svg")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit and compare reduced-mechanistic MSR corrosion models.")
    parser.add_argument("--targets", type=Path, default=ROOT / "data" / "targets.csv", help="Normalized validation targets CSV")
    parser.add_argument("--case-context", type=Path, default=ROOT / "results" / "advanced" / "advanced_case_context.csv")
    parser.add_argument("--parameters", type=Path, default=ROOT / "data" / "parameters.json", help="Calibrated effective-BV parameter JSON")
    parser.add_argument("--outdir", type=Path, default=ROOT / "generated" / "advanced", help="Output project directory")
    parser.add_argument("--prior-strength", type=float, default=0.25, help="Regularization strength for mechanistic physical priors")
    parser.add_argument("--max-nfev", type=int, default=5000, help="Maximum optimizer evaluations")
    args = parser.parse_args()

    targets = merge_advanced_case_context(pd.read_csv(args.targets), args.case_context)
    effective_model = load_model_from_json(args.parameters)
    fit = fit_mechanistic_model(
        targets,
        effective_model,
        prior_strength=args.prior_strength,
        max_nfev=args.max_nfev,
    )
    comparison, metrics = compare_mechanistic_to_effective(targets, fit["model"], effective_model)
    result = fit["optimizer_result"]
    metrics.update(
        {
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "optimizer_cost": float(result.cost),
            "optimizer_nfev": int(result.nfev),
            "prior_strength": float(args.prior_strength),
        }
    )
    outputs = {"fit": fit, "comparison": comparison, "metrics": metrics}
    save_mechanistic_outputs(outputs, args.outdir)
    make_comparison_plot(comparison, args.outdir / "figures")

    print("Mechanistic validation complete")
    print(f"Optimizer success: {metrics['optimizer_success']}")
    print(f"Mechanistic constraints: {metrics['n_mechanistic_constraints']}")
    print(f"Direct/range targets: {metrics['n_direct_or_range_targets']}")
    print(f"Effective log RMSE: {metrics['effective_log_rmse_direct']:.6f}")
    print(f"Mechanistic log RMSE: {metrics['mechanistic_log_rmse_direct']:.6f}")
    print(f"Effective median factor error: {metrics['effective_median_factor_error_direct']:.6f}")
    print(f"Mechanistic median factor error: {metrics['mechanistic_median_factor_error_direct']:.6f}")
    print(f"Effective within factor 2: {metrics['effective_within_factor_2_direct']:.1%}")
    print(f"Mechanistic within factor 2: {metrics['mechanistic_within_factor_2_direct']:.1%}")
    print(f"Effective constraint pass: {metrics['effective_constraint_pass_fraction']:.1%}")
    print(f"Mechanistic constraint pass: {metrics['mechanistic_constraint_pass_fraction']:.1%}")
    print(f"Outputs written to: {args.outdir}")


if __name__ == "__main__":
    main()
