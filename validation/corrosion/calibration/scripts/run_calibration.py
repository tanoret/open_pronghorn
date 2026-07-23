#!/usr/bin/env python3
"""Rerun the corrosion/plating calibration using only files in this repository.

The script is deliberately location-independent: all default paths are resolved
from this file, not from the caller's current working directory.  It writes to
an ignored generated directory and never modifies the frozen reference or the
OpenPronghorn production database.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Any


CALIBRATION_ROOT = Path(__file__).resolve().parents[1]
SRC = CALIBRATION_ROOT / "src"
DEFAULT_WORKBOOK = CALIBRATION_ROOT / "data" / "msr_corrosion_plating_validation_data.xlsx"
DEFAULT_REFERENCE = CALIBRATION_ROOT / "reference"
DEFAULT_OUTPUT = CALIBRATION_ROOT / "generated" / "latest"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from msr_corrosion_bv.calibrate import fit_model  # noqa: E402


CORE_METRICS = (
    "n_measurements_total",
    "n_active_constraints",
    "n_direct_or_range_targets",
    "n_upper_lower_ranking_constraints",
    "median_factor_error_direct",
    "within_factor_2_direct",
    "within_factor_5_direct",
    "constraint_pass_fraction",
)


def sha256(path: Path) -> str:
    """Return a content hash used to identify the exact local workbook."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def compare_with_reference(
    output_dir: Path,
    reference_dir: Path,
    parameter_atol: float,
    metric_atol: float,
) -> dict[str, Any]:
    """Compare a fresh fit with the frozen result and return an audit record."""
    fitted = load_json(output_dir / "results" / "parameters.json")
    expected = load_json(reference_dir / "parameters.json")

    missing = sorted(set(expected) - set(fitted))
    extra = sorted(set(fitted) - set(expected))
    common = sorted(set(expected) & set(fitted))
    differences = {
        name: abs(float(fitted[name]) - float(expected[name]))
        for name in common
    }
    worst_name = max(differences, key=differences.get)
    max_parameter_abs_difference = differences[worst_name]

    fitted_metrics = load_json(output_dir / "results" / "metrics.json")
    expected_metrics = load_json(reference_dir / "metrics.json")
    metric_differences: dict[str, float] = {}
    metric_failures: list[str] = []
    for name in CORE_METRICS:
        observed = fitted_metrics[name]
        reference = expected_metrics[name]
        if isinstance(reference, int):
            difference = 0.0 if observed == reference else float("inf")
        else:
            difference = abs(float(observed) - float(reference))
        metric_differences[name] = difference
        if difference > metric_atol:
            metric_failures.append(name)

    passed = (
        not missing
        and not extra
        and max_parameter_abs_difference <= parameter_atol
        and not metric_failures
    )
    return {
        "passed": passed,
        "parameter_count": len(fitted),
        "parameter_absolute_tolerance": parameter_atol,
        "max_parameter_absolute_difference": max_parameter_abs_difference,
        "worst_parameter": worst_name,
        "missing_parameters": missing,
        "extra_parameters": extra,
        "metric_absolute_tolerance": metric_atol,
        "metric_differences": metric_differences,
        "metric_failures": metric_failures,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workbook",
        type=Path,
        default=DEFAULT_WORKBOOK,
        help="local calibration workbook (default: calibration/data workbook)",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="generated-output directory",
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=DEFAULT_REFERENCE,
        help="frozen reference used for the reproducibility check",
    )
    parser.add_argument("--prior-strength", type=float, default=0.35)
    parser.add_argument("--max-nfev", type=int, default=4000)
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="skip Matplotlib diagnostics; useful for fast/headless checks",
    )
    parser.add_argument(
        "--no-reference-check",
        action="store_true",
        help="run the optimizer without comparing with the frozen result",
    )
    parser.add_argument(
        "--parameter-atol",
        type=float,
        default=1.0e-3,
        help="maximum absolute coefficient drift allowed by the reference check",
    )
    parser.add_argument(
        "--metric-atol",
        type=float,
        default=1.0e-6,
        help="maximum absolute drift allowed for scalar fit metrics",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workbook = args.workbook.expanduser().resolve()
    output_dir = args.outdir.expanduser().resolve()
    reference_dir = args.reference_dir.expanduser().resolve()

    if not workbook.is_file():
        raise FileNotFoundError(f"Calibration workbook does not exist: {workbook}")
    if not args.no_reference_check and not (reference_dir / "parameters.json").is_file():
        raise FileNotFoundError(f"Frozen parameter reference does not exist: {reference_dir}")
    if output_dir in {CALIBRATION_ROOT.resolve(), (CALIBRATION_ROOT / "data").resolve(), reference_dir}:
        raise ValueError(
            "Refusing to write generated output over calibration inputs or the frozen reference."
        )

    outputs = fit_model(
        workbook,
        output_dir=output_dir,
        prior_strength=args.prior_strength,
        max_nfev=args.max_nfev,
    )

    if not args.no_plots:
        # Import lazily so headless/no-plot runs do not initialize Matplotlib.
        from msr_corrosion_bv.plotting import make_all_plots

        make_all_plots(outputs, output_dir / "figures")

    comparison = None
    if not args.no_reference_check:
        comparison = compare_with_reference(
            output_dir,
            reference_dir,
            parameter_atol=args.parameter_atol,
            metric_atol=args.metric_atol,
        )

    manifest = {
        "workbook": str(workbook),
        "workbook_sha256": sha256(workbook),
        "output_directory": str(output_dir),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "prior_strength": args.prior_strength,
        "max_nfev": args.max_nfev,
        "plots_generated": not args.no_plots,
        "optimizer_success": outputs["metrics"]["optimizer_success"],
        "optimizer_message": outputs["metrics"]["optimizer_message"],
        "optimizer_nfev": outputs["metrics"]["optimizer_nfev"],
        "reference_comparison": comparison,
    }
    with (output_dir / "run_manifest.json").open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2)

    metrics = outputs["metrics"]
    print("Calibration complete")
    print(f"  workbook: {workbook}")
    print(f"  workbook SHA-256: {manifest['workbook_sha256']}")
    print(f"  optimizer success: {metrics['optimizer_success']}")
    print(f"  optimizer evaluations: {metrics['optimizer_nfev']}")
    print(f"  active constraints: {metrics['n_active_constraints']}")
    print(f"  median direct/range factor error: {metrics['median_factor_error_direct']:.9g}")
    print(f"  output: {output_dir}")

    if comparison is not None:
        print(
            "  max coefficient difference: "
            f"{comparison['max_parameter_absolute_difference']:.6g} "
            f"({comparison['worst_parameter']})"
        )
        print(f"  frozen-reference check: {'PASS' if comparison['passed'] else 'FAIL'}")
        if not comparison["passed"]:
            return 2
    return 0 if metrics["optimizer_success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

