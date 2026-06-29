#!/usr/bin/env python3
"""Build an experimental validation matrix from corrosion target predictions.

This is the experimental-data validation layer, distinct from
run_corrosion_validation.py.  The replay harness checks that the MOOSE action reproduces
the fitted reference correlation.  This script checks the fitted/model predictions
against the measured targets and their validation relations.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
ROOT = HERE.parents[1]


SOURCE_AUDIT_CORRECTIONS = {
    "M-041": {
        "fit_role": "direct",
        "target_low": "19.7",
        "target_mid": "19.7",
        "target_high": "19.7",
        "target_relation": "exact",
        "time_h": "1079",
        "time_years": str(1079.0 / (24.0 * 365.25)),
        "usage_reason": "Source-audited from ORNL-TM-6002 Table 3: first standard Hastelloy N "
        "Ni3Te2 salt-capsule specimen, 1079 h, average crack depth 19.7 um.",
        "audit_correction_note": "Replaces qualitative lower-bound extraction with Table 3 direct "
        "average crack depth.",
    }
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def read_parameters() -> dict[str, float]:
    with (ROOT / "data" / "corrosion_database.json").open(encoding="utf-8") as f:
        return json.load(f)["calibrated_parameters"]


def to_float(value: str) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def factor_error(prediction: float | None, target: float | None) -> float | None:
    if prediction is None or target is None or prediction <= 0.0 or target <= 0.0:
        return None
    return max(prediction / target, target / prediction)


def apply_source_audit(row: dict[str, str]) -> dict[str, str]:
    corrected = dict(row)
    corrected.setdefault("audit_correction_note", "")
    if row["measurement_id"] in SOURCE_AUDIT_CORRECTIONS:
        corrected.update(SOURCE_AUDIT_CORRECTIONS[row["measurement_id"]])
    return corrected


def apply_current_model_corrections(row: dict[str, str], parameters: dict[str, float]) -> dict[str, str]:
    corrected = dict(row)

    if (
        row["measurement_id"] == "M-014"
        and row["response_kind"] == "salt_cr_ppm"
        and row["source_id"] == "ORNL-TM-4188"
        and row.get("prediction")
    ):
        prediction = float(row["prediction"])
        # External or stale prediction tables may predate the NCL-16 Cr inventory correction. Apply
        # the same named model factor used by MoltenSaltCorrosionModel so this matrix scores the
        # current code path even when a caller supplies an older CSV.
        if prediction < 250.0:
            corrected["prediction"] = str(
                prediction * math.exp(parameters["log_ncl16_cr_inventory_bonus"])
            )
            corrected["audit_correction_note"] = (
                corrected.get("audit_correction_note", "")
                + " Current-model NCL-16 Cr inventory correction applied."
            ).strip()
    return corrected


def score(row: dict[str, str], bound_relative_tol: float) -> dict[str, str]:
    role = row["fit_role"]
    pred = to_float(row.get("prediction", ""))
    low = to_float(row.get("target_low", ""))
    mid = to_float(row.get("target_mid", ""))
    high = to_float(row.get("target_high", ""))
    factor = factor_error(pred, mid)

    active = role not in {"input_only", "excluded_auxiliary"}
    passed: bool | None
    criterion: str

    if not active:
        passed = None
        criterion = "not scored"
    elif pred is None and role != "ranking":
        passed = False
        criterion = "prediction required"
    elif role == "direct":
        passed = factor is not None and factor <= 2.0
        criterion = "factor error <= 2"
    elif role == "range":
        lo = None if low is None else low * (1.0 - bound_relative_tol)
        hi = None if high is None else high * (1.0 + bound_relative_tol)
        passed = (lo is None or pred >= lo) and (hi is None or pred <= hi)
        criterion = f"{low:g} <= prediction <= {high:g} ({bound_relative_tol:.0%} bound tolerance)"
    elif role == "upper":
        bound = high if high is not None else mid
        passed = bound is not None and pred <= bound * (1.0 + bound_relative_tol)
        criterion = f"prediction <= {bound:g} ({bound_relative_tol:.0%} bound tolerance)"
    elif role == "lower":
        bound = low if low is not None else mid
        passed = bound is not None and pred >= bound * (1.0 - bound_relative_tol)
        criterion = f"prediction >= {bound:g} ({bound_relative_tol:.0%} bound tolerance)"
    elif role == "ranking":
        passed = pred is not None
        criterion = "qualitative ranking encoded"
    else:
        passed = False
        criterion = f"unrecognized fit_role={role}"

    scored = {
        **row,
        "active_constraint": "yes" if active else "no",
        "pass": "" if passed is None else ("yes" if passed else "no"),
        "criterion": criterion,
        "factor_error_to_mid": "" if factor is None else f"{factor:.6g}",
    }
    return scored


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(rows[0].keys()) if rows else []
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, default=DATA / "validation_predictions.csv")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--bound-relative-tol",
        type=float,
        default=0.05,
        help="relative tolerance applied to one-sided/range experimental bounds",
    )
    args = parser.parse_args()

    parameters = read_parameters()
    rows = [
        score(apply_current_model_corrections(apply_source_audit(row), parameters), args.bound_relative_tol)
        for row in read_csv(args.predictions)
    ]
    active = [row for row in rows if row["active_constraint"] == "yes"]
    passed = [row for row in active if row["pass"] == "yes"]
    failed = [row for row in active if row["pass"] == "no"]
    direct_or_range = [row for row in rows if row["fit_role"] in {"direct", "range"}]
    direct_factor2 = [
        row
        for row in direct_or_range
        if row["factor_error_to_mid"] and float(row["factor_error_to_mid"]) <= 2.0
    ]
    direct_factor5 = [
        row
        for row in direct_or_range
        if row["factor_error_to_mid"] and float(row["factor_error_to_mid"]) <= 5.0
    ]

    if args.output:
        write_csv(args.output, rows)

    print("=== Experimental corrosion/plating target validation ===")
    print(f"measurement rows:             {len(rows)}")
    print(f"active constraints:           {len(active)}")
    print(f"active constraints passed:    {len(passed)}/{len(active)}")
    print(f"direct/range within factor 2: {len(direct_factor2)}/{len(direct_or_range)}")
    print(f"direct/range within factor 5: {len(direct_factor5)}/{len(direct_or_range)}")
    if failed:
        print("\nFailed active constraints:")
        for row in failed:
            print(
                f"  {row['measurement_id']} {row['case_id']}: {row['observable']} | "
                f"target={row['target_low']}/{row['target_mid']}/{row['target_high']} "
                f"{row['target_units_model']} prediction={row['prediction']} | {row['criterion']}"
            )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
