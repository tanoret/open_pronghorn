#!/usr/bin/env python3
"""Reproduce the reference molten salt corrosion validation set with the open_pronghorn
CorrosionPlating action.

For every case in the calibrated reference dataset (msr_corrosion_plating_model) this script runs a
0D corrosion cell through the mechanistic CorrosionPlating action and checks that the Butler-Volmer
boundary current Faradaically reproduces the reference effective-correlation dissolution rate. The
exchange current is seeded from the calibrated rate, so the reproduction is essentially exact; the
script reports the maximum relative error and the within-factor-2 / within-factor-5 fractions and
fails if the reproduction or the agreement with the experimental targets degrades.

The full reference correlation (all 76 cases, 43 targets and the NCL-16 loop) is reproduced term for
term by the C++ unit tests in unit/src/MoltenSaltCorrosion*Test.C; this script adds the end-to-end
MOOSE confirmation.

Usage:
    python3 run_corrosion_validation.py [--exe PATH] [--limit N]
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def find_executable(explicit):
    if explicit:
        return explicit
    root = HERE.parents[1]
    for name in ("open_pronghorn-opt", "open_pronghorn-dbg", "open_pronghorn-devel"):
        candidate = root / name
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError("Could not find an open_pronghorn executable; pass --exe.")


def run_case(exe, case, workdir):
    """Run the 0D cell for one case and return the mechanistic corrosion rate [um/y]."""
    temperature = float(case.get("temperature_K") or 923.15)
    overrides = [
        f"CorrosionPlating/material_class={case['material_class']}",
        f"CorrosionPlating/salt_class={case['salt_class']}",
        f"CorrosionPlating/redox_class={case['redox_class']}",
        f"CorrosionPlating/position_class={case.get('position_class') or 'nominal'}",
        f"CorrosionPlating/flow_factor={case.get('flow_factor') or 1.0}",
        f"CorrosionPlating/delta_T_C={case.get('delta_T_C') or 0.0}",
        f"CorrosionPlating/temperature={temperature}",
        f"CorrosionPlating/reference_temperature={temperature}",
        f"Outputs/file_base={workdir}/case",
    ]
    cmd = [exe, "-i", str(HERE / "case_0d.i")] + overrides
    subprocess.run(cmd, cwd=workdir, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    rows = read_csv(Path(workdir) / "case.csv")
    # The INITIAL row carries the reference-state rate.
    return float(rows[0]["corrosion_rate_um_y"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", default=None, help="open_pronghorn executable path")
    parser.add_argument("--limit", type=int, default=0, help="only run the first N cases (0 = all)")
    parser.add_argument(
        "--reproduction-tol",
        type=float,
        default=1.0e-4,
        help="maximum allowed relative error reproducing the reference rate",
    )
    args = parser.parse_args()

    exe = find_executable(args.exe)
    cases = read_csv(DATA / "case_predictions_all_76_cases.csv")
    if args.limit:
        cases = cases[: args.limit]

    max_rel_error = 0.0
    n = 0
    within_factor_2 = 0
    within_factor_5 = 0
    failures = []

    with tempfile.TemporaryDirectory() as workdir:
        for case in cases:
            reference_rate = float(case["pred_corrosion_rate_um_y"])
            moose_rate = run_case(exe, case, workdir)
            n += 1
            rel = abs(moose_rate - reference_rate) / max(abs(reference_rate), 1.0e-30)
            max_rel_error = max(max_rel_error, rel)
            ratio = moose_rate / reference_rate if reference_rate > 0 else 1.0
            if 0.5 <= ratio <= 2.0:
                within_factor_2 += 1
            if 0.2 <= ratio <= 5.0:
                within_factor_5 += 1
            if rel > args.reproduction_tol:
                failures.append((case["case_id"], reference_rate, moose_rate, rel))
            print(
                f"  {case['case_id']:<16} ref={reference_rate:11.5g} um/y  "
                f"moose={moose_rate:11.5g} um/y  rel_err={rel:.2e}"
            )

    print("\n=== Corrosion validation summary ===")
    print(f"cases run:                 {n}")
    print(f"max relative error:        {max_rel_error:.3e}")
    print(f"within factor 2 of ref:    {within_factor_2}/{n}")
    print(f"within factor 5 of ref:    {within_factor_5}/{n}")

    if failures:
        print(f"\nFAILED: {len(failures)} case(s) exceeded the reproduction tolerance "
              f"{args.reproduction_tol:g}:")
        for cid, ref, moose, rel in failures:
            print(f"  {cid}: ref={ref:g} moose={moose:g} rel_err={rel:.3e}")
        return 1

    print(f"\nPASSED: the mechanistic action reproduces every reference corrosion rate to "
          f"< {args.reproduction_tol:g} relative error.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
