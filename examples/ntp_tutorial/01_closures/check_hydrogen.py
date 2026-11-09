#!/usr/bin/env python3
"""Stage 1: compare the MOOSE hydrogen properties with the Python model of common/hydrogen.py.

Runs hydrogen_properties.i (a temperature sweep from 20 to 6000 K at 7 MPa through
NTPHydrogenFluidProperties) and compares density, enthalpy, isobaric specific heat, viscosity and
conductivity with common/hydrogen.py. Both implement the same equations, so the differences must be
at round-off. It also prints the reference values the tutorial uses: the element flow rate
mdot = Q / (h(3000 K) - h(30 K)) / 37.

Usage (inside the MOOSE environment):
    python3 check_hydrogen.py [--exe ../../../open_pronghorn-opt]
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import hydrogen as H  # noqa: E402

P = 7e6


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--exe", default=os.path.join(HERE, "..", "..", "..", "open_pronghorn-opt"))
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        base = os.path.join(tmp, "h2")
        out = subprocess.run(
            [args.exe, "-i", os.path.join(HERE, "hydrogen_properties.i"), f"p={P}",
             f"Outputs/file_base={base}", "--color", "off"],
            capture_output=True, text=True)
        if out.returncode != 0:
            sys.exit(out.stdout[-3000:] + out.stderr[-3000:])
        path = sorted(glob.glob(base + "_props_*.csv"))[-1]
        with open(path) as f:
            rows = list(csv.DictReader(f))

    T = np.array([float(r["T"]) for r in rows])
    moose = {q: np.array([float(r[q]) for r in rows]) for q in ("rho", "h", "cp", "mu", "k")}
    python = {"rho": H.rho(T, P), "h": H.h(T), "cp": H.cp(T), "mu": H.mu(T), "k": H.k(T)}
    print(f"{len(T)} temperatures from {T.min():g} to {T.max():g} K at {P / 1e6:g} MPa")
    worst = 0.0
    for q in moose:
        # enthalpy changes sign near 298 K: compare it relative to its range
        scale = np.ptp(python[q]) if q == "h" else np.abs(python[q])
        diff = np.max(np.abs(moose[q] - python[q]) / scale)
        worst = max(worst, diff)
        print(f"  {q:3s}  max relative difference {diff:.2e}")
    T_back = H.T_from_h(H.h(T))
    print(f"  T_from_h(h(T)) - T: max {np.max(np.abs(T_back - T)):.2e} K")

    Q, n_elements = 1e9, 37
    dh = float(H.h(3000.0) - H.h(30.0))
    print(f"\nh(3000 K) - h(30 K) = {dh:.6e} J/kg")
    print(f"bundle flow Q / dh = {Q / dh:.6f} kg/s, element flow = {Q / dh / n_elements:.6f} kg/s")
    if worst > 1e-10:
        sys.exit(f"FAILED: MOOSE and Python properties differ by {worst:.2e}")
    print("PASSED: MOOSE and Python properties agree to round-off")


if __name__ == "__main__":
    main()
