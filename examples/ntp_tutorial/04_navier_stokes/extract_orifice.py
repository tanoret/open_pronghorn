#!/usr/bin/env python3
"""Write the orifice of the Newton anchor (stage 3) as an input file for the stage-4 models.

Reads the final log_chi_0 ... log_chi_<n-1> of the anchor's CSV output and writes
orifice_profile.i: the cold-frit multiplier chi_cf(z) = exp(lambda_j), constant on each of the n
equal axial slices.

Usage: python3 extract_orifice.py [../03_ludewig_anchor/anchor_out.csv]
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
L = 0.7294


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "03_ludewig_anchor", "anchor_out.csv")
    with open(path) as f:
        row = list(csv.DictReader(f))[-1]
    n = sum(1 for k in row if k.startswith("log_chi_") and k[8:].isdigit())
    chi = [math.exp(float(row["log_chi_%d" % j])) for j in range(n)]
    starts = " ".join("%.10g" % (L * j / n) for j in range(n))
    values = " ".join("%.8g" % c for c in chi)
    text = f"""# Orifice of the Newton anchor (stage 3, orifice_implicit.i) at the design power: the
# cold-frit drag multiplier chi_cf = exp(lambda_j) of each of the {n} axial slices.
# Written by extract_orifice.py from {os.path.relpath(path, HERE)}; s_power = {float(row['s_power']):g}.

[Functions]
  [chi_cf]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '{starts}'
    y = '{values}'
  []
[]
"""
    with open(os.path.join(HERE, "orifice_profile.i"), "w") as f:
        f.write(text)
    print(f"chi_cf from {min(chi):.3f} to {max(chi):.3f} over {n} slices -> orifice_profile.i")


if __name__ == "__main__":
    main()
