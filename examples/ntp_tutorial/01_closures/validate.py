#!/usr/bin/env python3
"""Stage 1: validate the packed-bed pressure-drop closure of the fuel-element model.

Runs the 1-D packed bed of packed_bed.i (Newton Darcy-Forchheimer path) and packed_bed_simple.i
(segregated SIMPLE path) at every point of the Texas A&M air pressure-drop data set and compares the
predicted modified friction factor

    f Re_m = (dp / L) d^2 eps^3 / (mu u (1 - eps)^2),      Re_m = rho u d / (mu (1 - eps)),

with the measurements and with the KTA 3102.3 and Ergun correlations evaluated here in closed form.
f Re_m depends on Re_m only, so the air properties below cancel; they just set the velocity of
each run.

Data: data/tamu_air_horizontal_friction.csv, Table A-5 of G. Zou and R. Hu, "SAM Code Validation
on Frictional Pressure Drop through Pebble Beds", ANL/NSE-20/5, Argonne National Laboratory (2020);
original measurements by Hassan and Kang (Texas A&M University).

Usage (inside the MOOSE environment):
    python3 validate.py [--exe ../../../open_pronghorn-opt] [--model kta|ergun] [--no-simple]
Outputs: validation_points.csv, validation_summary.txt, validation.png
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RHO_AIR = 1.2041  # air at 1 atm, 20 C [kg/m^3]
MU_AIR = 1.8134e-5  # [Pa s]
BED_LENGTH = 1.0  # [m]; f Re_m does not depend on it


def kta(Re_m):
    """KTA 3102.3: psi = 160/Re_m + 3/Re_m^0.1, so f Re_m = 160 + 3 Re_m^0.9."""
    return 160.0 + 3.0 * Re_m**0.9


def ergun(Re_m):
    """Ergun: f Re_m = 150 + 1.75 Re_m."""
    return 150.0 + 1.75 * Re_m


def friction_from_dp(dp, eps, d, u):
    return dp / BED_LENGTH * d**2 * eps**3 / (MU_AIR * u * (1 - eps) ** 2)


def run_moose(exe, input_file, eps, d, u, model, workdir, tag):
    base = os.path.join(workdir, tag)
    cmd = [
        exe,
        "-i",
        os.path.join(HERE, input_file),
        f"eps={eps}",
        f"d={d}",
        f"L={BED_LENGTH}",
        f"u={u}",
        f"rho={RHO_AIR}",
        f"mu={MU_AIR}",
        f"model={model}",
        f"Outputs/file_base={base}",
        "--color",
        "off",
    ]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"MOOSE failed for {input_file} ({tag}):\n{out.stdout[-3000:]}\n{out.stderr[-3000:]}")
    with open(base + ".csv") as f:
        rows = list(csv.DictReader(f))
    return abs(float(rows[-1]["dp"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--exe", default=os.path.join(HERE, "..", "..", "..", "open_pronghorn-opt"))
    parser.add_argument("--model", default="kta", choices=["kta", "ergun"])
    parser.add_argument("--no-simple", action="store_true", help="skip the SIMPLE path")
    args = parser.parse_args()

    with open(os.path.join(HERE, "data", "tamu_air_horizontal_friction.csv")) as f:
        data = list(csv.DictReader(f))

    correlation = kta if args.model == "kta" else ergun
    points = []
    with tempfile.TemporaryDirectory() as tmp:
        for i, row in enumerate(data):
            d = float(row["particle_diameter_m"])
            eps = float(row["voidage"])
            Re_m = float(row["modified_reynolds"])
            measured = float(row["modified_friction_factor"])
            u = Re_m * MU_AIR * (1 - eps) / (RHO_AIR * d)
            newton = friction_from_dp(run_moose(args.exe, "packed_bed.i", eps, d, u, args.model, tmp, f"n{i}"), eps, d, u)
            simple = (
                np.nan
                if args.no_simple
                else friction_from_dp(
                    run_moose(args.exe, "packed_bed_simple.i", eps, d, u, args.model, tmp, f"s{i}"), eps, d, u
                )
            )
            points.append(
                dict(
                    d=d,
                    D_over_d=float(row["diameter_ratio"]),
                    eps=eps,
                    Re_m=Re_m,
                    measured=measured,
                    correlation=correlation(Re_m),
                    newton=newton,
                    simple=simple,
                )
            )
            print(f"{i + 1:3d}/{len(data)}  d={d:.5f}  Re_m={Re_m:7.0f}  measured={measured:7.0f}  "
                  f"Newton={newton:8.1f}  SIMPLE={simple:8.1f}")

    keys = list(points[0])
    with open(os.path.join(HERE, "validation_points.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(points)

    measured = np.array([p["measured"] for p in points])
    corr = np.array([p["correlation"] for p in points])
    newton = np.array([p["newton"] for p in points])
    simple = np.array([p["simple"] for p in points])
    err = newton / measured - 1
    lines = [
        f"Texas A&M air data, {len(points)} points, {args.model.upper()} closure",
        f"  MOOSE (Newton) vs data:        mean |error| {100 * np.mean(np.abs(err)):.2f} %, "
        f"max |error| {100 * np.max(np.abs(err)):.2f} %",
        f"  MOOSE (Newton) vs closed form: max |relative difference| {np.max(np.abs(newton / corr - 1)):.2e}",
    ]
    if not args.no_simple:
        lines.append(
            f"  MOOSE (SIMPLE) vs closed form: max |relative difference| {np.max(np.abs(simple / corr - 1)):.2e}"
        )
    lines.append("  KTA 3102.3 states a pressure-drop uncertainty of 15 %.")
    summary = "\n".join(lines)
    print(summary)
    with open(os.path.join(HERE, "validation_summary.txt"), "w") as f:
        f.write(summary + "\n")

    plot(points, args.model)


def replot():
    """Redraw the figure from validation_points.csv without rerunning MOOSE."""
    with open(os.path.join(HERE, "validation_points.csv")) as f:
        points = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
    plot(points, "kta")


def plot(points, model):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # categorical slots 1-3 of the reference palette (validated all-pairs), one per particle size
    colors = ["#2a78d6", "#eb6834", "#1baf7a"]
    markers = ["o", "s", "^"]
    groups = sorted({p["d"] for p in points})
    ink, muted, grid = "#0b0b0b", "#52514e", "#e4e3df"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for ax in (ax1, ax2):
        ax.set_facecolor("#fcfcfb")
        ax.grid(color=grid, lw=0.6)
        ax.set_axisbelow(True)
        for s in ax.spines.values():
            s.set_color(grid)
        ax.tick_params(colors=muted, labelsize=9)

    Re_all = [p["Re_m"] for p in points]
    Re = np.logspace(np.log10(0.8 * min(Re_all)), np.log10(1.2 * max(Re_all)), 100)
    curve = kta(Re) if model == "kta" else ergun(Re)
    ax1.plot(Re, curve, color=muted, lw=2, label=f"{model.upper()} correlation")
    ax1.fill_between(Re, 0.85 * curve, 1.15 * curve, color=muted, alpha=0.12, lw=0, label="±15 %")
    for g, c, m in zip(groups, colors, markers):
        sel = [p for p in points if p["d"] == g]
        D_over_d = sel[0]["D_over_d"]
        Re_g = [p["Re_m"] for p in sel]
        ax1.plot(Re_g, [p["measured"] for p in sel], m, ms=8, mfc="none", mec=c, mew=1.5,
                 label=f"data, d = {1000 * g:.2f} mm (D/d = {D_over_d:.1f})")
        ax1.plot(Re_g, [p["newton"] for p in sel], m, ms=4, color=c)
        ax2.plot(Re_g, [100 * (p["newton"] / p["measured"] - 1) for p in sel], m, ms=8, color=c,
                 mec="#fcfcfb", mew=1, label=f"d = {1000 * g:.2f} mm")
    ax1.plot([], [], "o", ms=4, color=muted, label="MOOSE (filled)")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("modified Reynolds number Re$_m$", color=ink)
    ax1.set_ylabel("modified friction factor f Re$_m$", color=ink)
    ax1.set_title("Texas A&M air packed beds", color=ink, fontsize=11, loc="left")
    ax1.legend(fontsize=8, frameon=False, labelcolor=ink)

    ax2.axhspan(-15, 15, color=muted, alpha=0.12, lw=0)
    ax2.set_ylim(-20, 20)
    ax2.axhline(0, color=muted, lw=1)
    ax2.set_xscale("log")
    ax2.set_xlabel("modified Reynolds number Re$_m$", color=ink)
    ax2.set_ylabel("MOOSE / measured − 1 [%]", color=ink)
    ax2.set_title("Prediction error (band: KTA uncertainty)", color=ink, fontsize=11, loc="left")
    ax2.legend(fontsize=8, frameon=False, labelcolor=ink)
    from matplotlib.ticker import LogLocator, NullFormatter, ScalarFormatter

    for ax in (ax1, ax2):
        ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
    ax1.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    ax1.yaxis.set_minor_formatter(NullFormatter())
    fig.savefig(os.path.join(HERE, "validation.png"), dpi=150, facecolor="#fcfcfb")


if __name__ == "__main__":
    main()
