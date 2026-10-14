#!/usr/bin/env python3
"""Training case 12: run coupon_rate.i for seven measured corrosion rates and compare.

Usage:  python run_scan.py            (uses ../../../../open_pronghorn-opt or $OPEN_PRONGHORN_EXE)

Writes scan_results.csv (case, features, measured, model, ratio) and the parity figure
case12_parity.png.
"""

import csv
import math
import os
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
import plotting  # noqa: E402

# case, alloy, salt, redox, position, T [K], flow factor, delta T [C], measured [um/y] (low, mid, high)
CASES = [
    ("ORNL-FL-04", "hastelloy_n", "fluoride_fuel", "purified_baseline", "nominal", 977.15, 1.0,
     166.0, (1.016, 1.016, 1.016)),
    ("FB-03", "hastelloy_n", "fluoroborate", "purified_baseline", "nominal", 960.15, 1.0, 249.0,
     (7.366, 7.366, 7.366)),
    ("ORNL-FL-01", "stainless_304l", "fluoride_fuel", "msre_or_fuel_baseline", "nominal", 961.15,
     1.0, 100.0, (21.844, 21.844, 21.844)),
    ("ORNL-FL-12", "stainless_316", "fluoride_fuel", "msre_or_fuel_baseline", "nominal", 923.15,
     1.0, 110.0, (25.4, 25.4, 25.4)),
    ("MOD-F-07", "gh3535", "flinak", "purified_baseline", "nominal", 973.15, 0.35, 0.0,
     (5.910, 23.266, 91.600)),
    ("ORNL-FL-15", "stainless_316", "flibe", "impure_moisture", "hot_leg", 923.15, 1.0, 0.0,
     (213.36, 213.36, 213.36)),
    ("MOD-F-07s", "hastelloy_n", "flinak", "stressed", "nominal", 973.15, 0.35, 0.0,
     (216.93, 216.93, 216.93)),
]


# label offsets [points] so that neighbouring markers stay readable
LABEL_OFFSET = {"ORNL-FL-12": (6, 4), "MOD-F-07": (-58, 6), "ORNL-FL-01": (6, -12),
                "ORNL-FL-15": (-72, 4), "MOD-F-07s": (6, -12)}


def run_case(exe, case):
    name, alloy, salt, redox, position, T, flow, dT, _ = case
    base = f"scan_{name}"
    args = [exe, "-i", "coupon_rate.i", f"T={T}", f"material={alloy}", f"salt={salt}",
            f"redox={redox}", f"position={position}", f"flow={flow}", f"dT={dT}",
            f"Outputs/file_base={base}"]
    subprocess.run(args, cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    return plotting.read_csv(HERE / f"{base}.csv")["corrosion_rate_um_y"][0]


def main():
    exe = os.environ.get("OPEN_PRONGHORN_EXE", str(HERE.parents[3] / "open_pronghorn-opt"))
    rows = []
    for case in CASES:
        model = run_case(exe, case)
        low, mid, high = case[-1]
        rows.append((case, model, mid, model / mid))
        print(f"{case[0]:11s} {case[1]:15s} {case[2]:14s} {case[3]:22s} "
              f"measured {mid:8.3f}  model {model:8.3f}  ratio {model / mid:5.2f}")

    with open(HERE / "scan_results.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "alloy", "salt", "redox", "T_K", "measured_um_y", "model_um_y",
                    "ratio"])
        for case, model, mid, ratio in rows:
            w.writerow([case[0], case[1], case[2], case[3], case[5], mid, f"{model:.6g}",
                        f"{ratio:.4f}"])

    log_err = [abs(math.log(r[3])) for r in rows]
    print(f"geometric-mean |ratio| = {math.exp(sum(log_err) / len(log_err)):.3f}; "
          f"all within a factor of 2: {all(e < math.log(2) for e in log_err)}")

    fig, ax = plotting.plt.subplots(figsize=(5.6, 5.0))
    lims = (0.5, 500)
    ax.plot(lims, lims, "k-", lw=1)
    ax.fill_between(lims, [l / 2 for l in lims], [l * 2 for l in lims], color="0.85",
                    label="factor of 2")
    for case, model, mid, _ in rows:
        low, _, high = case[-1]
        color = "C0" if case[1].startswith(("hastelloy", "gh3535")) else "C3"
        marker = "s" if case[3] in ("impure_moisture", "stressed") else "o"
        ax.errorbar(mid, model, xerr=[[mid - low], [high - mid]], fmt=marker, color=color,
                    ms=8, capsize=3)
        ax.annotate(case[0], (mid, model), textcoords="offset points",
                    xytext=LABEL_OFFSET.get(case[0], (6, -10)),
                    fontsize=8)
    ax.set(xscale="log", yscale="log", xlim=lims, ylim=lims,
           xlabel="measured corrosion rate [um/y]", ylabel="model corrosion rate [um/y]",
           title="Ni alloys (blue) vs stainless (red)")
    ax.legend(loc="upper left")
    plotting.save(fig, HERE, "case12_parity")


if __name__ == "__main__":
    main()
