#!/usr/bin/env python3
"""Stage 2: compare the cold flow of one element between the two solution methods.

Reads the Exodus output of the Newton Darcy-Forchheimer run (newton_cold.e) and of the segregated
SIMPLE run (simple_cold.e), both at zero power with uniform frits, and plots along the element:
  - the radial mass flux entering the fuel bed through the cold frit (the slice flow split);
  - the plenum pressure.
The Darcy model has wall friction but no inertia in the plenum; the Navier-Stokes model has both.

Usage: python3 compare_cold_flow.py [newton_cold.e simple_cold.e]  ->  cold_flow.png
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np


def axial_profile(path, block, var):
    d = netCDF4.Dataset(path)
    names = [b"".join(n).decode().strip("\x00") for n in d["name_elem_var"][:]]
    blocks = [b"".join(n).decode().strip("\x00") for n in d["eb_names"][:]]
    b = blocks.index(block) + 1
    c = d["connect%d" % b][:] - 1
    y = d["coordy"][:][c].mean(1)
    v = d["vals_elem_var%deb%d" % (names.index(var) + 1, b)][-1]
    # average over the radial cells of each axial row
    zs = np.unique(np.round(y, 6))
    return zs, np.array([v[np.isclose(np.round(y, 6), z)].mean() for z in zs])


def main():
    newton = sys.argv[1] if len(sys.argv) > 1 else "newton_cold.e"
    simple = sys.argv[2] if len(sys.argv) > 2 else "simple_cold.e"
    ink, muted, grid = "#0b0b0b", "#52514e", "#e4e3df"
    colors = {"Darcy-Forchheimer (Newton)": "#2a78d6", "Navier-Stokes (SIMPLE)": "#eb6834"}
    runs = {"Darcy-Forchheimer (Newton)": (newton, "G_r", "o"), "Navier-Stokes (SIMPLE)": (simple, "G_x", "s")}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for label, (path, gvar, marker) in runs.items():
        z, G = axial_profile(path, "cold_frit", gvar)
        zp, p = axial_profile(path, "plenum", "pressure")
        ax1.plot(1000 * z, -G, "-" + marker, color=colors[label], lw=2, ms=6, label=label)
        ax2.plot(1000 * zp, (p - 7e6) / 1e3, "-" + marker, color=colors[label], lw=2, ms=6, label=label)
        print(f"{label}: cold-frit inflow {(-G).min():.2f} .. {(-G).max():.2f} kg/(m2 s), "
              f"plenum p - p_out {(p.min() - 7e6) / 1e3:.1f} .. {(p.max() - 7e6) / 1e3:.1f} kPa")
    for ax in (ax1, ax2):
        ax.set_facecolor("#fcfcfb")
        ax.grid(color=grid, lw=0.6)
        ax.set_axisbelow(True)
        for s in ax.spines.values():
            s.set_color(grid)
        ax.tick_params(colors=muted, labelsize=9)
        ax.set_xlabel("axial position z [mm] (inlet at 0)", color=ink)
        ax.legend(fontsize=8, frameon=False, labelcolor=ink)
    ax1.set_ylabel("mass flux into the cold frit [kg/(m$^2$ s)]", color=ink)
    ax1.set_title("Flow split, cold, uniform frits", color=ink, fontsize=11, loc="left")
    ax2.set_ylabel("plenum pressure - 7 MPa [kPa]", color=ink)
    ax2.set_title("Plenum pressure", color=ink, fontsize=11, loc="left")
    fig.savefig("cold_flow.png", dpi=150, facecolor="#fcfcfb")


if __name__ == "__main__":
    main()
