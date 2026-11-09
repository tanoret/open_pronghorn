#!/usr/bin/env python3
"""Stage 5: transverse distribution of the 3-D anchor.

Reads the Exodus output of anchor_3d.i and plots, against the azimuth theta (0 = the hot side
x > 0 of a tilted power, 180 = the cold side), relative to the mean temperature rise of the gas,
so that runs at different powers compare:
  - the gas temperature rise leaving the hot frit, averaged over the length of the element;
  - the largest particle temperature rise of the fuel bed.
Several runs (e.g. tilt = 0 and tilt = 0.05) can be compared.

Usage: python3 plot_3d.py anchor_3d_out.e[:label] [more.e[:label] ...]  ->  transverse.png
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np

T_IN = 30.0  # inlet gas temperature [K]


def azimuthal(path, block, var, reduce):
    d = netCDF4.Dataset(path)
    names = [b"".join(n).decode().strip("\x00") for n in d["name_elem_var"][:]]
    blocks = [b"".join(n).decode().strip("\x00") for n in d["eb_names"][:]]
    b = blocks.index(block) + 1
    c = d["connect%d" % b][:] - 1
    x, z = d["coordx"][:][c].mean(1), d["coordz"][:][c].mean(1)
    theta = np.degrees(np.arctan2(z, x))
    v = d["vals_elem_var%deb%d" % (names.index(var) + 1, b)][-1]
    th = np.unique(np.round(theta, 3))
    return th, np.array([reduce(v[np.isclose(np.round(theta, 3), t)]) for t in th])


def main():
    runs = sys.argv[1:] or ["anchor_3d_out.e"]
    colors = ["#2a78d6", "#eb6834", "#1baf7a"]
    ink, muted, grid = "#0b0b0b", "#52514e", "#e4e3df"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for i, run in enumerate(runs):
        path, _, label = run.partition(":")
        label = label or path
        th, T_gas = azimuthal(path, "hot_frit", "T_fluid", np.mean)
        th2, T_fuel = azimuthal(path, "fuel_bed", "T_solid", np.max)
        rise = T_gas.mean() - T_IN
        ax1.plot(th, (T_gas - T_IN) / rise, "-o", color=colors[i % 3], lw=2, ms=6, label=label)
        ax2.plot(th2, (T_fuel - T_IN) / rise, "-o", color=colors[i % 3], lw=2, ms=6, label=label)
        print(f"{label}: hot-frit gas {T_gas.min():.0f}..{T_gas.max():.0f} K, "
              f"fuel max {T_fuel.max():.0f} K")
    for ax in (ax1, ax2):
        ax.set_facecolor("#fcfcfb")
        ax.grid(color=grid, lw=0.6)
        ax.set_axisbelow(True)
        for s in ax.spines.values():
            s.set_color(grid)
        ax.tick_params(colors=muted, labelsize=9)
        ax.set_xlabel("azimuth θ [deg] (0: x > 0, the hot side of a tilt)", color=ink)
        ax.set_xlim(0, 180)
        ax.legend(fontsize=8, frameon=False, labelcolor=ink)
    ax1.set_ylabel("gas temperature rise / mean rise [-]", color=ink)
    ax1.set_title("Gas leaving the bed", color=ink, fontsize=11, loc="left")
    ax2.set_ylabel("largest fuel temperature rise / mean gas rise [-]", color=ink)
    ax2.set_title("Fuel", color=ink, fontsize=11, loc="left")
    fig.savefig("transverse.png", dpi=150, facecolor="#fcfcfb")


if __name__ == "__main__":
    main()
