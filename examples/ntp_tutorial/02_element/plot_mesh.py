#!/usr/bin/env python3
"""Draw the element mesh (mesh_in.e from mesh.i) in the r-z plane, one color per region.

The radius is stretched: the element is 34 mm in radius and 729 mm long.
Usage: python3 plot_mesh.py [mesh_in.e]  ->  mesh.png
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import netCDF4

path = sys.argv[1] if len(sys.argv) > 1 else "mesh_in.e"
d = netCDF4.Dataset(path)
x, y = d["coordx"][:] * 1000, d["coordy"][:] * 1000
names = [b"".join(n).decode().strip("\x00") for n in d["eb_names"][:]]
# categorical slots 1-5 of the reference palette, in radial order from the axis
colors = {"chamber": "#2a78d6", "hot_frit": "#eb6834", "fuel_bed": "#1baf7a", "cold_frit": "#eda100",
          "plenum": "#e87ba4"}
labels = {"chamber": "chamber", "hot_frit": "hot frit", "fuel_bed": "fuel bed", "cold_frit": "cold frit",
          "plenum": "plenum"}

fig, ax = plt.subplots(figsize=(5.4, 5.6), constrained_layout=True)
for b, name in enumerate(names):
    c = d["connect%d" % (b + 1)][:] - 1
    polys = [list(zip(x[e], y[e])) for e in c]
    ax.add_collection(PolyCollection(polys, facecolor=colors[name], edgecolor="#fcfcfb", lw=0.6,
                                     label=labels[name]))
ax.set_xlim(0, x.max() * 1.02)
ax.set_ylim(-60, y.max() + 30)
# the inlet is the plenum annulus at z = 0, between the cold frit and the outer wall
y0 = y < 1e-6
r_in_mid = 0.5 * (x[y0].max() + sorted(set(x[y0]))[-3])
ax.annotate("inlet (H$_2$, 30 K)", xy=(r_in_mid, 0), xytext=(r_in_mid + 3, -45), fontsize=8,
            color="#0b0b0b", arrowprops=dict(arrowstyle="->", color="#52514e"))
ax.annotate("outlet (to nozzle)", xy=(6, y.max()), xytext=(2, y.max() + 12), fontsize=8, color="#0b0b0b",
            arrowprops=dict(arrowstyle="->", color="#52514e"))
ax.set_xlabel("r [mm] (stretched)", color="#0b0b0b")
ax.set_ylabel("z [mm]", color="#0b0b0b")
ax.set_aspect(0.05)
ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=9)
for s in ax.spines.values():
    s.set_color("#e4e3df")
ax.tick_params(colors="#52514e", labelsize=8)
fig.savefig("mesh.png", dpi=150, facecolor="#fcfcfb")
