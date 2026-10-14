#!/usr/bin/env python3
"""Case 08: where the noble metals go -- salt, plated surfaces, off-gas."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "noble_metal_partition.i", "noble_metal_partition_out.csv")

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
t = d["time_d"]
ax.plot(t, d["Mo99_salt"], label="in salt")
ax.plot(t, d["Mo99_pl"], label="plated")
ax.plot(t, d["Mo99_og"], label="off-gas")
plotting.reference(ax, t, d["Mo99_salt_analytic"], "salt, analytic")
ax.set(xlabel="time [d]", ylabel="Mo-99 [mol/m$^3$]", title="Mo-99 partition")
ax.legend()

names = ["Mo99", "Ru103", "Te132"]
salt = [d[n + "_salt"][-1] for n in names]
plated = [d[n + "_pl"][-1] for n in names]
offgas = [d[n + "_og"][-1] for n in names]
tot = [s + p + o for s, p, o in zip(salt, plated, offgas)]
x = range(len(names))
bottom = [0.0] * len(names)
for i, (vals, label) in enumerate([(salt, "salt"), (plated, "plated"), (offgas, "off-gas")]):
    frac = [100 * v / t_ for v, t_ in zip(vals, tot)]
    ax2.bar(x, frac, bottom=bottom, width=0.55, color=plotting.SERIES[i], label=label,
            edgecolor=plotting.SURFACE, linewidth=2)
    bottom = [b + f for b, f in zip(bottom, frac)]
ax2.set_xticks(list(x), ["Mo-99", "Ru-103", "Te-132"])
ax2.set(ylabel="share of the inventory after 30 d [%]", title="removal r$_p$ : r$_o$ = 5 : 1")
ax2.legend(loc="lower right")
ax2.grid(axis="x", visible=False)
plotting.save(fig, here, "case08_noble_metal_partition")
