#!/usr/bin/env python3
"""Case 04: MSRE fuel loop -- xenon poisoning and noble-metal partition against the MSRE record."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "msre_loop.i", "msre_loop_out.csv")
plotting.ensure_output(here, "msre_loop.i", "msre_loop_bubbles.csv",
                       ["r_bubble=2.1414e-3", "Outputs/file_base=msre_loop_bubbles"])
b = plotting.read_csv(here / "msre_loop_bubbles.csv")

# ---- figure 1: xenon poisoning -------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 4.3))
th = [t * 24 for t in d["time_d"]]
ax.axhspan(0.26, 0.28, color=plotting.SERIES[2], alpha=0.25, lw=0,
           label="measured, U-235 runs: 0.26-0.28 % (TM-3464)")
ax.axhline(1.08, color=plotting.INK_SECONDARY, ls=":", lw=1.5,
           label="ORNL design-model prediction: 1.08 % (TM-3464)")
ax.plot(th, d["xe_poison_pct"], color=plotting.SERIES[0], label="a priori (spray + graphite, no bubbles)")
ax.plot([t * 24 for t in b["time_d"]], b["xe_poison_pct"], color=plotting.SERIES[1],
        label="+ bubble stripping, $r_b$ = 2.14e-3 1/s (calibrated)")
plotting.reference(ax, th, d["xe_poison_analytic_pct"], "analytic steady state")
plotting.reference(ax, th, b["xe_poison_analytic_pct"], None)
ax.set_xlim(0, 120)
ax.set_ylim(0, 1.25)
ax.set(xlabel="time after the step to 7.4 MW [h]", ylabel="Xe-135 poisoning [% $\\delta k/k$]",
       title="MSRE xenon poisoning")
ax.legend(loc="center right", fontsize=8.5, framealpha=0.9)
plotting.save(fig, here, "case04_msre_xenon")

# ---- figure 2: noble-metal partition vs ORNL-4865 Table 12.3 -------------------------------------
fig, ax = plt.subplots(figsize=(8.5, 4.3))
regions = [("salt", "salt"), ("metal", "Hastelloy N"), ("graphite", "graphite"), ("offgas", "off-gas")]
nuclides = [("Mo", "Mo-99"), ("Ru", "Ru-103"), ("Te", "Te-132")]
# ORNL-4865 Table 12.3 ranges [% of the inventory]: noble metals, and Te/Sb (Te-132)
ornl = {"Mo": [(1, 20), (5, 30), (5, 30), None], "Ru": [(1, 20), (5, 30), (5, 30), None],
        "Te": [(1, 20), (20, 90), (5, 30), None]}
for k in ornl:
    ornl[k][3] = (0.01, 0.9)  # off-gas: "tenths to hundredths of one percent"
w = 0.2
for j, (key, label) in enumerate(nuclides):
    for i, (reg, rlabel) in enumerate(regions):
        x = i + (j - 1) * w
        val = 100 * d[f"{key}_{reg}_fraction"][-1]
        ax.bar(x, max(val, 1e-3), width=w * 0.9, color=plotting.SERIES[j],
               edgecolor=plotting.SURFACE, label=label if i == 0 else None)
        rng = ornl[key][i]
        if rng:
            ax.plot([x, x], rng, color=plotting.INK, lw=2.2, solid_capstyle="butt",
                    label="ORNL-4865 range (MSRE)" if (i == 0 and j == 0) else None)
            ax.plot([x - w * 0.3, x + w * 0.3], [rng[0]] * 2, color=plotting.INK, lw=1.5)
            ax.plot([x - w * 0.3, x + w * 0.3], [rng[1]] * 2, color=plotting.INK, lw=1.5)
ax.set_yscale("log")
ax.set_ylim(0.005, 1000)
ax.set_xticks(range(len(regions)), [r[1] for r in regions])
ax.set(ylabel="share of the inventory after 42 d [%]",
       title="Noble metals: mass-transfer model (bars) vs MSRE (whiskers)")
ax.legend(loc="upper right", fontsize=8.5, framealpha=0.9, ncol=2)
ax.grid(axis="x", visible=False)
plotting.save(fig, here, "case04_msre_noble_metals")
