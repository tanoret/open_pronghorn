#!/usr/bin/env python3
"""Case 03: the same NCL-16 loop with four corrosion models, against the measurements.

Runs (if needed) the three inputs of this case (BV capsule, flow wall law, DRIDN) and reads the
MSTDB-TC endpoint from the calibration results (it cannot be run without the licensed database)."""
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, cap = plotting.case_setup(__file__, "ncl16_capsule_bv.i", "ncl16_capsule_bv_out.csv")
_, loop = plotting.case_setup(__file__, "ncl16_loop_flow.i", "ncl16_loop_flow_out.csv")
_, dridn = plotting.case_setup(__file__, "ncl16_dridn.i", "ncl16_dridn_out.csv")
repo = here.parents[3]
summary = repo / "validation/corrosion/calibration/results/advanced/advanced_validation_summary.csv"
with open(summary, newline="") as f:
    mstdb = {r["measurement_id"]: float(r["thermochemical_prediction"]) for r in csv.DictReader(f)}

metrics = [  # (title, unit, measured, capsule, loop, DRIDN, MSTDB, which models are calibrated)
    ("hot-leg rate (M-013)", "um/y", 1.016, cap["mean_rate_um_y"][-1], loop["mean_rate_um_y"][-1],
     dridn["mean_rate_um_y"][-1], mstdb["M-013"], (False, False, True, True)),
    ("hot-leg loss (M-011)", "mg/cm$^2$", 2.9, cap["mass_loss_mg_cm2"][-1],
     loop["mass_loss_mg_cm2"][-1], dridn["mass_loss_mg_cm2"][-1], mstdb["M-011"],
     (False, False, True, True)),
    ("cold-leg gain (M-012)", "mg/cm$^2$", 1.7, 0.0, loop["mass_gain_mg_cm2"][-1],
     dridn["mass_gain_mg_cm2"][-1], mstdb["M-012"], (False, True, True, True)),
    ("salt Cr increase (M-014)", "ppm", 500.0, cap["salt_Cr_ppm"][-1], loop["salt_Cr_ppm"][-1],
     dridn["salt_Cr_ppm"][-1], mstdb["M-014"], (False, True, True, True)),
]
models = ["BV\ncapsule", "flow\nwall law", "DRIDN", "MSTDB-TC\n(not run)"]

fig = plt.figure(figsize=(14, 8.2))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.9])
for k, (title, unit, meas, *vals, calib) in enumerate(metrics):
    ax = fig.add_subplot(gs[0, k])
    xs = range(len(models))
    for i, (v, c) in enumerate(zip(vals, calib)):
        ax.bar(i, v, width=0.7, color=plotting.SERIES[i], edgecolor=plotting.SURFACE,
               hatch="//" if c else None, linewidth=0)
        ax.text(i, v, f"{v:.3g}", ha="center", va="bottom", fontsize=9, color=plotting.INK_SECONDARY)
    ax.axhline(meas, color=plotting.INK, ls="--", lw=1.4)
    ax.set_xticks(list(xs))
    ax.set_xticklabels(models, fontsize=8.5)
    ax.set(title=f"{title}\nmeasured {meas:g} (dashed)", ylabel=unit)
    ax.set_ylim(0, 1.18 * max(max(vals), meas))
    ax.grid(axis="x", visible=False)

# legend for the hatching (calibration vs prediction)
from matplotlib.patches import Patch  # noqa: E402

handles = [Patch(facecolor=plotting.INK_MUTED, label="prediction (not fitted to this value)"),
           Patch(facecolor=plotting.INK_MUTED, hatch="//", edgecolor=plotting.SURFACE,
                 label="calibration (in-sample: fitted to NCL-16)")]

# salt Cr history of the three models that are run here
ax = fig.add_subplot(gs[1, :2])
ax.plot(cap["time_y"], cap["salt_Cr_ppm"], color=plotting.SERIES[0], label="BV capsule (saturates)")
ax.plot(loop["time_y"], loop["salt_Cr_ppm"], color=plotting.SERIES[1], label="flow wall law (loop)")
ax.plot(dridn["time_y"], dridn["salt_Cr_ppm"], color=plotting.SERIES[2], label="DRIDN (calibrated)")
ax.plot([3.365], [500], "D", ms=8, color=plotting.INK, label="measured +500 ppm")
ax.plot([3.365], [mstdb["M-014"]], "s", ms=7, color=plotting.SERIES[3], label="MSTDB-TC endpoint")
ax.set(xlabel="time [y]", ylabel="salt Cr [ppm]", title="salt chromium")
leg = ax.legend(fontsize=9, loc="lower right")
ax.add_artist(leg)
ax.legend(handles=handles, fontsize=9, loc="lower center", bbox_to_anchor=(0.42, 0.0),
          title="bars above", title_fontsize=9)

# hot and cold wall mass change of the flow wall law
ax = fig.add_subplot(gs[1, 2:])
ax.plot(loop["time_y"], loop["mass_loss_mg_cm2"], color=plotting.SERIES[1], label="hot wall: loss")
ax.plot(loop["time_y"], loop["mass_gain_mg_cm2"], color=plotting.SERIES[1], ls=":",
        label="cold wall: gain")
ax.plot(cap["time_y"], cap["mass_loss_mg_cm2"], color=plotting.SERIES[0], label="capsule: loss")
ax.plot([3.365], [2.9], "D", ms=8, color=plotting.INK, label="measured loss 2.9")
ax.plot([3.365], [1.7], "D", ms=8, color=plotting.INK, mfc=plotting.SURFACE,
        label="measured gain 1.7")
ax.set(xlabel="time [y]", ylabel="[mg/cm$^2$]", title="wall mass change (flow wall law)")
ax.legend(fontsize=9, loc="upper left")
plotting.save(fig, here, "case03_corrosion_models")
