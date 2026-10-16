#!/usr/bin/env python3
"""Case 05: MSFR (EVOL / VTB) capstone -- flow, xenon around the loop, wall corrosion and plating."""
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here = pathlib.Path(__file__).resolve().parent
plotting.ensure_output(here, "msfr_flow.i", "msfr_flow_out.csv", ranks=8)
d = plotting.case_setup(__file__, "msfr_chemistry.i", "msfr_chemistry_out.csv", ranks=8)[1]
mid = plotting.read_csv(plotting.latest_vpp(here, "msfr_chemistry_out_midplane"))
loop = plotting.read_csv(plotting.latest_vpp(here, "msfr_chemistry_out_loop"))
axis = plotting.read_csv(plotting.latest_vpp(here, "msfr_flow_out_axis"))


def chain(prof, centre):
    """Order wall-cell samples by the angle around a point inside the wall's loop (clockwise from
    the core side, over the top, down the external loop and back) and return the arc length."""
    idx = sorted(range(len(prof["x"])),
                 key=lambda i: -math.atan2(prof["y"][i] - centre[1], prof["x"][i] - centre[0]))
    w = {k: [v[i] for i in idx] for k, v in prof.items()}
    s = [0.0]
    for i in range(1, len(idx)):
        s.append(s[-1] + math.hypot(w["x"][i] - w["x"][i - 1], w["y"][i] - w["y"][i - 1]))
    return s, w


T_hot, T_cold = d["T_bulk_hot"][-1], d["T_bulk_cold"][-1]

# --- Figure 1: flow and temperature ---------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
ax = axes[0]
ax.plot(mid["x"], mid["T_fluid"], color=plotting.SERIES[1], label="mid-plane (y = 0)")
ax.plot(axis["y"], axis["T_fluid"], color=plotting.SERIES[0], label="near axis (r = 0.05 m), vs y")
ax.axhline(T_hot, ls="--", lw=1.4, color=plotting.INK_MUTED)
ax.annotate(f"bulk hot leg {T_hot:.0f} K", (0.02, T_hot), textcoords="offset points",
            xytext=(0, 4), fontsize=9, color=plotting.INK_SECONDARY)
ax.set(xlabel="r (mid-plane) or y (axis) [m]", ylabel="salt temperature [K]",
       title="core cavity temperature")
ax.legend(fontsize=9, loc="lower left")
ax = axes[1]
ax.plot(mid["x"], mid["vel_y"], color=plotting.SERIES[0], label="mid-plane, axial velocity")
ax.plot(axis["y"], axis["vel_y"], color=plotting.SERIES[2], label="near axis, axial velocity vs y")
ax.axhline(0, lw=0.8, color=plotting.INK_MUTED)
ax.set(xlabel="r (mid-plane) or y (axis) [m]", ylabel="v_y [m/s]", title="core cavity velocity")
ax.legend(fontsize=9)
ax = axes[2]
ax.plot(loop["T_fluid"], loop["y"], color=plotting.SERIES[1])
ax.axhspan(0.2725, 0.6775, color=plotting.SERIES[2], alpha=0.12, lw=0)
ax.axhspan(-0.673, 0.023, color=plotting.SERIES[0], alpha=0.12, lw=0)
ax.text(T_cold + 3, 0.47, "pump", fontsize=9, color=plotting.INK_SECONDARY)
ax.text(T_cold + 3, -0.33, "heat exchanger", fontsize=9, color=plotting.INK_SECONDARY)
ax.set(xlabel="salt temperature [K]", ylabel="height in the loop (r = 2.13 m) [m]",
       title=f"external loop: {T_hot:.0f} -> {T_cold:.0f} K")
plotting.save(fig, here, "case05_msfr_flow")

# --- Figure 2: xenon around the loop, Cr wall exchange ---------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
ax = axes[0]
xe0 = d["Xe135_avg"][-1]
ax.plot([x / xe0 for x in loop["Xe135"]], loop["y"], color=plotting.SERIES[0], label="Xe-135")
ax.plot([m / d["Mo99_avg"][-1] for m in loop["Mo99"]], loop["y"], color=plotting.SERIES[3],
        label="Mo-99")
ax.axhspan(0.2725, 0.6775, color=plotting.SERIES[2], alpha=0.12, lw=0)
ax.axhspan(-0.673, 0.023, color=plotting.SERIES[0], alpha=0.12, lw=0)
ax.text(0.99, 0.40, "pump:\nHe bubbling", fontsize=9, color=plotting.INK_SECONDARY)
ax.text(1.03, -0.45, "HX: noble-metal\nplating", fontsize=9, color=plotting.INK_SECONDARY)
ax.set(xlabel="inventory / salt average", ylabel="height in the loop (r = 2.13 m) [m]",
       title="fission products around the loop")
ax.legend(fontsize=9, loc="lower right")

walls = [("shield_profile", "inner (shield) wall", plotting.SERIES[0], (1.6, 0.0)),
         ("reflector_profile", "outer (reflector) wall", plotting.SERIES[1], (1.0, 0.0))]
for name, label, col, centre in walls:
    s, w = chain(plotting.read_csv(plotting.latest_vpp(here, f"msfr_chemistry_out_{name}")), centre)
    hx = [si for si, x, y in zip(s, w["x"], w["y"]) if x > 1.93 and -0.68 < y < 0.03]
    if hx:
        for a in axes[1:]:
            a.axvspan(min(hx), max(hx), color=col, alpha=0.08, lw=0)
    axes[1].plot(s, w["T_fluid"], color=col, label=label)
    mg = [j * 51.9961 * 3.15576e7 * 0.1 for j in w["J_Cr_wall"]]  # mol/m^2/s -> mg Cr/cm^2/y
    axes[2].plot(s, mg, color=col, label=label)
axes[1].axhline(T_hot, ls="--", lw=1.2, color=plotting.INK_MUTED)
axes[1].axhline(T_cold, ls="--", lw=1.2, color=plotting.INK_MUTED)
axes[1].set(xlabel="distance along the wall [m] (core -> top -> loop -> bottom)",
            ylabel="wall-cell salt temperature [K]", title="wall temperature")
axes[1].legend(fontsize=9)
axes[2].axhline(0, lw=0.8, color=plotting.INK_MUTED)
axes[2].set(xlabel="distance along the wall [m] (shaded: HX ducts)",
            ylabel="Cr wall flux [mg/cm$^2$/y]",
            title="Cr: + dissolves (hot), - plates (cold)")
axes[2].legend(fontsize=9)
plotting.save(fig, here, "case05_msfr_chemistry")
