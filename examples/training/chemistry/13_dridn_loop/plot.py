#!/usr/bin/env python3
"""Case 13: three years of a Hastelloy N thermal-convection loop (DRIDN)."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "dridn_loop.i", "dridn_loop_out.csv")
t = d["time_y"]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax.plot(t, d["salt_Cr_ppm"], label="Cr")
ax.plot(t, d["salt_Ni_ppm"], label="Ni")
ax.plot(t, d["salt_Fe_ppm"], label="Fe")
ax.plot([t[-1]], [500], marker="D", ls="none", color=plotting.INK, label="NCL-16 measured Cr")
ax.set(xlabel="time [y]", ylabel="dissolved in salt [ppm]", title="salt inventory")
ax.legend()
ax2.plot(t, d["mass_loss_mg_cm2"], label="hot leg: loss")
ax2.plot(t, d["mass_gain_mg_cm2"], label="cold leg: plating")
ax2.set(xlabel="time [y]", ylabel="mass change [mg/cm$^2$]", title="hot-to-cold mass transfer")
ax2.legend()
plotting.save(fig, here, "case13_dridn_loop")
