#!/usr/bin/env python3
"""Case 14: radiolysis switches on wall corrosion; the salt remembers it after shutdown."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "irradiated_pool.i", "irradiated_pool_out.csv")
t = d["time_d"]

fig, axes = plt.subplots(3, 1, figsize=(8, 7.5), sharex=True)
axes[0].plot(t, [100 * p for p in d["power"]], color=plotting.INK_SECONDARY)
axes[0].set(ylabel="power [%]", title="power history")
axes[1].plot(t, d["Cr_II_avg"], label="Cr(II)")
axes[1].plot(t, d["Cr_III_avg"], label="Cr(III)")
axes[1].plot(t, d["Cr_total"], label="total Cr")
axes[1].set(ylabel="[mol/m$^3$]", title="the oxidant converts Cr(II) to Cr(III)")
axes[1].legend(ncol=3)
axes[2].plot(t, d["corrosion_rate_um_y"])
axes[2].set(xlabel="time [d]", ylabel="wall loss [um/y]", title="wall corrosion rate")
plotting.save(fig, here, "case14_irradiated_pool")
