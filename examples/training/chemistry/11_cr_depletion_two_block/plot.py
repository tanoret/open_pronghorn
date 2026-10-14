#!/usr/bin/env python3
"""Case 11: the Cr-depletion layer in 316 stainless and the falling corrosion rate."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "cr_depletion.i", "cr_depletion_out.csv")
prof = plotting.read_csv(plotting.latest_vpp(here, "cr_depletion_out_alloy_profile"))
L_salt = 0.5
cs0 = 0.17 * 8.0e6 / 51.9961

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax.plot(d["time_d"], d["corrosion_rate_um_y"])
ax.set(xlabel="time [d]", ylabel="corrosion rate [um/y]", yscale="log",
       title="diffusion-limited dissolution")
ax.yaxis.set_major_formatter(plotting.matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(plotting.matplotlib.ticker.NullFormatter())
ax.set_yticks([3, 5, 10, 20])
depth = [(x - L_salt) * 1e6 for x in prof["x"]]
ax2.plot(depth, [c / cs0 for c in prof["cs_Cr"]])
ax2.axhline(0.95, color=plotting.INK_MUTED, ls="--", lw=1.2)
ax2.annotate("95 % of bulk", (150, 0.95), textcoords="offset points", xytext=(0, -14),
             color=plotting.INK_SECONDARY, fontsize=9)
ax2.set(xlabel="depth into the alloy [um]", ylabel="Cr / bulk Cr", title="Cr profile after 1 year")
plotting.save(fig, here, "case11_cr_depletion")
