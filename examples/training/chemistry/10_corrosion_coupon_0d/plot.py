#!/usr/bin/env python3
"""Case 10: a static capsule corrodes until the salt reaches the Nernst equilibrium."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "corrosion_coupon_0d.i", "corrosion_coupon_0d_out.csv")
t = d["time_d"]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax.plot(t, d["corrosion_rate_um_y"])
ax.axhline(16.82, color=plotting.INK_MUTED, ls="--", lw=1.6, label="reference rate 16.8 um/y")
ax.set(xlabel="time [d]", ylabel="corrosion rate [um/y]", title="the rate decays to zero")
ax.legend()
ax2.plot(t, d["recession_um"], label="recession")
ax2.axhline(0.662, color=plotting.INK_MUTED, ls="--", lw=1.6, label="Faraday limit 0.66 um")
ax2.set(xlabel="time [d]", ylabel="wall recession [um]", title="self-limiting attack")
ax2.legend(loc="lower right")
plotting.save(fig, here, "case10_corrosion_coupon")
