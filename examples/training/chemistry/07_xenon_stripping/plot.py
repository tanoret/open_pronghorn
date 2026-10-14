#!/usr/bin/env python3
"""Case 07: Xe-135 poisoning vs the sparging (stripping) rate."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "xenon_stripping.i", "xenon_stripping_out.csv")
scan = plotting.read_csv(plotting.latest_vpp(here, "xenon_stripping_out_scan"))

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax.plot(d["time_h"], [100 * p for p in d["poison_no_strip"]], label="no stripping")
ax.plot(d["time_h"], [100 * p for p in d["poison_strip_1e-3"]], label="r = 10$^{-3}$ 1/s")
ax.set(xlabel="time [h]", ylabel="neutrons absorbed by Xe-135 [%]", title="approach to equilibrium")
ax.legend()
r = [max(v, 1e-6) for v in scan["r_strip"]]
ax2.semilogx(r, [100 * p for p in scan["xe_poison_fraction"]], marker="o")
ax2.axvline(1.0e-4, color=plotting.INK_MUTED, lw=1, ls="--")
ax2.annotate("$\\lambda_{Xe} + \\sigma\\phi$", (1.0e-4, 1.8), textcoords="offset points",
             xytext=(4, 0), color=plotting.INK_SECONDARY, fontsize=9)
ax2.set(xlabel="stripping rate r [1/s]  (r = 0 drawn at 10$^{-6}$)",
        ylabel="equilibrium poisoning [%]", title="stripping pays off once r > decay + burnout")
plotting.save(fig, here, "case07_xenon_stripping")
