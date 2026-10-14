#!/usr/bin/env python3
"""Case 06: I-135/Xe-135 equilibrium and the post-scram xenon peak."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "xenon_transient.i", "xenon_transient_out.csv")
t = d["time_h"]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
plotting.reference(ax, t, d["N_Xe135_analytic"])
ax.plot(t, d["N_Xe135"], label="Xe-135")
ax.plot(t, d["N_I135"], label="I-135")
ax.axvline(96, color=plotting.INK_MUTED, lw=1)
ax.annotate("scram", (96, max(d["N_I135"])), textcoords="offset points", xytext=(4, -4),
            color=plotting.INK_SECONDARY, fontsize=9)
ax.set(xlabel="time [h]", ylabel="inventory [mol/m$^3$]", title="iodine feeds the xenon")
ax.legend()
ax2.plot(t, d["rho_Xe_pcm"], color=plotting.SERIES[0])
i = min(range(len(t)), key=lambda k: d["rho_Xe_pcm"][k])
ax2.annotate(f"peak {t[i] - 96:.1f} h after scram", (t[i], d["rho_Xe_pcm"][i]),
             textcoords="offset points", xytext=(-110, -6), color=plotting.INK_SECONDARY,
             fontsize=9)
ax2.set(xlabel="time [h]", ylabel="xenon worth [pcm]", title="the xenon pit")
plotting.save(fig, here, "case06_xenon_transient")
