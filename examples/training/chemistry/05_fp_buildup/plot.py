#!/usr/bin/env python3
"""Case 05: saturation (I-131) vs linear growth (Cs-137) of fission products."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "fp_buildup.i", "fp_buildup_out.csv")
t = d["time_d"]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
plotting.reference(ax, t, d["N_I131_analytic"])
ax.plot(t, d["N_I131"], marker="o", ls="none", markevery=20, label="I-131 (8 d)")
ax.set(xlabel="time [d]", ylabel="inventory [mol/m$^3$]", title="short-lived: saturates")
ax.legend()
plotting.reference(ax2, t, d["N_Cs137_analytic"])
ax2.plot(t, d["N_Cs137"], marker="o", ls="none", markevery=20, label="Cs-137 (30 y)")
ax2.plot(t, d["N_Sr90"], label="Sr-90 (29 y)")
ax2.set(xlabel="time [d]", ylabel="inventory [mol/m$^3$]", title="long-lived: grows with burnup")
ax2.legend()
plotting.save(fig, here, "case05_fp_buildup")
