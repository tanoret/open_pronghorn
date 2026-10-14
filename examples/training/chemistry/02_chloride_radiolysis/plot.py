#!/usr/bin/env python3
"""Case 02: chloride radical chain -- transient and dose scaling of the steady products."""
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "chloride_radiolysis.i", "chloride_radiolysis_out.csv")
scan = plotting.read_csv(plotting.latest_vpp(here, "chloride_radiolysis_out_dose_scan"))
dose = [10 ** (2 + math.floor(x)) for x in scan["x"]]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
t = d["time"][1:]
for key, label in [("c_Cl_rad", "Cl"), ("c_Cl2m_rad", "Cl$_2^-$"), ("c_Cl3_ion", "Cl$_3^-$"),
                   ("c_Cl2_diss", "Cl$_2$"), ("c_e_sol", "e$^-_{sol}$")]:
    ax.loglog(t, [max(v, 1e-30) for v in d[key][1:]], label=label)
ax.set(xlabel="time [s]", ylabel="concentration [mol/m$^3$]", ylim=(1e-20, 1e1),
       title="first cell (10$^2$ W/m$^3$)")
ax.legend(ncol=2)

for i, (key, label) in enumerate([("Cl2m_rad", "Cl$_2^-$"), ("Cl3_ion", "Cl$_3^-$")]):
    ax2.loglog(dose, scan[key], marker="o", color=plotting.SERIES[i + 1], label=label)
    slope = math.log(scan[key][-1] / scan[key][0]) / math.log(dose[-1] / dose[0])
    ax2.annotate(f"slope {slope:.2f}", (dose[-2], scan[key][-2]), textcoords="offset points",
                 xytext=(6, -14), color=plotting.INK_SECONDARY, fontsize=9)
ax2.set(xlabel="dose rate [W/m$^3$]", ylabel="concentration at 1000 s [mol/m$^3$]",
        title="dose scaling: Cl$_2^-$ ~ D$^{1/2}$, Cl$_3^-$ ~ D")
ax2.legend()
plotting.save(fig, here, "case02_chloride_radiolysis")
