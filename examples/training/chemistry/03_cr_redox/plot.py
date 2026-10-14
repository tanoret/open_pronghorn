#!/usr/bin/env python3
"""Case 03: Cr(II)/Cr(III) redox under gamma irradiation and the two invariants."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "cr_redox.i", "cr_redox_out.csv")
t = d["time"][1:]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
for key, label in [("c_Cr_II", "Cr(II)"), ("c_Cr_III", "Cr(III)"), ("c_Cr_I", "Cr(I)"),
                   ("c_Cl2_diss", "Cl$_2$")]:
    ax.loglog(t, [max(v, 1e-30) for v in d[key][1:]], label=label)
ax.set(xlabel="time [s]", ylabel="concentration [mol/m$^3$]", ylim=(1e-14, 10),
       title="chromium oxidation states")
ax.legend()
ax2.loglog(t, [max(abs(v), 1e-18) for v in d["cr_error"][1:]], label="Cr conservation")
ax2.loglog(t, [max(abs(v), 1e-18) for v in d["redox_balance"][1:]], label="redox balance")
ax2.set(xlabel="time [s]", ylabel="relative error", title="invariants")
ax2.legend()
plotting.save(fig, here, "case03_cr_redox")
