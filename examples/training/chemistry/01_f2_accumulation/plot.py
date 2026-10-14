#!/usr/bin/env python3
"""Case 01: radiolytic F2 grows linearly with the analytic G-value source."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "f2_accumulation.i", "f2_accumulation_out.csv")
t_min = [t / 60 for t in d["time"]]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
plotting.reference(ax, t_min, d["analytic"], "analytic  S t")
ax.plot(t_min, d["c_F2_diss"], marker="o", ls="none", markevery=5, label="OpenPronghorn")
ax.set(xlabel="time [min]", ylabel="dissolved F$_2$ [mol/m$^3$]", title="F$_2$ accumulation")
ax.legend()
ax2.semilogy(t_min[1:], [max(abs(e), 1e-16) for e in d["rel_error"][1:]])
ax2.set(xlabel="time [min]", ylabel="relative error", title="model vs analytic")
plotting.save(fig, here, "case01_f2_accumulation")
