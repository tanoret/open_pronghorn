#!/usr/bin/env python3
"""Case 01: fission oxidation spends the U(III) buffer in ~340 days, then F2 reaches the cover gas."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "fuel_salt_redox.i", "fuel_salt_redox_out.csv")
t_d = [t / 86400 for t in d["time"]]

fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4))
plotting.reference(ax, t_d, [100 * u / 490.0 for u in d["U_III_analytic"]],
                   "analytic  U(III)$_0$ - R$_{ox}$ t")
ax.plot(t_d, [100 * f for f in d["U_III_fraction"]], marker="o", ls="none", markevery=0.08,
        label="OpenPronghorn")
ax.set(xlabel="time [days]", ylabel="U(III) / U [%]",
       title="redox buffer under fission oxidation (4 MW/m$^3$)")
ax.legend()

ax2.semilogy(t_d[1:], [max(v, 1e-16) for v in d["c_F2_diss"][1:]], label="F$_2$ in salt")
ax2.semilogy(t_d[1:], [max(v, 1e-16) for v in d["c_F2_gas"][1:]], label="F$_2$ in cover gas")
ax2.set(xlabel="time [days]", ylabel="concentration [mol/m$^3$]", ylim=(1e-14, 1e2),
        title="F$_2$ appears only after the buffer is spent")
ax2.legend()

for i, (key, label) in enumerate([("uranium_error", "uranium"), ("fluorine_error", "fluorine"),
                                  ("redox_error", "redox balance")]):
    ax3.semilogy(t_d[1:], [max(abs(v), 1e-16) for v in d[key][1:]], label=label,
                 color=plotting.SERIES[i])
ax3.set(xlabel="time [days]", ylabel="relative error", ylim=(1e-16, 1e-1),
        title="conservation checks (verification)")
ax3.legend()
fig.tight_layout()
plotting.save(fig, here, "case01_fuel_salt_redox")
