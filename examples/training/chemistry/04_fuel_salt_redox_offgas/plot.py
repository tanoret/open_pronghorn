#!/usr/bin/env python3
"""Case 04: the U(III) buffer is used up, then F2 builds up in the salt and the cover gas."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "fuel_salt_redox.i", "fuel_salt_redox_out.csv")
t_h = [t / 3600 for t in d["time"]]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax.plot(t_h, [100 * f for f in d["U_III_fraction"]], label="U(III) / U")
ax.set(xlabel="time [h]", ylabel="U(III) fraction [%]", title="redox buffer")
ax.legend()
ax2.semilogy(t_h[1:], [max(v, 1e-14) for v in d["c_F2_diss"][1:]], label="F$_2$ in salt")
ax2.semilogy(t_h[1:], [max(v, 1e-14) for v in d["c_F2_gas"][1:]], label="F$_2$ in cover gas")
ax2.set(xlabel="time [h]", ylabel="concentration [mol/m$^3$]", ylim=(1e-10, 10),
        title="fluorine release after the buffer is spent")
ax2.legend()
plotting.save(fig, here, "case04_fuel_salt_redox")
