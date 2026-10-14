#!/usr/bin/env python3
"""Case 09: a month of operation -- the different memories of the inventory."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "power_history.i", "power_history_out.csv")
t = d["time_d"]

fig, axes = plt.subplots(3, 1, figsize=(8, 7.5), sharex=True)
axes[0].plot(t, [100 * p for p in d["power"]], color=plotting.INK_SECONDARY)
axes[0].set(ylabel="power [%]", title="power history")
axes[1].plot(t, d["rho_Xe_pcm"])
axes[1].set(ylabel="Xe worth [pcm]", title="Xe-135 follows the power within hours")
i131 = max(d["N_I131"])
cs = max(d["N_Cs137"])
axes[2].plot(t, [v / i131 for v in d["N_I131"]], label="I-131 (weeks)")
axes[2].plot(t, [v / cs for v in d["N_Cs137"]], label="Cs-137 (integrates the energy)")
axes[2].set(xlabel="time [d]", ylabel="inventory / max", title="memory of the history")
axes[2].legend()
plotting.save(fig, here, "case09_power_history")
