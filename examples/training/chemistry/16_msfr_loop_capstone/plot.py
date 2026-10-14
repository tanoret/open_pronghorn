#!/usr/bin/env python3
"""Case 16: MSFR loop -- top-leg profiles and the comparison with the 0D cases."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, d = plotting.case_setup(__file__, "msfr_loop_capstone.i", "msfr_loop_capstone_out.csv")
line = plotting.read_csv(plotting.latest_vpp(here, "msfr_loop_capstone_out_loop_line"))

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax.plot(line["x"], line["T_fluid"], color=plotting.SERIES[1])
ax.set(xlabel="x along the top leg [m]", ylabel="salt temperature [K]",
       title="hot leg: core outlet to HX inlet")
labels = ["Xe-135\nvs 0D", "Mo-99\nvs 0D", "Xe-135\nstripped", "Mo-99\nplated"]
vals = [d["Xe135_vs_0d"][-1], d["Mo99_vs_0d"][-1], d["Xe135_stripped_fraction"][-1],
        d["Mo99_plated_fraction"][-1]]
bars = ax2.bar(range(4), vals, width=0.55, color=[plotting.SERIES[0]] * 2 + [plotting.SERIES[2]] * 2)
for b, v in zip(bars, vals):
    ax2.annotate(f"{v:.4f}" if v > 0.99 else f"{100 * v:.0f} %", (b.get_x() + b.get_width() / 2, v),
                 textcoords="offset points", xytext=(0, 3), ha="center", fontsize=9,
                 color=plotting.INK_SECONDARY)
ax2.set_xticks(range(4), labels)
ax2.set(ylim=(0, 1.15), title="loop inventory = 0D prediction")
ax2.grid(False)
ax2.set_axisbelow(True)
plotting.save(fig, here, "case16_msfr_loop")
