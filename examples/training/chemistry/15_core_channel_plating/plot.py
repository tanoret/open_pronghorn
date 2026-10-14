#!/usr/bin/env python3
"""Case 15: axial wall profiles and the radial Cr plume / Mo-99 boundary layer."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import plotting  # noqa: E402
from plotting import plt  # noqa: E402

here, _ = plotting.case_setup(__file__, "core_channel.i", "core_channel_out.csv")
wall = plotting.read_csv(plotting.latest_vpp(here, "core_channel_out_wall_profile"))
rad = plotting.read_csv(plotting.latest_vpp(here, "core_channel_out_radial_profile"))

fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
axes[0].plot(wall["T_fluid"], wall["y"], color=plotting.SERIES[1])
axes[0].set(xlabel="wall-cell temperature [K]", ylabel="height [m]", title="heated channel")
axes[1].plot(wall["corrosion_rate_um_y"], wall["y"])
axes[1].set(xlabel="corrosion rate [um/y]", title="Arrhenius wall attack")
cr0 = max(rad["c_Cr"])
mo0 = max(rad["Mo99"])
axes[2].plot([x * 100 for x in rad["x"]], [c / cr0 for c in rad["c_Cr"]], label="Cr (from wall)")
axes[2].plot([x * 100 for x in rad["x"]], [m / mo0 for m in rad["Mo99"]],
             label="Mo-99 (to wall)")
axes[2].set(xlabel="distance from the left wall [cm]", ylabel="value / max",
            title="mid-height profiles")
axes[2].legend()
plotting.save(fig, here, "case15_core_channel")
