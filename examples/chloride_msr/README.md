# Complete molten chloride fast reactor (MCFR) fuel-channel model

`chloride_msr.i` is a single, fully coupled, publication-quality model of a molten chloride fast
reactor fuel channel. It integrates every physics of an operating MSR in one segregated (SIMPLE)
linear finite-volume solve and verifies the result.

## Physics

| Physics | Model |
| :------ | :---- |
| Fluid flow | incompressible turbulent Navier-Stokes (`LinearWCNSFVMomentumFlux`, Rhie-Chow, pressure projection) |
| Nuclear power | simplified separable cosine neutron-flux shape, peaked at the core mid-plane |
| Energy deposition | volumetric fission heat, advected enthalpy, molecular + turbulent conduction |
| Turbulence | standard k-epsilon with wall functions (turbulent viscosity `mu_t`) |
| Radiolysis | chloride radical chain producing the oxidant Cl2.- proportional to power, advected and diffused (`[MoltenSaltRadiolysis]`, `metals = Cr`) |
| Corrosion | temperature-dependent Butler-Volmer dissolution of the structural-alloy walls, releasing chromium as Cr(II) (`[CorrosionPlatingFlow]`, `release_variables = 'Cr_II'`) |

## Geometry and operating point

A 2D vertical fuel-salt channel of half-pitch `W = 0.5 m` between two 316-stainless walls, over an
active core height `H = 2 m`. The fuel salt is a NaCl-UCl3 fast-reactor salt (rho 3300 kg/m^3,
mu 4e-3 Pa s, cp 1000 J/kg/K, k 0.4 W/m/K, ~650 C). The salt flows upward at 1.5 m/s
(Reynolds number ~1.2e6, fully turbulent), enters at 600 C, and is heated by ~200 MW/m^3 peak fission
power.

## The radiolysis-corrosion coupling

The walls corrode and release chromium into the salt **as the radiolysis-tracked Cr(II) species**, so
the corrosion product is exactly what the radiolysis network acts on. The radiolytic oxidant Cl2.-
then oxidizes it (`Cl2m_rad + Cr_II -> Cr_III`), setting the salt's Cr(III)/Cr(II) oxidation state.
This couples the two physics: radiolysis makes the salt oxidizing (driving corrosion), and the
dissolved corrosion product buffers the radiolytic oxidant. The wall corrosion is Arrhenius in the
local temperature, so it is fastest where the salt is hottest.

## Verification and validation

The example reports verification postprocessors:

- **Energy balance** — the integrated fission power equals the enthalpy rise of the flow:
  `energy_balance_error = |delta_T - P/(rho cp Q)| / (P/(rho cp Q))` is below 0.5 %.
- **Turbulent regime** — Reynolds number ~1.2e6; `mu_t_max` reports the turbulent viscosity.
- **Corrosion rate** — `corrosion_rate_wall_max` (~150 um/y) reproduces the calibrated effective
  Butler-Volmer correlation for this aggressive oxidizing chloride case.
- **Redox coupling** — `redox_ratio_CrIII_CrII` is the radiolytic oxidation of the corrosion-product
  chromium per pass (it accumulates over recirculation).

Each physics component is independently validated in this application:

| Component | Validation |
| :-------- | :--------- |
| Turbulent flow | `validation/free_flow` (ERCOFTAC channel flow and backward-facing step) |
| Radiolysis kinetics | `validation/msr` (pulse-radiolysis transient-absorption traces) |
| Corrosion kinetics | `validation/corrosion` (76 cases / 43 targets of the effective Butler-Volmer correlation) |

```
open_pronghorn-opt -i chloride_msr.i
```

The Exodus output shows the velocity and temperature fields, the turbulent viscosity, the radiolytic
oxidant Cl2.-, the dissolved chromium Cr(II)/Cr(III) and its oxidation state, and the wall corrosion
rate field.

## Notes and extensions

- The fast electron/atom radiolysis kinetics are lumped into a net oxidant yield (Cl2.- produced
  directly proportional to dose) to keep the steady solve affordable; the explicit radical network can
  be integrated by removing the `g_value_overrides`.
- The inlet Cr(II) represents the bulk corrosion inventory of the recirculating salt; over many passes
  the loop accumulates the inventory and oxidation state. A recirculating-loop closure (mapping the
  outlet back to the inlet) and a two-region (core + heat exchanger) model are natural extensions.
- For the same coupling over a power transient see `examples/operating_msr_corrosion`; for the
  detailed corrosion/plating electrochemistry (electromigration, potential, two-block interface) see
  the finite-element `[CorrosionPlating]` framework and `examples/corrosion_plating`.
