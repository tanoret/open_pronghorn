# MSFR core-cavity model (EVOL-benchmark-inspired)

`msfr_cavity.i` is a complete, fully coupled, publication-quality model of a molten salt fast reactor
(MSFR) core in a **realistic core-cavity geometry** rather than a channel. It is the spirit of the
EVOL/MSFR CFD benchmark — a heated core cavity with internal recirculation — carrying the full
multiphysics in one segregated (SIMPLE) linear finite-volume solve.

## Geometry

An open 2D core cavity (2 m x 2 m) bounded by solid structural-alloy walls. The fuel salt (a
NaCl-UCl3 fast-reactor salt) is pumped in through a central inlet at the bottom and collected at a
central outlet plenum at the top. The central jet rises through the core and entrains the surrounding
salt, so **two large recirculation zones** form on either side, with the salt flowing back down along
the side walls. This is the realistic MSR concern: those slow, recirculating regions

- run hotter (poor heat removal -> a hot spot),
- trap the radiolytic oxidant, and
- see altered wall corrosion.

The cavity boundaries are carved with `ParsedGenerateSideset` into the central `inlet`, the central
`outlet`, and the surrounding solid walls (`left`, `right`, `top_wall`, `bottom_wall`); all four solid
walls are corroding structural alloy.

## Physics (all coupled in one solve)

| Physics | Model |
| :------ | :---- |
| Fluid flow | turbulent incompressible Navier-Stokes (Rhie-Chow), with the central jet and side recirculation |
| Nuclear power | simplified separable cosine neutron-flux shape |
| Energy deposition | volumetric fission heat, advected enthalpy, molecular + turbulent conduction |
| Turbulence | standard k-epsilon with wall functions |
| Radiolysis | chloride radical chain producing the oxidant Cl2.- proportional to power (`metals = Cr`) |
| Corrosion | temperature-dependent Butler-Volmer dissolution of the walls, releasing chromium as Cr(II) |

The walls release chromium into the salt as the radiolysis-tracked `Cr_II`, which the radiolytic
oxidant then oxidizes to `Cr_III` — coupling corrosion and radiolysis through the dissolved-chromium
redox, with the temperature-dependent corrosion fastest where the salt is hottest.

## Verification and validation

Reported postprocessors:

- **Recirculation** — `vel_y_min` is negative (downward, reverse flow along the side walls), the
  signature of the recirculation zones; `vel_y_max` is the central jet.
- **Hot spot** — `T_hotspot` (the peak temperature, in the slow recirculation zone) exceeds the mean
  `T_outlet`, the realistic poor-heat-removal concern.
- **Energy balance** — `energy_balance_error` compares the fission power to the net enthalpy convected
  out, using the Rhie-Chow **mass-flux-weighted** outlet/inlet temperature (so it is a true global
  conservation check even with the recirculation); it closes to well under 1%.
- **Turbulence / corrosion** — `mu_t_max`, `corrosion_rate_wall_max`, and the `redox_ratio_CrIII_CrII`
  coupling observable.

Each physics component is independently validated in this application:

| Component | Validation |
| :-------- | :--------- |
| Turbulent / recirculating flow | `validation/free_flow` (ERCOFTAC channel flow and backward-facing step) |
| Radiolysis kinetics | `validation/msr` (pulse-radiolysis transient-absorption traces) |
| Corrosion kinetics | `validation/corrosion` (76 cases / 43 targets of the effective Butler-Volmer correlation) |

```
open_pronghorn-opt -i msfr_cavity.i
```

The Exodus output shows the velocity field with the central jet and corner recirculation, the
temperature field with the recirculation-zone hot spot, the turbulent viscosity, the radiolytic
oxidant Cl2.-, the dissolved chromium Cr(II)/Cr(III) and its oxidation state, and the wall corrosion
rate.

## Notes and extensions

- This is a heated recirculating cavity inspired by the EVOL/MSFR benchmark; the inlet/outlet
  placement and the cavity aspect are simplified. The dead-zone hot spot illustrates why the inlet
  geometry (which sweeps the floor corners here) and the recirculation pattern are central MSFR design
  concerns -- they limit the allowable power density.
- A recirculating-loop closure (mapping the outlet back to the inlet through an external heat
  exchanger), buoyancy (Boussinesq) coupling, and delayed-neutron-precursor drift are natural next
  steps toward a full MSFR plant model.
- The straight-channel variant (`examples/chloride_msr`) and the transient power-cycle coupling
  (`examples/operating_msr_corrosion`) complement this cavity model; all share the same validated
  physics and database.
