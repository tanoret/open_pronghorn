# Molten salt fast reactor (MSFR) primary loop

`msfr_loop.i` is a complete, fully coupled model of an MSFR **primary fuel-salt loop** — a closed,
**pump-driven** circulation around a central solid reflector, enclosed by a solid reflector/vessel,
with conjugate heat transfer. It is the corrected, realistic MSR concept: the fuel salt is not pushed
through a channel by inlet/outlet boundary conditions, it *circulates* around a loop driven by a pump,
and the surrounding reactor structure is modeled as real solid subdomains.

## Geometry

A 2D rectangular salt annulus (the primary loop) around a central solid reflector block, enclosed by
an outer solid reflector/vessel frame — built with `SubdomainBoundingBoxGenerator` /
`ParsedSubdomainMeshGenerator` and `SideSetsBetweenSubdomainsGenerator` (the `inner_wall` and
`outer_wall` fluid-solid interfaces). The four legs of the loop are:

- **core leg** (left) — the active core: volumetric fission heat + radiolytic dose;
- **top leg** — the hot leg;
- **heat-exchanger leg** (right) — heat removed to the secondary circuit;
- **bottom leg** — the pump.

## Physics (all coupled in one SIMPLE solve)

| Physics | Model |
| :------ | :---- |
| Flow | turbulent Navier-Stokes (k-epsilon, Rhie-Chow), **driven by a pump momentum source** (`LinearFVSource` body force in the bottom leg) and buoyancy (`LinearFVMomentumBuoyancy`); closed loop, pressure pinned |
| Nuclear power | cosine-shaped volumetric fission heat in the core leg |
| Heat exchanger | Newton cooling toward the cold-leg temperature (`LinearFVReaction` + `LinearFVSource`), which also sets the loop temperature level |
| Conjugate heat transfer | separate `T_fluid` (salt) and `T_solid` (reflector + vessel) temperatures, coupled at the interfaces by `LinearFVConvectiveHeatTransferBC`; the solids conduct heat |
| Radiolysis | chloride radical chain producing the oxidant Cl2.- in the core leg (`[MoltenSaltRadiolysis]`, `metals = Cr`) |
| Corrosion | temperature-dependent Butler-Volmer dissolution **at the true fluid-solid interfaces** (`[CorrosionPlatingFlow]`, `reaction_boundary = 'inner_wall outer_wall'`, `release_variables = 'Cr_II'`) |

The walls release chromium into the salt as the radiolysis-tracked `Cr_II`; as the salt recirculates,
the radiolytic oxidant produced in the core leg oxidizes it to `Cr_III`. Because the loop is closed,
this redox conversion **accumulates over many passes** (`redox_ratio_CrIII_CrII` is substantial,
unlike the single-pass channel examples).

The corrosion tracks chromium, iron and nickel (`elements = 'Cr Fe Ni'`). Only chromium has chloride
radiolysis chemistry, so it is released into the radiolysis `Cr_II`
(`release_variables = 'Cr_II none none'`); iron and nickel have no radiolysis reactions and are
tracked as their own corrosion-product scalars `c_Fe` and `c_Ni` (the `none` entries), whose solver
systems `c_Fe_sys` / `c_Ni_sys` are listed in `[Problem] linear_sys_names` and the executioner's
`passive_scalar_systems`. (Use `time_derivative = false` so the steady SIMPLE solve adds no time
derivative to these created variables.)

## Verified operating point

| Quantity | Value |
| :------- | :---- |
| Loop velocity (pump-driven) | ~1.5 m/s up the core leg, ~1.5 m/s down the HX leg |
| Salt temperature | ~700 C cold leg to ~765 C hot leg |
| Solid reflector temperature | tracks the salt (conjugate heat transfer) |
| **Energy balance** | fission power = heat-exchanger removal to within 0.01% (`energy_balance_error`) |
| Radiolysis-corrosion redox | Cr(III)/Cr(II) ~ 0.5 (accumulated over recirculation) |

The fuel salt is a NaCl-UCl3 fast-reactor salt (rho 3300 kg/m^3, mu 4e-3 Pa s, cp 1000 J/kg/K); the
reflector/vessel are a nickel alloy (k 23 W/m/K). The geometry, the closed pump-driven loop, the
solid structures and the operating conditions follow the spirit of the EVOL/MSFR benchmark.

```
open_pronghorn-opt -i msfr_loop.i
```

The Exodus output shows the recirculating velocity field, the fluid and solid temperatures (conjugate
heat transfer), the turbulent viscosity, the radiolytic oxidant Cl2.-, and the dissolved chromium
Cr(II)/Cr(III) released from the walls.

## Validation provenance and limitations

Each physics component is independently validated in this application: turbulent flow against
`validation/free_flow` (ERCOFTAC), radiolysis kinetics against `validation/msr` (pulse radiolysis),
and corrosion kinetics against `validation/corrosion` (the 76-case effective Butler-Volmer
correlation). The integrated loop conserves energy to 0.01%.

This is a coarse (60x60) 2D Cartesian representation of the loop: the true MSFR core is axisymmetric
with a curved ("hourglass") cavity, 16 inlet/outlet nozzles and external pump+heat-exchanger loops
(Brovchenko et al., EVOL; Rouch et al., Ann. Nucl. Energy 2014). The rectangular loop here captures
the pump-driven circulation, the conjugate heat transfer with the surrounding solids, and the
coupled radiolysis/corrosion, but a mesh-refinement study and the true axisymmetric EVOL geometry
(with the curved core and a proper pump/heat-exchanger loop closure) are the next steps before the
numbers are used quantitatively.
