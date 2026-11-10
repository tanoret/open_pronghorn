# Thermal-hydraulics of a particle-bed nuclear thermal propulsion fuel element

This tutorial builds a resolved model of one fuel element of a particle-bed reactor (PBR) for
nuclear thermal propulsion. The element is Ludewig et al.'s Design 1, 1 GW in 37 elements,
heating hydrogen from 30 K to 3000 K at 7 MPa. The tutorial takes it from its closures to the
design point, in two and three dimensions.

It reproduces and extends a one-dimensional channel model of the element (the reference model):
the same hydrogen, the same packed-bed closures and the same anchor state. Every flow path is
resolved: the plenum, the frits, the bed and the chamber, in `r` and `z` (and `theta`). The flow
split between the axial slices, and the orifice that corrects it, are part of the solution
instead of an input.

Each stage builds on the previous one and states what it verifies. The equations and closures
are in [THEORY.md](THEORY.md).

| stage | directory | what it does |
|---|---|---|
| 1 | [01_closures](01_closures/README.md) | hydrogen properties; KTA packed-bed pressure drop vs Texas A&M air data, through both solvers |
| 2 | [02_element](02_element/README.md) | the element mesh; the Newton Darcy-Forchheimer and the SIMPLE Navier-Stokes models; cold flow and 5 % power |
| 3 | [03_ludewig_anchor](03_ludewig_anchor/README.md) | the 1 GW anchor: steady power continuation with the orifice designed together with the flow |
| 4 | [04_navier_stokes](04_navier_stokes/README.md) | the segregated Navier-Stokes model (SIMPLE and PIMPLE) at 1 GW with an orifice designed on the slice flows |
| 5 | [05_element_3d](05_element_3d/README.md) | the element in 3-D, Newton and segregated: verification against 2-D, a power tilt across the element, and the azimuthal flow excursion of a ring orifice |

## Main results

- **Closures.** The hydrogen properties of MOOSE and of the Python post-processing agree to
  7e-14. The KTA pressure drop agrees with 46 air measurements within the KTA uncertainty (mean
  5.3 %, max 13.7 %), and both solver paths reproduce the correlation to round-off.
- **Anchor.** Mixed outlet gas temperature 2989.4 K (reference model 2988.6 K), largest fuel
  temperature 3193 K, inlet pressure 7.84 MPa. Energy and mass are conserved to 3e-7 and 1e-13.
- **Flow excursion.** With uniform frits, the steady branch of the element has a turning point at
  8.5 % of the design power: the heated element is unstable without orificing. The orifice
  designed with the flow removes it. The design needs a mean reference: referencing a single
  slice can require a starved slice's drag to vanish, where the design is singular.
- **Channel momentum.** The resolved momentum of the plenum and chamber turns the cold flow split
  from inlet-end-overfed (Darcy model) into a U shape (Navier-Stokes). The decelerating plenum
  flow recovers static pressure towards the closed end. The Darcy model, like the reference
  model, omits this.
- **Navier-Stokes at power.** With an orifice designed on the slice mass flows, the segregated
  model reaches the design point in RZ with SIMPLE and with PIMPLE: 2989.8 K mixed outlet, slice
  outlet temperatures within 0.1 % of the rise, fuel at most 3185 K, energy balance to 1e-8.
  A design on the slice temperatures has no stable fixed point (a slice-alternating mode), and the
  Newton orifice does not fit the Navier-Stokes split. Frozen, the designed orifice drifts away in
  the segregated iterations (stage 4), but in physical time it is stable: kicked, the Newton
  element returns to the design point with a 2 s time constant (stage 3, `stability.i`).
- **3-D.** The 3-D model reproduces the 2-D anchor (fuel temperature within 0.12 %). A transverse
  power tilt of ±5 % across the element amplifies to a +47 % / −24 % spread of the gas
  temperature rise and limits the element to 18 % of the design power. The axially designed
  orifice cannot correct an azimuthal maldistribution.
- **Segregated model in 3-D.** The segregated SIMPLE and PIMPLE models reach the design point in
  3-D as well, but only with the orifice designed per column (axial slice and azimuthal sector).
  With rings, a perturbation around the element at round-off grows into a flow excursion between
  the columns of a slice past about 44 % power. SIMPLE needs a flowing start in 3-D (a plug-flow
  estimate); PIMPLE starts from rest.

## Running

All inputs run with the OpenPronghorn executable from their own directory, for example

```bash
cd 03_ludewig_anchor
mpirun -n 4 ../../../open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i
```

inside the MOOSE environment. The Python scripts need numpy, matplotlib and, for reading Exodus
files, netCDF4.

## Code

The tutorial uses these objects of OpenPronghorn (`src/ntp`, documented under `doc/content/source/ntp`):

| object | role |
|---|---|
| `NTPHydrogenFluidProperties` | ideal-gas hydrogen with NASA polynomials, Sutherland and Eucken |
| `NTPPackedBedFunctorMaterial`, `NTPDuctFlowFunctorMaterial` | closed-form Darcy-Forchheimer mobilities, heat transfer, conductivity (Newton) |
| `NTPDarcyMassFlux`, `NTPDarcyAdvection`, `NTPGasMassTimeDerivative` | mass and energy of the Darcy model, with series face mobility |
| `NTPDarcyMassFluxFunctorMaterial`, `NTPDarcyFluxIntegral`, `NTPMixedTemperature` | cell mass flux (AD), boundary flows, mixed temperature |
| `NTPSliceUniformityConstraint` | implicit orifice design (Newton) |
| `NTPPackedBedSegregatedFunctorMaterial`, `NTPDuctFrictionSegregatedFunctorMaterial`, `NTPGasPropertiesFunctorMaterial` | closures of the segregated model |
| `NTPRhieChowMassFlux`, `NTPRhieChowMassFluxAux`, `NTPRhieChowFluxIntegral` | series Rhie-Chow coupling, cell mass flux, boundary flows (segregated) |
| `NTPLinearFVEnthalpyExchange` | gas-particle exchange for the enthalpy, with partial elimination |
| `NTPOrificeFlowDesignAux`, `NTPRhieChowSliceFlow`, `NTPElementColumnsDivision` | orifice designed on the column mass flows, column flows from the Rhie-Chow fluxes, columns (slices x sectors) of the element (segregated) |
| `NTPFVCellIntegral` | volume integrals with the finite-volume cell volumes (the power the source kernels deposit) |
| `NTPOrificeDesignAux` | lagged orifice design on the slice temperatures (segregated) |

Regression tests are in `test/tests/ntp` and unit tests in `unit/src/NTP*`.
