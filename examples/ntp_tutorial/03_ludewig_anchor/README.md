# Stage 3: the 1 GW anchor

The anchor is Ludewig et al.'s Design 1 at 1 GW, 7 MPa chamber pressure and a 3000 K outlet,
with 37 elements. This stage reaches it with the Newton Darcy-Forchheimer model by a steady
power continuation. The cold-frit orifice is designed together with the flow.

```bash
open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i     # about 10 s, 4 ranks
```

The inputs are:

- `anchor.i`: the element, its mesh and the 1 GW power distribution. It includes
  `anchor_body.i`, which also serves the 3-D stage. The fuel-bed power density averages
  4e10 W/m³, with the inlet-peaked axial shape of the published anchor (12 slices, from 1.28 to
  0.82).
- `continuation.i`: the time variable becomes a continuation parameter (THEORY.md 6.1). Each step
  is a steady Newton solve started from the previous one. The power fraction runs as follows:
  - 0 (cold flow), then
  - 5 % up to `t_design`, then
  - a linear ramp to 1 at `t = 1`.
- `orifice_implicit.i`: the orifice `chi_j = exp(lambda_j)` of each axial slice is a scalar
  unknown of the Newton system (`NTPSliceUniformityConstraint`). Its equations make every slice's
  gas outlet temperature equal to the mean over the slices, with the mean of `lambda` fixed at
  `orifice_level = 2`. The design is switched on gradually at 5 % power, between `design_on` and
  `design_full`.

## Results

| quantity | this model | one-dimensional reference model |
|---|---|---|
| element power | 26.92 MW | 1 GW / 37 (design) |
| mixed outlet gas temperature | 2989.4 K | 2988.6 K |
| largest fuel temperature | 3193 K | |
| inlet pressure | 7.84 MPa | |
| orifice multiplier, inlet slice → exit slice | 0.15 → 20.6 | |
| energy balance (H_out − H_in − Q)/Q | 3e-7 | |
| mass balance | 1e-13 | |

The mixed outlet temperature follows from the first law with the hydrogen model of stage 1. It
is within 0.03 % of the one-dimensional channel model this tutorial reproduces, which uses the
same hydrogen model and anchor state. The orifice is the result specific to the 2-D model: it opens
the inlet slice, where the power peaks, and throttles the cold frit more and more towards the exit
end.

## What the continuation shows

Without the orifice (`anchor.i continuation.i` with `s_design=0.3 t_design=0.6`), the
uniform-frit element has a turning point of the steady branch at 8.5 % of the design power. The
Newton steps converge quadratically up to it, and no steady solution exists above it. This is the
flow excursion of a heated parallel-channel system: a hotter slice has a lower gas density and a
larger pressure drop, so it receives less flow and heats further.

The orifice removes the excursion, but only if it is designed with the flow. The design phase
taught three lessons, each reflected in the inputs:

- **Reference the mean, not one slice.** With one slice fixed at `chi = 1`, the starved slices
  needed their cold-frit drag removed altogether (`lambda → −∞`). There the slice temperature no
  longer depends on the orifice, and Newton stalls. The mean reference with a level keeps every
  slice controllable.
- **The level is a margin, not a fix.** A higher uniform throttling before the design is active
  moves the uniform-frit fold to *lower* power: 3.3 % at `chi = 30`, 2.0 % at `chi = 100`. The
  orifice must shape the flow, not just add resistance.
- **The closures must be differentiated.** The interphase heat transfer depends on the cell mass
  flux. Computed in an auxiliary variable, that flux is lagged by one residual evaluation. The
  Jacobian is then wrong in the exchange terms, by 93 % (relative Frobenius norm) in one of our
  tests. `NTPDarcyMassFluxFunctorMaterial`
  computes it as an AD functor, which makes the Jacobian exact to finite-difference precision.

## Mesh

`nz` and the radial cells can be changed on the command line. The orifice has one component per
axial cell, so its order must follow `nz`:

```bash
open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i nz=24 log_chi_order=TWENTYFOURTH
```

A coarse version of this run is the regression test `test/tests/ntp/element_anchor`.
