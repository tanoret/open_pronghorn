# Stage 5: the element in three dimensions

Stages 2–4 assume the element is axisymmetric. In a core it is not: an element next to the edge
of the active core sees a power gradient across it. This stage solves the anchors of stages 3 and
4 in three dimensions, on half of the element (`0 <= theta <= pi`, symmetry planes at `theta = 0`
and `pi`): the Newton Darcy-Forchheimer model with the implicit orifice design, and the segregated
Navier-Stokes model (SIMPLE and PIMPLE) with the orifice designed on the flows.

```bash
open_pronghorn-opt -i mesh_3d.i --mesh-only --allow-unused      # the mesh only
mpirun -n 12 open_pronghorn-opt -i anchor_3d.i ../03_ludewig_anchor/continuation.i \
    ../03_ludewig_anchor/orifice_implicit.i tolerance_3d.i                  # uniform, ~40 s
mpirun -n 8 open_pronghorn-opt -i anchor_3d.i ../03_ludewig_anchor/continuation.i \
    ../03_ludewig_anchor/orifice_implicit.i tolerance_3d.i tilt=0.05 \
    Outputs/file_base=anchor_3d_tilt_out                                   # tilted power
python3 plot_3d.py anchor_3d_out.e:"no tilt" anchor_3d_tilt_out.e:"tilt 0.05"   # transverse.png

# segregated Navier-Stokes model (linear FV)
mpirun -n 8 open_pronghorn-opt -i element_simple_3d.i                       # cold flow, ~20 s
mpirun -n 8 open_pronghorn-opt -i anchor_simple_3d.i \
    ../04_navier_stokes/orifice_flow.i orifice_sectors=8                     # SIMPLE, 1 GW, ~5 min
mpirun -n 8 open_pronghorn-opt -i anchor_pimple_3d.i \
    ../04_navier_stokes/orifice_flow.i orifice_sectors=8                     # PIMPLE, 1 GW, ~8 min
mpirun -n 8 open_pronghorn-opt -i anchor_pimple_3d.i \
    ../04_navier_stokes/orifice_flow.i t_end=10 pimple_steady_state_detection=false \
    Outputs/file_base=anchor_pimple_3d_rings_out                           # ring orifice: excursion
```

## Model

- **Mesh.** `element_mesh_3d.i` extrudes the reference strip of the axisymmetric mesh in a third
  coordinate and maps it onto the cones in cylindrical coordinates about the element axis `y`
  (`n_theta` cells over the half element). The axis, where the cells would collapse, is a
  0.4 mm symmetry cylinder.
- **Faceting.** The faceted annuli are 2.5 % smaller than the circular ones at `n_theta = 8`. The
  inlet mass flux therefore uses the mesh inlet area (`AreaPostprocessor`), and the power density
  is scaled by the exact over the meshed bed volume. The flow and power of the half element are
  then exact on any mesh.
- **Physics.** The physics, closures, continuation and orifice design are the stage-3 inputs.
  `anchor_3d.i` includes `../03_ludewig_anchor/anchor_body.i`. Two changes make the materials
  dimension-independent: the duct mobility treats both `x` and `z` as transverse, and the cell
  mass flux has three components.
- **Tolerance.** `tolerance_3d.i` raises the absolute tolerance of the continuation steps to
  1e-5. The 3-D residual has a higher round-off floor, and the balances stay at 3e-7.
- **Power tilt.** The power density is the axial anchor shape times `1 + tilt x / r_bed(z)`. This
  is a linear power gradient across the element, symmetric about `z = 0`. It does not change the
  element power.

## Results

| | axisymmetric (stage 3) | 3-D, no tilt |
|---|---|---|
| mixed outlet gas temperature | 2989.4 K | 2989.4 K |
| largest fuel temperature | 3193.2 K | 3197.1 K |
| inlet pressure | 7.841 MPa | 7.836 MPa |
| orifice multiplier, exit slice | 20.63 | 20.63 |
| energy balance error | 3e-7 | 3e-7 |

Without a tilt, the 3-D solution is axisymmetric to plotting accuracy and reproduces the
axisymmetric anchor. The largest fuel temperature differs by 0.12 % and the inlet pressure by
0.06 %. This verifies the 3-D model against the 2-D one.

With a tilt, the element no longer reaches the design power. The orifice of the anchor is a set
of rings: it shapes the flow along the element, not around it. A transverse power gradient
drives the same flow excursion as the axial one, now around the circumference.

| power tilt (±, across the element) | largest power reached (fraction of design) |
|---|---|
| 0 | 1.00 |
| 0.01 | 0.513 |
| 0.05 | 0.178 |
| 0.15 | 0.050 |

![transverse distribution](transverse.png)

At the turning point of the `tilt = 0.05` element:

- **Gas.** The ±5 % power difference has become a +47 % / −24 % spread of the gas temperature
  rise around the element.
- **Fuel.** The largest fuel temperature rise is 2.2 times the mean gas rise on the hot side,
  against 1.07 times without a tilt.

An element at the edge of a core with a power gradient across it needs an azimuthally shaped
orifice, a flattened power, or both. Neither the axisymmetric model of this tutorial nor the
one-dimensional channel model it reproduces can show this.

## Segregated Navier-Stokes model in 3-D

The segregated model of stages 2 and 4 (THEORY.md 6.2) runs on the same mesh with three velocity
components (`../common/element_segregated_3d.i`).

- `element_simple_3d.i`: the stage-2 element, cold or at a fixed power, with SIMPLE.
- `anchor_simple_3d.i`, `anchor_pimple_3d.i`: the stage-4 anchors, SIMPLE with a fixed-point
  continuation in the power and PIMPLE with a pseudo-time ramp, with `../04_navier_stokes/orifice_flow.i`.

**Start of the flow.** Started from rest, the first SIMPLE iteration sees only the gas viscosity
in the free channels: the face mass fluxes are zero, so the momentum diagonal of the plenum and
the chamber has no advection part. In 3-D the first pressure correction then gives velocities of
order 1e5 m/s, and the iterations diverge at any relaxation. Any start with a flow converges, so
SIMPLE starts from a plug-flow estimate of the cold velocity, `segregated_plug_velocity` (the inlet
speed in the plenum falling linearly to zero at the closed end, the outlet speed in the chamber
rising from zero). PIMPLE starts from rest: its time derivative bounds the momentum mobility by
`dt / rho`. The two give the same cold flow, an inlet pressure 11.97 kPa above the outlet
(11.55 kPa axisymmetric), with mass conserved to 1e-11.

**Rings are not enough.** With the orifice of stage 4, one multiplier per axial slice (rings), the
3-D element loses its symmetry on the way up in power. A perturbation around the element, at
round-off (1e-7) at the start, grows by a factor of 5 to 12 per second of pseudo-time, faster as
the power rises. It is a `cos 2 theta` mode: the columns next to the symmetry planes lose their
flow, the column at 90° takes it. Past about 44 % power it becomes visible, and at 50 % the
starved columns exceed 5000 K, their cold frit carries almost no flow, and the slice outlet
temperatures spread by more than the mean rise. The growth starts from round-off, so the
run-to-run detail of the runaway differs (the suite checks the threshold). Within a slice the
columns share one multiplier, so the design cannot hold their flows apart; the growth itself is
the drift of a fixed orifice in the segregated iterations (stage 4, finding (d)), here between the
columns of a ring. The Newton model, whose transient is physical, shows the limit as a turning
point instead, at a power that falls with the size of the azimuthal imperfection (the tilt table
above: 51 % for a ±1 % tilt). Either way, a ring orifice cannot take the element to the design
power once the power varies around it.

**A column orifice reaches the design point.** With `orifice_sectors = 8` the design acts on every
column of the cold frit, 12 slices times 8 sectors (`NTPElementColumnsDivision`), each designed
for its share of the power. The perturbation then decays, and both executioners reach the design
power:

| | RZ, SIMPLE (stage 4) | 3-D, SIMPLE | 3-D, PIMPLE |
|---|---|---|---|
| mixed outlet gas temperature | 2989.81 K | 2989.42 K | 2989.42 K |
| largest fuel temperature | 3185 K | 3188.2 K | 3188.2 K |
| inlet pressure | 8.875 MPa | 8.902 MPa | 8.901 MPa |
| slice outlet spread / mean rise | 1.2e-3 | 1.1e-3 | 1.1e-3 |
| asymmetry around the element (fuel) | — | 2.1e-7 | 2.1e-7 |
| orifice range | 1.18 – 130 | 1.01 – 132 | 1.02 – 131 |
| energy balance error | -2e-8 | -1e-8 | -4e-10 |
| mass balance error | 2e-13 | 5e-13 | -2e-12 |
| path | 50 power steps, 200 solves | 50 power steps, 200 solves | ramp over 20 s, steady at t = 34.1 s |
| run time | 95 s, 1 rank | 250 s, 8 ranks | 450 s, 8 ranks |

The 3-D element reproduces the axisymmetric one: the fuel temperature to 0.1 %, the inlet
pressure to 0.3 %. The mixed outlet temperature is 0.4 K lower because the 3-D element receives
26.916 MW, the exact power of the axial shape over the conical bed (the power density is
normalized with the exact bed volume), against 26.920 MW on the axisymmetric mesh, whose
trapezoidal bed cells hold 1.5e-4 more volume than the cones.

**Energy balance in 3-D.** The cells of the 3-D mesh are mapped onto cones, so their faces are not
planar. The finite-volume kernels use the exact cell volume; the first-order element quadrature of
`ElementIntegralFunctorPostprocessor` and `VolumePostprocessor` overestimates it, by 6.0e-4 at
`nz = 6` and 1.5e-4 at `nz = 12`. The segregated postprocessors therefore integrate the power and
the bed volume with `NTPFVCellIntegral`, the sum over the finite-volume cells. With the quadrature
integrals the energy balance error was -1.5e-4 at every power, that volume difference exactly.
