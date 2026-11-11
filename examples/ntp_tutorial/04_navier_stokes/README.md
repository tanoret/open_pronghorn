# Stage 4: the anchor with the segregated Navier-Stokes model

Stage 3 reaches the anchor with the Darcy-Forchheimer model, which, like the one-dimensional
reference model, reduces the plenum and the chamber to channels with wall friction. Stage 2
showed that resolving their momentum changes the cold flow split qualitatively: the decelerating
plenum flow recovers pressure towards the closed end. This stage takes the segregated
Navier-Stokes model (THEORY.md 6.2) up to the design power with the SIMPLE and PIMPLE algorithms,
which needs an orifice designed for this model.

**Status.** With an orifice designed on the slice mass flows (`orifice_flow.i`), the segregated
model reaches the design point in RZ with both SIMPLE and PIMPLE. Frozen, the designed orifice is
not held by the segregated iterations, although it is physically stable: see (d) below.

## Inputs

```bash
open_pronghorn-opt -i anchor_simple.i orifice_flow.i          # (c) SIMPLE, flow-designed orifice
open_pronghorn-opt -i anchor_pimple.i orifice_flow.i          # (c) PIMPLE, flow-designed orifice
open_pronghorn-opt -i anchor_simple.i orifice_flow.i simple_fixed_point=320 orifice_max_updates=200
                                                              # (d) the orifice frozen after 200 solves
python3 extract_orifice.py ../03_ludewig_anchor/anchor_out.csv    # orifice of the Newton anchor
open_pronghorn-opt -i anchor_pimple.i                             # (a) that orifice, ramp to 1 GW
open_pronghorn-opt -i anchor_pimple.i orifice_lagged.i s_final=0.5 t_ramp=30 t_end=60 \
                   pimple_steady_state_detection=false            # (b) lagged temperature design

# (e) turbulent free channels, resolved by 8 cells across, and the energy closures of the beds
open_pronghorn-opt -i anchor_simple.i orifice_flow.i ../common/channel_turbulence.i \
                   ../common/channel_turbulence_rz.i n_ch=8 n_pl=8 simple_iterations=600
open_pronghorn-opt -i anchor_simple.i orifice_flow.i ../common/channel_turbulence.i \
                   ../common/channel_turbulence_rz.i n_ch=8 n_pl=8 simple_iterations=600 \
                   bed_emissivity=0.8 bed_dispersion=0.1 segregated_static_outputs=csv
```

- `anchor_simple.i`: the stage-2 segregated model with the SIMPLE executioner and a fixed-point
  continuation in the power: one SIMPLE solve per fixed-point iteration, the power stepped from
  `s_start` to `s_final` over the first `n_ramp` iterations, then held. The density, the friction
  speed and the orifice are updated between solves. The CSV has one line per solve.
- `anchor_pimple.i`: the same model with the PIMPLE executioner (`../common/element_pimple.i`). The
  power ramps in pseudo-time from `s_start` to `s_final` over `t_ramp`. The density, the friction
  speed and the orifice are lagged by one step.
- `orifice_flow.i`: the orifice designed on the slice mass flows (`NTPOrificeFlowDesignAux`,
  THEORY.md 4.5). After every solve each slice's multiplier moves so that the slice carries its
  power share of the flow, at the geometric mean `exp(orifice_level)`. The slice flows are the exact
  Rhie-Chow fluxes through the plenum / cold-frit interface (`NTPRhieChowSliceFlow`). With
  `orifice_share_exponent > 0` the shares are also corrected from the slice outlet temperatures.
  The design columns are the axial slices (`NTPElementColumnsDivision`); in 3-D,
  `orifice_sectors = n_theta` splits each slice into azimuthal sectors (stage 5).
- `orifice_profile.i`: the cold-frit orifice of the Newton anchor, written by
  `extract_orifice.py`.
- `orifice_lagged.i`: an orifice designed on the slice outlet temperatures (`NTPOrificeDesignAux`).
  After every step each slice's multiplier moves towards a uniform outlet temperature, at the
  geometric mean `exp(orifice_level)`. It is frozen at `orifice_freeze_time`.
- `../common/channel_turbulence.i` with `../common/channel_turbulence_rz.i`: the plenum and the
  chamber resolved radially (`n_pl`, `n_ch`), with the `k`-`epsilon` model of OpenPronghorn, wall
  functions on their solid walls, the full compressible viscous stress, and no body-force duct
  friction (THEORY.md 4.4 and 6.3).
- `bed_emissivity`, `bed_dispersion`: the radiative conductivity and the thermal dispersion of the
  packed beds (THEORY.md 4.3); `segregated_static_outputs = csv` writes the static outlet
  temperature `T_mixed_static` next to the stagnation one `T_mixed`.

## Findings

**(a) The Newton orifice does not fit the Navier-Stokes element.** At the design power the mixed
outlet temperature stays near the first-law value (it oscillates around 2980 – 3000 K), but the
middle slices are starved and the fuel there overheats without bound, to several 1e5 K. The
Darcy-model orifice was designed for an inlet-overfed split. The Navier-Stokes split is U-shaped (stage 2), so the two models need
different orifices.

**(b) A lagged design tracks the ramp but has no stable fixed point.** With `orifice_level = 4`
and a damping exponent of 0.3, the march tracks the flow along the ramp. The fuel temperature
stays close to the gas temperature, and the energy balance converges to 4e-5 at the end of the
ramp:

| power | mixed gas temperature | largest fuel temperature | inlet pressure | orifice range |
|---|---|---|---|---|
| 10 % | 384 K | 432 K | 7.20 MPa | 2.1 – 152 |
| 25 % | 872 K | 948 K | 7.49 MPa | 0.90 – 147 |
| 50 % | 1634 K | 1751 K | 7.97 MPa | 0.86 – 140 |

The command freezes the orifice at the end of the ramp (`orifice_freeze_time = t_ramp`). Held at
50 %, the frozen march then loses one slice about 3 s later, and its fuel runs away (above
10 000 K at t = 34 s): the drift of a fixed orifice in the segregated iterations, finding (d).
With the previous face value of the particle conductivity, two more observations were made:

- **With the design kept active**, the orifice ran away after about 10 more steps: its largest
  multiplier grew from 77 to over 300, and the fuel overheated.
- **With the orifice frozen**, the march converged, with the energy balance at 2e-7, to a state
  in which every other slice was starved (slices 0, 2, 3, 6 and 9 carried the flow): a steady
  solution of the frozen orifice, not the orificed element.

A smaller damping exponent (0.1) loses the flow earlier in the ramp, at 20 %. The orifice then
cannot follow the power.

The alternating pattern shows that the slice-by-slice update, lagged by one step, excites a
mode with one slice of wavelength. The axial slices of the element are hydraulically coupled
only through the plenum and the chamber, and nothing in the update couples neighbouring slices.
The Newton design of stage 3 solves the orifice together with the flow and has no such mode.

**(c) An orifice designed on the slice flows reaches the design point.** The temperature design
(b) closes its loop through the density of the hot gas: a slice that heats up expands, its flow
falls, and it heats up further before the next update. The flow design acts on a hydraulic
quantity instead: the flow through a slice falls monotonically with its cold-frit drag, at the
cold density of the frit, and the targets, the slice power shares, do not depend on the solution.
At the defaults (`orifice_level = 4`, `orifice_flow_exponent = 0.7`, shares fixed to the power
shares) both executioners converge to the same state:

| | SIMPLE (`anchor_simple.i`) | PIMPLE (`anchor_pimple.i`) |
|---|---|---|
| path | 50 power steps, 150 solves at 1 GW | ramp over 20 s, steady at t = 33.9 s (685 steps) |
| run time (1 rank) | 95 s | 157 s |
| mixed gas temperature | 2989.81 K | 2989.81 K |
| largest fuel temperature | 3185 K | 3185 K |
| slice outlet temperatures | 2989.0 – 2992.4 K | 2989.0 – 2992.3 K |
| spread / mean rise | 1.2e-3 | 1.1e-3 |
| inlet pressure | 8.875 MPa | 8.874 MPa |
| orifice range | 1.18 – 130 | 1.2 – 129 |
| energy balance error | -2e-8 | -5e-10 |
| mass balance error | 2e-13 | 6e-14 |
| slice flow / target | within 4e-5 | within 4e-8 |

The mixed outlet temperature is the first-law value: with the power-share targets every slice
gets the same enthalpy rise. The Newton anchor (stage 3) gives 2989.4 K and 3193 K with the
Darcy-Forchheimer model and its own orifice. The Navier-Stokes orifice rises monotonically from
the inlet, 1.2 at the inlet slice to 130 at the closed end. The plenum loses static pressure
along its length by only 37 kPa, so the cold frit must throttle the downstream slices against a
chamber whose pressure falls by 800 kPa: the slice pressure drop grows from 1.06 MPa at the inlet
to 1.83 MPa at the closed end, most of it across the cold frit there.

The level matters. The inlet slice needs the least cold-frit resistance, so the level must leave
it room to be opened: at `orifice_level = 3` its multiplier falls to 1e-5, at the bound of the
design, and the fuel overshoots to 3840 K during the ramp. Level 4 keeps every multiplier between
1 and 130 at an inlet pressure only 2 kPa higher. Above about 4.5 the uniform starting orifice
throttles the cold start too hard and the first solves diverge. Neither the fixed-point
iterations nor the time steps show the slice-alternating mode of (b) at any power. The ramp
length matters little: SIMPLE with all the power in the first solve (`n_ramp=1`) reaches the same
state after 150 solves, through a large fuel temperature excursion in the first few.

Correcting the shares further from the slice outlet temperatures
(`orifice_share_exponent > 0`) brings the density feedback back. With `gamma = 0.5` (PIMPLE)
the march settles into a small limit cycle at the design power: spread 0.3 – 2 %, largest fuel
temperature 3145 – 3193 K, the inlet slice oscillating with a period of about 50 updates. With
`gamma = 0.15` and the correction applied at every solve (SIMPLE), the inlet-slice multiplier
oscillates between 0.09 and 1.6 and the fuel reaches 3170 – 3790 K. `gamma = 0.05`, limited to 1 % per solve,
converges slowly to the same state as the power shares. The power shares are the design, and the
temperature correction is not needed on this element.

**(d) The designed orifice, frozen, is a steady solution the iterations do not hold.** After 200
designed solves (`orifice_max_updates = 200`) the state is converged, and freezing the orifice
leaves it unchanged at first: the spread stays at its converged value for 40 – 80 solves. But a
perturbation grows, by a factor of 1.20 per solve at the design power, whatever the level. A
slice near the closed end (9 or 10) loses flow, heats up and runs away. Its neighbour takes the
flow and cools. The same happens at partial power, more slowly at 25 %. PIMPLE behaves the same
way: frozen at t = 34.7 s, the fuel exceeds 8000 K 2.5 s later.

| frozen at | inlet pressure | growth per solve | after 120 frozen solves: spread, largest fuel temperature |
|---|---|---|---|
| level 3, 1 GW | 8.87 MPa | 1.20 | 1.3, 5500 K |
| level 3.5, 1 GW | 9.10 MPa | 1.20 | 11.8, 36 000 K |
| level 4, 1 GW | 9.65 MPa | 1.20 | 0.23, 3650 K |
| level 4.5, 1 GW | 10.56 MPa | 1.20 | 0.014, 3164 K (onset later) |
| level 3, 75 % / 50 % / 25 % | | 1.19 / 1.18 / 1.13 | runaway |
| level 3, relaxation 0.15 | 8.87 MPa | 1.06 | runaway after 150 solves |
| level 4, harmonic particle conduction (the defaults) | 8.88 MPa | 1.18 | 9.5, 30 000 K |

The rows above the last were computed before the face conductivity of the particle conduction was
made harmonic (THEORY.md 6.2); the last row, with the current inputs, shows the same drift.

A higher level only delays the onset. Halving the lagged-state relaxation slows the growth per
solve without removing it: the lagged density map amplifies the mode. Is it the flow-excursion
instability of heated parallel channels, or a mode of the iteration? The segregated pressure
equation has no compressibility, so its march is not the physical transient. Stage 3 answers it
with the Newton model, which has the gas mass storage: frozen at the design point and kicked, the
element returns to its steady state with a 2.0 s time constant (`../03_ludewig_anchor/stability.i`).
The growing mode here is therefore a mode of the lagged-density iteration, not of the element.
With the design active, the flow loop acts as a flow controller of every slice, and the converged
state is stable for SIMPLE and PIMPLE.

**(e) Turbulent free channels.** With one cell across the plenum (2 mm) and the chamber, every
channel cell touches a wall and a `k`-`epsilon` model has no bulk to act on; their wall shear is
the duct-friction closure. Resolved by 8 cells across, the channels carry the `k`-`epsilon` model of
OpenPronghorn with wall functions (`common/channel_turbulence.i`, THEORY.md 6.3). Three things had
to be right for it to converge:

- the strain rate of the channel flow, without the jump of the superficial velocity into the
  porous cells (`NTPChannelVelocity`): with it, `mu_t` ran to its limiter;
- turbulent diffusion of momentum and enthalpy on every block with harmonic faces: restricted to
  the channels, it injected four to eight times the fission power into the gas;
- the flowing start (`segregated_plug_velocity`), in RZ too.

Cold, the inlet pressure converges with the radial mesh: 3.60, 3.82 and 4.03 kPa above the outlet
for 4, 8 and 16 cells across (one cell with the duct friction: 11.6 kPa), mass conserved to 1e-9.
At the design power the solves settle from 550 iterations during the ramp to 180, and the state is
steady over the last 20 of 200 solves (fuel within 0.1 K):

| | laminar, 1 cell, duct friction | turbulent, 8 cells | turbulent + bed radiation and dispersion |
|---|---|---|---|
| mixed outlet (stagnation) | 2989.81 K | 2989.81 K | 2989.81 K |
| mixed outlet (static) | 2976.3 K | — | 2976.1 K |
| largest fuel temperature | 3184.8 K | 3185 K | 3182 K |
| inlet pressure | 8.875 MPa | 8.811 MPa | 8.851 MPa |
| slice outlet spread / rise | 1.2e-3 | 1.1e-3 | 8e-4 |
| orifice range | 1.18 – 130 | 1.21 – 129 | 4.0 – 112 |
| energy balance error | -2e-8 | -4e-9 | -1e-8 |
| run time (1 rank) | 95 s | 28 min | 34 min |

The turbulent channels leave the fuel temperature unchanged and lower the inlet pressure by
64 kPa: the turbulent wall shear of the resolved plenum is below the duct friction of the one-cell
closure (`f = 0.02`). The bed radiation and dispersion lower the fuel temperature by 2.9 K and make
the slice outlets more uniform (THEORY.md 4.3); the orifice then needs less range (4.0 – 112), at
40 kPa more inlet pressure. The kinetic energy of the gas leaving at 700 m/s is 13.5 K of the outlet
temperature.

## Next steps

- **A frozen orifice in the segregated iterations.** The frozen design is physically stable
  (stage 3, `stability.i`), but the lagged-density iteration does not converge to it. A pressure
  equation with the compressibility of the gas (a weakly compressible PIMPLE) would make the
  march a physical transient.
- **3-D.** Done in stage 5: the 3-D element needs the design per column (slice and sector), as
  rings it loses its azimuthal symmetry in a flow excursion between the columns of a slice.
- **Cost of the turbulent model.** At power each SIMPLE solve stops on the pressure tolerance
  (1e-8), after 180 – 550 iterations; the fixed-point loop does not need that much. A looser
  tolerance during the ramp would shorten the 30-minute run.
