# Stage 2: one fuel element, two solution methods

This stage sets up the axisymmetric model of one Ludewig Design 1 fuel element. It solves the
element with both solution methods of the tutorial, cold and at 5 % of the design power, before
stage 3 takes it to 1 GW. The equations are in [THEORY.md](../THEORY.md).

## 2.1 The mesh

```bash
open_pronghorn-opt -i mesh.i --mesh-only --allow-unused
python3 plot_mesh.py               # mesh.png
```

`mesh.i` includes `../common/ludewig_design1.i` (radii, porosities, particle sizes) and
`../common/element_mesh.i`. That file builds a reference strip with one column per region and maps
its nodes onto the cones of the element. All five regions taper linearly from the inlet end
(`z = 0`) to the chamber exit (`z = L = 729.4 mm`).

The resolution is set by `nz` (axial cells) and `n_ch n_hf n_bed n_cf n_pl` (radial cells per
region). The default is 12 axial cells, one cell per region. The boundaries are:

- `inlet` (the plenum annulus at `z = 0`);
- `outlet` (the chamber at `z = L`);
- `closed_ends`;
- `axis`;
- `outer_wall`.

![mesh](mesh.png)

## 2.2 Newton Darcy-Forchheimer model

The model is `../common/element_physics.i` (THEORY.md 6.1). The stage-3 anchor input runs it
with a steady power continuation:

```bash
open_pronghorn-opt -i ../03_ludewig_anchor/anchor.i ../03_ludewig_anchor/continuation.i \
                   Executioner/end_time=0   Outputs/file_base=newton_cold        # cold flow
open_pronghorn-opt -i ../03_ludewig_anchor/anchor.i ../03_ludewig_anchor/continuation.i \
                   Executioner/end_time=0.1 Outputs/file_base=newton_5pc         # 5 % power
```

## 2.3 Segregated Navier-Stokes model (SIMPLE)

`element_simple_rz.i` assembles the segregated model from these pieces:

- `../common/element_physics_segregated.i`: physics and closures;
- `element_segregated_rz.i`: the axisymmetric part;
- `element_simple.i`: the SIMPLE executioner;
- `element_segregated_postprocessors.i`: exact balances from the Rhie-Chow face fluxes.

The model is THEORY.md 6.2. It solves the superficial-velocity momentum equation in all regions,
keeping the inertia and wall friction of the plenum and the chamber.

```bash
open_pronghorn-opt -i element_simple_rz.i Outputs/file_base=simple_cold          # cold flow
open_pronghorn-opt -i element_simple_rz.i s_power=0.05 simple_fixed_point=40     # 5 % power
```

Each SIMPLE solve converges the flow, the gas energy and the particles with the lagged density
and friction state held fixed. The fixed-point iterations then update that state.

## 2.4 Results

```bash
python3 compare_cold_flow.py       # cold_flow.png, from newton_cold.e and simple_cold.e
```

| | Newton (Darcy-Forchheimer) | SIMPLE (Navier-Stokes) |
|---|---|---|
| cold: inlet pressure − 7 MPa | 16.2 kPa | 11.6 kPa |
| cold: mass balance (out − in)/in | −9e-7 | 5e-13 |
| 5 %: mixed outlet temperature | 218.9 K | 218.9 K |
| 5 %: energy balance error | 8e-8 | −4e-8 |
| 5 %: largest fuel temperature (uniform frits) | 325 K | 1003 K |

![cold flow](cold_flow.png)

The two models share their closures and the plenum wall friction (`f_ann = 0.02`). They differ
only in the plenum and chamber momentum:

- **Darcy-Forchheimer model.** The plenum pressure falls monotonically along the element, from
  friction. The inlet end, where the frits see the largest pressure difference, is overfed.
- **Navier-Stokes model.** The decelerating plenum flow recovers static pressure towards the
  closed end. The flow split becomes U-shaped, and the middle of the element is the least fed.
  With uniform frits that region overheats first, which is why the fuel temperature at 5 % is
  much higher than in the Darcy model.

Neither split is acceptable for a 3000 K outlet. Stage 3 adds the orifice that corrects it.

## Numerical notes

Three treatments make the segregated model work across the porosity jumps:

- **Series face coupling.** The packed beds and frits are one cell thick next to free channels
  whose momentum diagonal is ten orders of magnitude smaller. `NTPRhieChowMassFlux` gives every
  face the series mobility of its two half cells and corrects the cell pressure gradient
  accordingly (`test/tests/ntp/rhie_chow_series`). With MOOSE's linear interpolation, a one-cell
  frit loses almost all of its resistance.
- **Lagged density.** The density is lagged (THEORY.md 6.2). Updated inside the pressure–velocity
  iterations, it makes them diverge or cycle.
- **Partial elimination.** The gas-to-particle exchange is solved with partial elimination
  (`NTPLinearFVEnthalpyExchange`). At 5 % power with the flow frozen, the gas and particle
  energies reach a residual of 1e-8 in 166 iterations with it; without it they are still at 5e-6
  after 3000.
