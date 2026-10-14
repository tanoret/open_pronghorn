# Case 15: a heated turbulent core channel -- corrosion plumes and noble-metal plating

**Module 4 (coupled) | ~30 min | 2D, 24 ranks, 12 s (4 ranks in the test)**

## Objective

The first flow-resolved case: run all three chemistry actions on a turbulent, heated flow solved by
SIMPLE, and close the energy and fission-product balances.

## Physics

A 2 m fuel-salt channel (FLiBe-like properties, 0.5 m/s upward, k-epsilon) between two Hastelloy N
walls, heated by fission with a cosine flux shape. One segregated SIMPLE solve gives the flow,
turbulence and temperature, and on that field:

- **radiolysis**: F2 produced with the flux shape and partly evolved as gas (case 04),
- **corrosion**: Arrhenius Butler-Volmer release of Cr, Fe and Ni at the local wall temperature
  (cases 10–12), carried up the channel as plumes,
- **fission products**: Mo-99, Ru-103 and Te-132 born with the fission rate and plating on the
  walls with $J = k_m N$. Case 08 treated the same process as a 0D rate $r_p = k_m A/V$.

Turbulence matters twice. $\mu_t$ sets $D_{eff} = D + \mu_t/(\rho\,Sc_t)$, which carries the
corrosion products from the wall into the core and the fission products from the core to the wall.

$k_m$ is not a free parameter. It comes from the Chilton-Colburn correlation
$Sh = 0.023\,Re^{0.8}Sc^{1/3}$: with $Re = 6\times10^5$ and $Sc = 500$, $k_m = 1.3\times10^{-5}$
m/s.

## The input, step by step

[core_channel.i](core_channel.i) builds on `examples/flowing_msr_corrosion`:

- The flow, k-epsilon and energy blocks are standard linear-FV Navier-Stokes. `[GlobalParams]` sets
  `rhie_chow_user_object` and `advected_interp_method_name` once for every block.
- `[MoltenSaltRadiolysis]`, `[CorrosionPlatingFlow]` and `[FissionProductInventory]` all take
  `rhie_chow_user_object = 'rc'` and `diffusivity = D_eff`, so they are advected by the solved flow.
  `time_derivative = false` because SIMPLE is steady.
- `plating_boundaries = ${walls}` with `plating_mass_transfer_coefficient = ${k_plate}` adds the
  Robin plating BC. Inlet (clean salt) and outlet BCs are given for each nuclide.
- Every species system appears in `passive_scalar_systems`, each with a relaxation and a tolerance.
- Balances: `energy_balance_error` (fission power against $\rho c_p \int v(T - T_{in})\,dA$ at the
  outlet) and `Mo99_balance_error` (production = outlet + plated + decayed).

## Run

```bash
mpirun -n 24 open_pronghorn-opt -i core_channel.i
python plot.py
# mesh sensitivity
mpirun -n 24 open_pronghorn-opt -i core_channel.i Mesh/core/nx=12 Mesh/core/ny=48 Outputs/file_base=fine
```

## Expected results

![](../figures/case15_core_channel.png)

| quantity | 6x24 | 12x48 |
|---|---|---|
| T outlet [K] | 932.3 | 932.1 |
| energy balance error | 2.8e-4 | 3.4e-4 |
| Mo-99 balance error | 4.3e-4 | 4.9e-4 |
| max corrosion rate [µm/y] | 6.97 | 6.98 |
| Cr at outlet [mol/m^3] | 9.01e-7 | 9.28e-7 |
| Mo-99 plated fraction | 1.41e-4 | 1.26e-4 |

- The wall attack is fastest at the hot top (Arrhenius).
- Cr is high at the walls and diluted in the core. Mo-99 is depleted near the walls.
- Only ~$10^{-4}$ of the Mo-99 produced plates during the 4 s transit: plating is a slow,
  many-pass process. Case 16 closes the loop.
- The outlet quantities are mesh-converged to ~3%. The plated fraction, a wall quantity, changes
  by 10%.

## Exercises

1. Halve the flow (`u_in=0.25`). Predict the outlet temperature and the plated fraction, then run.
2. Plating needs Mo-99 at the wall. What limits it here, $k_m$ or the turbulent transport to the
   wall? Try `FissionProductInventory/plating_mass_transfer_coefficient=1e-4`.
3. Estimate the plated fraction of case 08 ($r_p = 5\times10^{-4}$ 1/s) from this channel: what
   $A/V$ and residence time would give $r_p$?

<details><summary>Solutions</summary>

1. The temperature rise doubles: $T_{out}$ = 964.6 K (from 932.3 K). Both the transit time and the
   plating fraction grow, but $k_m$ falls as $u^{0.8}$ (Sherwood), so the plated fraction only
   rises to $1.6\times10^{-4}$. The hotter wall corrodes at 8.4 µm/y.
2. With $k_m$ 8x larger the plated fraction barely moves ($1.41\times10^{-4}$), and the Mo-99
   balance error grows to $1.4\times10^{-3}$ as the wall gradient steepens. The turbulent
   transport to the wall is the bottleneck, not the surface kinetics.
3. $A/V = 1/W$ = 3.3 1/m per wall, so $r_p = k_m A/V \approx 2 \cdot 1.3\times10^{-5} / 0.3$
   = $8.5\times10^{-5}$ 1/s for the channel. Case 08's $5\times10^{-4}$ 1/s corresponds to a heat
   exchanger with ~6x the area per volume (tube bundles).

</details>
