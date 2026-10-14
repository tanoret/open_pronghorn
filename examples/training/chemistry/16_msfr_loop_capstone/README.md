# Case 16 (capstone): the MSFR primary loop

**Module 4 (coupled) | ~45 min | 2D conjugate, 24 ranks, 1 min 35 s (heavy test, 8 ranks)**

## Objective

Put the whole day together in one reactor-grade calculation. A closed, pump-driven fuel-salt loop
with conjugate heat transfer carries radiolysis, wall corrosion and the fission-product inventory
with xenon stripping and heat-exchanger plating. Validate the 2D result against the 0D cases of the
morning.

## Physics

A NaCl-UCl3 fuel-salt loop around a central reflector, inside a reflector and vessel
(`examples/msfr_loop`):

- k-epsilon flow driven by a **pump** (a momentum source in the bottom leg) and buoyancy,
- **conjugate heat transfer** with the solid reflector and vessel,
- fission heat in the core leg (left) and heat removal in the heat-exchanger leg (right),
- **chloride radiolysis** with Cr(II)/Cr(III) (module 1, case 14),
- **Butler-Volmer corrosion** of the 316 walls releasing Cr, Fe and Ni (module 3),
- **fission products** (module 2) with fast-spectrum yields: I-135 -> Xe-135 **stripped only in
  the top leg**, and Mo-99, Ru-103 and Te-132 **plating only on the heat-exchanger walls**, with
  the Chilton-Colburn $k_m$ of case 15.

**Why steady state is enough.** The loop time (~5 s) is much shorter than the half-lives (Xe-135
9.1 h, Mo-99 66 h). The inventories are nearly uniform around the loop, and their equilibrium is
the 0D prediction of cases 07 and 08 with the loop-integrated removal rates:

$$\bar N_{Xe} = \frac{P_{135}/V}{\lambda_{Xe} + \int r\,dV/V}, \qquad
  \bar N_{Mo} = \frac{P_{Mo}/V}{\lambda_{Mo} + \int k_m\,dA/V}$$

**Corrosion in a closed loop.** There is no outlet, so the dissolved metals settle at the wall
(Nernst) equilibrium, $c \sim c_{ref}\,e^{zf\eta} \approx 30$ mol/m^3 per metal. The radiolytic
oxidant holds Cr(III)/Cr(II) at ~5. Hot-wall dissolution then balances cold-wall plating: the
temperature-gradient mass transfer of case 13.

## The input, step by step

[msfr_loop_capstone.i](msfr_loop_capstone.i):

- `[Mesh]`: a 60x60 square with `reflector` and `vessel` subdomains, and `inner_wall`/`outer_wall`
  sidesets on the fluid-solid interfaces.
- The flow and heat blocks are as in `examples/msfr_loop`. `energy_balance_error` compares the
  fission power with the heat-exchanger removal.
- `[FissionProductInventory]`:

  ```
  block = 'salt'
  nuclides = 'I135 Xe135 Mo99 Tc99m Ru103 Te132'
  yield_set = U235_fast
  power_density = fission_shape          # the same functor as the energy equation's heat source
  removal_groups = 'noble_gas'
  removal_rates = 'strip_rate'           # a function: r_sparge in the top leg, 0 elsewhere
  removal_destinations = 'offgas'
  plating_boundaries = 'inner_wall outer_wall'
  plating_mass_transfer_coefficient = plate_coefficient   # k_m on the HX walls, 0 elsewhere
  ```

  Removal is placed with **spatial functors**, not with blocks.
- The Xe-135 burnout is left out (`neutron_flux` = 0): the database's 2.6e6 b is a thermal value
  and does not apply to a fast spectrum.
- **Initial conditions of the metals.** The segregated iterations approach the closed-loop
  equilibrium very slowly, because the wall exchange is the only coupling. The salt therefore starts
  from the converged loop averages of a previous run (`initial_condition_values`,
  `initial_concentration`), and the metals and fission products use relaxation 1.

## Run

```bash
mpirun -n 24 open_pronghorn-opt -i msfr_loop_capstone.i                 # ~1.5 min
python plot.py
# mesh sensitivity (~3.5 min)
mpirun -n 24 open_pronghorn-opt -i msfr_loop_capstone.i Mesh/gen/nx=90 Mesh/gen/ny=90 Outputs/file_base=fine
```

## Expected results

![](../figures/case16_msfr_loop.png)

| quantity | 60x60 | 90x90 |
|---|---|---|
| T min / max [K] | 970.8 / 1036.1 | 970.9 / 1035.4 |
| energy balance error | 5.9e-5 | 2.3e-4 |
| core-leg velocity [m/s] | 1.56 | 1.68 |
| Xe-135 / 0D prediction | 1.00002 | 1.00002 |
| Mo-99 / 0D prediction | 1.0004 | 1.0005 |
| Xe-135 stripped fraction | 0.826 | 0.826 |
| Mo-99 plated fraction | 0.939 | 0.939 |
| Xe-135 max/min - 1 | 6.6e-4 | 6.6e-4 |
| Cr(III)/Cr(II) | 5.14 | 5.14 |
| Fe, Ni [mol/m^3] | 30.38 | 30.45 |

- The 2D inventories reproduce the 0D predictions to $4\times10^{-4}$ and are uniform around the
  loop to $7\times10^{-4}$.
- 83% of the Xe-135 is stripped (the rest decays), and 94% of the Mo-99 ends up on the
  heat-exchanger walls.
- The velocity is the most mesh-sensitive quantity (8%, the k-epsilon wall treatment). The
  inventories and balances are insensitive.

## Exercises

1. Weaker sparging (`r_sparge=3e-5`). Predict the Xe-135 inventory with the 0D formula before
   running.
2. Let the noble metals plate on every wall of the loop:
   `"Functions/plate_coefficient/expression=\${k_plate}"`. What happens to the plated fraction and
   to the salt inventory, and why does the HX still matter for the source term?
3. Why can the loop-average Xe-135 equal a 0D model to $2\times10^{-5}$ even though the stripping
   happens only in the top leg?
4. (Open) What would be needed to follow a load-follow transient (case 09) in this loop?

<details><summary>Solutions</summary>

1. From the base stripped fraction, $\int r\,dV/V = 0.826/0.174\,\lambda_{Xe} = 4.75\,\lambda_{Xe}$.
   A 10x weaker sparge gives $0.475\,\lambda_{Xe}$, so $\bar N_{Xe}$ rises by
   $(1 + 4.75)/(1 + 0.475)$ = **3.90x**. The run gives 3.90x ($6.85\times10^{-4}$ ->
   $2.67\times10^{-3}$), still `Xe135_vs_0d` = 1, with 32% stripped instead of 83%.
2. The plated fraction rises only from 94% to **98%** (the removal was already much faster than
   decay), but the salt inventory falls 2.9x ($1.57\times10^{-3}$ -> $5.5\times10^{-4}$, still
   `Mo99_vs_0d` = 1). The location changes: the plated activity (and its decay heat) is spread
   over the core and piping instead of concentrated in the heat exchanger.
3. The loop time (~5 s) is $10^{-4}$ of the Xe-135 lifetime. A xenon atom makes thousands of
   passes through the sparged leg before it decays, so it only sees the volume-averaged rate.
4. A transient executioner (PIMPLE), or a frozen-flow restart of the scalars with
   `time_derivative = true`: the inventory time scale (hours) is much longer than the flow time
   scale (seconds).

</details>
