# Case 14: radiolysis-driven corrosion -- an irradiated chloride pool between 316 walls

**Module 4 (coupled) | ~25 min | 1D, serial, 51 s**

## Objective

Couple the two actions: wall corrosion (`[CorrosionPlatingFlow]`) releases Cr straight into the
radiolysis Cr(II) of `[MoltenSaltRadiolysis]`, and both evolve through a power cycle.

## Physics

A molten chloride pool (600 C) between two 316 walls, 20 cm apart:

```
wall:  Cr(alloy) <-> Cr(II) + 2 e-         net rate ~ 1 - c_Cr(II),wall / c_eq
salt:  Cl2- + Cr(II)  -> Cr(III) + 2 Cl-   (k = 7.2e6)
       Cl2- + Cr(III) -> Cr(II) + Cl2      (k = 1.4e6)
```

The salt starts **pre-equilibrated** with the walls: Cr(II) = $c_{eq}$ = 14.3 mol/m^3 (~470 wppm).
Without irradiation the net current is zero and the salt is passive (case 10).

Under irradiation the radiolytic oxidant converts Cr(II) to Cr(III), reaching 84% Cr(III)
(= 7.2/(7.2 + 1.4)). The wall then sees Cr(II) far below $c_{eq}$ and corrodes again. **The alloy is
the ultimate sink of the radiolytic oxidant.** After shutdown nothing reduces the Cr(III), so the
walls keep corroding until Cr(II) is back at $c_{eq}$: the salt remembers the irradiation.

## The input, step by step

[irradiated_pool.i](irradiated_pool.i):

- `Cr_init = ${fparse c_ref * exp(2 F / (R T) * eta)}`: the initial condition is the Nernst value,
  computed in the input.
- `[MoltenSaltRadiolysis]` with `dose_rate = dose_rate`, a function following the power history.
  The fast e_sol/Cl recombination is lumped into a direct Cl2- yield
  (`g_value_overrides = '0 0 0.30'`), which keeps the 11-day transient affordable.
- `[CorrosionPlatingFlow]` with `release_variables = 'Cr_II'`: the corrosion flux is a boundary
  source of the radiolysis variable. `reaction_boundary = 'left right'` makes both walls corrode.
- `diffusivity = ${D_mix}` = $10^{-5}$ m^2/s is an eddy diffusivity for the natural-convection
  mixing (mixing time $L^2/D$ ~ 1 h). Case 15 resolves the flow instead.
- `corrosion_rate_um_y` is the wall loss inferred from the growth of the dissolved Cr:
  $(V/A)\,\dot c_{Cr}\,M/\rho$.

## Run

```bash
open_pronghorn-opt -i irradiated_pool.i
python plot.py
```

## Expected results

![](../figures/case14_irradiated_pool.png)

- Before start-up: `c_eq_ratio` = 1 and no corrosion.
- At power: Cr(III) = 84% within a day, and the walls corrode at **35 µm/y**, 0.9x the 39.6 µm/y
  fresh-salt rate.
- After shutdown: 33 µm/y for days. Cr(III) is frozen and Cr(II) slowly recovers.
- Cl2 accumulates (64 mol/m^3 after 5 days) because the pool has no off-gas; case 04 shows the
  Henry-law sink.

## Exercises

1. `dose_full=0`: what happens, and why is this the right control case?
2. Poor mixing: `D_mix=1e-7` (`Executioner/end_time=172800`). Does the corrosion rate change? Why?
3. Where does the Cl2 come from, and what removes it in a real reactor?

<details><summary>Solutions</summary>

1. Nothing happens: Cr(II) stays at 14.3, `c_eq_ratio` = 1 and the corrosion rate is 0. The walls
   and the salt are in equilibrium, so any corrosion in the base case is caused by the radiation.
2. After 2 days the rate is 35.4 µm/y, the same as with good mixing. The Cl2- oxidant is produced
   everywhere, so the wall boundary layer sees the same Cr(II) depletion: the coupling is local,
   not transport-limited.
3. From Cl2- + Cr(III) -> Cr(II) + Cl2: Cr(III) catalyses Cl2- -> Cl2 + Cl-. In a reactor the Cl2
   is stripped to the off-gas (case 04) or reduced by a redox buffer.

</details>
