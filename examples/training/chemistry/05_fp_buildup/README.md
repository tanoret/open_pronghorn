# Case 05: fission-product build-up -- saturation versus growth

**Module 2 (fission products) | ~15 min | 1 cell, serial, 3 s**

## Objective

Meet the `[FissionProductInventory]` action: database-driven nuclides, fission yields and decay.
See the two limits of the Bateman equation.

## Physics

A well-mixed fuel salt at a loop-averaged 20 MW/m^3 (a thermal MSR: 8 MW in ~0.4 m^3 of salt) for
one year. For a nuclide with no tracked precursor, no capture and no removal:

$$\frac{dN}{dt} = \frac{Y F}{N_A} - \lambda N \quad\Rightarrow\quad
  N(t) = \frac{Y F}{N_A \lambda}\left(1 - e^{-\lambda t}\right), \qquad F = \frac{P}{200\,\mathrm{MeV}}$$

- **I-131** (8 d) saturates within weeks. Its activity $\lambda N \to Y F$ measures the *power*.
- **Cs-137** and **Sr-90** (~30 y) grow linearly. They measure the *burnup*.

## The input, step by step

[fp_buildup.i](fp_buildup.i):

```
[FissionProductInventory]
  nuclides = 'I131 Cs137 Sr90 Kr85'
  power_density = ${power_density}   # converted to fissions with 200 MeV/fission
  verbose = true                     # prints yields, decay constants and the system order
[]
```

- The action creates one linear-FV variable and system per nuclide, named `<nuclide>_sys`. Copy the
  `linear_sys_names` line that `verbose` prints into `[Problem]` and the executioner.
- The yields are ENDF/B-VIII.0 U-235 thermal cumulative yields, and half-lives are from NUBASE2020
  (`data/fission_product_database.json`). Untracked precursors are lumped in automatically (the
  `effective` yield).
- The `*_analytic` postprocessors hold the formula above.

## Run

```bash
open_pronghorn-opt -i fp_buildup.i
python plot.py
```

## Expected results

![](../figures/case05_fp_buildup.png)

- I-131 matches the analytic curve to $10^{-13}$ once saturated, and its activity equals
  `I131_saturation_activity` = $Y F$.
- After one year Cs-137 = 2.000 mol/m^3, within $2\times10^{-7}$ of the analytic value.
- Kr-85 (10.7 y) is in between: 0.090 mol/m^3 after a year.

## Exercises

1. Run with `Executioner/dt=21600` (6 h steps instead of 1 day). How does the I-131 error during
   the first week change? What does that say about the implicit-Euler time step?
2. The salt is drained after one year. How long until the I-131 activity falls by 1000x? And the
   Cs-137 activity?
3. Read the `verbose` output. Why is the I-131 yield used here (0.0289) larger than its
   independent yield?

<details><summary>Solutions</summary>

1. With 1-day steps the I-131 error at day 7 is $5\times10^{-3}$; with 6-h steps it is
   $3\times10^{-4}$. Implicit Euler is first order, but the error decays away once the
   inventory saturates (both runs are within $10^{-12}$ at one year).
2. $\ln(1000)/\lambda$ = 10 half-lives: 80 days for I-131, 300 years for Cs-137.
3. The `effective` yield is the *cumulative* yield minus the part fed by tracked parents. Here no
   parent (Te-131) is tracked, so the yield is cumulative and includes all of the short-lived
   precursor chain.

</details>
