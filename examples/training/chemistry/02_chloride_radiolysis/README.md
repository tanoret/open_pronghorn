# Case 02: the chloride radical chain in irradiated LiCl-KCl

**Module 1 (radiolysis) | ~20 min | 4 cells, serial, 6 s**

## Objective

Solve a stiff radical network, see how the segregated linear-FV systems are coupled by a fixed-point
iteration, and use a 1D mesh as a parameter scan.

## Physics

Gamma radiolysis of a chloride melt (LiCl-KCl eutectic, 400 C) produces solvated electrons and
chlorine atoms, with G = 0.3 for each:

```
Cl + Cl-       -> Cl2-       (fast, pseudo-first order in 2e4 mol/m^3 chloride)
e_sol + Cl     -> Cl-
2 Cl           -> Cl2
2 Cl2-         -> Cl3- + Cl- (second order: sets the radical steady state)
Cl3-           -> Cl2 + Cl-  (k = 1 1/s)
```

With $S = G D/(100\,\mathrm{eV}\,N_A)$ the radical source, the analytic limits are:

| species | limit | scales as |
|---|---|---|
| Cl2- | $\sqrt{S/(2k)}$ | $D^{1/2}$ |
| Cl3- | $S/2$ | $D$ |
| e_sol | $S\,t$ (no scavenger, accumulates) | $D$ |
| Cl2 | $S\,t/2$ (the stable product) | $D$ |

## The input, step by step

[chloride_radiolysis.i](chloride_radiolysis.i):

- `[Mesh]` has 4 cells and **no transport**: each cell is an independent 0D sample.
- `[Functions] dose_rate = 10^(2 + floor(x))` gives cell $i$ a dose rate of $10^{2+i}$ W/m^3, so
  one run is a dose scan from $10^2$ to $10^5$ W/m^3.
- `[MoltenSaltRadiolysis] salt_type = chloride` loads the network above.
- `[Executioner]` has `multi_system_fixed_point = true` with an `IterationCountConvergence` of 20
  sweeps. Each species is its own linear system, and the reaction coupling between them is
  converged by sweeping all the systems several times per time step.
- The time step starts at $10^{-9}$ s and grows geometrically. The chemistry has time scales from
  nanoseconds (Cl + Cl-) to seconds (Cl3- decay).
- `[VectorPostprocessors] dose_scan` writes one row per cell at the end of the run.

## Run

```bash
open_pronghorn-opt -i chloride_radiolysis.i
python plot.py
```

## Expected results

![](../figures/case02_chloride_radiolysis.png)

- Left (cell 0): Cl reaches its steady state within microseconds, then Cl2- and Cl3-. e_sol and Cl2
  grow linearly.
- Right: fitted slopes of **0.50** (Cl2-) and **1.00** (Cl3-). Cl3- is $S/2$ to 4 digits in every
  cell. Back-computing $S/c_{Cl_2^-}^2$ gives $2k = 4.4\times10^6$ m^3/(mol s), the
  Cl2- recombination constant of the database.

## Exercises

1. Why does Cl2- scale as $D^{1/2}$ but Cl3- as $D$?
2. e_sol keeps growing. Which species would remove it in a real salt? (Hint: case 03.)
3. Reduce the fixed-point sweeps (`Convergence/fp/max_iterations=2`) and compare. What goes wrong?

<details><summary>Solutions</summary>

1. Cl2- is removed by a second-order reaction: $2k c^2 = S$ gives $c = \sqrt{S/2k}$. Cl3- is
   produced by that reaction at the rate $S/2$ and decays first-order at $k = 1$ 1/s, so
   $c_{Cl_3^-} = S/2$, linear in $D$.
2. Any reducible species: dissolved metals such as Cr(III) (case 03) or U(IV) (case 04) scavenge
   the electrons.
3. With too few sweeps per step the coupled species lag each other by a time step. The radical
   concentrations overshoot and the invariants (case 03) are no longer conserved.

</details>
