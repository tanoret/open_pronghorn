# Case 03: chromium redox in an irradiated chloride salt

**Module 1 (radiolysis) | ~20 min | 1 cell, serial, 14 s**

## Objective

Add a dissolved corrosion product to the radical chain of case 02 and verify the two exact
invariants of the chemistry: chromium conservation and the redox (electron) balance.

## Physics

Chromium is the main corrosion product of Ni-Cr alloys and stainless steels. The radiolysis
products cycle its oxidation state:

```
e_sol + Cr(III) -> Cr(II)           e_sol + Cr(II) -> Cr(I)
Cl2-  + Cr(II)  -> Cr(III) + 2 Cl-  Cl2-  + Cr(III) -> Cr(II) + Cl2
Cr(III) + Cr(I) -> 2 Cr(II)
```

The Cr(III)/Cr(II) ratio sets the redox potential of the salt, which controls how corrosive it is.

Invariants:

- **Chromium**: Cr(I) + Cr(II) + Cr(III) = Cr0 (`cr_error`).
- **Redox balance**: radiolysis makes electrons and holes in equal numbers. The reducing equivalents
  (e_sol, Cr(I)) must equal the oxidizing ones (Cl, Cl2-, 2 Cl3-, 2 Cl2, Cr(III)).
  `redox_balance` is the difference, relative to the electrons produced.

## The input, step by step

[cr_redox.i](cr_redox.i):

- `metals = 'Cr'` in `[MoltenSaltRadiolysis]` adds the three Cr species and their reactions.
- `initial_condition_species = 'Cl_ion Cr_II'`: 1 mol/m^3 of Cr(II) (~50 wppm), all reduced.
- `n_fp = 30` fixed-point sweeps per step. The Cr reactions couple 9 linear systems, so more sweeps
  are needed than in case 02.
- The invariants are `ParsedPostprocessor`s, written out term by term so you can read the
  bookkeeping.

## Run

```bash
open_pronghorn-opt -i cr_redox.i
python plot.py
```

## Expected results

![](../figures/case03_cr_redox.png)

- The salt **stays reduced**: Cr(III)/Cr = $4.3\times10^{-6}$ after 1000 s. e_sol reduces Cr(III),
  and Cr(III) + Cr(I) comproportionate, as fast as Cl2- oxidizes Cr(II). Meanwhile Cl2 slowly
  accumulates.
- `cr_error` stays below $5\times10^{-7}$ and `redox_balance` below $2\times10^{-6}$. Both are
  limited by the fixed-point tolerance, not by the model.

## Exercises

1. Increase the dose 100x (`dose_rate=1e6`). Does the salt become oxidizing?
2. Reduce the sweeps (`n_fp=5`). Which invariant degrades first, and why?
3. Most electrons recombine with Cl atoms. Use the outputs to estimate what fraction escape
   recombination and end up as reduced Cr or Cl2.

<details><summary>Solutions</summary>

1. Only slightly. Cr(III)/Cr rises 10x, to $4.2\times10^{-5}$, while Cl2 grows 1000x, to
   $2.5\times10^{-4}$ mol/m^3. The electrons keep the chromium reduced, so the oxidizing
   equivalents pile up as Cl2. `cr_error` grows to $5\times10^{-4}$ at the larger dose with
   30 sweeps; more sweeps or smaller steps restore it.
2. With `n_fp=5`, `cr_error` hardly changes ($5.2\times10^{-7}$), but `redox_balance` grows 8x, to
   $1.3\times10^{-5}$. Every Cr reaction moves chromium between its own three species, so the Cr
   sum closes regardless. The redox balance sums radicals with very different time scales, so any
   lag between the segregated systems shows up there first.
3. At 1000 s the electrons produced are 0.311 mol/m^3, against Cl2 at $2.6\times10^{-7}$ and
   Cr(I) and Cr(III) at $4\times10^{-6}$. Only about $10^{-5}$ of the electrons survive; the rest
   recombine.

</details>
