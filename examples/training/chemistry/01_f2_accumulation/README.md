# Case 01: radiolytic F2 in a fluoride salt -- the G value and its units

**Module 1 (radiolysis) | ~15 min | 1 cell, serial, 2 s**

## Objective

Write the smallest possible chemistry input, connect a G value to a production rate by hand, and
check the result against the analytic solution.

## Physics

One well-mixed cell of fluoride salt under a constant gamma dose rate $D$. Only the direct F2 yield
is kept (the solvated-electron and F-atom G values are overridden to zero), so dissolved F2 has a
constant source and no sink:

$$\frac{dC_{F_2}}{dt} = S, \qquad S = \frac{G\,D}{100\,\mathrm{eV}\cdot N_A}\quad[\mathrm{mol/(m^3\,s)}]$$

with $G$ in molecules per 100 eV, $D$ in W/m^3 (= J/(m^3 s)) and $N_A$ Avogadro's number.
$D = 10^6$ W/m^3 is the order of the gamma heating in an MSR core, a few percent of a
$\sim 10^8$ W/m^3 fission power density.

## The input, step by step

[f2_accumulation.i](f2_accumulation.i):

- `[Mesh]`: a 1D mesh with one element. Every chemistry input is a mesh, even a 0D one.
- `[Problem] linear_sys_names`: the radiolysis action solves each species as its own **linear
  finite-volume system**, so every species needs a system name. They are listed in `[Problem]` and
  again in `[Executioner] system_names`.
- `[MoltenSaltRadiolysis]`: `salt_type = fluoride` loads the fluoride reaction network and its G
  values from the database. `g_value_species`/`g_value_overrides` switch the radical yields off.
  `time_derivative = true` makes it a transient problem.
- `[Postprocessors]`: `analytic` is $S\,t$, written with a `ParsedPostprocessor`, and `rel_error`
  compares it with the model.

## Run

```bash
cd examples/training/chemistry/01_f2_accumulation
open_pronghorn-opt -i f2_accumulation.i
python plot.py            # writes ../figures/case01_f2_accumulation.png
```

## Expected results

![](../figures/case01_f2_accumulation.png)

After one hour $C_{F_2} = 1.866$ mol/m^3, and `rel_error` is $3\times10^{-15}$ (round-off).

## Exercises

1. Compute $S$ by hand from the G value in the input and $D = 10^6$ W/m^3. Is it consistent with
   1.866 mol/m^3 after 3600 s?
2. Double the dose rate (`open_pronghorn-opt -i f2_accumulation.i dose_rate=2e6`). What do you
   expect, and what do you get?
3. Restore the radical yields (delete the `g_value_species`/`g_value_overrides` lines). Why is F2
   production no longer linear? (Keep this question for case 02.)

<details><summary>Solutions</summary>

1. $S = 1.866/3600 = 5.18\times10^{-4}$ mol/(m^3 s). Inverting, $G = S \cdot 100\,\mathrm{eV}
   \cdot N_A / D = 5.18\times10^{-4} \times 1.602\times10^{-17} \times 6.022\times10^{23} / 10^6
   = 0.005$ molecules/100 eV, the direct F2 yield of the database.
2. $C_{F_2} = 3.731$ mol/m^3, exactly twice the base case: the source is linear in the dose.
3. F atoms and electrons now recombine (e + F -> F-) and F atoms combine (2F -> F2), so the F2 rate
   depends on the radical concentrations. These reach a quasi-steady state quickly, but the net
   yield is no longer the bare $G_{F_2}$.

</details>
