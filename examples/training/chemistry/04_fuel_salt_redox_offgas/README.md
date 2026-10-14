# Case 04: radiolytic F2, the U(IV)/U(III) redox buffer and the cover gas

**Module 1 (radiolysis) | ~20 min | 1 cell, serial, ~5 min (heavy test)**

## Objective

Model the redox control of a fluoride fuel salt and the transfer of a radiolytic gas to the
headspace with Henry's law.

## Physics

In a fluoride fuel salt (FLiBe-UF4, 600 C) the uranium couple controls the chemistry. MSRE kept
U(III)/U at ~1% to protect the Hastelloy N:

```
e_sol + U(IV)    -> U(III)            (reduces the salt)
F2    + 2 U(III) -> 2 U(IV) + 2 F-    (U(III) scavenges radiolytic fluorine)
2 F -> F2,  e_sol + F -> F-,  plus a direct F2 yield G = 0.005
```

Dissolved F2 exchanges with the cover gas:

$$\frac{dC_l}{dt} = -k_La\,(C_l - k_H R T\,C_g), \qquad
  \frac{dC_g}{dt} = \frac{V_l}{V_g}\,k_La\,(C_l - k_H R T\,C_g)$$

The direct F2 yield has no matching electron. It slowly oxidizes the salt at
$2G_{F_2}D/(100\,\mathrm{eV}\,N_A) \approx 10^{-3}$ mol/(m^3 s), so the 1% buffer is spent in
about 1.6 h. Only then does F2 build up.

## The input, step by step

[fuel_salt_redox.i](fuel_salt_redox.i):

- `metals = 'U'` adds U(III)/U(IV). `gas_species = 'F2'` with `kLa` and `volume_ratio` adds the
  well-mixed headspace variable `F2_gas` and the Henry-law exchange.
- `initial_condition_species = 'F_ion U_IV U_III'` starts at 1% U(III).
- `dtmax = 10` s. The comment in the input records the time-step study: the redox error is
  $2\times10^{-4}$ at 10 s and $1.4\times10^{-2}$ at 60 s.
- Checks: `uranium_error`, `redox_error` and `henry_ratio`, the distance from gas-liquid
  equilibrium.

## Run

```bash
open_pronghorn-opt -i fuel_salt_redox.i          # 6 h of operation, ~5 min
python plot.py
```

## Expected results

![](../figures/case04_fuel_salt_redox.png)

- U(III)/U falls from 1% to below 0.01% in 1.6 h. From then on F2 accumulates in the salt and
  the cover gas.
- At 6 h, `henry_ratio` is 1.44: the headspace is still filling and has not reached equilibrium.
- `uranium_error` stays below $3\times10^{-6}$ and the maximum `redox_error` is $2\times10^{-4}$.

## Exercises

1. Remove the buffer (`U_III_0=0`, and `Executioner/end_time=3600` to keep it short). What happens
   to F2?
2. MSRE replenished U(III) by adding Be metal. Estimate the Be addition rate [mol/(m^3 h)] that
   would hold the buffer at 1% at this dose rate.
3. Why is `henry_ratio` above 1 while the gas is still filling?

<details><summary>Solutions</summary>

1. With nothing to scavenge it, F2 appears immediately. After 1 h the salt holds
   $7.9\times10^{-3}$ mol/m^3 and the gas 0.37 mol/m^3, with `henry_ratio` = 2.9.
2. Net oxidation is ~$10^{-3}$ mol/(m^3 s) of equivalents. Be -> Be2+ supplies 2 equivalents per
   atom, so ~$5\times10^{-4}$ mol/(m^3 s) ≈ 1.8 mol/(m^3 h) of Be.
3. The ratio compares $C_l$ with $k_H R T\,C_g$. While the gas is filling, $C_l$ is above its
   equilibrium with the gas, and the excess is what drives the flux into the headspace.

</details>
