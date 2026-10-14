# Case 09: a month of operation -- start-up, load follow, shutdown and decay

**Module 2 (fission products) | ~20 min | 1 cell, serial, 27 s**

## Objective

Combine cases 05–07 in one realistic power history, and read the different "memories" of the
inventory.

## Physics

The relative power $p(t)$ (the flux follows the power):

| days | operation |
|---|---|
| 0–1 | start-up ramp, 0 -> 100 % |
| 1–10 | full power |
| 10–13 | load follow at 50 % |
| 13–20 | full power |
| 20 | shutdown, then 10 days of decay |

- **Xe-135** follows the power with a lag of hours and swings the reactivity on every power change.
- **I-131** (8 d) remembers the last few weeks.
- **Cs-137** (30 y) integrates the energy: $N_{Cs} = Y E/(E_f N_A)$ to within its decay. It is
  the classic burnup monitor.

## The input, step by step

[power_history.i](power_history.i):

- `power_density` and `flux` are `ParsedFunction`s of the same `PiecewiseLinear` power shape.
- Xe stripping at $r = 10^{-4}$ 1/s, from case 07.
- `energy` is a `TimeIntegratedPostprocessor` of the power. `Cs137_burnup_check` compares
  $N_{Cs}$ with $Y E/(E_f N_A)$.
- `scheme = bdf2` with a `FunctionDT` time stepper: 30-min steps, cut to 1 min around every fast
  power change.

## Run

```bash
open_pronghorn-opt -i power_history.i
python plot.py
```

## Expected results

![](../figures/case09_power_history.png)

- Equilibrium xenon worth is **-1066 pcm**. At 50% load it relaxes to -666 pcm.
- After the shutdown the stripped salt peaks at only **-1325 pcm** (2.5 h after the shutdown).
- After 10 days of decay I-131 is at 42% of its maximum. Cs-137 is 0.9987 of the integrated
  energy (the deficit is decay during the month).
- Burnup: 360 MWd/m^3.

## Exercises

1. Turn the stripping off (`FissionProductInventory/removal_rates=1e-12`). What is the
   post-shutdown peak?
2. The operator wants to restart 6 h after the shutdown. With and without stripping, how much
   excess reactivity is needed?
3. Why is Cs-137 still a good burnup monitor even though the salt circulates and is sparged?

<details><summary>Solutions</summary>

1. **-4030 pcm, 8.5 h** after the shutdown (compared with -1325 pcm at 2.5 h with stripping).
2. Read `rho_Xe_pcm` at day 20.25 from the two runs: **1100 pcm** with stripping (already past
   its peak) and **3930 pcm** without. Stripping turns a restart problem into a routine one.
3. Caesium is salt-seeking (it forms a stable fluoride or chloride), is not removed by the
   sparging, and decays negligibly over the operating time.

</details>
