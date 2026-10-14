# Case 12: which alloy, which salt, which redox state? -- comparison with loop data

**Module 3 (corrosion) | ~15 min | 7 runs, serial, ~5 s total**

## Objective

Use command-line overrides to drive a parameter study from Python, and compare the calibrated rate
with seven measured loop and capsule corrosion rates.

## Physics

The initial Butler-Volmer rate of a coupon (case 10) is the calibrated correlation evaluated for the
case features: alloy, salt, redox state, temperature, loop position, flow and thermal gradient. The
cases come from `validation/corrosion` (`data/validation_predictions.csv`):

| case | alloy | salt | redox | T [C] | measured [µm/y] |
|---|---|---|---|---|---|
| ORNL-FL-04 | Hastelloy N | fuel salt | purified | 704 | 1.0 |
| FB-03 | Hastelloy N | fluoroborate | purified | 687 | 7.4 |
| ORNL-FL-01 | 304L | fuel salt | MSRE fuel | 688 | 21.8 |
| ORNL-FL-12 | 316 | fuel salt | MSRE fuel | 650 | 25.4 |
| MOD-F-07 | GH3535 | FLiNaK | purified | 700 | 23.3 (5.9–91.6) |
| ORNL-FL-15 | 316 | FLiBe | moisture-contaminated | 650 | 213 |
| MOD-F-07s | Hastelloy N | FLiNaK | stressed | 700 | 217 |

## The input, step by step

[coupon_rate.i](coupon_rate.i) declares the features as **top-level HIT variables**:

```
T = 977.15
material = hastelloy_n
salt = fluoride_fuel
...
[CorrosionPlating]
  material_class = ${material}
  ...
```

so any of them can be overridden from the command line:

```bash
open_pronghorn-opt -i coupon_rate.i material=stainless_316 salt=flibe redox=impure_moisture T=923.15
```

[run_scan.py](run_scan.py) loops over the table, writes `scan_results.csv` and draws the parity
plot.

## Run

```bash
python run_scan.py        # or: python plot.py
```

## Expected results

![](../figures/case12_parity.png)

| case | model [µm/y] | model / measured |
|---|---|---|
| ORNL-FL-04 | 1.38 | 1.36 |
| FB-03 | 6.96 | 0.94 |
| ORNL-FL-01 | 16.8 | 0.77 |
| ORNL-FL-12 | 21.6 | 0.85 |
| MOD-F-07 | 20.4 | 0.87 |
| ORNL-FL-15 | 191 | 0.90 |
| MOD-F-07s | 193 | 0.89 |

All seven are within a factor of 2 (geometric-mean error 1.18). Stainless steels corrode an order of
magnitude faster than Hastelloy N in fuel salt. Impurities (moisture) or stress raise the rate by
another order of magnitude.

## Exercises

1. Run ORNL-FL-12 but with Hastelloy N. How much does the alloy choice buy?
2. Run ORNL-FL-04 with `redox=oxidizing_fef2`. What does a loss of redox control cost?
3. Add an eighth case of your own to `CASES` in `run_scan.py`.

<details><summary>Solutions</summary>

1. `open_pronghorn-opt -i coupon_rate.i T=923.15 material=hastelloy_n redox=msre_or_fuel_baseline
   dT=110` gives **2.8 µm/y**, against 21.6 µm/y for 316: the nickel alloy is 7.7x more
   resistant in the same salt.
2. **9.1 µm/y**, 6.6x the purified-salt rate of 1.38 µm/y. The `oxidizing_fef2` class represents
   a salt whose redox buffer is lost (case 04), which is exactly what MSRE operated to avoid.
3. Append a tuple `(name, alloy, salt, redox, position, T_K, flow, dT, (low, mid, high))`.

</details>
