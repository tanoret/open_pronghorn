# Case 07: removing the Xe-135 poison by sparging -- an MSR advantage

**Module 2 (fission products) | ~15 min | 5 cells, serial, 8 s**

## Objective

Add a removal pathway to the inventory (off-gas stripping), follow the removed inventory as it
keeps decaying, and scan the removal rate in one run.

## Physics

Xenon is insoluble in the salt. An MSR can strip it continuously by sparging helium bubbles and
venting the cover gas. At MSRE this was a spray ring in the pump bowl. Stripping adds a first-order
removal $r$ to case 06:

$$N_{Xe} = \frac{Y_{Xe}F/N_A + \lambda_I N_I}{\lambda_{Xe} + \sigma_{Xe}\phi + r}$$

The stripped xenon goes to an off-gas inventory, `Xe135_offgas`, where it decays to Cs-135. The
off-gas system therefore holds a decay-heat and Cs-135 source term.

## The input, step by step

[xenon_stripping.i](xenon_stripping.i):

```
[FissionProductInventory]
  nuclides = 'I135 Xe135 Cs135'
  removal_groups = 'noble_gas'            # which chemical group is removed
  removal_rates = 'stripping_rate'        # a functor [1/s]: here it varies from cell to cell
  removal_destinations = 'offgas'         # creates Xe135_offgas (and Cs135_offgas from its decay)
[]
```

The 5 cells have $r$ = 0, 1e-5, 1e-4, 1e-3 and 1e-2 1/s, so one run is the scan.

## Run

```bash
open_pronghorn-opt -i xenon_stripping.i
python plot.py
```

## Expected results

![](../figures/case07_xenon_stripping.png)

| r [1/s] | 0 | 1e-5 | 1e-4 | 1e-3 | 1e-2 |
|---|---|---|---|---|---|
| neutrons absorbed by Xe-135 | 2.13 % | 1.93 % | 1.07 % | 0.19 % | 0.02 % |

Stripping only pays off once $r$ is comparable with $\lambda_{Xe} + \sigma_{Xe}\phi = 10^{-4}$
1/s. MSRE reached a Xe-135 poisoning of a few tenths of a percent (ORNL-4865), the $r \sim 10^{-3}$
column.

## Exercises

1. At $\phi_0 = 10^{18}$ (`phi0=1e18`), what stripping rate halves the poisoning?
2. Where does the stripped xenon's decay heat end up, and why must the off-gas system be designed
   for it?
3. Combine with case 06: would a stripped MSR still have a xenon pit after shutdown?

<details><summary>Solutions</summary>

1. Burnout is now $\sigma\phi = 2.6\times10^{-4}$ 1/s, so the denominator without stripping is
   $3.3\times10^{-4}$ 1/s. Halving needs $r \approx 3.3\times10^{-4}$ 1/s. The scan gives 2.49%
   at r = 0, 1.85% at 1e-4 and 0.55% at 1e-3.
2. In the off-gas holdup (charcoal beds at MSRE): `Xe135_offgas` and its Cs-135 daughter decay
   there, outside the core.
3. Much less of one. The stripping continues after shutdown and removes the xenon fed by I-135
   decay (case 09 shows -1330 pcm at the peak instead of -4030 pcm).

</details>
