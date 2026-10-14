# Case 10: a stainless-steel coupon in a static fluoride-salt capsule

**Module 3 (corrosion) | ~20 min | 1D, serial, 3 s**

## Objective

Meet the `[CorrosionPlating]` action: a Butler-Volmer wall seeded from the calibrated corrosion
correlation. Then see why corrosion in a closed capsule is self-limiting.

## Physics

The wall (left boundary) is a Butler-Volmer electrode. Chromium, the least noble component, dissolves
as Cr2+:

$$i = i_0\left[e^{\alpha z f\eta} - \frac{c_{Cr}}{c_{ref}}e^{-\alpha z f\eta}\right], \qquad
  f = \frac{F}{RT}, \qquad J_{Cr} = \frac{i}{zF}$$

The action **seeds $i_0$** so that the initial rate reproduces the reference correlation for the
alloy/salt/redox class. Here that is case ORNL-FL-01 of `validation/corrosion`: 304L in fuel salt at
688 C, model 16.82 µm/y, measured 21.8 µm/y.

In a closed capsule the dissolved Cr builds up until the net current vanishes at the **Nernst
equilibrium**:

$$c_{eq} = c_{ref}\,e^{z f \eta} = 11.2\ \mathrm{mol/m^3}\ (\approx 180\ \mathrm{wppm}), \qquad
  d_{final} = \frac{(c_{eq} - c_0)\,L\,M_{Cr}}{\rho_{alloy}} = 0.66\ \mu\mathrm{m}$$

## The input, step by step

[corrosion_coupon_0d.i](corrosion_coupon_0d.i):

```
[CorrosionPlating]
  topology = salt_only               # salt domain + a Butler-Volmer boundary
  reaction_boundary = left
  elements = 'Cr'
  material_class = stainless_304l    # alloy / salt / redox / position classes of the correlation
  salt_class = fluoride_fuel
  redox_class = msre_or_fuel_baseline
  flow_factor = 1  delta_T_C = 100   # loop features of the reference case
  applied_overpotential = ${eta}
  verbose = true                     # prints the reference rate and the seeded i0
[]
```

- The salt depth $L$ = volume / wetted area = 1 cm, on a graded 20-cell 1D mesh.
- Checks: `rate_vs_reference` (= 1 at t = 0), `c_eq_ratio` (wall concentration over the Nernst
  value) and `faraday_ratio` (recession over the dissolved Cr inventory, = 1).

## Run

```bash
open_pronghorn-opt -i corrosion_coupon_0d.i
python plot.py
```

## Expected results

![](../figures/case10_corrosion_coupon.png)

- The initial rate is 16.82 µm/y (`rate_vs_reference` = 1.000).
- The rate decays with a time constant of about 2 weeks. After 120 days `c_eq_ratio` is 0.9997 and
  the recession has stopped at **0.662 µm**, the analytic Faraday limit. Had the initial rate been
  sustained, it would have been 5.5 µm.
- `faraday_ratio` = 1.000.

**Lesson:** a static capsule corrodes until the salt is saturated. Sustained attack needs a sink:
the cold leg of a loop (case 13) or an oxidant that keeps removing Cr(II) (case 14).

## Exercises

1. Make the salt more oxidizing (`eta=0.15`). What happens to the initial rate, $c_{eq}$ and the
   final recession? Why is the initial rate unchanged?
2. A 10x deeper capsule (`L=0.1`, `Executioner/end_time=3.1536e7`). How do the time constant and the
   final recession scale?
3. Check $c_{eq}$ and $d_{final}$ by hand for the base case.

<details><summary>Solutions</summary>

1. The initial rate is still **16.82 µm/y**: the action re-seeds $i_0$ so that the reference
   correlation is reproduced at the given $\eta$. The equilibrium moves to
   $c_{eq} = e^{2f \cdot 0.15}$ = 37.4 mol/m^3, and the final recession to 2.37 µm. After 120 days
   it is 2.13 µm, with `c_eq_ratio` = 0.90: a higher equilibrium takes longer to approach.
2. Both scale with $L$. After one year the recession is 5.7 µm (`c_eq_ratio` = 0.90), approaching
   the 6.62 µm limit; the time constant is ~10x longer (months).
3. $f = F/RT$ = 24.15 1/V at 961 K, so $c_{eq} = e^{2 \cdot 24.15 \cdot 0.1}$ = 11.19 mol/m^3.
   $d = (11.19 - 1) \cdot 0.01 \cdot 51.996 / 8.0\times10^6$ m = 0.662 µm.

</details>
