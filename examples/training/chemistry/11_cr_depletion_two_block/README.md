# Case 11: the chromium-depletion layer in 316 stainless steel

**Module 3 (corrosion) | ~20 min | 1D two-block, serial, 3 s**

## Objective

Resolve the alloy as a second block, couple it to the salt with a Butler-Volmer interface, and see
dissolution become limited by solid-state diffusion.

## Physics

Case 10 treated the wall as an infinite Cr source. In reality Cr must diffuse through the alloy to
the surface, and the near-surface alloy becomes depleted. In molten fluorides this is the voided,
Cr-depleted zone seen in metallography of loop specimens.

$$\text{salt: } \frac{\partial c}{\partial t} = \nabla\cdot(D_{salt}\nabla c), \qquad
  \text{alloy: } \frac{\partial c_s}{\partial t} = \nabla\cdot(D_{alloy}\nabla c_s)$$

$$i = i_0\left[\frac{c_s}{c_{s,0}}e^{\alpha z f\eta} - \frac{c}{c_{ref}}e^{-\alpha z f\eta}\right]
  \qquad (\texttt{solid\_activity = true})$$

The anodic branch carries the **alloy-surface activity** $c_s/c_{s,0}$, so dissolution slows as the
surface depletes. Without it the surface concentration can go negative when $D_{alloy}$ is small.

## The input, step by step

[cr_depletion.i](cr_depletion.i):

- `[Mesh]` is a `CartesianMeshGenerator` with a 0.5 m salt reservoir (25 cells) and a 300 µm alloy
  (300 cells of 1 µm), named `salt` and `alloy`, with an `interface` sideset.
- `[CorrosionPlating] topology = two_block`, `solid_diffusivity = ${D_alloy}`,
  `solid_activity = true`, and the alloy initial Cr is $c_{s,0}$ = 0.17 · 8 g/cm^3 / M_Cr.
- $D_{alloy} = 10^{-16}$ m^2/s is an effective value lumping lattice ($\sim10^{-19}$ at 650 C) and
  grain-boundary transport.
- The `alloy_profile` VectorPostprocessor samples the depletion profile at the end.

## Run

```bash
open_pronghorn-opt -i cr_depletion.i
python plot.py
```

## Expected results

![](../figures/case11_cr_depletion.png)

- The rate drops from **23.4 to 3.2 µm/y** within the year. The surface Cr falls to **48%** of the
  bulk, and the zone below 95% of the bulk reaches ~125 µm ($\sim\sqrt{4Dt}$ = 112 µm).
- `cr_balance_error` = $10^{-13}$, `faraday_ratio` = 1, and `depth_over_sqrt_Dt` = 1.06
  ($2/\sqrt{\pi}$ = 1.13 for a fixed surface concentration).

## Exercises

1. A small capsule (`L_salt=0.05`). What limits the corrosion now?
2. Slower alloy diffusion (`D_alloy=1e-17`). Compare the surface Cr and the dissolved Cr.
3. Run with `CorrosionPlating/solid_activity=false` and `D_alloy=1e-18`. Look at the minimum of
   `cs_Cr` and explain.

<details><summary>Solutions</summary>

1. The salt reaches the Nernst limit of case 10 within weeks (143.7 wppm), dissolution stops, and
   Cr diffusing from the bulk **refills** the surface: it falls to 64% and then recovers to 77%.
   The rate at one year is only 0.6 µm/y.
2. The surface depletes deeper (to 23%), and only 12.7 wppm dissolves, compared with 25.6 in the
   base case. Dissolution is limited by the supply of Cr from the alloy.
3. Without the surface activity the anodic branch keeps pulling Cr at the fresh-alloy rate. The
   surface cell reaches **-9.2x the bulk** (unphysical) because the diffusion cannot resupply it,
   while `cr_balance_error` is still $10^{-13}$: a conservative scheme can still be wrong. This is
   why the option exists.

</details>
