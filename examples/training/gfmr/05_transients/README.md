# Stage 5: reactivity-insertion transient

`gfmr_rod_step.i` is an unprotected **control-rod step of 0.3 $** (204 pcm, about 4 % of the
total rod worth) at t = 1 s, from the nominal steady state, without scram. It is a lumped model
of the whole core: point kinetics with circulating-fuel precursors and two-temperature
feedback, solved as MOOSE scalar ODEs (BDF2, one rank, a few seconds).

Only this one scenario is done. Loss of flow, overcooling and loss of heat sink need a loop
thermal model, which the specification does not define (see the end of this README).

## Model

| Part | Equation | Data |
|---|---|---|
| Power | $\dot n = \frac{\rho - \beta + \rho_0}{\Lambda} n + \sum_i \lambda_i c_i$ | β_eff, Λ from OpenMC IFP |
| Core precursors | $\dot c_i = \frac{\beta_i}{\Lambda} n - \lambda_i c_i - \frac{c_i - d_i}{t_\text{core}}$ | 6 groups, ENDF/B-VIII.0 |
| Loop precursors | $\dot d_i = \frac{c_i - d_i}{t_\text{ex}} - \lambda_i d_i$ | t_core 3.02 s (CFD), t_ex 6.04 s ASSUMED |
| Salt | $M_s c_s \dot T_s = (1-f_g) P - 2\dot m c_s (T_s - T_\text{in}) + UA (T_g - T_s)$ | stage 3 |
| Graphite | $M_g c_g \dot T_g = f_g P - UA (T_g - T_s)$ | stage 3 |
| Feedback | $\rho = \rho_\text{step} + \alpha_f (T_s - T_{s0}) + \alpha_g (T_g - T_{g0})$ | α_f −2.7, α_g −4.1 pcm/K (OpenMC) |

The model's constants:

| Constant | Value | Origin |
|---|---|---|
| Graphite heat fraction f_g | 7.87 % | OpenMC coupled n-γ heating |
| Graphite-to-salt conductance UA | 1.66×10⁵ W/K | 7.87 % of 150 MW across the 71 K graphite–salt difference of the CFD |
| Graphite heat capacity M_g c_g | 7.28×10⁷ J/K | lattice graphite mass × 1900 J/kg/K (assumed c_g) |
| Graphite time constant | 438 s | follows |

`T_s` is the core-average salt temperature, with $T_\text{out} = 2T_s - T_\text{in}$. The
inlet temperature is held at 613.8 °C: a perfect heat exchanger, the most conservative boundary
for a fast insertion because it gives no loop feedback.

### Precursor drift

Precursors born in the core are carried out with the salt, and part of them decay in the loop.
Their neutrons are lost. At steady state this costs the **drift reactivity**:

$$\rho_0 = \beta - \Lambda \sum_i \lambda_i c_{i0} = 295\ \text{pcm} = 0.43\ \$,$$

which the core must supply to stay critical. With t_ex = 6.04 s, 43 % of β_eff is lost. The
OpenMC k-eff of a static core does not include this loss: it makes the 490 pcm subcriticality of
the specified core (stage 2) worse, not better (exercise 1).

## Verification

| Check | Result |
|---|---|
| Null transient (ρ_step = 0, 300 s) | n, T_s, T_g constant to 8 digits |
| Final state vs analytic feedback balance | n = 1.3789 at 600 s vs 1.3797 analytic (still settling); ΔT_s = 13.8 K, ΔT_g = 40.8 K |

The analytic final state solves $\rho_\text{step} + \alpha_f \Delta T_s + \alpha_g \Delta T_g = 0$
with the steady $\Delta T_s$ and $\Delta T_g$ of the thermal model.

## Results

| Time | Power | T_salt avg | T_graphite | T_outlet | Net reactivity |
|---|---|---|---|---|---|
| 1.0 s (step) | 1.00 | 650.1 °C | 721.1 °C | 686.4 °C | 0.30 $ |
| 2.0 s | 2.32 | 664.0 | 721.3 | 714.2 | 0.24 $ |
| 10.6 s | 2.76 | 706.1 | 724.4 | 798.5 | 0.06 $ |
| 15.2 s (peak) | **2.78 (417 MW)** | — | — | — | — |
| 61 s | 2.34 | 694.6 | 742.2 | 775.4 | 0.00 $ |
| 301 s | 1.40 | 664.7 | 761.8 | 715.5 | 0.00 $ |
| 600 s | 1.38 | 663.9 | 761.9 | 714.0 | 0.00 $ |

**Reading the transient.**

- **Prompt jump, then a fuel-feedback turnover.** Drift leaves only β − ρ₀ = 386 pcm of delayed
  neutrons in the core, so the 204 pcm step is 0.53 $ of the *effective* β: still below prompt
  critical. The power first jumps by about (β − ρ₀)/(β − ρ₀ − ρ) = 2.1 (2.3 simulated at 2 s). The salt then heats within
  seconds (its time constant is ~3 s) and the fuel coefficient turns the power over at 2.8 ×
  nominal at 15 s.
- **Graphite takes over slowly.** The graphite heats on its 438 s time constant. Its larger
  coefficient (−4.1 pcm/K) brings the power down to the final 1.38 × nominal over ~10 minutes.
- **Outlet temperature.** The peak is 798 °C, reached at about 11 s.

The temperatures are core-average values; the hot-channel peaks are higher (stage 3: the
hottest salt cell is 33 K above the mixed outlet).

## Limits

- Point kinetics: the power shape is frozen at the nominal shape. A real rod withdrawal tilts
  the flux toward the top of the core.
- Lumped temperatures: one salt and one graphite node.
- The inlet temperature is fixed (no loop).
- t_ex and c_g are assumed.

## Not done

Loss of flow, overcooling and loss of heat sink. Each needs a loop model with a heat exchanger
and a pump coast-down, which the GFMR specification does not define. With such a loop, the same
input extends directly:

- **Loss of flow:** make ṁ(t) a function. A lower flow slows the precursor exit (t_core, t_ex grow)
  and returns reactivity while the salt heats.
- **Overcooling:** make T_in(t) fall.
- **Loss of heat sink:** add an ex-core salt node in place of the fixed T_in.

## How to run

```bash
open_pronghorn-opt -i gfmr_rod_step.i                         # 0.3 $ step, 600 s
open_pronghorn-opt -i gfmr_rod_step.i rho_step=0              # null transient (verification)
open_pronghorn-opt -i gfmr_rod_step.i "rho_step=${fparse 0.6 * 0.00681035}"   # a 0.6 $ step
```

## Exercises

1. With the drift reactivity ρ₀ = 295 pcm, how much excess reactivity does the GFMR need at
   operating conditions? Compare with the 490 pcm subcriticality of the static OpenMC core.
2. Run a 0.6 $ step. How does the peak power scale? The circulating core is prompt critical at
   ρ = β − ρ₀ = 386 pcm, i.e. at 0.57 $ of the static β: what happens just above it?
3. Set t_ex to 1e6 s (no loss in the loop). How do ρ₀ and the peak power change?
