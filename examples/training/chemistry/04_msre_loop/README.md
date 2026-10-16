# Case 04: the MSRE fuel loop, a lumped multi-region model checked against the reactor record

**Module 4 (coupled, reactor-grade) | ~50 min | 0D multi-region, serial, 2 runs of ~8 s**

## Objective

Model the fuel loop of a real reactor, the Molten Salt Reactor Experiment (MSRE, ORNL, 1965-1969),
and compare the result with what was measured there: the xenon poisoning, where the noble metals
went, and how much chromium the salt picked up. Learn to tell a code-to-model check, a calibration
and a validation apart, on data that the models were not built from.

## What you learn

- How to combine actions in **one problem**: `[FissionProductInventory]` (linear systems) and
  `[CorrosionPlating]` (a nonlinear system, `nl_sys_names`) side by side.
- How to represent **regions** (core, pump-bowl spray, graphite, Hastelloy N, off-gas) in a 0D
  model through their removal rates, and add region inventories by hand (`FissionProductFeed`,
  `FissionProductLoss`).
- How to derive every rate from the plant data in the input (`fparse`), and run a calibrated
  variant from the command line (`r_bubble=...`, `Outputs/file_base=...`).
- How a model that reproduces the design calculations can still fail against the reactor.

## What is new compared with case 03

| cases 01-03 | case 04 |
|---|---|
| one action per input | two actions, a nonlinear and 16 linear systems in one problem |
| one salt volume, generic or loop-experiment data | a **multi-region** reactor model: spray, graphite, metal walls, off-gas, each from MSRE design data |
| removal to `offgas` only | several pathways per group and destination (`offgas`, `plated`, `processed`) and **hand-built graphite inventories** |
| verification and calibration | **VALIDATION** against independent reactor measurements (ORNL-TM-3464, ORNL-4865), plus a code-to-model check against ORNL's own model |
| one run per input | a second, **calibrated run** of the same input written to its own file (`Outputs/file_base`) |

The radiolysis action of case 01 is not used here: at MSRE's measured redox state its net effect is
zero (case 01, exercise 1), and the corrosion correlation is selected for the MSRE redox class
(`msre_or_fuel_baseline`). All three actions run together, in 2D, in case 05.

## Physics

The derivations are in [THEORY.md, section 3](../THEORY.md#3-fission-product-inventory) (xenon,
noble metals; one rate per (group, destination) in
[3.4](../THEORY.md#34-removal-pathways-one-rate-per-group-destination)) and
[section 4](../THEORY.md#4-corrosion-and-plating) (corrosion).

**Why 0D.** The salt goes round the loop in 26 s, far faster than any chemistry time here (Xe-135
9.1 h, Mo-99 66 h, corrosion years), so it is well mixed; ORNL's own xenon model (Engel and Steffy)
was lumped too. The regions enter through their **rates**; all inventories are per m^3 of salt.

| Region | Role | MSRE data (source) |
|---|---|---|
| core | fission, flux | 0.65 of 2.01 m^3 of salt, 88 % of the fissions, loop-averaged flux 2.9e12 n/cm^2/s (ORNL-4865) |
| pump-bowl spray ring | Xe stripping | 50 of the 1200 gpm, loop transit 26 s (ORNL-TM-3464) |
| graphite | Xe and noble-metal sink, Xe burnout | 225 m^2 (ORNL-4865); salt-to-graphite Xe mass transfer 0.01-0.18 ft/h (TM-3464) |
| Hastelloy N | noble-metal plating, Cr corrosion | 79.2 m^2 (pump 30 + piping 45 + heat exchanger 346 + vessel 431 ft^2, ORNL-4865) |

Power 7.4 MW (isotopic; 8 MW nominal): loop-averaged $P$ = 3.68 MW/m^3.

**1. Xenon.** The salt Xe-135 is lost by decay, burnout, and three removal paths:

$$r_s = \frac{50}{1200}\,\frac{\varepsilon_s}{26\ \mathrm{s}} = 1.60\times10^{-4}\ \mathrm{s^{-1}}\ (\varepsilon_s = 0.10),
\qquad r_g = \frac{k_g A_g}{V} = 5.69\times10^{-4}\ \mathrm{s^{-1}}\ (k_g = 0.06\ \mathrm{ft/h}),
\qquad r_b\ \text{(bubbles)}.$$

The graphite xenon ($r_g N$) decays and burns out in the core flux $\phi_c = \phi\,(V/V_c)\,0.88$ =
7.9e16 n/m^2/s. Poisoning, without importance weighting:

$$\rho = \frac{\sigma_X\,(\phi\,N_{salt} + \phi_c\,N_{graphite})}{\nu \dot F},\qquad \nu = 2.43.$$

Circulating helium bubbles were not in ORNL's design model, so the a priori run has $r_b$ = 0.

**2. Noble metals** (Mo-99, Ru-103, Te-132) deposit at the mass-transfer limit with perfect
sticking: on the metal $r_m = k_m A_{metal}/V$ = 3.9e-4 1/s ($k_m$ = 1e-5 m/s, turbulent pipe
flow), on the graphite $r_{mg}$ = 2.2e-4 1/s ($k_{mg}$ = 2e-6 m/s, laminar core channels), and to
the off-gas 5e-7 1/s. Salt fraction at steady state: $\lambda/(\lambda + r_m + r_{mg} + r_o)$.

**3. Corrosion.** `[CorrosionPlating]` salt_only on a salt column of depth $V/A_{metal}$ = 2.54 cm,
with the MSRE features of targets M-003/M-004 (hastelloy_n, fluoride_fuel, msre_or_fuel_baseline,
core_or_reactor, flow factor 2, 650 C): 3.00 µm/y. From 40 ppm Cr (the U-233 start) the salt Cr
follows from Faraday's law; McCoy's "uniform Cr removal" depth is the recession over the 7 wt% Cr.

## The input, step by step

[msre_loop.i](msre_loop.i) is one input in three parts on a single-cell mesh.

**1. MSRE data and rates**, all derived in the file, so that each pathway can be changed:

```
eps_spray = 0.10
r_spray = ${fparse 50 / 1200 * eps_spray / 26}       # 1.60e-4 1/s
k_g = ${fparse 0.06 * 0.3048 / 3600}                 # 0.06 ft/h in m/s
r_graphite = ${fparse k_g * A_graphite / V}          # 5.69e-4 1/s
r_bubble = 0                                         # calibrated run: 2.1414e-3
r_xe_offgas = ${fparse r_spray + r_bubble}
r_metal = ${fparse k_metal * A_metal / V}            # 3.9e-4 1/s
```

**2. The fission-product action** with five pathways:

```
[FissionProductInventory]
  nuclides = 'I135 Xe135 Mo99 Ru103 Te132'
  yield_set = U235_thermal
  power_density = ${P}                    # 7.4 MW / 2.01 m^3
  neutron_flux = ${phi}                   # loop-averaged 2.9e16 n/m^2/s
  removal_groups = 'noble_gas noble_gas noble_metal noble_metal noble_metal'
  removal_rates = '${r_xe_offgas} ${r_graphite} ${r_metal} ${r_mgraph} ${r_nm_offgas}'
  removal_destinations = 'offgas processed plated processed offgas'
[]
```

The action takes **one rate per (group, destination)** and rejects a second `noble_gas offgas`
entry, so spray and bubbles are summed in `r_xe_offgas`. The graphite pathways go to `processed`
(removed from the salt, not tracked by the action).

**3. Hand-built graphite inventories**, as case 01 added kernels to action variables: four
`<X>_graphite` linear variables, each with a time derivative, a feed from the salt and a loss:

```
[Xe135_graphite_uptake]
  type = FissionProductFeed
  variable = Xe135_graphite
  parent = Xe135                          # the action's salt inventory
  transfer_rate = ${r_graphite}
[]
[Xe135_graphite_loss]                     # decay + burnout in the core flux
  type = FissionProductLoss
  variable = Xe135_graphite
  decay_constant = ${fparse log(2) / 32904}
  capture_cross_section = 2647600
  neutron_flux = ${phi_core}
[]
```

**4. Corrosion** in the same problem (`[Problem] nl_sys_names = 'nl0'` next to the 16 linear
systems). `salt_diffusivity = k_metal L / 2` makes the wall-to-cell resistance $1/k_m$, and
`reference_concentration = c0` keeps the seeded rate at the 40 ppm start.

**5. Postprocessors**: the poisoning and its analytic steady state (`xe_poison_ratio`), the
fractions `<X>_<region>_fraction`, `Mo_salt_ratio`, `Ru_salt_fraction_equilibrium`, the salt Cr
rise, `rate_vs_correlation` and `faraday_ratio`. 42 days at full power, adaptive steps from 60 s to
4 h (Newton, since the corrosion system is nonlinear).

## Run

```bash
cd examples/training/chemistry/04_msre_loop
$OPEN_PRONGHORN_EXE -i msre_loop.i                                    # a priori, serial, ~8 s
$OPEN_PRONGHORN_EXE -i msre_loop.i r_bubble=2.1414e-3 Outputs/file_base=msre_loop_bubbles
python3 plot.py          # writes ../figures/case04_msre_xenon.png and case04_msre_noble_metals.png
```

Tests: `./run_tests --test-root examples/training --re=msre_loop` (two CSVDiff tests).

## Expected results

![](../figures/case04_msre_xenon.png)

| quantity (gold, 42 d) | a priori | with bubbles ($r_b$ = 2.14e-3 1/s) | MSRE |
|---|---|---|---|
| Xe-135 poisoning [% δk/k] | 1.033 | 0.270 | 0.26-0.28 measured; 1.08 ORNL model |
| share of the worth in the graphite | 97 % | 97 % | |
| time to 50 % / 90 % of equilibrium | 14 h / 38 h | 14 h / 38 h | |
| `xe_poison_ratio` | 1.0000 | 1.0000 | |

![](../figures/case04_msre_noble_metals.png)

| fraction after 42 d | Mo-99 | Ru-103 | Te-132 | ORNL-4865 Table 12.3 |
|---|---|---|---|---|
| salt | 0.47 % | 0.063 % | 0.40 % | 1-20 % |
| Hastelloy N | 63.4 % | 63.7 % | 63.5 % | 5-30 % (Te 20-90 %) |
| graphite | 36.0 % | 36.2 % | 36.1 % | 5-30 % |
| off-gas | 0.08 % | 0.08 % | 0.08 % | tenths to hundredths of 1 % |

| corrosion | model | MSRE (McCoy, TM-3063) |
|---|---|---|
| initial rate (`rate_vs_correlation` = 1.000) | 3.00 µm/y | |
| salt Cr rise after 3 weeks / 6 weeks | +26 ppm / +50 ppm | +60 ppm "over a few weeks" |
| Cr-only removal depth after 3 / 6 weeks | 2.4 / 4.7 µm | 7.6 µm (0.3 mil) |
| rate at 42 d | 2.67 µm/y (89 %) | |

| check | kind | result |
|---|---|---|
| `xe_poison_ratio`, `Mo_salt_ratio`, `faraday_ratio`, `rate_vs_correlation` | VERIFICATION | 1.0000 |
| a priori Xe 1.033 % vs ORNL's no-bubble model 1.08 % | CODE-TO-CODE (code-to-model) | 4 % apart |
| a priori Xe vs measured 0.26-0.28 % | VALIDATION | **fails, 4x too high** |
| $r_b$ fitted to 0.27 % | CALIBRATION (one parameter) | 0.270 % |
| noble-metal fractions vs ORNL-4865 | VALIDATION | **fails**: salt 2-15x too low, metal 2x too high |
| salt Cr rise vs McCoy | CALIBRATION (M-003/M-004 are correlation fit targets) | right order |

- **Xenon.** The a priori model reproduces ORNL's design calculation (our 225 m^2 of graphite, one
  flux and no importance weighting explain the 4 %), but it is 4x the measurement. This was the
  historical MSRE finding: 0.1-1 vol% of circulating helium bubbles strip xenon before it reaches
  the graphite. Is the fitted $r_b$ plausible? Xe partitions ~440:1 into helium (Henry constant
  3e-8 mol/cm^3/atm, TM-3464), so 0.15 vol% voids hold 0.66x the dissolved xenon; stripping half
  the bubbles every 26 s would remove up to ~1.3e-2 1/s, and the fit is ~17 % of that. Plausible,
  but one number cannot separate the void fraction from the stripping efficiency (Engel and Steffy:
  0.15 vol% at 50-100 % or 0.5-1 vol% at 10-15 %). The build-up times (14 h, 38 h) are a
  prediction; TM-3464 gives no clean measured value.
- **Noble metals.** Deposition at the mass-transfer limit with perfect sticking does not describe
  the MSRE, which is also what ORNL concluded: they pointed to incomplete sticking, colloidal metal
  carried in the salt and gas, and sampling artefacts. Ru-103 (39 d) is still building up at 42 d;
  its equilibrium salt fraction is 0.033 % (`Ru_salt_fraction_equilibrium`).
- **Corrosion.** Same magnitude as McCoy, but not independent: the targets table assigns McCoy's
  7.6 µm to the whole 22,533 h exposure, while McCoy relates it to the first weeks after the U-233
  start (0.4 mil for the whole operation).

## Limitations

- 0D: no spatial flux or temperature, no importance weighting of the xenon worth.
- The bubble stripping is one fitted rate, not a bubble model (void fraction, size, interfacial
  area).
- Noble metals: perfect, irreversible sticking; no colloids, no bubble attachment.
- Corrosion: the correlation rate with the Cr build-up; the redox state is not coupled to the
  fission oxidation of case 01, and there is no hot/cold transfer (case 03).

## Exercises

1. *(easy)* Estimate by hand the a priori Xe poisoning without the graphite:
   $N = Y\dot F/(\lambda + \sigma\phi + r_s)$, $\rho = \sigma\phi N/(\nu\dot F)$.
2. *(easy)* Take the graphite out of the xenon problem (`k_g=0`). What is the poisoning, and why
   does the graphite matter so much although it *removes* xenon from the salt?
3. *(medium)* Half power with the bubbles: `P_th=3.7e6 phi=1.45e16 r_bubble=2.1414e-3`. How does
   the poisoning change, and why is it not halved?
4. *(hard)* Assume the noble-metal deposition is 10x slower than the mass-transfer limit (sticking
   probability ~0.1): `k_metal=1e-6 k_mg=2e-7`. Does the salt fraction reach the MSRE range? What
   else would you need to match the metal and graphite shares?

<details><summary>Solutions</summary>

1. $\dot F = P/E_f$ = 3.68e6 / 3.204e-11 = 1.15e17 fissions/(m^3 s); with $Y$ = 0.0654,
   $\lambda + \sigma\phi + r_s$ = 2.11e-5 + 7.7e-6 + 1.60e-4 = 1.89e-4 1/s, $N$ = 3.98e19 atoms/m^3,
   and $\rho$ = 2.65e-22 x 2.9e16 x 3.98e19 / (2.43 x 1.15e17) = **0.109 %**.
2. **0.109 %** (gold run with `k_g=0`: 0.1093 %), the hand estimate of exercise 1, and 10x lower
   than with graphite. Without graphite the spray removes the salt xenon at 5.5x its decay
   constant. With graphite, the salt xenon is taken up 3.5x faster than the spray removes it, and
   in the graphite it can leave only by decay or burnout, in the core at 2.7x the loop-averaged
   flux. The graphite is the main xenon poison of a graphite-moderated MSR, which is why ORNL wanted
   sealed graphite for the MSBR.
3. **0.179 %**: 66 % of the full-power value, not 50 %. At lower flux the burnout is smaller
   ($\sigma\phi_c$ = 2.1e-5 1/s at full power, about $\lambda_X$), so each fission leaves more
   xenon behind. The rise to 50 % of equilibrium takes 18 h instead of 14 h. (A priori, without
   bubbles: 0.687 %.) Engel and Steffy did not report a clean half-power value: a prediction.
4. Salt 4.5 % (Mo-99), 0.63 % (Ru-103), 3.9 % (Te-132): now inside the 1-20 % range for Mo and Te.
   The metal and graphite shares hardly change (60 % and 34-36 %), because the removal is still
   ~100x faster than the decay; the split between them is the ratio $r_m/r_{mg}$, which a common
   sticking factor does not change. Matching ORNL-4865 also needs a different metal/graphite split
   (lower $k_m$ relative to $k_{mg}$, or a colloid that is carried to the off-gas and the
   graphite), which is what ORNL concluded.

</details>

## Further reading

- [THEORY.md, section 3](../THEORY.md#3-fission-product-inventory) and
  [section 4](../THEORY.md#4-corrosion-and-plating).
- J. R. Engel, R. C. Steffy, "Xenon behavior in the Molten Salt Reactor Experiment", ORNL-TM-3464
  (1971).
- E. L. Compere et al., "Fission product behavior in the Molten Salt Reactor Experiment",
  ORNL-4865 (1975).
- H. E. McCoy, ORNL-TM-3063 (1970): Hastelloy N corrosion and the salt Cr in the MSRE.
- † P. N. Haubenreich, J. R. Engel, "Experience with the Molten-Salt Reactor Experiment",
  Nucl. Appl. Technol. 8 (1970) 118.

*References marked † were added from general literature knowledge and have not been checked against
the repository's sources; verify the bibliographic details before citing them.*
