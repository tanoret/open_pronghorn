# Case 03: one loop, four corrosion models -- Hastelloy N in ORNL loop NCL-16

**Module 3 (corrosion and plating) | ~50 min | 0D and 1D, serial, 3 inputs, ~20 s in total**

## Objective

Put the corrosion models of OpenPronghorn side by side on one well-documented experiment, the ORNL
thermal-convection loop NCL-16, and learn what each one can and cannot tell you. Be strict about
what is a calibration (the value was fitted) and what is a prediction.

## What you learn

- How the third action, `[CorrosionPlating]` / `[CorrosionPlatingFlow]`, seeds a Butler-Volmer
  wall law from the calibrated corrosion correlation (`data/corrosion_database.json`).
- Why a closed isothermal capsule stops corroding (the Nernst limit) and why a loop does not
  (the temperature dependence of the equilibrium).
- How to build a **1D model with two walls** at different temperatures, declare your own
  variables, and let the action release the dissolved metals into them.
- How a reduced-order long-term model (DRIDN, a user object) reproduces its calibration data, and
  how to read such a fit critically.

## What is new compared with case 02

| case 02 | case 03 |
|---|---|
| `[FissionProductInventory]` | `[CorrosionPlating]` and `[CorrosionPlatingFlow]` (the third action), plus the `DRIDNCorrosionUserObject` |
| 0D, one cell | a **1D salt gap** with a hot and a cold wall (`reaction_boundary = 'left right'`, a temperature function) |
| the action creates its variables | **you declare** `c_Cr c_Fe c_Ni` and the action releases into them (`release_variables`); your own time-derivative and diffusion kernels |
| one model | **several models** of the same experiment, compared in one plot; a `solve = false` user-object model |
| verification only | CALIBRATION against measured loop data, with honest labels |

## Physics

The derivations are in [THEORY.md, section 4](../THEORY.md#4-corrosion-and-plating): the Nernst
limit in [4.3](../THEORY.md#43-c-the-nernst-limit-and-why-a-constant-eta-is-wrong-for-a-loop), the
flow wall law in [4.4](../THEORY.md#44-d-the-element-resolved-temperature-dependent-equilibrium-law-corrosionplatingflow),
MSTDB-TC in [4.5](../THEORY.md#45-e-the-mstdb-tc-standard-state-model), DRIDN in
[4.6](../THEORY.md#46-f-dridn-the-dynamic-redox-inventory-depletion-network) and the comparison in
[4.8](../THEORY.md#48-comparison-of-the-corrosion-models).

**The experiment.** NCL-16 (ORNL-TM-4188; validation case ORNL-FL-04): Hastelloy N in purified
LiF-BeF2-ThF4-UF4 fuel salt, 704 C hot leg, 538 C cold leg, 29500 h (3.37 y) before the FeF2
additions. Measured (`validation/corrosion/data/targets.csv`):

| id | quantity | measured |
|---|---|---|
| M-011 | hot-leg weight loss | 2.9 mg/cm^2 |
| M-012 | cold-leg weight gain | 1.7 mg/cm^2 |
| M-013 | uniform corrosion rate | 1.016 µm/y (0.04 mil/y) |
| M-014 | salt Cr increase | +500 ppm |
| M-015 | salt Fe change | -100 ppm |

**The wall law.** Both actions use a Butler-Volmer law whose exchange current is seeded from the
calibrated rate of the NCL-16 features (Hastelloy N, fluoride fuel, purified salt, ΔT = 166 C):
1.382 µm/y on fresh salt at 704 C. The element-resolved form drives each element $k$ toward its
own temperature-dependent equilibrium:

$$c_{eq,k}(T) = c_{eq,k}\,e^{-\Delta H_k/R\,(1/T - 1/T_{ref})},\qquad
J_k = J_{k,ref}\,e^{E_a/R\,(1/T_{ref}-1/T)}\sqrt{\frac{c_{eq,k}(T)}{c_{eq,k}}}\,\Big(1 - \frac{c_k}{c_{eq,k}(T)}\Big).$$

With $\Delta H > 0$, a hot wall has the higher equilibrium: it dissolves, the cold wall plates,
and the salt settles between the two. $J_{k,ref}$ is the `seed_fractions` share of the calibrated
rate (Cr/Fe/Ni = 0.90/0.07/0.03: selective Cr leaching).

**The four models.**

1. **Butler-Volmer capsule** ([ncl16_capsule_bv.i](ncl16_capsule_bv.i)): the hot-leg wall alone,
   in a closed isothermal salt volume of depth V/A = 3.9e-4 m. Its Nernst limit
   $c_{ref}e^{zf\eta}$ is set to 900 ppm (the hot-wall equilibrium of model 2).
2. **Thermodynamic flow wall law** ([ncl16_loop_flow.i](ncl16_loop_flow.i)): a closed salt gap
   between a hot wall (977 K) and a cold wall (811 K) of equal area, V/A_total = 3.9e-4 m.
3. **DRIDN** ([ncl16_dridn.i](ncl16_dridn.i)): the Dynamic Redox-Inventory-Depletion Network, a 0D
   19-state model of hot-leg dissolution, salt inventory, cold-leg capture, an empirical
   **bulk-capture sink** and alloy depletion, integrated over years. Its parameters were fitted
   (in-sample) to NCL-16 and other loops.
4. **MSTDB-TC standard-state model**: thermodynamic affinities from MSTDB Gibbs data, with a
   physical cold capture $1 - a_{sat}(T_c)/a_{sat}(T_h)$. **It cannot be run here**: it needs the
   licensed MSTDB-TC V3.1 files. Its NCL-16 endpoint is read from
   `validation/corrosion/calibration/results/advanced/advanced_validation_summary.csv`.

**The equilibria of model 2 are a two-parameter calibration.** $c_{eq,Cr}$ = 900 ppm at 704 C
puts the salt Cr near +500 ppm (M-014), and $\Delta H_{Cr}$ = 100 kJ/mol (72 ppm at 538 C) puts
the cold gain near 1.7 mg/cm^2 (M-012). The hot-leg loss and the rate are then **predictions**. Fe
starts at an assumed 150 ppm with $c_{eq,Fe}$ = 40 ppm, so Cr from the alloy reduces the salt
Fe(II) and Fe plates out. Ni gets 1 ppm (the noble component).

## The input, step by step

**Model 1, [ncl16_capsule_bv.i](ncl16_capsule_bv.i).** A salt column of depth V/A with the wall on
the left. The overpotential is chosen so that the Nernst limit is 900 ppm:

```
eta = ${fparse log(c_eq) * 8.314462618 * T_hot / (2 * 96485.33212)}   # 0.170 V
[CorrosionPlating]
  topology = salt_only                  # the action creates the salt Cr variable and its kernels
  reaction_boundary = left
  elements = 'Cr'
  material_class = hastelloy_n          # these four classes select the calibrated correlation
  salt_class = fluoride_fuel
  redox_class = purified_baseline
  position_class = nominal
  delta_T_C = 166.0
  applied_overpotential = ${eta}
  default_salt_concentration = 0.0      # fresh (purified) salt
[]
```

`faraday_ratio` (recession over the dissolved inventory) is a VERIFICATION of the wall flux.

**Model 2, [ncl16_loop_flow.i](ncl16_loop_flow.i).** A 1D gap of 4 cells, `gap = 7.8e-4` m, with a
linear temperature function (hot on the left, cold on the right). You declare the three dissolved
metals yourself, so that Fe can start at 150 ppm, and give them a time derivative and a mixing
diffusion ($D$ = 1e-8 m^2/s, about 1 min across the gap: the salt is well mixed, as in a loop):

```
[Variables]
  [c_Fe]
    type = MooseLinearVariableFVReal
    solver_sys = c_Fe_sys
    initial_condition = ${fparse Fe0_ppm * ppm / M_Fe}
  []
  ...                                   # c_Cr, c_Ni start at 0
[]
[LinearFVKernels]                       # Cr_time, Cr_diffusion, ... for each metal
[]
[CorrosionPlatingFlow]
  elements = 'Cr Fe Ni'
  release_variables = 'c_Cr c_Fe c_Ni'  # release into your variables
  reaction_boundary = 'left right'      # both walls
  temperature = T                       # a function: each wall has its own T
  reference_temperature = ${T_hot}
  temperature_dependent_kinetics = true
  equilibrium_concentrations = '${c_eq_Cr} ${c_eq_Fe} ${c_eq_Ni}'
  dissolution_enthalpies = '${dH_Cr} ${dH_Fe} ${dH_Ni}'
  seed_fractions = '${seed_Cr} ${seed_Fe} ${seed_Ni}'
  ...                                   # the same four correlation classes as model 1
  verbose = true                        # prints i0, the planned rate and c_eq per element
[]
```

`J_*_hot` and `J_*_cold` evaluate the wall law from the wall-cell concentrations, and
`TimeIntegratedPostprocessor`s turn them into the hot-wall loss and the cold-wall gain (mg/cm^2).
`inventory_balance_error` compares loss minus gain with the dissolved inventory (VERIFICATION).
`scheme = bdf2` and `dtmax` = 0.005 y keep its time-integration error below $10^{-3}$.

**Model 3, [ncl16_dridn.i](ncl16_dridn.i).** No PDE: `[Problem] solve = false`, and the
`DRIDNCorrosionUserObject` integrates its ODE network once per time step (40 steps over 3.37 y).
The parameters are grouped and labelled in the input: design data (alloy composition and density)
and CALIBRATED values (kinetics, mass transfer, cold capture, inventory, depletion):

```
[UserObjects]
  [dridn]
    type = DRIDNCorrosionUserObject
    mass_fractions = '0.07 0.05 0.71'        # design data: Hastelloy N Cr / Fe / Ni
    mass_transfer_rate = 1004.9168083798145  # um/y, CALIBRATED
    cold_capture_fraction = '0.314 0.666 0.765'   # CALIBRATED (full precision in the input)
    area_to_salt_mass = 7.767901069021746    # cm^2/g, CALIBRATED (shared by all loops)
    inventory_capacity = 839.1043065935307   # ppm, CALIBRATED
    ...
  []
[]
```

`AdvancedCorrosionModelPostprocessor`s expose the state: dissolved ppm, cumulative release
(`released_*_ppm`), bulk capture (`bulk_captured_*_ppm`), `mean_rate_um_y` (the time-averaged
recession, the measured quantity), the attack depths and `mass_balance_error`.

## Run

```bash
cd examples/training/chemistry/03_corrosion_models
$OPEN_PRONGHORN_EXE -i ncl16_capsule_bv.i      # serial, seconds
$OPEN_PRONGHORN_EXE -i ncl16_loop_flow.i       # serial, ~15 s (the slowest)
$OPEN_PRONGHORN_EXE -i ncl16_dridn.i           # serial, < 1 s, writes ncl16_dridn_out.csv
python3 plot.py                                # runs missing inputs; writes ../figures/case03_corrosion_models.png
```

Tests: `./run_tests --test-root examples/training --re=ncl16` (three CSVDiff tests).

## Expected results

![](../figures/case03_corrosion_models.png)

Top: the four measured quantities (dashed) and the four models (hatched bars = fitted to that
value). Bottom left: salt Cr against time. Bottom right: the hot-wall loss and cold-wall gain of
the flow wall law.

| quantity (3.37 y) | measured | (1) BV capsule | (2) flow wall law | (3) DRIDN | (4) MSTDB-TC (not run) |
|---|---|---|---|---|---|
| hot-leg rate [µm/y] (M-013) | 1.016 | 0.039 | 0.614 (P) | 1.045 (C) | 1.004 (C) |
| hot-leg loss [mg/cm^2] (M-011) | 2.9 | 0.116 | 1.84 (P) | 3.13 (C) | 3.00 (C) |
| cold-leg gain [mg/cm^2] (M-012) | 1.7 | 0 | 1.74 (C) | 1.61 (C) | 1.71 (C) |
| salt Cr [ppm] (M-014) | +500 | 900 | 511 (C) | 533 (C) | 491 (C) |
| salt Fe [ppm] (M-015) | -100 | -- | -127 (P, from an assumed 150 ppm) | **+106 (wrong sign)** | -100 (C, post-hoc formula) |

C = in-sample calibration; P = prediction of a calibrated model (that value was not fitted).

| check | kind | result |
|---|---|---|
| capsule `faraday_ratio` (recession = dissolved inventory) | VERIFICATION | 1.0008 |
| capsule salt Cr = Nernst limit (`c_eq_ratio`) | VERIFICATION | 0.99997 |
| flow law `inventory_balance_error` (loss - gain = dissolved) | VERIFICATION | 7.4e-4 |
| DRIDN `mass_balance_error` | VERIFICATION | 4e-16 |
| flow law: salt Cr and cold gain | CALIBRATION (2 parameters) | 511 ppm, 1.74 mg/cm^2 |
| flow law: hot loss and rate | prediction of a calibrated model, in-sample loop | 0.63x and 0.60x measured |
| DRIDN and MSTDB-TC, all five targets | CALIBRATION (fitted to NCL-16) | within 13 % except DRIDN Fe |

**None of this is a VALIDATION**: NCL-16 was used to build or tune every loop model here.

**(1) The capsule stalls at its Nernst limit**: half saturated after 31 d, 99 % after 192 d, then
the wall stops (0.13 µm of recession in total) and there is no cold leg. A closed isothermal system
cannot sustain corrosion for years.

**(2) The thermodynamic wall law makes a real loop.** After ~0.3 y the salt sits at 511 ppm,
between the wall equilibria (900 and 72 ppm), and from then on Cr is carried hot to cold at a steady
0.48 mg/cm^2/y (0.54 µm/y of hot-wall Cr recession). Loss minus gain (0.10 mg/cm^2) is the
dissolved inventory: no bulk sink is needed. The predicted loss and rate are 0.6x the measurement;
with the gain calibrated, only an extra sink or a hotter, more depleted surface could raise them,
and the single gap has no temperature profile along the legs.

**(3) DRIDN** matches the targets it was fitted to, but read it critically:

- **Fe has the wrong sign.** The salt Fe *decreased* by 100 ppm in the loop (Fe2+ reduced by Cr
  from the alloy). DRIDN starts the salt at ~1 ppm of each metal, so it can only add Fe (+106 ppm).
  The calibration matched M-015 through a separate post-hoc stoichiometric formula, not this state.
- **Ni is an artefact.** 378 ppm of dissolved Ni (as much as the Cr) is unconstrained by any
  measurement. Ni is the noblest component; purified fluoride loops keep it at tens of ppm at most.
- **Most of the metal goes to the bulk-capture sink.** The hot leg releases 6864 ppm of Cr: 72 %
  (4958 ppm) is bulk-captured, 20 % plates on the cold leg and only 8 % (533 ppm) stays in the salt.
  For all three metals, 64 % is bulk-captured. Loss minus gain (1.52 mg/cm^2) is therefore **not**
  the dissolved inventory (0.13 mg/cm^2). The sink stands for deposition and precipitation that the
  model does not resolve.
- `area_to_salt_mass` = 7.77 cm^2/g is a shared fitted value; a ~1.7 cm tube full of salt has
  $4/(D\rho)$ ≈ 1.1 cm^2/g. The instantaneous `front_rate_um_y` (1.54 µm/y) is not the measured
  time-averaged rate (`mean_rate_um_y` = 1.045). Intergranular depth: 37 µm.

**(4) MSTDB-TC** gives the closest endpoints (in-sample too), from physical saturation activities;
it cannot be reproduced without the licensed database.

## Limitations

- The flow wall law has no film mass transfer in series, no Cr(III)/Cr(II) chemistry and no
  coupling of $c_{eq}$ to the salt redox state (case 01): the equilibria are inputs.
- One hot and one cold wall of equal area stand for a loop with a temperature profile.
- DRIDN is a fitted closure; its bulk sink, Ni and Fe states are not physical predictions.
- Every model here is calibrated on the loop it is compared with.

## Exercises

1. *(easy)* Close the loop thermally: run `ncl16_loop_flow.i` with `T_cold=977.15` (both walls hot).
   Compare with model 1.
2. *(easy)* DRIDN: switch the cold-leg capture off
   (`"UserObjects/dridn/cold_capture_fraction=0 0 0"`), then double the inventory capacity
   (`UserObjects/dridn/inventory_capacity=1678.2`). What sets the hot-leg loss, and what sets the
   salt Cr plateau?
3. *(medium)* The MSTDB-TC CrF2 saturation activities (0.704 at 704 C, 0.483 at 538 C) imply
   $\Delta H = -R\ln(a_c/a_h)/(1/T_c - 1/T_h) \approx 15$ kJ/mol. Scan
   `dH_Cr=60 dH_Fe=60 dH_Ni=60`, `=150` and `=15` in model 2. How do the loss, the gain and the salt
   Cr respond, and can one $\Delta H$ fit all three?
4. *(hard)* Replace the selective split by congruent dissolution at the alloy mass fractions
   (`seed_Cr=0.07 seed_Fe=0.05 seed_Ni=0.71`). What happens to the Cr flux, the total loss and the
   salt Ni? Is congruent dissolution consistent with NCL-16? Which model would you use to predict a
   new loop, and what would you measure first?

<details><summary>Solutions</summary>

1. With both walls at 977 K the gap is the capsule of model 1: the salt Cr rises to the wall
   equilibrium (900 ppm), each wall loses 0.10 mg/cm^2 (the "gain" column reads -0.10: the second
   wall dissolves too), and then all corrosion stops (`hot_Cr_rate_um_y` ~ 0). Sustained loop
   corrosion needs the *difference* between the wall equilibria.
2. Without cold capture the hot-leg loss is unchanged (3.13 mg/cm^2, 1.045 µm/y): the fitted
   hot-leg kinetics and mass transfer set it. Nothing plates (gain 0), so the bulk sink takes more
   Cr (6340 instead of 4958 ppm) and the salt holds more Fe (106 -> 133 ppm) and Ni (378 -> 488);
   the Cr plateau barely moves (529 ppm). With twice the capacity the plateau rises to 831 ppm and
   the gain to 2.26 mg/cm^2, again with the same loss. The plateau is set by `inventory_capacity`
   and the bulk-capture closure, not by a hot/cold balance.
3. | $\Delta H$ [kJ/mol] | loss [mg/cm^2] | gain [mg/cm^2] | salt Cr [ppm] | cold $c_{eq}$ [ppm] |
   |---|---|---|---|---|
   | 15 | 0.45 | 0.26 | 834 | 617 |
   | 60 | 1.21 | 1.08 | 656 | 199 |
   | 100 (gold) | 1.84 | 1.74 | 511 | 72 |
   | 150 | 2.54 | 2.49 | 351 | 21 |

   A larger $\Delta H$ widens the gap between the hot and cold equilibria; the salt sits further
   from both, both walls are driven harder, and the transfer rises while the salt Cr falls. No
   single $\Delta H$ fits the loss (2.9), the gain (1.7) and the salt Cr (500) together. The ideal
   15 kJ/mol gives far too little transfer: the real driving force also needs the activity
   coefficients and the U(IV)/U(III) redox control that MSTDB-TC includes.
4. The steady hot-wall Cr flux falls 10x (0.048 against 0.48 mg/cm^2/y; Cr recession
   0.054 µm/y), but the total loss (1.83 mg/cm^2) and gain (1.74) barely change: they are now mostly
   **Ni** carried hot to cold, with the salt Ni at 0.57 ppm between its wall equilibria (1 and
   0.08 ppm). The flux is set by the seeded i0 and the distance from each $c_{eq}$, not by how much
   Ni the salt holds. NCL-16 deposits were Cr-rich and the attack Cr-selective, so the 90/7/3 split
   is the physical choice. For a new loop none of these models is validated; DRIDN and MSTDB-TC
   are calibrated on the same 76 cases (median factor error 1.1 in-sample,
   `validation/corrosion/calibration/MODEL.md`). The flow wall law with $c_{eq,k}(T)$ from
   thermochemistry, the measured salt redox state and the Fe/Ni impurity levels is the
   transferable choice, checked against an independent loop that was not used in any fit.

</details>

## Further reading

- [THEORY.md, section 4](../THEORY.md#4-corrosion-and-plating): Butler-Volmer seeding, the Nernst
  limit, thermal-gradient mass transfer, DRIDN.
- `CorrosionPlatingAction.md`, `CorrosionPlatingFlowAction.md`, `DRIDNCorrosionUserObject.md`
  (OpenPronghorn documentation); `validation/corrosion/calibration/MODEL.md`.
- ORNL-TM-4188 (NCL-16 loop data, https://www.osti.gov/servlets/purl/4597158) and the case and
  target tables `validation/corrosion/data/case_features.csv`, `targets.csv`.
- J. W. Koger, ORNL-TM-4189 (1972): Hastelloy N after nine years in a fluoride-salt loop.
- † J. R. Keiser, "Compatibility studies of potential molten-salt breeder reactor materials in
  molten fluoride salts", ORNL/TM-5783 (1977).

*References marked † were added from general literature knowledge and have not been checked against
the repository's sources; verify the bibliographic details before citing them.*
