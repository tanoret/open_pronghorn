# Case 01: fission oxidation, the U(IV)/U(III) redox buffer and the cover gas

**Module 1 (radiolysis and redox) | ~45 min | 0D (1 cell), serial, ~1.5 min**

## Objective

Model the redox state of the MSRE fuel salt over 450 days of operation. Fission slowly oxidizes
the salt, the small U(III) inventory absorbs the oxidant, and only when this buffer is spent does
F2 survive and pass to the cover gas by Henry's law. Check the result against an analytic line and
three conservation laws.

## What you learn

- How one action, `[MoltenSaltRadiolysis]`, builds a whole species network (variables, reactions,
  G values, initial conditions, gas exchange) from `data/msr_database.json`.
- How to **add your own kernels** (`[LinearFVKernels]`) to the variables that an action creates,
  here for the fission oxidant and the beryllium addition.
- How to verify a chemistry model with **invariants** (uranium, fluorine, redox equivalents) and an
  analytic solution, and how to do a time-step study.
- Why radiolysis does not change the redox state of a hot fuel salt, and what does.

## What is new

This is the first case: one action, one cell, hand-added kernels, a `Transient` with an adaptive
time step, and `ParsedPostprocessor` checks. Every later case reuses these building blocks.

## Physics

The full derivation is in [THEORY.md, section 2](../THEORY.md#2-radiolysis-and-redox) (fission
oxidation: [2.8](../THEORY.md#28-fission-oxidation-and-redox-control); Henry exchange:
[2.6](../THEORY.md#26-henry-law-gas-exchange-and-degassing); invariants:
[2.7](../THEORY.md#27-invariants-element-and-charge-redox-conservation)). In short:

**Salt.** MSRE fuel salt 65LiF-29.1BeF2-5ZrF4-0.9UF4 at 600 C (2.27 g/cm^3): $8.0\times10^4$
mol/m^3 of fluoride and 490 mol/m^3 of uranium. MSRE kept U(III)/U near 1 % (4.9 mol/m^3) to
protect the Hastelloy N.

**Radiolysis moves charge but does not oxidize.** Each F atom made by radiation has a matching
solvated electron; in the hot salt they recombine (net F2 yield about zero; MSRE saw F2 only above
frozen, cold salt). The direct F2 yield is therefore set to zero (`g_value_overrides = '0'`).

**Fission oxidizes.** A fissioned UF4 frees four fluorides, but the fission products bind only
about 3.2 of them, so each fission leaves about 0.8 oxidizing equivalents. With
$\dot F = P/(200\ \mathrm{MeV}\,N_A)$ mol of fissions per m^3 and s, at the MSRE loop-averaged
4 MW/m^3 (8 MW in about 2 m^3):

$$R_{ox} = 0.8\,\dot F = 1.66\times10^{-7}\ \mathrm{eq\,m^{-3}\,s^{-1}},\qquad
t_{buffer} = \frac{[\mathrm{U(III)}]_0}{R_{ox}} = \frac{4.9}{1.66\times10^{-7}}\ \mathrm{s} = 342\ \mathrm{d}.$$

The oxidant is represented as nascent F2 (a fluoride sink $-R_{ox}$ and an F2 source
$+R_{ox}/2$, so fluorine and charge are conserved), which the database chemistry scavenges:

```
F2 + 2 U(III) -> 2 U(IV) + 2 F-      k = 5e4 m^6/(mol^2 s)      (database)
Be + 2 U(IV)  -> Be2+ + 2 U(III)     Be_rate                    (input kernels, exercise 3)
e_sol + U(IV) -> U(III),  2 F -> F2,  e_sol + F -> F-           (radiolysis, exercise 1)
```

**Cover gas.** Once the buffer is gone, F2 exchanges with a closed headspace:

$$\frac{dC_l}{dt} = -k_La\,(C_l - k_H R T\,C_g),\qquad
\frac{dC_g}{dt} = \frac{V_l}{V_g}\,k_La\,(C_l - k_H R T\,C_g)$$

with $k_La = 10^{-3}$ 1/s and $V_l/V_g = 20$ (MSRE-like), and an **assumed**
$k_H(\mathrm{F_2}) = 10^{-6}$ mol/(m^3 Pa) ($k_H R T = 7.26\times10^{-3}$): F2 is too reactive for
a measured Henry constant.

**Invariants** (all VERIFICATION):

| check | expression | postprocessor |
|---|---|---|
| uranium | U(III) + U(IV) = 490 mol/m^3 | `uranium_error` |
| fluorine | F- + F + 2 (F2,l + F2,g / (V_l/V_g)) = F0 | `fluorine_error` |
| redox | (U(III) + e) - (F + 2 F2,l + 2 F2,g/(V_l/V_g)) = U(III)_0 + (2 Be_rate - R_ox) t | `redox_error` |
| analytic | U(III) = max(U(III)_0 + (2 Be_rate - R_ox) t, 0) | `U_III_analytic` |

## The input, step by step

[fuel_salt_redox.i](fuel_salt_redox.i). **1. Parameters.** The salt, the power and the derived
oxidation rate are top-level variables, so every one can be changed on the command line:

```
P_fission = 4.0e6        # W/m^3, MSRE loop average
eq_per_fission = 0.8
R_ox = ${fparse eq_per_fission * P_fission / (200e6 * 1.602176634e-19) / 6.02214076e23}
Be_rate = 0              # mol Be/m^3/s (exercise 3)
dose_rate = 0            # W/m^3 (exercise 1)
```

**2. The action.** One block creates the seven species (`e_sol F_ion F_rad F2_diss U_IV U_III`
and the headspace `F2_gas`), their reactions and the gas exchange. The linear systems it builds
must be listed in `[Problem] linear_sys_names` and `Executioner/system_names`.

```
[MoltenSaltRadiolysis]
  salt_type = fluoride
  metals = 'U'                       # adds U(III)/U(IV) and their reactions
  gas_species = 'F2'                 # adds F2_gas and the Henry-law exchange
  kLa = ${kLa}
  volume_ratio = ${volume_ratio}
  dose_rate = ${dose_rate}
  g_value_species = 'F2_diss'
  g_value_overrides = '0'            # charge-balanced radiolysis: no direct F2
  initial_condition_species = 'F_ion U_IV U_III'
  initial_condition_values = '${F0} ${U_IV_0} ${U_III_0}'
[]
```

**3. Your own kernels.** Plain `LinearFVSource` kernels act on the action's variables:

```
[LinearFVKernels]
  [fission_F_ion_sink]
    type = LinearFVSource
    variable = F_ion
    source_density = ${fparse -R_ox}
  []
  [fission_F2_source]
    type = LinearFVSource
    variable = F2_diss
    source_density = ${fparse R_ox / 2}
  []
  ...                                # radiolysis_F_ion_sink, Be_U_III, Be_U_IV
[]
```

`radiolysis_F_ion_sink` removes the fluoride each radiolytic F atom comes from (the database
source does not), so the fluorine balance closes when the dose is on.

**4. Checks.** `ParsedPostprocessor`s compute `U_III_fraction`, `U_III_analytic`, `henry_ratio`
($C_l/(k_H R T\,C_g)$, 1 at equilibrium) and the three errors of the table above.

**5. Time integration** ([THEORY.md, section 1.2](../THEORY.md#12-segregated-fixed-point-coupling-transient)).
The segregated systems are coupled by a multi-system fixed point
(`multi_system_fixed_point = true`, 30 iterations). `IterationAdaptiveDT` grows the step from 1 s
to `dtmax = 12 h`. Time-step study (max `redox_error`): $2.3\times10^{-2}$ at 1 day,
$1.9\times10^{-3}$ at 12 h, $2\times10^{-5}$ at 6 h (twice the runtime). The error appears only
after the buffer is spent, from the lagged F2/U(III) coupling.

## Run

```bash
cd examples/training/chemistry/01_fuel_salt_redox
$OPEN_PRONGHORN_EXE -i fuel_salt_redox.i       # serial, ~1.5 min (931 time steps)
python3 plot.py                                # writes ../figures/case01_fuel_salt_redox.png
```

Test: `./run_tests --test-root examples/training --re=fuel_salt_redox`.

## Expected results

![](../figures/case01_fuel_salt_redox.png)

Left: U(III)/U falls on the analytic line. Middle: F2 is negligible until the buffer is spent,
then fills the salt and the cover gas. Right: the three invariants.

| quantity (gold) | value |
|---|---|
| U(III)/U at day 100 / 300 | 0.708 % / 0.123 % |
| U(III)/U < 0.01 % from | day 339 (analytic exhaustion 342) |
| F2 in salt / cover gas at day 450 | 0.0976 / 13.4 mol/m^3 |
| `henry_ratio` at day 450 | 1.0007 |
| max `uranium_error` / `fluorine_error` / `redox_error` | 1.1e-8 / 1.2e-7 / 1.9e-3 |

The headspace follows the salt within $1/(k_La(1 + (V_l/V_g)k_HRT)) \approx 15$ min, so the gas
and the salt stay at Henry equilibrium.

| check | kind | result |
|---|---|---|
| U(III) against the analytic line | VERIFICATION | agreement to 1e-12 relative before exhaustion |
| uranium, fluorine, redox invariants | VERIFICATION | 1.1e-8, 1.2e-7, 1.9e-3 (12 h steps) |
| buffer life (342 d) vs MSRE Be additions every months to a year | plausibility, not VALIDATION | same order |

No CALIBRATION or VALIDATION is claimed: $k_H(\mathrm{F_2})$ and the 0.8 eq/fission are
literature estimates, and no measured F2 history exists for a hot fuel salt.

## Limitations

- 0D: one well-mixed salt volume at a single temperature; no core/loop distinction.
- The fission products are not tracked as species; their fluorine deficit enters as a fixed
  0.8 eq/fission (it depends on the fuel and on the fission-product valences, so on the redox
  state itself).
- The uranium inventory does not burn (about 1.6 % of it would fission in 450 days).
- No corrosion: in a real loop, Cr from the alloy would also reduce the oxidant (case 03).
- The F2 Henry constant is an assumed value; the cover gas is a closed volume (no helium purge).

## Exercises

1. *(easy)* Switch the radiolysis on for one hour at the fission-fragment dose:
   `dose_rate=3.5e6 Executioner/end_time=3600 Executioner/dtmax=10`. Does it consume the buffer?
2. *(easy)* Raise the power density to an MSBR-like 40 MW/m^3:
   `P_fission=4e7 Executioner/end_time=5.184e6` (60 days). When is the buffer spent?
3. *(medium)* MSRE replenished U(III) with Be metal. Compute the Be rate that holds U(III)/U at
   1 % at 4 MW/m^3, check it with `Be_rate=...`, and convert it to grams per year for 2 m^3 of salt.
4. *(hard)* Remove the buffer (`U_III_0=0 Executioner/end_time=2.592e6`, 30 days). What happens to
   F2, and why does `henry_ratio` stay slightly above 1?

<details><summary>Solutions</summary>

1. No. The radiolytic source $S = G D/(100\,\mathrm{eV}\,N_A) = 7.3\times10^{-2}$ mol/(m^3 s) of
   e and F is about $4\times10^5\,R_{ox}$, but balanced: the electrons reduce U(IV) and the F atoms
   re-oxidize U(III). After 1 h F sits at a quasi-steady $5.7\times10^{-5}$ mol/m^3 and U(III) is
   4.89948 mol/m^3, within $8\times10^{-5}$ of the no-radiolysis line (4.89940). `redox_error` is
   $4.2\times10^{-6}$, `fluorine_error` $1.9\times10^{-10}$, and dissolved F2 only
   $3.0\times10^{-8}$ mol/m^3.
2. $R_{ox}$ scales with the power, so the buffer lasts 10x less: U(III)/U drops below 0.01 % at
   day 33.9 (analytic 34.2 d). At day 60: 0.233 mol/m^3 of F2 in the salt and 32.1 in the gas.
3. Each Be supplies two equivalents: $\dot n_{Be} = R_{ox}/2 = 8.29\times10^{-8}$ mol/(m^3 s).
   With `Be_rate=8.2914e-8` U(III)/U stays at 1.0000 % for 450 days (minimum 0.0099999750). For
   2 m^3: $8.29\times10^{-8}\times3.156\times10^7\times2 = 5.2$ mol, i.e. **47 g of Be per
   full-power year**, a few small Be-rod additions, which is what MSRE did.
4. With nothing to scavenge it, the fission F2 appears at once: after 1 day
   $9.5\times10^{-4}$ mol/m^3 in the salt and 0.121 in the gas (`henry_ratio` 1.08); after 30 days
   0.0272 and 3.73 mol/m^3 (`henry_ratio` 1.003). $C_l$ must stay above its equilibrium with the
   gas to drive the flux into the headspace; the excess, $\dot C_g/(k_La\,V_l/V_g\,k_HRT\,C_g)$,
   falls as the gas inventory grows.

</details>

## Further reading

- [THEORY.md, section 2](../THEORY.md#2-radiolysis-and-redox): radiolysis, G values, the
  fission-oxidation balance and Henry's law.
- `MoltenSaltRadiolysisAction.md` (OpenPronghorn documentation) and `data/msr_database.json`.
- † C. F. Baes, "The chemistry and thermodynamics of molten salt reactor fuels", J. Nucl. Mater.
  51 (1974) 149: the U(IV)/U(III) redox control and the oxidation per fission.
- † R. C. Robertson, ORNL-TM-0728 (1965) and P. N. Haubenreich, J. R. Engel, "Experience with the
  MSRE", Nucl. Appl. Technol. 8 (1970) 118: MSRE Be additions.

*References marked † were added from general literature knowledge and have not been checked against
the repository's sources; verify the bibliographic details before citing them.*
