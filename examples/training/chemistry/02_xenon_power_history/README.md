# Case 02: xenon and the fission-product inventory over a month of operation

**Module 2 (fission products) | ~40 min | 0D (1 cell), serial, ~30 s**

## Objective

Follow the fission-product inventory of a fuel salt through a realistic month: start-up, full
power, a load-follow dip, shutdown and ten days of decay. Read the different "memories" of the
inventory: Xe-135 (hours), I-131 (weeks) and Cs-137 (the integrated energy), and see how
sparging the noble gases changes the xenon reactivity.

## What you learn

- How the second action, `[FissionProductInventory]`, builds the Bateman chains (yields, decay,
  neutron capture and removal) of the nuclides you list, from `data/fission_product_database.json`.
- How to drive an input with **time-dependent functions**: the power and the flux follow one
  `PiecewiseLinear` power history.
- How a removal pathway (sparging to the off-gas) changes the xenon equilibrium and the
  post-shutdown xenon peak.
- How to convert an inventory into a reactivity worth and a burnup, and why the Cs-137 burnup
  monitor of a circulating, sparged fuel has a hidden bias.

## What is new compared with case 01

| case 01 | case 02 |
|---|---|
| `[MoltenSaltRadiolysis]` | `[FissionProductInventory]`: a different database and different physics |
| constant source terms | `power_density` and `neutron_flux` given as **functions** of time |
| steady operation | a **transient** with fast power changes: `FunctionDT` refines the time step around them, `scheme = bdf2` |
| hand-added kernels | the action's own removal pathway (`removal_groups`, `removal_rates`, `removal_destinations`) |
| conservation checks | integrated quantities: `TimeIntegratedPostprocessor` for the energy and the burnup |

## Physics

The full derivation is in [THEORY.md, section 3](../THEORY.md#3-fission-product-inventory) (xenon:
[3.7](../THEORY.md#37-xenon-equilibrium-stripping-and-worth); burnup monitoring:
[3.8](../THEORY.md#38-burnup-monitoring)). For
each nuclide $i$ (atoms in mol/m^3 of salt),

$$\frac{dN_i}{dt} = y_i\,\dot F(t) + \sum_j b_{j\to i}\lambda_j N_j - (\lambda_i + \sigma_i\phi(t) + r_i)\,N_i,
\qquad \dot F = \frac{P(t)}{E_f N_A}, \quad E_f = 200\ \mathrm{MeV}.$$

**Xe-135.** Fed by its direct yield and by I-135 (6.6 h); lost by decay (9.1 h,
$\lambda_X = 2.11\times10^{-5}$ 1/s), burnout ($\sigma_X\phi$ = 2.65 Mb x
$3\times10^{17}$ n/m^2/s = $7.9\times10^{-5}$ 1/s) and sparging $r = 10^{-4}$ 1/s. At equilibrium

$$N_X = \frac{y_{X,cum}\,\dot F}{\lambda_X + \sigma_X\phi + r},\qquad
\rho_X = -\frac{\sigma_X N_X\,\phi}{\nu\,\dot F}\quad(\nu = 2.43).$$

After a shutdown the burnout stops, but I-135 keeps feeding Xe: the xenon peak. Sparging removes
both the equilibrium xenon and the peak.

**I-131** (8.0 d) remembers the last weeks. **Cs-137** (30.1 y) integrates the energy:
$N_{Cs} \approx y_{Cs}E/(E_f N_A)$, the classic burnup monitor.

**Power history** (relative power $p(t)$; the flux follows the power, $\phi = \phi_0 p$):

| day | 0-1 | 1-10 | 10-13 | 13-20 | 20-30 |
|---|---|---|---|---|---|
| power | ramp 0 -> 100 % | 100 % | 50 % | 100 % | 0 (decay) |

$P_0$ = 20 MW/m^3 (a power reactor core average) and $\phi_0 = 3\times10^{17}$ n/m^2/s, the
salt-residence-weighted (loop-averaged) flux, consistent with the loop-averaged power density.

## The input, step by step

[power_history.i](power_history.i). **1. The power history** is one function, reused for the
power and the flux:

```
[Functions]
  [relative_power]
    type = PiecewiseLinear
    x = '0 ${day} ${fparse 10 * day} ${fparse 10 * day + 3600} ...'
    y = '0 1 1 0.5 0.5 1 1 0 0'
  []
  [power_density]
    type = ParsedFunction
    expression = '${P0} * p'
    symbol_names = 'p'
    symbol_values = 'relative_power'
  []
  [flux]  ...  expression = '${phi0} * p'
[]
```

**2. The action.** Four nuclides; the parents and daughters needed by the chain are linked from
the database. One removal pathway sends the `noble_gas` group (here Xe-135) to an off-gas
inventory `Xe135_offgas`:

```
[FissionProductInventory]
  nuclides = 'I135 Xe135 I131 Cs137'
  power_density = power_density        # a function name, not a number
  neutron_flux = flux
  removal_groups = 'noble_gas'
  removal_rates = '1e-4'               # 1/s sparging
  removal_destinations = 'offgas'
[]
```

The action creates one linear system per inventory; they are listed in `[Problem]
linear_sys_names = 'I135_sys Xe135_sys I131_sys Cs137_sys Xe135_offgas_sys'`.

**3. Derived quantities.** `rho_Xe_pcm` is the worth formula above with the full-power flux.
`energy` is a `TimeIntegratedPostprocessor` of the power density; `burnup_MWd_per_m3` and
`Cs137_burnup_check` = $N_{Cs}/(y_{Cs}E/(E_fN_A))$ follow from it.

**4. Time stepping.** `FunctionDT` gives 30 min steps, cut to 1 min for two hours after every fast
power change, with second-order `bdf2`.

## Run

```bash
cd examples/training/chemistry/02_xenon_power_history
$OPEN_PRONGHORN_EXE -i power_history.i         # serial, ~30 s
python3 plot.py                                # writes ../figures/case02_xenon_power_history.png
```

Test: `./run_tests --test-root examples/training --re=power_history`.

## Expected results

![](../figures/case02_xenon_power_history.png)

Top: the power history. Middle: the xenon worth follows every power change within hours, with a
small dip after each power reduction and a peak after the shutdown. Bottom: I-131 decays after
the shutdown while Cs-137 keeps the memory of the whole energy.

| quantity (gold) | value |
|---|---|
| equilibrium xenon worth at full power (day 10) | -1066 pcm |
| worth at the end of the 50 % plateau (day 13) | -665 pcm |
| transient dip after the step to 50 % | -1176 pcm, 2.5 h after it |
| post-shutdown xenon peak | **-1325 pcm, 2.5 h** after the shutdown |
| worth 6 h after the shutdown (day 20.25) | -1104 pcm |
| I-131 after 10 days of decay | 42 % of its maximum |
| burnup at day 30 | 360 MWd/m^3 |
| `Cs137_burnup_check` at day 10 / day 30 | 0.9997 / 0.9987 |

| check | kind | result |
|---|---|---|
| equilibrium worth against $\sigma_XN_X\phi/(\nu\dot F)$ with $N_X = y\dot F/(\lambda+\sigma\phi+r)$ = 1065 pcm | VERIFICATION | -1066 pcm |
| `Cs137_burnup_check` against $y_{Cs}E/(E_fN_A)$ | VERIFICATION | 0.9987: the 0.13 % deficit is the Cs-137 decay over the month |
| I-131 after 10 d of decay against $e^{-10\ln2/8.02}$ = 42 % | VERIFICATION | 42 % |

No CALIBRATION or VALIDATION is claimed: the power history and the flux are generic.

**The hidden Xe-137 loss.** Most of the mass-137 chain passes through Xe-137 (3.8 min), a noble
gas. Xe-137 is not in the database, so the action lumps it into the Cs-137 cumulative yield and it
is never sparged. With $r = 10^{-4}$ 1/s a fraction
$r/(\lambda_{Xe137}+r) = 10^{-4}/3.1\times10^{-3} \approx 3\,\%$ of the Cs-137 precursors would
really go to the off-gas, 20x the decay deficit (MSRE found Cs-137 in its off-gas system). A scratch
run with Xe-137 added to a copy of the database gives `Cs137_burnup_check` = 0.967, not 0.9987.

## Limitations

- 0D: the salt is one well-mixed volume. The core and the loop are lumped into a loop-averaged
  power and flux (case 05 shows when this is accurate).
- One-group worth, without importance weighting or spatial flux shape.
- Only the listed nuclides are tracked (no Xe-137, see above). U-235 thermal yields.
- The sparging rate is a single constant; in a real pump bowl it depends on the gas flow, the bubble
  size and the xenon solubility (case 04).

## Exercises

1. *(easy)* Check the equilibrium worth by hand from the formulas above
   ($y_{Xe,cum}$ = 0.0654, $\sigma_X$ = 2.65 Mb).
2. *(easy)* Turn the stripping off (`FissionProductInventory/removal_rates=1e-12`). What are the
   equilibrium worth and the post-shutdown peak?
3. *(medium)* The operator wants to restart 6 h after the shutdown. With and without stripping,
   how much excess reactivity is needed to override the xenon?
4. *(hard)* Increase the sparging 10x (`FissionProductInventory/removal_rates=1e-3`). What happens
   to the equilibrium and to the peak, and what limits further gains? Why is Cs-137 still a good
   burnup monitor in a circulating, sparged salt, and what would you correct?

<details><summary>Solutions</summary>

1. $\dot F = 2\times10^7/3.204\times10^{-11} = 6.24\times10^{17}$ fissions/(m^3 s);
   $\lambda_X+\sigma_X\phi+r = 2.11\times10^{-5}+7.94\times10^{-5}+10^{-4} = 2.01\times10^{-4}$ 1/s;
   $N_X = 0.0654\times6.24\times10^{17}/2.01\times10^{-4} = 2.04\times10^{20}$ atoms/m^3
   ($3.38\times10^{-4}$ mol/m^3, as `N_Xe135` at day 10);
   $\rho = 2.65\times10^{-22}\times2.04\times10^{20}\times3\times10^{17}/(2.43\times6.24\times10^{17}) = 1065$ pcm.
2. Without stripping the equilibrium is **-2127 pcm** (the denominator halves), the 50 % plateau
   -1760 pcm, and the post-shutdown peak **-4032 pcm at 8.5 h** (against -1325 pcm at 2.5 h with
   stripping): the I-135 inventory decays into xenon that can only leave by its own decay.
3. Read `rho_Xe_pcm` at day 20.25: **1104 pcm** with stripping (already past its peak) and
   **3928 pcm** without. Stripping turns a xenon-override problem into a routine restart.
4. With $r = 10^{-3}$ 1/s the equilibrium worth falls to **-194 pcm** and the peak almost vanishes
   (-196 pcm, 13 min after the shutdown). Once $r \gg \lambda_X + \sigma_X\phi$ the xenon is removed
   before it can absorb, and the worth scales as $1/r$; what remains is set by the I-135 that decays
   in the core. In practice the gas-liquid mass transfer limits $r$ (case 04: MSRE reached about
   $2\times10^{-3}$ 1/s with its circulating bubbles). Cs-137 is a good monitor because caesium is
   salt-seeking and barely decays over the operating time, but about 3 % of its chain is stripped as
   Xe-137 at $r = 10^{-4}$ 1/s, and **25 %** ($10^{-3}/(3.03\times10^{-3}+10^{-3})$) at
   $r = 10^{-3}$ 1/s. The monitor must count the Cs-137 in the off-gas system, or be corrected by
   $r/(\lambda_{Xe137}+r)$.

</details>

## Further reading

- [THEORY.md, section 3](../THEORY.md#3-fission-product-inventory): Bateman chains, xenon worth,
  removal pathways.
- `FissionProductInventoryAction.md` (OpenPronghorn documentation) and
  `data/fission_product_database.json` (yields from ENDF/B-VIII.0).
- J. R. Engel, R. C. Steffy, "Xenon behavior in the Molten Salt Reactor Experiment", ORNL-TM-3464
  (1971): the MSRE xenon poisoning and its stripping (case 04).
- † D. L. Hetrick, *Dynamics of Nuclear Reactors* (1971), chapter on xenon transients.

*References marked † were added from general literature knowledge and have not been checked against
the repository's sources; verify the bibliographic details before citing them.*
