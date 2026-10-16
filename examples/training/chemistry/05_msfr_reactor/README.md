# Case 05 (capstone): the MSFR primary circuit -- CFD with all three chemistry actions

**Module 4 (coupled, reactor-grade) | ~60 min | 2D-RZ, two steps on 8 ranks: flow ~1 min + chemistry ~6 min (heavy tests)**

## Objective

Put the whole day together on a reactor design: the 3 GWth **Molten Salt Fast Reactor** of the
Euratom EVOL project, in the geometry, salt, power shape and pump / heat-exchanger model of the INL
Virtual Test Bed (VTB). On a turbulent flow and energy solution, the fuel salt carries its redox
chemistry, element-resolved wall corrosion and plating, and the fission-product inventory with
helium bubbling and noble-metal plating. Check every chemistry result against a conservation law
or a 0D loop balance of the earlier cases, and learn what such a calculation can and cannot claim.

## What you learn

- How to run the three actions on a **CFD** solution: the advection by the Rhie-Chow velocity
  (`rhie_chow_user_object`), the turbulent diffusivity, and `time_derivative = false` for a steady
  passive-scalar solve inside `SIMPLE`.
- How to split a heavy calculation into two steps with **restart**: the flow is solved once and
  written to an Exodus file; the chemistry reads it (`initial_from_file_var`).
- How to make spatially varying source and removal terms from functor materials (fission power,
  bubbling in the pump, plating in the heat exchanger) and boundary-restricted plating.
- How to check a 2D field against the **0D balance** of cases 01-04, and why the hot/cold wall map
  follows from the temperature dependence of the equilibria (case 03).

## What is new compared with case 04

| case 04 | case 05 |
|---|---|
| 0D, regions as rates | **2D-RZ CFD** (8272 cells): k-epsilon, Boussinesq buoyancy, energy, pump and HX models |
| two actions | **all three actions** (`[MoltenSaltRadiolysis]`, `[CorrosionPlatingFlow]`, `[FissionProductInventory]`) in one solve |
| serial, seconds | **8 MPI ranks**, two steps, **restart** from `msfr_flow_restart.e`; heavy tests |
| constant rates | **functor** rates localized to the pump (bubbling) and HX (plating), and wall plating on boundaries |
| transient, adaptive steps | steady state: passive scalars converged inside the `SIMPLE` executioner with a direct (MUMPS) solve |
| measured MSRE data | a design without measured chemistry: VERIFICATION and CODE-TO-CODE only |

## Physics

The theory is in [THEORY.md, section 1.4](../THEORY.md#14-coupling-to-the-flow) (coupling to the
flow), [section 2.8](../THEORY.md#28-fission-oxidation-and-redox-control) (redox control),
[section 4.4](../THEORY.md#44-d-the-element-resolved-temperature-dependent-equilibrium-law-corrosionplatingflow)
(the wall law) and [section 3.5](../THEORY.md#35-removed-inventories-wall-plating) (plating).

**Plant.** The VTB MSFR (`virtual_test_bed/msr/msfr`, `steady/run_ns.i`): a 2D-RZ toroidal core
cavity (height ~2.6 m, outer radius 2.3 m) with one lumped external loop in revolution around the
axis: hot leg -> **pump** (body force) -> **heat exchanger** (porous block, Forchheimer friction,
volumetric conductance 600 m^2/m^3 x 20 kW/m^2/K to an 873 K secondary) -> cold leg. Fuel salt
LiF-ThF4-UF4 (77.5/20/2.5 mol%), rho 4284 kg/m^3, cp 1594 J/kg/K, mu 0.0166 Pa s, k 1.7 W/m/K.
The mesh and the Griffin fission power density (3 GW) come from the VTB coupled steady state,
stored as `msfr_vtb_mesh_power.e` (a copy of `steady/restart/run_ns_coupled_restart.e`, 1.3 MB).

**Flow (step 1).** Linear-FV SIMPLE with k-epsilon and wall functions. Differences from the VTB
(nonlinear INSFV, capped mixing length): the linear-FV solver (needed by the chemistry actions),
k-epsilon, the HX Forchheimer friction linearised with a lagged speed, and a pump force of
24 kN/m^3 instead of 20 kN/m^3 so that the mass flow rate matches the VTB one.

**Redox (case 01, in 2D).** U(III)/U ~ 1 % buffers the salt.

- *Radiolysis is charge consistent.* The dose is ~0.87 x the fission power (fission fragments).
  Every F atom is made with an electron, and in the hot, buffered salt the pair recombines within
  ~1e-10 s, so the net radiolytic oxidation is zero. The full radical set is too stiff for the
  segregated solve, so this zero net effect is imposed directly: G_e = G_F = G_F2 = 0 (case 01,
  exercise 1, shows the zero net effect with the radicals resolved).
- *Fission oxidizes the salt* by 0.8 equivalents per fission, U(III) -> U(IV), with the local
  Griffin fission rate.
- *Redox control.* The MSFR fuel-processing unit is modelled in the pump region as a relaxation of
  U(III) and U(IV) toward the target composition with a 1-day loop-averaged time constant. At steady
  state its net reduction equals the fission oxidation.

**Corrosion and plating (case 03, in 2D).** The `[CorrosionPlatingFlow]` walls (shield and
reflector, Hastelloy N class) use the thermodynamic wall law of case 03:

$$c_{eq,k}(T) = c_{eq,k}\,e^{-\Delta H/R\,(1/T-1/T_{ref})},\qquad
J_k = \frac{i_{0,k}(T)}{zF\sqrt{c_{eq,k}(T)}}\,\bigl(c_{eq,k}(T) - c_k\bigr)$$

with $c_{eq}$ = 4 (Cr), 0.05 (Fe), 1e-4 (Ni) mol/m^3 at 973 K, $\Delta H$ = 50 kJ/mol and selective
Cr leaching (`seed_fractions = '0.90 0.08 0.02'`). **These equilibria are illustrative
assumptions** (Cr ~50 wppm is the order of MSRE salt Cr at U(III)/U ~ 1 %); there are no MSFR data.

**Fission products (cases 02 and 04, in 2D).** I-135 -> Xe-135 and Mo-99, Ru-103, Te-132, with
**U235_fast** yields as a stand-in for U-233 (no U-233 set in the database). Xe-135 is stripped by
**helium bubbling in the pump** with the EVOL 30 s gas-extraction time (local rate
V_salt / V_pump / 30 s). The noble metals plate on the **HX plates** (a = 600 m^2/m^3,
Chilton-Colburn k_m from the VTB h = 20 kW/m^2/K) and on the **cavity walls** (turbulent Sherwood
correlation). The loop time (5.3 s) is far below the half-lives, so the loop inventory must match
the 0D balance of case 04:

$$\bar N = \frac{P/V}{\lambda + \int r\,dV/V + k_m A_{wall}/V}.$$

## The input, step by step

**Step 1, [msfr_flow.i](msfr_flow.i).**

- `[Mesh]` reads `msfr_vtb_mesh_power.e` and adds the `hx_top` / `hx_bot` sidesets, where the bulk
  temperatures are taken; `power_density` is read from the same file (`initial_from_file_var`).
- `[LinearFVKernels]` hold the momentum, pressure, k-epsilon and energy kernels, the pump
  (`LinearFVSource` in block `pump`), the HX friction (`LinearFVMomentumFriction`) and the HX sink.
- `[Postprocessors]`: `mdot`, `T_hx_in` / `T_hx_out` (mass-flow weighted bulk hot / cold leg),
  `T_max`, `energy_balance_error`. `[VectorPostprocessors]` sample the mid-plane, a line near the
  axis and the external loop.
- An extra Exodus output writes the converged fields for step 2:

```
[Outputs]
  [restart]
    type = Exodus
    file_base = msfr_flow_restart
  []
[]
```

**Step 2, [msfr_chemistry.i](msfr_chemistry.i).** `[Mesh]` reads `msfr_flow_restart.e`, and every
flow variable has `initial_from_file_var`, so SIMPLE starts converged; loose flow tolerances stop it
after its first iteration, and the passive-scalar (chemistry) loop runs to 2e-9 with MUMPS. Each
action is told it lives on a flow:

```
[MoltenSaltRadiolysis]
  block = ${fluid}
  salt_type = fluoride
  metals = 'U'
  temperature = T_fluid                   # the CFD temperature field
  dose_rate = dose_rate                   # 0.87 x the Griffin power density
  g_value_species = 'e_sol F_rad F2_diss'
  g_value_overrides = '0 0 0'             # charge-consistent: zero net radiolysis
  rhie_chow_user_object = 'rc'            # advection by the SIMPLE face velocities
  diffusivity = D_eff                     # molecular + turbulent (mu_t / Sc_t)
  time_derivative = false                 # steady passive scalar
[]
[CorrosionPlatingFlow]
  reaction_boundary = ${walls}            # shield_wall reflector_wall
  temperature = T_fluid
  equilibrium_concentrations = '${c_eq_Cr} ${c_eq_Fe} ${c_eq_Ni}'
  dissolution_enthalpies = '${dH_diss}'
  seed_fractions = '${f_Cr} ${f_Fe} ${f_Ni}'
  ...
[]
[FissionProductInventory]
  nuclides = 'I135 Xe135 Mo99 Ru103 Te132'
  yield_set = U235_fast
  power_density = power_density
  removal_groups = 'noble_gas noble_metal'
  removal_rates = 'bubbling_rate hx_plating_rate'     # functor materials, nonzero in pump / HX
  removal_destinations = 'offgas plated'
  plating_boundaries = ${walls}
  plating_mass_transfer_coefficient = ${k_m_wall}
  ...
[]
```

- Hand-added kernels, as in case 01: `fission_oxidation_UIII/UIV` (`LinearFVSource` of the power
  density x 0.8 eq/fission) and `control_*` (`LinearFVReaction` + `LinearFVSource` with the
  pump-region rate `ctrl_rate`).
- The functor materials `J_Cr`, `J_Fe`, `J_Ni` evaluate the wall law to report the wall fluxes
  (the exchange currents `i0_*` are copied from the action's verbose output).
- The postprocessors hold the checks below; `shield_profile` and `reflector_profile` sample the wall
  cells, `loop` and `midplane` the salt.

## Run

```bash
cd examples/training/chemistry/05_msfr_reactor
mpiexec -n 8 $OPEN_PRONGHORN_EXE -i msfr_flow.i        # ~1 min, writes msfr_flow_restart.e
mpiexec -n 8 $OPEN_PRONGHORN_EXE -i msfr_chemistry.i   # ~6 min
python3 plot.py        # writes ../figures/case05_msfr_flow.png and case05_msfr_chemistry.png
```

Tests: `./run_tests --heavy --test-root examples/training --re=msfr_` (three tests on 8 ranks,
~7 min in total; the profile test compares the VectorPostprocessor CSVs of the chemistry run).

## Expected results

**Flow and energy** (CODE-TO-CODE benchmark against the VTB, not a validation):

| quantity | this case | VTB (nonlinear INSFV, mixing length, rerun in this MOOSE) |
|---|---|---|
| mass flow rate | 19,270 kg/s | 19,220 kg/s |
| bulk hot / cold leg (HX in / out) | 997.4 / 902.0 K | dT 97.6 K |
| salt-average temperature | 956.8 K | 959.3 K |
| T_max (cell) | 1095 K | 1304 K |
| max speed | 2.07 m/s | 2.03 m/s |
| `energy_balance_error` | 2.2e-3 | - |

The EVOL nominal is 923 / 1023 K; the VTB secondary at 873 K puts both legs ~25 K lower.

*Why T_max is 210 K lower.* The hottest salt is in the slow, recirculating core centre, where the
power peaks. The VTB's capped mixing length gives a small, fixed eddy viscosity there; k-epsilon
transports the turbulence of the shear layers into the centre (mu_t up to 63 Pa s, ~4000 x mu) and
mixes the hot spot. The bulk temperatures and the flow rate agree because the pump/HX balance sets
them; the peak is a turbulence-model property. Neither model is validated for this flow (the VTB
tuned its mixing length on Nek5000 LES of an isothermal cavity and reports T_max ~1280 K): treat
T_max as uncertain by ~+/-200 K. *Mesh:* one uniform refinement (33,088 cells, 12 ranks, not
tested) gives mdot +6.7 %, legs 992.3 / 903.2 K and T_max 1102 K.

![](../figures/case05_msfr_flow.png)

**Chemistry checks:**

| check | kind | value (gold) |
|---|---|---|
| `redox_balance_error` (control = fission oxidation, 1.24e-4 mol eq/s) | VERIFICATION | 5e-10 |
| `U_total_avg` (uranium conserved) | VERIFICATION | 1200.0 mol/m^3 |
| `F2_max` (no F2 with charge-consistent radiolysis) | VERIFICATION | 0 |
| `Cr_wall_balance` (Cr dissolved hot = Cr plated cold) | VERIFICATION | 5e-6 |
| `Cr_vs_0d` (salt Cr = 0D wall balance) | VERIFICATION | 1.000 |
| `Xe135_vs_0d` (0D loop equilibrium) | VERIFICATION | 1.023 (2 % from the pump-local stripping) |
| `Mo99_vs_0d` | VERIFICATION | 1.008 |
| `Mo99_balance_error` (production = HX + walls + decay) | VERIFICATION | 3e-11 |

There is no measured MSFR chemistry, so no CALIBRATION or VALIDATION is possible here; case 04 is
the validation of the same models against the MSRE.

| model result (gold) | value |
|---|---|
| U(III)/U, uniformity | 0.962 %, uniform to 5e-5 (redox time days, loop 5.3 s) |
| fission oxidation = reductant demand | 10.7 mol eq/day (~2.7 mol, 0.6 kg, of Th per day) |
| salt Cr | 4.02 mol/m^3 = 48.8 wppm, between the hot / cold leg equilibria 4.65 / 2.46 |
| dissolving wall area | 40 % (hot upper reflector wall, hot leg and pump ducts) |
| mean Cr loss there / peak | 1.12 / ~4.7 mg/cm^2/y |
| Cr moved hot -> cold | 0.52 kg/y |
| Fe, Ni | at their equilibria (0.050, 1.0e-4 mol/m^3) |
| Xe-135 stripped by the bubbling | 99.94 % of its production |
| Mo-99 on the HX plates / cavity walls | 99.5 % / 0.5 % |

![](../figures/case05_msfr_chemistry.png)

Left: Xe-135 drops by 15 % across the pump (1.08 -> 0.93 x the salt average) and the fission source
restores it in the core; Mo-99 is removed in the HX. Middle and right: walls hotter than the
~960-970 K balance temperature dissolve, colder ones (HX ducts, cold leg, lower cavity, most of the
inner shield wall) plate at -0.7 to -1.0 mg/cm^2/y. The sign is the one of case 03: hot dissolves.

## Limitations

- The corrosion equilibria and $\Delta H$ are illustrative and are not coupled to U(III)/U. The
  wall law has no film mass transfer in series and no Cr(III)/Cr(II) chemistry.
- The noble metals are a perfect, irreversible sink (case 04 shows this fails for the MSRE), with
  no bubble attachment or particulates. U-235 fast yields stand in for U-233.
- 2D-RZ lumps the 16 loops into one annulus; k-epsilon sets T_max (see above).
- Steady state only: no start-up or power transients (case 02 is the transient counterpart).

## Exercises

1. *(easy)* The database default for a fluoride salt has a lumped G(F2) = 0.005 with no matching
   electrons (`data/msr_database.json`). Estimate by hand, from `dose_power`, the oxidation it
   would add here and compare it with `fission_oxidation_rate`. Would `redox_balance_error`
   notice it, and what does that tell you about the check?
2. *(medium)* Weaken the redox control 10x (`tau_redox=8.64e5`, 8 ranks, ~6 min). What happens to
   U(III)/U and its non-uniformity, and why does the reductant demand not change?
3. *(medium)* Make the corrosion congruent (`f_Cr=0.07 f_Fe=0.04 f_Ni=0.70`, the alloy's
   Cr/Fe/Ni mass fractions). Which results change, and which do not?
4. *(hard)* Turn the temperature dependence of the equilibria off (`dH_diss=1e-6`). Where does the
   Cr go now? What would the legacy constant-overpotential law of case 03, model 1, do in this loop?

<details><summary>Solutions</summary>

1. By hand: 2 x 0.005 x 2.61e9 W / (100 eV x N_A) = **2.7 mol eq/s**, ~22,000 x the fission
   oxidation (1.24e-4 mol eq/s). That is 2.3e5 mol/day: the whole ~28,000 mol uranium inventory in
   hours, clearly unphysical; the buffer would be exhausted and F2 would build up. The
   `redox_balance_error` compares the control with the *fission* oxidation only, so it would fail.
   A conservation check verifies only what it counts: a lumped F2 yield without electrons is a
   charge source. (`surviving_fraction_equivalent` = 2.3e-6: the fission oxidation equals only
   2 ppm of the gross radiolytic pairs, 54 mol/s.)
2. U(III)/U drops from 0.962 % to **0.622 %**: the steady state needs a 10x larger deficit below the
   1 % target to deliver the same reduction. It stays uniform (U(III) 7.457-7.458 mol/m^3) because
   10 days is still far longer than the 5 s loop. The reductant demand stays at **10.75 mol eq/day**:
   at steady state the control always delivers exactly the fission oxidation (0.8 eq/fission). The
   corrosion does not change at all (Cr 4.024 mol/m^3, 0.518 kg/y), because the Cr equilibria are
   inputs that are not coupled to U(III)/U, a limitation: a more oxidized salt raises $c_{eq,Cr}$.
3. The dissolved concentrations hardly change (Cr 4.020, Fe 0.050, Ni 1.0e-4 mol/m^3): each
   element settles at its own wall balance, set by the equilibria. The exchange currents change, and
   with them the fluxes: Cr's i0 drops 0.90/0.07 = 13x, and the Cr redistribution falls from 0.52 to
   **0.040 kg/y** (mean Cr loss 0.087 instead of 1.12 mg/cm^2/y). Ni's i0 grows 35x, but its flux
   stays negligible at its 1e-4 equilibrium. The dissolving area (40 %) is unchanged: it is set by
   the temperature field. Selective leaching is the observed mode for Hastelloy N (case 03).
4. With $c_{eq}$ independent of T the salt equilibrates at $c_{eq}$ (Cr 4.000 mol/m^3) and every wall
   flux vanishes: Cr redistribution 1e-5 kg/y and loss 1e-5 mg/cm^2/y, numerical noise (the
   "dissolving area" of 73 % is then the sign of that noise). No temperature dependence, no
   thermal-gradient mass transfer, as in case 03, exercise 1. With a constant overpotential, the
   Nernst limit $c_{ref}e^{zF\eta/RT}$ *falls* with T, so the loop would plate hot and dissolve cold:
   the wrong direction, which is why the thermodynamic wall law is used.

</details>

## Further reading

- [THEORY.md, section 1](../THEORY.md#1-numerical-basis) (linear FV, SIMPLE, coupling to the flow),
  [section 5](../THEORY.md#5-usage) (action usage and pitfalls) and
  [section 6](../THEORY.md#6-verification-calibration-validation).
- INL Virtual Test Bed, MSFR model: `virtual_test_bed/msr/msfr` (https://mooseframework.inl.gov/virtual_test_bed/).
- M. Brovchenko et al., "Neutronic benchmark of the molten salt fast reactor in the frame of the
  EVOL and MARS collaborative projects", EPJ Nuclear Sci. Technol. 5 (2019) 2.
- † J. Serp et al., "The molten salt reactor (MSR) in generation IV: overview and perspectives",
  Prog. Nucl. Energy 77 (2014) 308.
- `MoltenSaltRadiolysisAction.md`, `CorrosionPlatingFlowAction.md`, `FissionProductInventoryAction.md`.

*References marked † were added from general literature knowledge and have not been checked against
the repository's sources; verify the bibliographic details before citing them.*
