# Chemistry models of the OpenPronghorn MSR training: theory and usage manual

This is the reference for the training day. It covers every chemistry model that the five training
cases use: radiolysis and redox, the fission-product inventory, and corrosion and plating. Each
equation was checked against the C++ source in `src/` and `include/` and against the bundled data
files in `data/`. Where the code and the existing documentation (`doc/content/theory/*.md`,
`doc/content/source/**`) disagree, **the code wins**, and the difference is marked
**[code vs doc]**.

## 0. Purpose, notation and model map

### 0.1 Notation

| symbol | meaning | unit |
|---|---|---|
| $C_i$, $c_k$ | concentration of a dissolved species or cation | mol/m$^3$ |
| $N_i$ | fission-product inventory per unit salt volume | mol/m$^3$ |
| $\mathbf u$ | salt velocity | m/s |
| $D$, $D_\text{eff}$ | molecular, effective (molecular + turbulent) diffusivity | m$^2$/s |
| $\mu_t$, $\rho$, $\mathrm{Sc}_t$ | turbulent viscosity, salt density, turbulent Schmidt number | Pa s, kg/m$^3$, - |
| $T$, $T_\text{ref}$ | temperature, reference temperature | K |
| $R$, $F$, $N_A$ | gas constant 8.314462618, Faraday constant 96485.33212, Avogadro 6.02214076e23 | J/(mol K), C/mol, 1/mol |
| $f = F/(RT)$ | inverse thermal voltage | 1/V |
| $\dot D$ | volumetric dose rate | W/m$^3$ = J/(m$^3$ s) |
| $G_i$ | radiolytic yield | molecules / 100 eV |
| $k(T)$, $A$, $E_a$ | rate constant, pre-exponential factor, activation energy | SI, SI, J/mol |
| $k_La$, $k_H$ | gas-liquid mass-transfer coefficient, Henry coefficient | 1/s, mol/(m$^3$ Pa) |
| $V_l/V_g$ | liquid-to-headspace volume ratio | - |
| $q'''$, $\dot F$ | fission power density, fission rate density | W/m$^3$, fissions/(m$^3$ s) |
| $E_f$ | recoverable energy per fission (default 200 MeV) | MeV |
| $Y_i$ | fission yield (independent IY, cumulative CY, effective) | per fission |
| $\lambda_i$, $b_{pi}$ | decay constant, branching ratio of parent $p$ to $i$ | 1/s, - |
| $\sigma_i$, $\phi$ | one-group (2200 m/s) capture cross section, one-group flux | barn (1e-28 m$^2$), 1/(m$^2$ s) |
| $r$ | first-order removal rate | 1/s |
| $k_m$, $A/V$ | wall mass-transfer coefficient, wetted area per salt volume | m/s, 1/m |
| $i$, $i_0$ | current density, exchange current density | A/m$^2$ |
| $z$, $\alpha_a$, $\alpha_c$ | charge number, anodic and cathodic transfer coefficients | - |
| $\eta$, $E_0$ | overpotential, standard potential | V |
| $c_\text{ref}$, $c_\text{eq}$ | reference concentration of $i_0$, equilibrium concentration | mol/m$^3$ |
| $\Delta H$ | dissolution enthalpy of the equilibrium | kJ/mol |
| $M$, $\rho_\text{alloy}$ | molar mass of the metal, alloy density | g/mol, g/cm$^3$ |
| $x_k$ | alloy mass fraction of element $k$ | - |
| $J$ | molar wall flux into the salt (positive = dissolution) | mol/(m$^2$ s) |

### 0.2 Model map

| physics | action / object | source | training case |
|---|---|---|---|
| radiolysis, mass-action chemistry, metal redox | `[MoltenSaltRadiolysis]` -> `MSRReaction`, `LinearFVSource` | `src/msr/actions/MoltenSaltRadiolysisAction.C`, `src/msr/linearfvkernels/MSRReaction.C` | 01, 05 |
| Henry-law gas exchange with a headspace | `MSRGasExchange` (built by the radiolysis action) | `src/msr/linearfvkernels/MSRGasExchange.C` | 01 |
| surface degassing, noble-metal wall plating | `MSRDegassingBC` | `src/msr/linearfvbcs/MSRDegassingBC.C` | 05 (plating) |
| fission oxidation, redox control | hand-added `LinearFVSource` kernels | the inputs | 01, 05 |
| Bateman chains, removal pathways | `[FissionProductInventory]` -> `FissionProductFeed`, `FissionProductLoss` | `src/msr/actions/FissionProductInventoryAction.C` | 02, 04, 05 |
| finite-element Butler-Volmer corrosion | `[CorrosionPlating]` -> `ButlerVolmerBC`, `ButlerVolmerInterface`, `CorrosionNernstPlanckFlux`, `CorrosionCurrentContinuity`, `CorrosionSolidDiffusion`, `CorrosionRateAux` | `src/corrosion/actions/CorrosionPlatingAction.C` | 03 (capsule), 04 |
| linear-FV Butler-Volmer wall law in a flow | `[CorrosionPlatingFlow]` -> `CorrosionLinearFVButlerVolmerBC` | `src/corrosion/actions/CorrosionPlatingFlowAction.C` | 03 (flow wall law), 05 |
| seed of $i_0$ (reduced correlation) | `MoltenSaltCorrosionModel` | `src/corrosion/base/MoltenSaltCorrosionModel.C`, `data/corrosion_database.json` | 03, 04, 05 |
| MSTDB-TC standard-state endpoint | `MSTDBTCStandardStateCorrosionModel`, `MSTDBTCCorrosionUserObject` | `src/corrosion/base/MSTDBTCStandardStateCorrosionModel.C` | 03 (read from CSV, not run) |
| DRIDN network | `DRIDNModel`, `DRIDNCorrosionUserObject`, `AdvancedCorrosionModelPostprocessor` | `src/corrosion/base/DRIDNModel.C` | 03 (DRIDN) |

The five cases:

| case | dimension | content |
|---|---|---|
| 01 `fuel_salt_redox.i` | 0D | MSRE fuel salt: fission oxidation, U(III)/U(IV) buffer, F2 to the cover gas |
| 02 `power_history.i` | 0D | I-135/Xe-135, I-131, Cs-137 through a month of power history, sparging, xenon worth |
| 03 `ncl16_*.i` | 1D / 0D | NCL-16 loop: BV capsule, element-resolved hot/cold wall law, DRIDN, MSTDB-TC endpoint |
| 04 `msre_loop.i` | 0D multi-region | MSRE: xenon (spray, graphite, bubbles), noble metals, corrosion |
| 05 `msfr_flow.i`, `msfr_chemistry.i` | 2D-RZ | VTB MSFR primary circuit: flow + all three actions |

## 1. Numerical basis

### 1.1 One linear finite-volume system per species

The three MSR actions (`[MoltenSaltRadiolysis]`, `[FissionProductInventory]`,
`[CorrosionPlatingFlow]`) create one `MooseLinearVariableFVReal` per species, each in its own linear
solver system named `<species>_sys` (`systemName()` in the action headers). Each species obeys

$$\frac{\partial C_i}{\partial t} + \nabla\cdot(\mathbf u\,C_i) - \nabla\cdot(D_\text{eff}\nabla C_i)
= a_i(\mathbf C^*) - b_i(\mathbf C^*)\,C_i ,$$

where the source $a_i$ and the implicit loss coefficient $b_i$ are evaluated with the latest
available values $\mathbf C^*$ of the *other* species. The kernels are:

| term | kernel | added when |
|---|---|---|
| $\partial C/\partial t$ | `LinearFVTimeDerivative` | `time_derivative = true` (default) |
| advection by a flow solution | `LinearFVScalarAdvection` | `rhie_chow_user_object` given |
| advection by a constant velocity | `LinearFVAdvection` | `velocity` given (and no Rhie-Chow) |
| diffusion | `LinearFVDiffusion` (`use_nonorthogonal_correction = false`) | `diffusivity` given |
| source density $\times$ factor | `LinearFVSource` | G values, fission yields, hand-added sources |
| reaction, decay, removal | `MSRReaction`, `FissionProductFeed`, `FissionProductLoss`, `MSRGasExchange` | from the databases |

A loss is always put on the matrix diagonal (implicit in the species itself). A production is always
put on the right-hand side, using the lagged value of its parent or reactant. Consequences:

- every system is linear and positive (a diagonal M-matrix for upwind advection), so a single
  linear solve per system per iteration is enough;
- a reaction couples systems only through the lagged values, so the coupled system is solved by a
  **segregated fixed-point (Gauss-Seidel) iteration** over the systems;
- the conservation of an element shared by several species (for example U between U(III) and U(IV))
  is exact only **at the converged fixed point**, when the lagged and current values coincide.

The finite-element `[CorrosionPlating]` action is different: its variables are nonlinear
finite-element variables solved by Newton with AD Jacobians (`solve_type = NEWTON`, cases 03 capsule
and 04).

### 1.2 Segregated fixed-point coupling (transient)

For a chemistry-only or prescribed-velocity problem, the executioner is a `Transient` with

```
[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 30
    converge_at_max_iterations = true
  []
[]
[Executioner]
  type = Transient
  system_names = 'e_sol_sys F_ion_sys ...'     # every system, in solve order
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
[]
```

`IterationCountConvergence` (MOOSE framework) does **not** test a residual: it reports convergence
after `max_iterations` sweeps when `converge_at_max_iterations = true` (and divergence otherwise). The
number of sweeps is therefore a numerical parameter of the model. `doc/content/theory/moltenSaltRadiolysis.md`
shows that the total-Cr conservation error of a closed irradiated cell falls roughly as
1/(number of sweeps) (1.5e-3 at 1 sweep, 4.7e-5 at 50 sweeps). Case 01 uses 30 sweeps; case 03 (flow
wall law) uses 10.

Without `multi_system_fixed_point` (case 02, and case 04 which also holds the nonlinear corrosion
system `nl0`), there are no fixed-point sweeps: each system is solved once per time step. For a
decay chain that is enough, because `[FissionProductInventory]` sorts its inventories parent before
daughter (section 3.6): the lagged parent value is already the new one.

### 1.3 Steady SIMPLE versus transient

| executioner | use | time derivative |
|---|---|---|
| `Transient` + `multi_system_fixed_point` | 0D/1D chemistry, prescribed flow (cases 01, 03 flow) | `time_derivative = true` |
| `Transient` (one pass per step) | parent-before-daughter chains (case 02), mixed with a Newton system (case 04) | `true` |
| `SIMPLE` | steady state on a segregated Navier-Stokes solution (case 05) | `time_derivative = false` |
| `PIMPLE` | transient on a segregated Navier-Stokes solution (not used in the training) | `true` |

In `SIMPLE` the species are `passive_scalar_systems`. Each SIMPLE iteration solves the momentum,
pressure, turbulence and energy systems and then each passive scalar once, so the fixed-point
coupling of the chemistry is iterated by the SIMPLE loop itself. `num_iterations`,
`passive_scalar_absolute_tolerance` and `passive_scalar_equation_relaxation` control it.

### 1.4 Coupling to the flow

- **Advection.** With `rhie_chow_user_object = 'rc'` the species are advected by the face mass flux
  of the `RhieChowMassFlux` user object of the flow solve (`LinearFVScalarAdvection`); the advected
  interpolation method is set with `advected_interp_method_name` (case 05: `upwind`).
  `MoltenSaltRadiolysisAction` can add a buoyant slip velocity to the gas-phase species
  (`gas_buoyancy_velocity` -> `u_slip, v_slip, w_slip`) and a gas dispersivity.
- **Turbulent diffusion.** The actions take a diffusivity *functor*. In case 05 it is

  $$D_\text{eff} = D + \frac{\mu_t}{\rho\,\mathrm{Sc}_t},\qquad D = 2\times10^{-9}\ \mathrm{m^2/s},\ \mathrm{Sc}_t = 1,$$

  built with a `ParsedFunctorMaterial` named `D_eff` from the k-epsilon turbulent viscosity `mu_t`.
- **Temperature.** Every rate constant, Henry exchange and Butler-Volmer law reads the temperature as
  a functor (case 05: `T_fluid`), so the chemistry follows the energy solution one-way.
- **Dose and fission rate.** These are functors too (case 05: the Griffin `power_density` read from
  the restart file).

### 1.5 System ordering

The order of `linear_sys_names` and of `system_names` / `passive_scalar_systems` is the Gauss-Seidel
order. Two rules:

1. **List every system in both places**: in `[Problem] linear_sys_names` (which creates the system)
   and in the executioner (which solves it). A system that is created but not listed in the executioner is not solved by it.
2. **Parents before daughters.** `[FissionProductInventory]` prints the order it needs with
   `verbose = true` (and the full `linear_sys_names` string). For the radiolysis species the order
   matters less, because the fixed-point iteration converges the coupling anyway.

### 1.6 When tolerances and relaxation matter

- **Lagged stiff coupling.** In case 01 the error of the redox balance comes from the lagged
  F2 / U(III) coupling when the buffer runs out: 1.9e-3 at `dtmax` = 12 h, 2e-5 at 6 h, 2.3e-2 at
  1 d (comment in the input). The radical chemistry is very stiff (bimolecular constants up to
  ~1e10 m^3/mol/s); `doc/content/theory/moltenSaltRadiolysis.md` notes that strongly driven
  problems converge slowly in the segregated scheme.
- **Closed loops (case 05).** A closed loop with only a weak wall exchange has, for each dissolved
  metal, an almost singular steady-state system: the steady level is set by the small balance
  $\sum_\text{walls} w\,(c_\text{eq}(T) - c) = 0$, not by a source/sink of order one. Case 05 therefore uses
  - a **direct** linear solve for the passive scalars (`lu` with `mumps`), because an iterative
    solver stalls on the nearly singular system;
  - **no relaxation** of the passive scalars (`passive_scalar_equation_relaxation = 1`): under-relaxing
    an almost singular equation slows the convergence of its floating level further;
  - a tight `passive_scalar_absolute_tolerance = 2e-9` and `num_iterations = 3000`, with loose flow
    tolerances so that the restarted (converged) flow stops after its first iterations and the rest
    of the loop is chemistry;
  - an initial Cr concentration at the reference equilibrium (`initial_concentration = c_eq_Cr`),
    which shortens the approach to the floating level.
- **Residual-level checks.** Balances of nearly cancelling sums (case 05: `Cr_wall_balance`,
  `redox_balance_error`, `Mo99_balance_error`) depend on the tolerance; their test only checks the
  order of magnitude (`override_rel_err` in `05_msfr_reactor/tests`).

## 2. Radiolysis and redox

### 2.1 The dose

The radiolytic source is driven by the volumetric dose rate $\dot D$ [W/m$^3$], given to
`[MoltenSaltRadiolysis]` as the functor `dose_rate` (default 0). In a **fuel** salt, the dose is not
gamma dose: the fission fragments deposit their kinetic energy (~84 % of the ~200 MeV recoverable per
fission) in the salt itself, and the betas and part of the gammas add a few %. The training therefore
takes

$$\dot D \approx (0.85\text{-}0.9)\;q''' ,$$

with `dose_fraction = 0.87` in case 05 (`dose_mat`: `dose_rate = 0.87 * power_density`) and
$\dot D \approx 0.875\,q'''$ in exercise 1 of case 01. The database yields, however, are **gamma**
yields (`radiation = gamma` is the only option of the action); fission-fragment yields in molten
fluorides are not in the database (limitation, section 2.9).

### 2.2 G values and the radiolytic source

$$S_i = \frac{G_i\,\dot D}{100\ \mathrm{eV}\cdot N_A}\quad[\mathrm{mol/(m^3\,s)}]$$

(`MSR::gToSource`: factor $1/(100 \cdot 1.602176634\times10^{-19}\cdot N_A)$). The action adds one
`LinearFVSource` per species with a nonzero G value, with `source_density = dose_rate` and
`scaling_factor` $= G_i/(100\,\mathrm{eV}\,N_A)$. Default gamma G values (`data/msr_database.json`):

| salt | $G(e^-_s)$ | $G(\mathrm{X}^\cdot)$ | $G(\mathrm X_2)$ |
|---|---|---|---|
| chloride | 0.3 | 0.3 (Cl$^\cdot$) | 0.0 |
| fluoride | 0.2 | 0.2 (F$^\cdot$) | 0.005 |

They are overridden with `g_value_species` / `g_value_overrides`. **Note:** the radiolytic source of
X$^\cdot$ is not balanced by a matching sink of X$^-$ in the database; case 01 adds it by hand
(`radiolysis_F_ion_sink`) to close the fluorine balance when the dose is on. A nonzero $G(\mathrm F_2)$
without the matching electrons breaks the charge balance; cases 01 and 05 set it to 0.

### 2.3 Mass-action reactions and the rate constant

Each reaction $\sum_r \nu_r R_r \to \sum_p \nu_p P_p$ has the rate

$$\dot r = k(T)\prod_r C_r^{\,o_r},\qquad o_r = \nu_r \text{ unless an 'orders' entry overrides it},$$

and the rate constant is chosen, in this order (`ReactionData::rate`):

$$k(T) = \begin{cases}
A\,e^{-E_a/(RT)} & A>0 \text{ and } E_a>0\\
k_\text{ref}\,\exp\!\left[\dfrac{E_a}{R}\left(\dfrac{1}{T_\text{ref}}-\dfrac{1}{T}\right)\right] & k_\text{ref},T_\text{ref},E_a>0\\
k_\text{ref} & \text{only } k_\text{ref}>0
\end{cases}$$

(a zero marks an unsupplied value). The action splits every reaction into `MSRReaction` kernels:

- **consumption** (one per reactant, matrix): coefficient $\nu_r\,k\,\prod_{s\ne r}C_s^{o_s}\,C_r^{\,o_r-1}$
  on the diagonal of reactant $r$ (the self factor $C_r^{o_r-1}$ is used only when $o_r>1$, with the
  lagged $C_r$);
- **production** (one per product, right-hand side): $\nu_p\,k\,\prod_r C_r^{o_r}$ with the lagged
  reactants.

Partner concentrations are floored at 0. Species that are not tracked in the problem are skipped.

**[code vs doc]** For the 2 X$^\cdot$ $\to$ X$_2$ reactions the database lists the reactant once
with coefficient 2 (`{"Cl_rad": 2}`), so the rate is $k\,C^2$ and the consumption of X$^\cdot$ is
$2k\,C^2$; the reaction names (`Cl_rad + Cl_rad -> Cl2_diss`) suggest two separate terms, but the code
uses one term with $\nu = o = 2$.

### 2.4 The reaction networks (`data/msr_database.json`)

$k_\text{ref}$ in m$^3$/(mol s) for bimolecular and 1/s for unimolecular steps;
$k_\text{ref}$ for $\mathrm{F_2} + 2\,\mathrm{M}$ is m$^6$/(mol$^2$ s) (third order, $o = 1, 2$).

**Chloride kernel** (species `e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss`):

| reaction | $k_\text{ref}$ | $T_\text{ref}$ [K] | $E_a$ [kJ/mol] |
|---|---|---|---|
| Cl$^\cdot$ + Cl$^-$ -> Cl$_2^{\cdot-}$ | 1e7 | 673.15 | 20 |
| e$^-_s$ + Cl$^\cdot$ -> Cl$^-$ | 1e7 | 673.15 | 25 |
| 2 Cl$^\cdot$ -> Cl$_2$ | 5e6 | 673.15 | 20 |
| 2 Cl$_2^{\cdot-}$ -> Cl$_3^-$ + Cl$^-$ | 2.2e6 | 673.15 | 26 |
| Cl$_3^-$ -> Cl$_2$ + Cl$^-$ | 1 | 673.15 | 10 |

**Fluoride kernel** (species `e_sol F_ion F_rad F2_diss`):

| reaction | $k_\text{ref}$ | $T_\text{ref}$ [K] | $E_a$ [kJ/mol] |
|---|---|---|---|
| e$^-_s$ + F$^\cdot$ -> F$^-$ | 1e7 | 673.15 | 25 |
| 2 F$^\cdot$ -> F$_2$ | 5e6 | 673.15 | 20 |

### 2.5 The metals (`metals = 'Zn Cr U'`)

| metal | species | salt | reaction | rate constant |
|---|---|---|---|---|
| Zn | `Zn_II Zn_I` | chloride | e$^-_s$ + Zn(II) -> Zn(I) | $A$ = 2.4e10, $E_a$ = 35.6 kJ/mol |
| | | chloride | Cl$_2^{\cdot-}$ + Zn(I) -> 2 Cl$^-$ + Zn(II) | 2e7 at 673.15 K, 25 kJ/mol |
| | | fluoride | e$^-_s$ + Zn(II) -> Zn(I) | 1e7 at 873.15 K, 30.31 kJ/mol |
| | | fluoride | F$_2$ + 2 Zn(I) -> 2 Zn(II) + 2 F$^-$ | 1e4 at 873.15 K, 30 kJ/mol |
| Cr | `Cr_II Cr_III Cr_I` | chloride | e$^-_s$ + Cr(II) -> Cr(I) | $A$ = 1.7e10, 33.5 kJ/mol |
| | | chloride | e$^-_s$ + Cr(III) -> Cr(II) | $A$ = 2.0e10, 31.8 kJ/mol |
| | | chloride | Cl$_2^{\cdot-}$ + Cr(II) -> Cr(III) + 2 Cl$^-$ | 7.2e6 (constant) |
| | | chloride | Cl$_2^{\cdot-}$ + Cr(III) -> Cr(II) + Cl$_2$ | 1.4e6 (constant) |
| | | chloride | Cr(III) + Cr(I) -> 2 Cr(II) | 1.7e7 (constant) |
| U | `U_IV U_III` | fluoride | e$^-_s$ + U(IV) -> U(III) | 1e7 at 873.15 K, 30 kJ/mol |
| | | fluoride | F$_2$ + 2 U(III) -> 2 U(IV) + 2 F$^-$ | 5e4 at 873.15 K, 30 kJ/mol |

Cr has no fluoride reactions and U no chloride reactions: in that salt the metal adds its species but
no chemistry (`metalReactions` returns an empty list). The Cr rates without $E_a$ are temperature
independent.

### 2.6 Henry-law gas exchange and degassing

With `gas_species = 'F2'` (or `Cl2`) the action adds the headspace species `F2_gas` and two
`MSRGasExchange` kernels:

$$\frac{dC_l}{dt} = -k_La\,\bigl(C_l - k_H R T\,C_g\bigr),\qquad
\frac{dC_g}{dt} = \frac{V_l}{V_g}\,k_La\,\bigl(C_l - k_H R T\,C_g\bigr).$$

In each equation the own concentration is implicit and the partner is lagged. Henry coefficients:
$k_H(\mathrm{Cl_2}) = 10^{-5}$, $k_H(\mathrm{F_2}) = 10^{-6}$ mol/(m$^3$ Pa); both are **assumed**
(no measurement; F2 is reactive). The gas is a well-mixed closed headspace (no purge). The invariant
is $C_l + C_g/(V_l/V_g)$ (moles per unit liquid volume). Equilibrium is $C_l = k_H R T\,C_g$, the
`henry_ratio` of case 01.

`MSRDegassingBC` is the boundary-flux form, for a free surface swept by cover gas:

$$J_\text{out} = k_\text{surf}\,\bigl(C - k_H\,p_\text{cover}\bigr),$$

with the face value taken as the cell value (no penetration of the melt). Defaults: $k_H = 0$,
$p_\text{cover} = 0$. It is not used for gases in the five cases, but it is the wall-plating BC of
section 3.5.

### 2.7 Invariants: element and charge (redox) conservation

Each database reaction conserves the elements and the charge. Mass-action kinetics therefore leave
linear combinations of the concentrations invariant, which are the **verification** checks of the
redox cases. For the fluoride kernel with U and F2 exchange (case 01):

| invariant | expression | postprocessor |
|---|---|---|
| uranium | $C_\text{U(III)} + C_\text{U(IV)} = U_0$ | `uranium_error` |
| fluorine | $C_\mathrm{F^-} + C_\mathrm{F^\cdot} + 2\,(C_\mathrm{F_2,l} + C_\mathrm{F_2,g}/(V_l/V_g)) = F_0$ | `fluorine_error` |
| redox | $(C_\text{U(III)} + C_{e^-_s}) - (C_\mathrm{F^\cdot} + 2C_\mathrm{F_2,l} + 2C_\mathrm{F_2,g}/(V_l/V_g)) = C_{\text{U(III)},0} + 2\dot n_\text{Be}\,t - R_\text{ox}\,t$ | `redox_error` |

The redox invariant counts reducing equivalents (U(III), solvated electrons) minus oxidizing
equivalents (F atoms, F2). Radiolysis creates e$^-_s$ and F$^\cdot$ in equal numbers when
$G(e^-_s) = G(\mathrm F^\cdot)$ and $G(\mathrm F_2) = 0$, so it does not change the invariant: **radiolysis
only moves charge around**. As section 1.1 explains, the invariants hold exactly only at the
converged fixed point; their residual measures the iteration and time-step error.

### 2.8 Fission oxidation and redox control

The net oxidant of a fuel salt is **fission**. A fissioned UF4 frees four fluorides, but the fission
products bind on average only ~3.2 of them, so each fission leaves ~0.8 oxidizing equivalents (ORNL
MSRE experience; the value used in cases 01 and 05):

$$R_\text{ox} = 0.8\,\dot F = 0.8\,\frac{q'''}{E_f N_A}\quad[\mathrm{eq/(m^3\,s)}].$$

This term is **not** part of the radiolysis action; the inputs add it with `LinearFVSource`:

- **case 01** (charge-consistent nascent F2): a fluoride sink $-R_\text{ox}$ on F$^-$ and an F2 source
  $+R_\text{ox}/2$, which the database reaction F2 + 2 U(III) -> 2 U(IV) + 2 F$^-$ scavenges. Net:
  U(III) -> U(IV) at $R_\text{ox}$. With $q''' = 4$ MW/m$^3$, $R_\text{ox} = 1.66\times10^{-7}$ eq/(m$^3$ s), and
  the 4.9 mol/m$^3$ U(III) buffer (U(III)/U = 1 %) lasts $4.9/R_\text{ox} = 342$ d; only then does F2
  survive and pass to the cover gas;
- **case 05** (direct): U(III) sink and U(IV) source $\pm R_\text{ox}$ with
  `scaling_factor` $= 0.8/(E_f N_A)$ [mol eq/J] times the `power_density`.

**Redox control.** A real plant needs a reductant. Case 01 provides a beryllium addition
Be + 2 U(IV) -> Be$^{2+}$ + 2 U(III) (`Be_rate`, 0 by default) as two hand-added sources.
Case 05 models the fuel-processing unit as a relaxation in the pump region,
$dC_\text{U(III)}/dt = k\,(C_{\text{U(III)},0} - C_\text{U(III)})$ and the same for U(IV), with
$k = (V_\text{salt}/V_\text{pump})/\tau$, $\tau$ = 1 d. At steady state its net reduction must equal
the fission oxidation (`redox_balance_error`, 5e-10).

In case 05 the radiolytic yields are set to zero (`g_value_overrides = '0 0 0'` for `e_sol F_rad F2_diss`):
the e$^-_s$/F$^\cdot$ pairs recombine within ~1e-10 s in the hot salt, so their net effect is zero,
and the full radical set is too stiff for the segregated iterations. Case 01 also runs with
`dose_rate = 0` in the base run: the radiolytic cycle is closed and turns over ~4e5 times faster than
the fission oxidation, which would force second-scale steps over a year.

### 2.9 Limitations of the radiolysis and redox model

- **Assumed G values.** The database holds gamma yields for generic chloride and fluoride melts. They
  are illustrative; fission-fragment (high-LET) yields are not available.
- **Recombination.** Molten fluorides recombine almost completely; MSRE saw F2 only from frozen salt
  (below ~100 C). The model recombines only through the listed reactions, so a nonzero dose with the
  default G values gives a transient radical population, not a measured net yield.
- **Missing back-reactions.** All reactions are irreversible forward steps (no equilibrium constant,
  no reverse rate). There is no thermodynamic floor: for example, F2 + 2 U(III) runs to completion,
  and the Cr chloride cycle has no reverse of Cl$_2^{\cdot-}$ + Cr(III) -> Cr(II) + Cl$_2$.
- **Fission oxidation and processing** are hand-added sources with a fixed 0.8 eq/fission; the
  actual value depends on the fission-product valences, and impurities (moisture, oxides) dominate the
  redox of real salts.
- **Uranium burnup** is not modelled (at 4 MW/m$^3$ it would be ~1.6 % of U in 450 d, case 01).
- **Stiffness.** The segregated scheme converges slowly for strongly driven radical chemistry; the
  theory page recommends a monolithic Newton solve for such problems.

## 3. Fission-product inventory

### 3.1 The transport Bateman equations

`[FissionProductInventory]` solves, for every tracked nuclide $i$ in the salt,

$$\frac{\partial N_i}{\partial t} + \nabla\cdot(\mathbf u N_i) - \nabla\cdot(D\nabla N_i)
= \frac{Y_i\,\dot F}{N_A} + \sum_{p\ \text{tracked}}\bigl(\lambda_p b_{pi} + \sigma_p\phi\,\delta_{p\to i}\bigr)N_p
- \Bigl(\lambda_i + \sigma_i\phi + \sum_{k\in\text{pathways}(i)} r_k\Bigr)N_i ,$$

with $\lambda = \ln 2/t_{1/2}$ (0 for a stable nuclide) and $\sigma$ in barn ($10^{-28}$ m$^2$).
$\delta_{p\to i} = 1$ when $i$ is the capture product of $p$. The fission rate density is either given
(`fission_rate_density`) or computed from the power density:

$$\dot F = \frac{q'''}{E_f\,(1.602176634\times10^{-13}\ \mathrm{J/MeV})},\qquad E_f = 200\ \mathrm{MeV}\ (\texttt{energy\_per\_fission}).$$

Exactly one of the two must be given. The kernels: `LinearFVSource` (direct yield, factor
$Y_i/(N_A E_f)$ or $Y_i/N_A$), `FissionProductFeed` (right-hand side,
$(\lambda_p b + \sigma_p\phi + r)\,N_p$ with the lagged parent) and `FissionProductLoss` (diagonal,
$\lambda + \sigma\phi + \sum r_k$).

### 3.2 Yields: effective, cumulative, independent; lumping

The database holds, per yield set (`U235_thermal` at 0.0253 eV, `U235_fast` at 500 keV), the
independent yield IY (ENDF MT454) and the cumulative yield CY (MT459). `yield_type` selects:

| `yield_type` | direct source $Y_i$ | use |
|---|---|---|
| `effective` (default) | $\max\bigl(\mathrm{CY}_i - \sum_{p\ \text{tracked}} b_{pi}\,\mathrm{CY}_p,\ \mathrm{IY}_i\bigr)$ | chains: no double counting, untracked precursors lumped |
| `cumulative` | $\mathrm{CY}_i$ | a single nuclide in saturation, no tracked parents |
| `independent` | $\mathrm{IY}_i$ | every precursor tracked explicitly |

Since $\mathrm{CY}_i = \mathrm{IY}_i + \sum_p b_{pi}\mathrm{CY}_p$, the effective yield is the part of
the cumulative yield that does **not** pass through a tracked parent. Short-lived untracked
precursors (e.g. Te-135 for I-135) are therefore **lumped** into their first tracked descendant and
produced instantly. Untracked isomers are lumped in the database itself (I-135 -> Xe-135m -> Xe-135
becomes I-135 -> Xe-135; `provenance.lumping`). Examples (U235_thermal): I-135 CY 0.0628, Xe-135 CY
0.0654, I-131 CY 0.0289, Cs-137 CY 0.0619.

### 3.3 Decay branchings and capture

Each nuclide entry lists its decay branches (`daughter`, `branching`) and a capture
(`cross_section`, `product`). Links between **tracked** nuclides are added automatically: a decay
feed $\lambda_p b_{pi}$ and, if $\sigma_p>0$ and the product is tracked, a capture feed
$\sigma_p\phi$. Untracked daughters are simply lost from the system. Examples: I-135 -> Xe-135
($b$ = 0.9995) and -> Cs-135 (0.0005); Mo-99 -> Tc-99m (0.879) and -> Tc-99 (0.121);
$\sigma(\mathrm{Xe\text{-}135}) = 2.6476\times10^6$ b, $\sigma(\mathrm{Sm\text{-}149}) = 4.01\times10^4$ b.
`neutron_flux` is a one-group functor (default 0).

### 3.4 Removal pathways: one rate per (group, destination)

Each nuclide belongs to one chemical group:

| group | nuclides in the database | behaviour |
|---|---|---|
| `noble_gas` | Kr-85, Xe-133, Xe-135 | insoluble, to bubbles and the cover gas |
| `noble_metal` | Nb-95, Mo-95, Mo-99, Tc-99m, Tc-99, Ru-103, Rh-103, Ru-106, Te-132 | insoluble metals: plate on walls, entrained in the off-gas |
| `halogen` | I-131, I-132, I-133, I-135 | stay in the salt as iodide |
| `salt_seeking` | Sr-90, Y-90, Zr-95, Cs-135, Cs-137, Pm-149, Sm-149 | stable fluorides/chlorides |

A removal pathway is a triple of equal-length lists entries:
`removal_groups` x `removal_rates` (functors, 1/s) x `removal_destinations`
(`offgas`, `plated`, `processed`). A group may appear several times with **different** destinations.
**A duplicate (group, destination) pair is rejected** with the error "The group ... is removed to ...
more than once. Combine the rates of such pathways into one functor." (constructor of
`FissionProductInventoryAction`). This is why case 04 sums spray and bubble stripping into one
`offgas` rate (`r_xe_offgas = r_spray + r_bubble`) and sends the graphite uptake to `processed`.

### 3.5 Removed inventories; wall plating

With `track_removed_inventories = true` (default), the action creates `<nuclide>_offgas` and
`<nuclide>_plated` inventories, per unit **salt** volume, for every nuclide removed to that
destination **and, recursively, its tracked daughters**:

$$\frac{dN_i^{d}}{dt} = r_{g(i),d}\,N_i + \sum_{p}\lambda_p b_{pi}\,N_p^{d} - \lambda_i N_i^{d},\qquad d\in\{\text{offgas},\text{plated}\}.$$

They are out of the flux (no capture), not transported, and they keep decaying, so a plated Te-132
feeds a plated I-132 if it is tracked. `processed` removes the inventory without tracking it; a user
can rebuild it by hand (case 04 does this for the graphite: `Xe135_graphite` with a
`FissionProductFeed` from the salt and a `FissionProductLoss` with decay and burnout in the core
flux).

**Wall plating.** With `plating_boundaries` and `plating_mass_transfer_coefficient` (and a
`diffusivity`, which is required because the flux is a boundary term of the diffusion operator), the
action adds an `MSRDegassingBC` with $k_H = 0$ on every `noble_metal` salt inventory:

$$J = k_m\,N\quad\text{(perfect sink at the wall)}.$$

Its 0D equivalent is a volumetric rate $r = k_m A/V$, which is how cases 04 (`r_metal = k_metal
A_metal/V`) and 05 (HX plates: `r_hx_plate = k_m,hx a_hx`) enter plating as a removal pathway. Note:
the boundary-plated inventory is **not** accumulated into the `_plated` variables, which only collect
the volumetric `plated` pathways (case 05 sets `track_removed_inventories = false` and reports the
wall flux with a postprocessor).

For a well-mixed loop the salt inventory then follows the 0D balance used by the verification of
cases 04 and 05:

$$\bar N = \frac{Y\dot F/N_A}{\lambda + \sigma\phi + \overline{r} + k_m A_\text{wall}/V}.$$

### 3.6 Ordering

The action sorts the inventories topologically, parent before daughter (stable with respect to the
input order). Because the loss is implicit and the feed uses the latest parent value, one sweep
in this order makes each implicit-Euler step exact for a decay chain. A cycle through capture
produces a warning and keeps the input order (several fixed-point sweeps are then needed).

### 3.7 Xenon: equilibrium, stripping and worth

With I-135 in equilibrium (Te-135 lumped into it),

$$N_I = \frac{Y_I\dot F}{N_A(\lambda_I + \sigma_I\phi)},\qquad
N_{Xe} = \frac{Y_{Xe}\dot F/N_A + b\,\lambda_I N_I}{\lambda_{Xe} + \sigma_{Xe}\phi + r_\text{strip}}.$$

The xenon worth used in the cases is the one-group poison estimate, **without importance
weighting**:

$$\rho_{Xe} \approx -\frac{\sigma_{Xe}\,\phi\,N_{Xe}N_A}{\nu\,\dot F}\qquad(\nu = 2.43).$$

- Case 02 (`rho_Xe_pcm`) evaluates it at the full-power flux $\phi_0$ and power $P_0$, i.e. as an
  inventory indicator during the transient, not a kinetics calculation.
- Case 04 (`xe_poison_pct`) adds the graphite xenon burning in the mean core flux
  $\phi_c = \phi\,(V/V_\text{core})\,0.88$: $\rho = \sigma(\phi N_\text{salt} + \phi_c N_\text{graphite})/(\nu\dot F)$.

**Caveats.** (i) The ratio $\sigma\phi N/(\nu\dot F)$ is the absorption by Xe per fission neutron;
the reactivity worth of a real core also depends on the spatial importance, the thermal utilization
and the spectrum. (ii) $\sigma$ is the 2200 m/s value; in an epithermal/fast spectrum (MSFR) the
effective one-group cross section is much smaller. (iii) The worth formula keeps the flux at its
full-power value in case 02, so it is not the reactivity during the shutdown itself. (iv) In the MSRE
the circulating helium bubbles dominate the stripping: case 04 overpredicts the measured
0.26-0.28 % by 4x without them and needs a **calibrated** bubble rate $r_b$ = 2.14e-3 1/s.

### 3.8 Burnup monitoring

Cs-137 ($t_{1/2}$ = 30 y) integrates the energy: without removal,

$$N_{Cs}(t) \approx \frac{\mathrm{CY}_{Cs}}{E_f N_A}\int_0^t q'''\,dt' = \frac{\mathrm{CY}_{Cs}\,E(t)}{E_f N_A},$$

to within its decay (0.13 % over a month). Case 02 checks it (`Cs137_burnup_check` = 0.9987,
**verification** of the energy integral and of the decay).

### 3.9 Database provenance and limitations

`data/fission_product_database.json` (23 nuclides, generated by
`scripts/fission_products/build_fission_product_database.py`):

| data | source (as recorded in `provenance`) |
|---|---|
| half-lives, branchings | ENDF/B-VIII.0 decay sub-library, MF8/MT457 |
| yields | ENDF/B-VIII.0 neutron fission yields, U-235, MT454 (independent), MT459 (cumulative), 0.0253 eV and 500 keV |
| capture | 2200 m/s capture cross sections from the MF1/MT451 descriptions of the ENDF/B-VIII.0 evaluations (Mughabghab, *Atlas of Neutron Resonances*) |
| molar masses | AWR of the MF8 head record times the neutron mass (generator script) |

**[code vs doc]** The database and its generator cite only ENDF/B-VIII.0 (and Mughabghab); NUBASE2020
is not referenced anywhere in the repository. If NUBASE2020 was intended as the source of the
half-lives, the database does not say so.

Limitations:

- **Xe-137 is not tracked.** Most of the mass-137 chain passes through Xe-137 (3.8 min), a noble gas.
  It is lumped into the Cs-137 cumulative yield and is never sparged. In a sparged salt a fraction
  $r/(\lambda_{Xe137} + r)$ (3 % at $r$ = 1e-4 1/s) of the Cs-137 precursors would go to the off-gas;
  case 02's burnup check would read ~0.967 instead of 0.9987.
- **One-group, 2200 m/s cross sections**, no spectrum, no resonance self-shielding; the flux is a
  given functor, not coupled to the inventory.
- **U-235 yields only** (thermal and 500 keV fast); case 05 uses `U235_fast` as a stand-in for U-233.
- **Noble metals are a perfect mass-transfer sink** (no sticking probability, colloids, or bubble
  attachment); case 04 shows that this overpredicts the metal deposits of the MSRE.
- **Chemistry is group-wise**: there is no redox dependence of the noble-metal behaviour (e.g.
  Nb, Te), and the halogens and salt seekers have no removal unless one is given.
- **23 nuclides**: chains not in the database cannot be tracked.

## 4. Corrosion and plating

### 4.1 (a) The Butler-Volmer law

All the mechanistic corrosion objects use the same kernel, `Corrosion::bvCurrent`
(`include/corrosion/base/ButlerVolmerKinetics.h`), for the reaction M(s) $\rightleftharpoons$ M$^{z+}$ + $z$ e$^-$:

$$\eta = (\phi_\text{metal} - \phi_\text{salt}) - E_0,\qquad
i_{BV} = i_0\Bigl[a_s\,e^{\alpha_a z f\eta} - \frac{\max(c,0)}{c_\text{ref}}\,e^{-\alpha_c z f\eta}\Bigr],\qquad
f = \frac{F}{RT}.$$

Each exponent is clipped to $[-60, 60]$. $a_s$ is the metal activity (default 1, floored at 0). A
positive current dissolves the metal (corrosion), a negative one plates it. The molar flux into the
salt is Faraday's law $J = i_{BV}/(zF)$, and the recession (or plating) rate of the metal is

$$v = \frac{M\,i_{BV}}{zF\rho_\text{alloy}}\quad(\texttt{corrosionCurrentToUmY}:\ \text{A/cm}^2 \to \mu\text{m/y}).$$

`CorrosionRateAux` outputs `current`, `penetration_rate` [µm/y] or `recession` [µm], the last
accumulated explicitly as $u^{n} = u^{n-1} + v\,\Delta t$. Database properties (Cr, Fe, Ni:
$z$ = 2, $\alpha_a = \alpha_c$ = 0.5, $E_0$ = 0, $c_\text{ref}$ = 1 mol/m$^3$, $D$ = 1e-9 m$^2$/s).

**Forms of the law, by topology** (`[CorrosionPlating] topology`):

| form | object | what is modelled | $\phi_\text{metal}$, $\phi_\text{salt}$, $c$, $a_s$ |
|---|---|---|---|
| `salt_only` | `ButlerVolmerBC` (species flux on the salt variable) | salt block, wall as a boundary | metal = `applied_potential` (default `applied_overpotential` + pin), salt = 0 or the solved potential; $c$ = the salt variable; $a_s$ = 1 |
| `solid_only` | `ButlerVolmerBC` with `metal_domain = true` (sink on the solid variable) | alloy block, salt as a boundary | salt = `applied_potential` (default $-\eta$); $c$ = `counter_concentration` (required) |
| `two_block` | `ButlerVolmerInterface` | salt and alloy blocks, interface sideset | metal = `applied_overpotential` (+ pin), salt = 0 or the solved potential; salt gains, solid loses $J$ |
| `two_block` + `solid_activity = true` | `ButlerVolmerInterface` with `solid_reference_concentration` | as above, with surface depletion | $a_s = c_s/c_{s,0}$ ($c_{s,0}$ = initial solid concentration) |
| potential-solved | + `CorrosionNernstPlanckFlux`, `CorrosionCurrentContinuity` | `solve_potential = true`, `pin_potential_boundary` | $\phi_\text{salt}$ is a variable |

In the finite-element objects the temperature of the BV exponent is the scalar
`reference_temperature` (the temperature functor is used by the transport kernels only).

**Transport in the FE action.** Salt cations obey the Nernst-Planck flux (`CorrosionNernstPlanckFlux`
plus `CorrosionAdvection` when velocities are given):

$$\frac{\partial c_k}{\partial t} + \nabla\cdot\mathbf N_k = 0,\qquad
\mathbf N_k = -D_k\nabla c_k - z_k f D_k c_k\nabla\phi + \mathbf u c_k .$$

The migration term is present only when the potential is solved. The alloy-side metal diffuses with
`CorrosionSolidDiffusion`, $\partial c_s/\partial t = \nabla\cdot(D_s\nabla c_s)$, with $D_s$ from the
calibrated correlation (section 4.2) unless `solid_diffusivity` is given.

**Potential-solved mode.** With `solve_potential = true`, `CorrosionCurrentContinuity` solves
$\nabla\cdot\mathbf i = 0$ for $\phi_\text{salt}$, with

$$\mathbf i = -F\sum_k z_k D_k\nabla c_k - \kappa\nabla\phi,\qquad \kappa = \frac{F^2}{RT}\sum_k z_k^2 D_k\max(c_k,0).$$

The constant mode is pinned with an `ADDirichletBC` (`pin_potential_boundary`, `pin_potential_value`);
the default metal potential carries the same offset so that $\phi_\text{metal} - \phi_\text{salt} = \eta$
at the reference state. With `supporting_electrolyte = true` (default) the BV current is *not* fed back
into the potential equation (the corroding cations are dilute tracers in an inert electrolyte); with
`false` and `two_block`, a charge-flux `ButlerVolmerBC` (`flux_type = charge`) adds $i_{BV}$ to the
salt charge balance. **[code vs doc]** `ButlerVolmerPotentialInterface` (documented as the
interface coupling of the current into both potential equations) exists but is not created by the
action. None of the five cases solves the potential.

**Linear-FV wall law** (`CorrosionLinearFVButlerVolmerBC`, built by `[CorrosionPlatingFlow]`). With a
prescribed $\eta$ the inward flux is linear in the wall-cell concentration $c_P$:

$$J_\text{in} = A - B\,c_P,\qquad A = \frac{i_0(T)}{zF}e^{\alpha_a z f\eta},\qquad
B = \frac{i_0(T)}{zF\,c_\text{ref}}e^{-\alpha_c z f\eta},$$

assembled without lagging ($B$ on the matrix, $A$ on the right-hand side). Here the temperature is the
local functor value at the face, and the exchange current is Arrhenius:

$$i_0(T) = i_0\exp\!\left[\frac{E_a}{R}\left(\frac{1}{T_\text{ref}} - \frac{1}{T}\right)\right],$$

with $E_a$ = `Ea_corr_kJ_mol` = 54.84 kJ/mol from the corrosion database when
`temperature_dependent_kinetics = true` (default), else 0. A diffusion kernel must be present on the
variable for the wall flux to be assembled.

**The limiting current.** There is **no explicit limiting-current term** in the BV objects
**[code vs doc]**. Transport limitation enters in three ways:
(i) the resolved diffusion layer in the FE salt block (the wall concentration rises until
$J$ equals what diffusion can carry away; case 04 sets `salt_diffusivity` $= k_m L/2$ on its
one-element salt column so that the wall-to-centre resistance is $1/k_m$); (ii) in the linear-FV BC,
the wall-cell value $c_P$ is coupled to the bulk only through the diffusion operator, so the
wall-to-cell resistance acts as a film resistance;
(iii) the reduced correlation used for the seed is itself a harmonic mean of a kinetic and a
mass-transfer-capacity branch (section 4.2), so the seeded $i_0$ already includes the calibrated
transport cap.

### 4.2 (b) Seeding $i_0$ from the reduced correlation

$i_0$ is not an input. The actions compute it so that the BV law reproduces, on fresh salt at the
reference state, the dissolution rate of the calibrated reduced correlation
`MoltenSaltCorrosionModel::corrosionRateUmY` (`data/corrosion_database.json`, `calibrated_parameters`;
a C++ port of the reference effective model, fitted to 43 measurement rows). Its structure:

$$\ln r_\text{kin} = \ln r_0 + \frac{E_\text{corr}}{R}\Bigl(\frac{1}{T_0} - \frac{1}{T}\Bigr)
+ \Delta_\text{mat} + \Delta_\text{salt} + \Delta_\text{redox} + \Delta_\text{pos}
+ \gamma_c\ln f_\text{flow} + \theta_c\ln\Bigl(1 + \frac{\Delta T}{100}\Bigr),$$

$$\ln r_\text{mt} = \ln r_\text{mt,0} + \frac{E_\text{mt}}{R}\Bigl(\frac{1}{T_0} - \frac{1}{T}\Bigr)
+ 0.35\,\Delta_\text{salt} + \gamma_m\ln f_\text{flow} + \theta_m\ln\Bigl(1 + \frac{\Delta T}{100}\Bigr),$$

$$r_\text{corr} = \Bigl(\frac{1}{r_\text{kin}} + \frac{1}{r_\text{mt}}\Bigr)^{-1}\quad[\mu\mathrm m/\mathrm y],\qquad T_0 = 923.15\ \mathrm K .$$

The offsets $\Delta$ are fitted class constants (`material_class`, `salt_class`, `redox_class`,
`position_class`; Hastelloy N and `fluoride_fuel` are the zero references); $f_\text{flow}$ =
`flow_factor` (floored at 1e-3, default 0.75 inside the model when 0) and $\Delta T$ = `delta_T_C`.
Selected values: $\ln r_0$ = 0.265, $E_\text{corr}$ = 54.84 kJ/mol, $\ln r_\text{mt,0}$ = 8.11,
$E_\text{mt}$ = 27.35 kJ/mol, $\gamma_c$ = 0.15, $\theta_c$ = 0.05, $\gamma_m$ = 0.603,
$\theta_m$ = 0.604; redox offsets `purified_baseline` -0.385, `msre_or_fuel_baseline` 0.731,
`oxidizing_fef2` 1.503. `no_salt` or $f_\text{flow} < 0.05$ returns 1e-6 µm/y.

The seed is

$$i_\text{ref} = \frac{r\,zF\rho_\text{alloy}}{M}\ (\text{unit-converted}),\qquad
i_0 = \frac{i_\text{ref}}{\mathcal B},\qquad
\mathcal B = e^{\alpha_a z f_\text{ref}(\eta_\text{app} - E_0)} - e^{-\alpha_c z f_\text{ref}(\eta_\text{app} - E_0)},$$

where $r$ is the planned rate of the element and $\mathcal B$ (`butlerVolmerSeedBracket`) is the BV
bracket at $c = c_\text{ref}$, unit activity and $T_\text{ref}$ = `reference_temperature`. The
bracket must be positive (`applied_overpotential` > $E_0$), else the action stops. Without
`seed_fractions`, **every** tracked element is seeded with the full rate (legacy behaviour; the total
then counts the rate once per element). The correlation also supplies the solid Cr diffusivity
$\ln D_s = \ln D_{s,0} + (E_D/R)(1/T_0 - 1/T)$ ($E_D$ = 240 kJ/mol).

**Interpretation.** The redox offset can be read as an effective overpotential,
$\Delta_\text{redox} = \alpha z F\eta_\text{eff}/(RT)$ (`bvOverpotentialEquivalentV`), but
`validation/corrosion/calibration/MODEL.md` warns that $\eta_\text{eff}$ is not a measured electrode
potential. The fit metrics (median factor error 1.11 on direct targets, 93 % within a factor 2) are
in-sample.

### 4.3 (c) The Nernst limit, and why a constant $\eta$ is wrong for a loop

At constant $\eta$ the current vanishes at

$$c_\text{eq} = c_\text{ref}\,\exp\bigl[(\alpha_a + \alpha_c)\,zf\eta\bigr] = c_\text{ref}\exp\!\Bigl[\frac{(\alpha_a+\alpha_c)zF\eta}{RT}\Bigr].$$

A closed isothermal system corrodes until $c = c_\text{eq}$, then stops. In case 03 (capsule) the
applied $\eta$ = 0.170 V gives $c_\text{eq}$ = 900 ppm Cr at 704 C; the capsule is half saturated after
31 d and stops after a total recession of $c_\text{eq}(V/A)M/\rho_\text{alloy}$ = 0.13 µm, against the
measured years-long loss of NCL-16 (model-form check).

With $\eta > 0$ constant, $c_\text{eq}$ **falls** as $T$ rises (the exponent scales with $1/T$). A
closed loop with one $\eta$ would therefore plate at its hot walls and dissolve at its cold walls,
the opposite of what is observed. The flow action's equilibrium law (4.4) fixes this.

### 4.4 (d) The element-resolved, temperature-dependent equilibrium law (`[CorrosionPlatingFlow]`)

Opt-in with `equilibrium_concentrations` (one positive value per tracked element, at
`reference_temperature`), optionally `dissolution_enthalpies` (one value for all elements, or one per
element, kJ/mol) and `seed_fractions`. The action passes `equilibrium_concentration` and
`dissolution_enthalpy` to each element's `CorrosionLinearFVButlerVolmerBC`, which then replaces the
potential functors by

$$c_{\text{eq},k}(T) = c_{\text{eq},k}\exp\!\Bigl[-\frac{\Delta H_k}{R}\Bigl(\frac{1}{T} - \frac{1}{T_\text{ref}}\Bigr)\Bigr],\qquad
\eta_k(T) = \frac{\ln\bigl(c_{\text{eq},k}(T)/c_\text{ref}\bigr)}{(\alpha_a + \alpha_c)\,z f}.$$

This $\eta_k$ is used directly in the exponents (no $E_0$ subtraction), so the net flux vanishes at
$c = c_{\text{eq},k}(T)$. A positive $\Delta H$ makes the equilibrium rise with $T$: **hot walls
dissolve, cold walls plate**, and each element can be given its own equilibrium (a noble Ni barely
dissolves). Giving the potential functors as well is an error; `dissolution_enthalpy` without a
positive equilibrium concentration is an error.

Inserting $\eta_k$ in $A$ and $B$:

$$J_k = \frac{i_{0,k}(T)}{zF}\Bigl(\frac{c_{\text{eq},k}(T)}{c_\text{ref}}\Bigr)^{\frac{\alpha_a}{\alpha_a+\alpha_c}}
\Bigl[1 - \frac{c_k}{c_{\text{eq},k}(T)}\Bigr]
\ \xrightarrow{\ \alpha_a=\alpha_c=1/2,\ c_\text{ref}=1\ }\
\frac{i_{0,k}(T)}{zF\sqrt{c_{\text{eq},k}(T)}}\bigl(c_{\text{eq},k}(T) - c_k\bigr).$$

**Seeding at $c = 0$.** In this mode the action seeds $i_0$ with the *fresh-salt anodic branch* at the
reference equilibrium, i.e. the bracket is $\mathcal B = (c_{\text{eq},k}/c_\text{ref})^{\alpha_a/(\alpha_a+\alpha_c)}$
(exponent clipped at 60) instead of the $c = c_\text{ref}$ bracket of 4.2. Hence, at
$T = T_\text{ref}$ and $c = 0$, $J_k$ reproduces exactly the planned (calibrated) rate of element $k$,
and in general (case 03 header)

$$J_k = J_{k,\text{ref}}\;e^{\frac{E_a}{R}\left(\frac{1}{T_\text{ref}} - \frac{1}{T}\right)}
\sqrt{\frac{c_{\text{eq},k}(T)}{c_{\text{eq},k}}}\Bigl(1 - \frac{c_k}{c_{\text{eq},k}(T)}\Bigr),\qquad
J_{k,\text{ref}} = \frac{r_k\,\rho_\text{alloy}}{M_k}.$$

**`seed_fractions`.** Nonnegative weights $w_k$, one per element, normalized to sum to 1, split the
calibrated rate: $r_k = r_\text{corr}\,w_k/\sum_j w_j$, so the total metal loss at the seed state
equals the calibrated rate. `'1 0 0'` is selective Cr dissolution; the alloy mass fractions give
congruent dissolution. Empty (default) seeds every element with the full rate. `seed_fractions` is
rejected in MSTDB-TC mode (which has its own source fractions).

**Training values.** Case 03 (flow wall law): $c_\text{eq}$ = 900 / 40 / 1 ppm (Cr / Fe / Ni) at
704 C, $\Delta H$ = 100 kJ/mol, seeds 0.90 / 0.07 / 0.03; $c_{\text{eq,Cr}}$ and $\Delta H_\text{Cr}$ are a
two-parameter **calibration** on M-014 and M-012. Case 05: 4 / 0.05 / 1e-4 mol/m$^3$ at 973.15 K,
$\Delta H$ = 50 kJ/mol, seeds 0.90 / 0.08 / 0.02; **illustrative assumptions**. In a closed loop the
salt settles where the area-weighted wall fluxes cancel, between the hot- and cold-wall equilibria,
and then carries metal from hot to cold at a steady rate. The MSTDB-TC CrF2 saturation activities
(0.704 at 704 C, 0.483 at 538 C) correspond to only $\Delta H \approx$ 15 kJ/mol (case 03, exercise 2).

**Limitations.** The equilibria and enthalpies are user assumptions, not computed from
thermochemistry; there is no coupling of $c_\text{eq}$ to the salt redox (U(III)/U(IV)) or to the
alloy surface depletion; the noble element equilibria are free parameters; there is no wall
comproportionation or film mass-transfer model beyond the cell resolution.

### 4.5 (e) The MSTDB-TC standard-state model

`MSTDBTCStandardStateCorrosionModel` (`kinetics_model = mstdb_tc_standard_state` in both corrosion
actions, or `MSTDBTCCorrosionUserObject`) is a **lumped endpoint correlation** for Cr/Fe/Ni whose
driving forces come from MSTDB-TC standard-state Gibbs energies. It is not a native
Thermochimica/SUBQ equilibrium calculation. Activities, charge transfer, transport, geometry and cold
capture are fitted closures (`thermochemical_parameters` of `data/advanced_corrosion_models.json`).

**Activities.** Dissolved product (mole fraction):
$a_{MX_2} = \max\bigl(10^{-6}\,\mathrm{ppm}\,\bar M_\text{salt}/M_{MX_2},\ 10^{-30}\bigr)$ (product
floor $e^{\texttt{log\_product\_floor\_ppm}}$ ppm when none is given). Metal: the mole fraction of the
element in the alloy, from the Cr/Fe/Ni mass fractions (remainder counted with M = 60 g/mol).

**Affinities** ($\mathcal A_k = \ln(K/Q)$, $\ln K = -\Delta_r G^\circ/(RT)$ from the MSTDB Gibbs functions):

| salt / redox | reaction | $\mathcal A_k$ |
|---|---|---|
| `fluoride_fuel` | M + 2 UF$_4$ -> MF$_2$ + 2 UF$_3$ | $\ln K - \ln a_{MF_2} + \ln a_M + 2\ln(\mathrm{U^{IV}/U^{III}})$, with U$^\text{IV}$/U$^\text{III}$ $= 100\,e^{0.5\,s}$ |
| `flibe` + `reducing_be` | M + BeF$_2$ -> MF$_2$ + Be | $\ln K - \ln a_{MF_2} + \ln a_M$ |
| other salts | M + FeX$_2$ -> MX$_2$ + Fe | $\ln K - \ln a_{MX_2} - \ln a_\text{Fe} + \ln a_M + \ln a_{\text{FeX}_2}$; **0 for Fe** (identity reaction) |

$s$ = `redoxLogShift` = $\Delta_\text{redox}(\text{class}) - \Delta_\text{redox}(\text{purified})$ from the
reduced correlation; the Fe(II) buffer is 500 ppm for `oxidizing_fef2`, otherwise
100 (30 FLiBe, 250 chloride) $\times e^{s}$ ppm.

**Source fractions** (how the front rate is shared among Cr, Fe, Ni):

$$\xi_k = \frac{\exp\bigl[\ln x_k + o_k + \sigma\,\mathrm{clip}(\mathcal A_k - \mathcal A_\text{Cr}, \pm 80)\bigr]}
{\sum_j \exp[\ldots]_j},\qquad o = (0,\ \texttt{log\_exchange\_fe\_relative},\ \texttt{log\_exchange\_ni\_relative}).$$

**Front rate:** $r_f = (1/r_\text{ct} + 1/r_\text{mt})^{-1}$ with

$$\ln r_\text{ct} = \ln r_{f,0} + \frac{E_\text{corr}}{R}\Bigl(\frac{1}{T_0} - \frac{1}{T_h}\Bigr)
+ \beta\ln\frac{x_\text{Cr}}{0.07} + \Delta_\text{salt}^\text{thermo} + s,\qquad
\ln r_\text{mt} = \ln r_{\text{mt},0} + \frac{E_\text{mt}}{R}\Bigl(\frac{1}{T_0} - \frac{1}{T_h}\Bigr) + \gamma\ln f_\text{flow},$$

where $\Delta_\text{salt}^\text{thermo}$ is the Cr affinity of the salt relative to `fluoride_fuel` minus a
fitted activity correction (0 for the fuel salt).

**Saturation activities and cold capture:**

$$a_{\text{sat},k}(T) = \mathrm{clip}\Bigl(\exp\Bigl[\mathrm{clip}\Bigl(\frac{G^\circ_{\text{solid}} - G^\circ_{\text{liquid}}}{RT}, -80, 0\Bigr)\Bigr], 10^{-30}, 1\Bigr),\qquad
\chi_k = \mathrm{clip}\Bigl(1 - \frac{a_{\text{sat},k}(T_c)}{a_{\text{sat},k}(T_h)}, 0, 1\Bigr).$$

Endpoint outputs: front depth $d = r_f t$; corrosion rate $= r_f \times$ (mass-loss fraction,
a fitted logistic clipped to [0.03, 1]); mass gain $= d\,\rho\,0.1 \times \kappa_\text{cap}\sum_k \xi_k\chi_k$
[mg/cm$^2$]; dissolved inventory
$\mathrm{ppm}_k \propto d\,\rho\,\xi_k\,(A/m)\,\varepsilon\,\max(0, 1 - \kappa_\text{cap}\chi_k)$ with the
`inventory_coupling_factor` $\varepsilon$ applied to the dissolved inventory only, then saturated as
$C_\text{cap}(1 - e^{-\text{raw}/C_\text{cap}})$; IGC depth = front + grain-boundary diffusion length
$2\sqrt{D_s t/\pi}$ times fitted multipliers.

**Use in the actions.** The total front rate $r_f$ is apportioned as $r_k = r_f\,\xi_k$ and seeded
as in 4.2; `elements` must be exactly `Cr Fe Ni`; `reference_cold_temperature`,
`reference_exposure_time`, `area_to_salt_mass` and `inventory_coupling_factor` are required.
The MSTDB-TC model is not re-evaluated locally: with `temperature_dependent_kinetics` the local
temperature enters only through the Arrhenius $i_0(T)$ of the base database.

**Licensed data required.** MSTDB-TC is not redistributed. The action needs the exact
`MSTDB-TC_V3.1_Fluorides_No_Func.dat` and `..._Chlorides_No_Func.dat` files (`fluoride_database`,
`chloride_database`), whose version (3.1) and SHA-256 are checked against `advanced_database`;
extrapolation beyond the last Gibbs interval is rejected unless `allow_extrapolation = true`. In the
training (case 03) the MSTDB-TC NCL-16 endpoint is therefore only **read** from the calibration
results, not run.

### 4.6 (f) DRIDN, the Dynamic Redox-Inventory-Depletion Network

`DRIDNModel` (`include/corrosion/base/DRIDNModel.h`), wrapped by `DRIDNCorrosionUserObject` and read by
`AdvancedCorrosionModelPostprocessor`, is a **well-mixed 0D** network of 19 scalar states integrated in
years (adaptive RK4 with step doubling, inside each MOOSE time step; MOOSE time is in seconds). The
user object executes exactly once per step at `TIMESTEP_END`; the input sets `[Problem] solve = false`.

**States** (Cr, Fe, Ni for each triplet):

| index | state | unit |
|---|---|---|
| 0-2 | surface availability $s_k \in [0.015, 1]$ | - |
| 3-5 | dissolved inventory $C_k$ | ppm |
| 6-8 | cumulative source | ppm |
| 9-11 | coupon (cold-leg) deposit | mg/cm$^2$ |
| 12 | reaction-front depth $d_f$ | µm |
| 13 | mass-equivalent recession | µm |
| 14 | grain-boundary depth squared | µm$^2$ |
| 15 | redox log shift $\varrho$ | - |
| 16-18 | bulk-captured inventory | ppm |

**Main equations** (`ratesImpl`, `rhsImpl`; $C_\text{floor}$ = `product_activity_floor`,
$\tilde C_k = \max(C_k, C_\text{floor})$):

- affinity: $\mathcal A_k = \mathcal A_k^0 + \varrho - \ln(\tilde C_k/C_\text{floor})$;
- species fractions: $\xi_k \propto x_k\,e^{o_k + \sigma\,\mathrm{clip}(\mathcal A_k - \mathcal A_\text{Cr}) + p\ln s_k}$ (plus a Cr bias), normalized;
- charge-transfer rate (clipped at $\pm 60$ in log):

  $$\ln r_\text{ct} = \ln r_\text{ct}^0 + \varrho + \ln k_\text{scale} + \ln\Phi_\text{int} + a_f\Bigl(-\sum_k\xi_k\ln\frac{\tilde C_k}{C_\text{floor}}\Bigr) + p\sum_k\xi_k\ln s_k ;$$

  ($r_\text{ct}^0$ = `log_charge_base_no_redox`, $\Phi_\text{int}$ = stress / fluoride-impurity
  interfacial factors, $a_f$ = `affinity_feedback_scale`, $p$ = `surface_availability_exponent`); in `legacy_irreversible` mode (default) this is the rate; in `affinity_gated_dissolution` it is
  multiplied by $1 - e^{-\bar{\mathcal A}}$ for $\bar{\mathcal A} = \sum\xi_k\mathcal A_k > 0$ and is 0 otherwise;
- mass transfer with inventory inhibition: $r_\text{mt} = r_\text{mt}^0\,\Phi_\text{int}/(1 + \beta\,\sum_k C_k/C_\text{cap})$;
- front rate $r_f = (1/r_\text{ct} + 1/r_\text{mt})^{-1}$;
- surface: $\dot s_k = \kappa\,\rho_k(1 - s_k) - r_f\,\xi_k/L_\text{res}$ ($L_\text{res}$ = reservoir depth, $\propto x_\text{Cr}$ ratio);
- dissolved inventory: $\dot C_k = S_k - (k_{d,k} + k_{b,k} + k_p)\,C_k$ with source
  $S_k = r_f\,\rho\,100\,(A/m)\,\varepsilon\,\xi_k$ [ppm/y];
- cold-leg deposition $k_{d,k} = k_d^0\,\chi_k\sqrt{\max(f_\text{flow}, 0.05)}$ (Fe multiplied by a fitted bias);
  coupon deposit $\dot m_k = k_{d,k}C_k\times$(ppm -> mg/cm$^2$ conversion $10^{-3}/(A/m)\times$ area factors);
- **bulk capture** $k_{b,k} = k_{d,k}\times$ (area multiplier) and **bulk precipitation**
  $k_p = \min(200,\ 0.35\,e^{\ln k_p'}\,(\sum C/C_\text{cap})^2)$ [1/y], both accumulated in the
  bulk-captured states;
- front depth $\dot d_f = r_f$; recession $= r_f\times$ mass-loss fraction; grain-boundary depth
  $\dot{(\delta^2)} = (4/\pi)\,D_s\,g^2\,e^{a_d\max(\varrho,0)}$;
- optional transient redox: $\dot\varrho = e^{\ln k_\varrho}(\theta\varrho_0 - \varrho) - \mathrm{sgn}(\varrho)\,c_\varrho r_f$.

The reported `mass_balance_relative_error` checks dissolved + deposits + bulk capture = cumulative
source per element; it is not a proof of alloy, charge or redox conservation.

**What is calibrated.** The 26 closure parameters (`dynamic_parameters` of
`data/advanced_corrosion_models.json`) were fitted, together with the thermochemical layer, to 32
supported constraints of the corrosion benchmark (relation-aware log RMSE 0.092, 96.9 % pass; in
sample). In case 03 the context is also calibrated: `affinity_baseline`, `log_exchange_offsets`,
`log_charge_base_no_redox`, `selectivity_scale`, `mass_transfer_rate`, `cold_capture_fraction`,
`deposit_area_factor`, `area_to_salt_mass` (7.77 cm$^2$/g, a shared fitted value, not the loop
geometry), `inventory_capacity`, `cr_diffusion`, `mass_loss_fraction`,
`grain_boundary_length_multiplier`. Only the alloy mass fractions, density and flow factor are design
data. NCL-16 targets M-011 to M-015 are fit targets, so case 03 is a **calibration**.

**Known shortcomings** (case 03): the salt Fe has the wrong sign (+106 vs -100 ppm measured, since the
salt starts with ~1 ppm Fe and can only gain), the 378 ppm dissolved Ni is a calibration artefact, and
most of the released metal goes to the empirical bulk-capture sink (72 % of the Cr), so loss minus
gain is not the dissolved inventory. The calibration and metrics are for `legacy_irreversible`, which
does not enforce detailed balance; the affinity-gated mode is not calibrated.

### 4.7 (g) Noble-metal plating versus corrosion-product plating

| | noble-metal fission products (Mo, Ru, Te, Nb, ...) | corrosion products (Cr, Fe, Ni) |
|---|---|---|
| state in the salt | reduced metal atoms / colloids, insoluble | dissolved cations (CrF$_2$, FeF$_2$, NiF$_2$) |
| driving force | none needed: any contact sticks (perfect sink) | supersaturation: $c > c_\text{eq}(T_\text{wall})$ |
| model | first-order removal $r$ or wall flux $J = k_m N$ (`MSRDegassingBC`, $k_H = 0$) | BV wall law $J = A - Bc$, cathodic when $c > c_\text{eq}$ |
| reversibility | irreversible (plated inventory only decays) | reversible: the same wall dissolves or plates depending on $T$ and $c$ |
| where | every wetted surface, rate set by mass transfer ($k_m$) | cold walls (with $\Delta H > 0$), hot walls dissolve |
| action | `[FissionProductInventory]` (`plated` destination, `plating_boundaries`) | `[CorrosionPlating]`, `[CorrosionPlatingFlow]` |
| cases | 04 (MSRE, against ORNL-4865), 05 (HX plates and walls) | 03, 04, 05 |

The reduced correlation also has a separate empirical deposition branch (`depositionRateUmY`, with
surface classes graphite / laminar / turbulent metal and a cold-leg bonus), but no training case and
no action use it to drive plating.

### 4.8 Comparison of the corrosion models

| model | equations | spatial? | plating | temperature dependence | calibrated on | limitations |
|---|---|---|---|---|---|---|
| reduced correlation (`MoltenSaltCorrosionModel`) | log-linear class offsets, harmonic mean of kinetic and mass-transfer branches | no (rate at one state) | separate empirical deposition branch | Arrhenius ($E_\text{corr}$, $E_\text{mt}$) | 43 measurement rows, 76 cases (in-sample) | not mechanistic; $\eta_\text{eff}$ is not a potential; class-bound |
| BV, constant $\eta$ (`[CorrosionPlating]`, `[CorrosionPlatingFlow]` default) | $i_0[a_s e^{\alpha_a zf\eta} - (c/c_\text{ref})e^{-\alpha_c zf\eta}]$ | yes (FE or linear FV) | cathodic branch when $c > c_\text{eq}$ | $i_0(T)$ Arrhenius (flow BC only); $c_\text{eq}$ falls with $T$ | $i_0$ seeded from the correlation | stalls at the Nernst limit in a closed system; wrong hot/cold sense in a loop |
| BV + potential solve | + Nernst-Planck and current continuity | yes (FE) | as above | as above | as above | potential pinned; supporting-electrolyte default; not used in the training |
| element-resolved equilibrium law (`[CorrosionPlatingFlow]` opt-in) | $J_k \propto (c_{\text{eq},k}(T) - c_k)/\sqrt{c_{\text{eq},k}(T)}$ | yes (linear FV, with flow) | cold walls ($\Delta H > 0$) | $c_\text{eq}(T)$ van 't Hoff + $i_0(T)$ | rate seeded at $c = 0$; $c_\text{eq}$, $\Delta H$ by hand (case 03: 2 parameters) | equilibria are assumptions; no redox or surface depletion coupling |
| MSTDB-TC standard state | Gibbs affinities, source fractions, saturation activities, cold capture | no (endpoint) | cold capture $1 - a_\text{sat}(T_c)/a_\text{sat}(T_h)$ | Gibbs functions + Arrhenius closures | 32 constraints (in-sample) | needs licensed MSTDB-TC V3.1; not an equilibrium solver; Cr/Fe/Ni only |
| DRIDN | 19-state ODE network | no (well mixed) | cold-leg capture + bulk sink | through the context (affinities, $D_s$) | 32 constraints + per-case context (in-sample) | Fe sign in NCL-16; bulk sink dominates; legacy law without detailed balance |

## 5. Usage

### 5.1 `[MoltenSaltRadiolysis]` (case 01)

```
[Problem]
  linear_sys_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys U_IV_sys U_III_sys F2_gas_sys'
[]
[MoltenSaltRadiolysis]
  salt_type = fluoride                 # kernel of data/msr_database.json: chloride | fluoride
  metals = 'U'                         # adds U_IV, U_III and their fluoride reactions
  gas_species = 'F2'                   # adds F2_gas and the Henry exchange with F2_diss
  temperature = ${T}                   # functor [K]
  dose_rate = ${dose_rate}             # functor [W/m^3]; 0 switches the radiolysis off
  g_value_species = 'F2_diss'          # override default G values ...
  g_value_overrides = '0'              # ... here: no direct F2 yield (charge balance)
  kLa = ${kLa}                         # 1/s
  volume_ratio = ${volume_ratio}       # V_liq / V_gas
  initial_condition_species = 'F_ion U_IV U_III'
  initial_condition_values = '${F0} ${U_IV_0} ${U_III_0}'   # mol/m^3; others start at 0
[]
```

In a flow (case 05) add `rhie_chow_user_object = 'rc'`, `diffusivity = D_eff`,
`time_derivative = false` and `block = ${fluid}`. Use `verbose = true` to print the species and the
system names to list. Fission oxidation and redox control are separate `LinearFVSource` kernels
(sections 2.8 and 5.6).

### 5.2 `[FissionProductInventory]` (case 02)

```
[Problem]
  linear_sys_names = 'I135_sys Xe135_sys I131_sys Cs137_sys Xe135_offgas_sys'
[]
[FissionProductInventory]
  nuclides = 'I135 Xe135 I131 Cs137'   # links between them are added automatically
  power_density = power_density        # W/m^3 functor (or fission_rate_density)
  neutron_flux = flux                  # one-group flux functor [1/(m^2 s)]
  removal_groups = 'noble_gas'         # one entry per pathway ...
  removal_rates = '1e-4'               # ... rate functor [1/s] ...
  removal_destinations = 'offgas'      # ... and destination: offgas | plated | processed
[]
```

Defaults: `yield_set = U235_thermal`, `yield_type = effective`, `energy_per_fission = 200`,
`track_removed_inventories = true`. Case 04 shows several pathways per group:

```
  removal_groups = 'noble_gas noble_gas noble_metal noble_metal noble_metal'
  removal_rates = '${r_xe_offgas} ${r_graphite} ${r_metal} ${r_mgraph} ${r_nm_offgas}'
  removal_destinations = 'offgas processed plated processed offgas'
```

Case 05 adds wall plating and transport:

```
  plating_boundaries = ${walls}
  plating_mass_transfer_coefficient = ${k_m_wall}   # m/s
  rhie_chow_user_object = 'rc'
  diffusivity = D_eff                               # required by plating_boundaries
  track_removed_inventories = false
  time_derivative = false
```

### 5.3 `[CorrosionPlating]` (finite element; case 03 capsule)

```
[CorrosionPlating]
  topology = salt_only                 # two_block | salt_only | solid_only
  reaction_boundary = left
  elements = 'Cr'                      # creates c_Cr (and cs_Cr for two_block/solid_only)
  temperature = ${T_hot}
  reference_temperature = ${T_hot}     # T of the seed and of the BV exponent
  material_class = hastelloy_n         # correlation features ...
  salt_class = fluoride_fuel
  redox_class = purified_baseline
  position_class = nominal
  flow_factor = 1.0
  delta_T_C = 166.0                    # ... used only to seed i0
  applied_overpotential = ${eta}       # V; also sets c_eq = c_ref exp(z f eta)
  default_salt_concentration = 0.0     # fresh salt
  reference_concentration = 1.0        # c_ref of i0
  transient = true
  verbose = true                       # prints the seeded i0
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
[]
```

Case 04 adds `applied_potential`, `salt_diffusivity` (film resistance) and starts from 40 ppm
(`default_salt_concentration = reference_concentration = c0`). The action creates the aux variables
`corrosion_rate_um_y` and `recession_um` used by the postprocessors.

### 5.4 `[CorrosionPlatingFlow]` (linear FV; case 03 flow wall law, case 05)

```
[CorrosionPlatingFlow]
  elements = 'Cr Fe Ni'
  release_variables = 'c_Cr c_Fe c_Ni'  # existing variables (declared in [Variables]); 'none' = create c_<El>
  reaction_boundary = 'left right'      # the walls
  temperature = T                       # functor: local wall temperature
  reference_temperature = ${T_hot}
  material_class = hastelloy_n
  salt_class = fluoride_fuel
  redox_class = purified_baseline
  position_class = nominal
  delta_T_C = 166
  temperature_dependent_kinetics = true # Arrhenius i0(T) with Ea_corr
  equilibrium_concentrations = '${c_eq_Cr} ${c_eq_Fe} ${c_eq_Ni}'   # mol/m^3 at T_ref
  dissolution_enthalpies = '${dH_Cr} ${dH_Fe} ${dH_Ni}'           # kJ/mol (one, or one per element)
  seed_fractions = '${seed_Cr} ${seed_Fe} ${seed_Ni}'             # split of the calibrated rate
  verbose = true
[]
```

When a variable is released into (`release_variables`), the action adds only the wall BC; the time,
advection and diffusion kernels of that variable are the user's (case 03) or the owning action's
(e.g. a radiolysis `Cr_II`). In case 05 the action owns `c_Cr`, `c_Fe`, `c_Ni` and adds
`rhie_chow_user_object`, `diffusivity = D_eff`, `time_derivative = false`, `initial_concentration`.

### 5.5 DRIDN (case 03)

```
[Problem]
  solve = false
[]
[UserObjects]
  [dridn]
    type = DRIDNCorrosionUserObject
    mass_fractions = '0.07 0.05 0.71'   # design: alloy Cr/Fe/Ni
    density = 8.89
    affinity_baseline = '18.11 8.84 0.34'   # calibrated context (full precision in the input)
    ...                                       # all context entries are required
  []
[]
[Postprocessors]
  [salt_Cr_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = dissolved_inventory_cr_ppm
  []
[]
```

### 5.6 Hand-added kernels (cases 01, 04, 05)

```
[LinearFVKernels]
  [fission_oxidation_UIII]             # 0.8 eq per fission: U(III) -> U(IV)
    type = LinearFVSource
    variable = U_III
    source_density = power_density
    scaling_factor = ${fparse -fission_eq_per_J}
  []
  [Xe135_graphite_uptake]              # case 04: a hand-built graphite inventory
    type = FissionProductFeed
    variable = Xe135_graphite
    parent = Xe135
    transfer_rate = ${r_graphite}
  []
[]
```

### 5.7 Key parameters

| action | parameter | default | meaning |
|---|---|---|---|
| radiolysis | `salt_type` | required | `chloride` or `fluoride` network |
| | `metals`, `gas_species` | none | `Zn Cr U`; `Cl2 F2` |
| | `dose_rate` | 0 | dose functor [W/m$^3$] |
| | `g_value_species` / `g_value_overrides` | database | G values [molecules/100 eV] |
| | `kLa`, `volume_ratio` | 0, 1 | gas exchange [1/s], $V_l/V_g$ |
| | `time_derivative` | true | false for SIMPLE |
| | `database` | `msr_database.json` | user database possible |
| fission products | `nuclides` | required | tracked nuclides (database names, e.g. `Xe135`) |
| | `power_density` or `fission_rate_density` | one required | source |
| | `energy_per_fission` | 200 MeV | |
| | `yield_set`, `yield_type` | `U235_thermal`, `effective` | |
| | `neutron_flux` | 0 | one-group flux functor |
| | `removal_groups/rates/destinations` | none | equal length, unique (group, destination) |
| | `track_removed_inventories` | true | `_offgas`, `_plated` variables |
| | `plating_boundaries`, `plating_mass_transfer_coefficient` | none | wall flux $k_m N$ (needs `diffusivity`) |
| both corrosion | `material_class`, `salt_class`, `redox_class`, `position_class`, `flow_factor`, `delta_T_C` | `hastelloy_n`, `fluoride_fuel`, (`purified_baseline` FE / `oxidizing_fef2` flow), `nominal`, 1, 0 | features of the seed correlation |
| | `reference_temperature` | 923.15 K | seed and Arrhenius anchor |
| | `applied_overpotential` | 0.1 V | metal-minus-salt potential |
| | `reference_concentration` | FE: initial salt conc.; flow: 1 | $c_\text{ref}$ [mol/m$^3$] |
| | `kinetics_model` | `reduced_empirical` | or `mstdb_tc_standard_state` (licensed data) |
| FE corrosion | `topology` | required | `two_block`, `salt_only`, `solid_only` |
| | `solve_potential`, `supporting_electrolyte` | false, true | potential solve |
| | `solid_activity` | false | $a_s = c_s/c_{s,0}$ (two_block) |
| | `salt_diffusivity`, `solid_diffusivity` | database, correlation $D_s$ | |
| flow corrosion | `equilibrium_concentrations`, `dissolution_enthalpies`, `seed_fractions` | empty | section 4.4 |
| | `temperature_dependent_kinetics` | true | Arrhenius $i_0(T)$ |
| | `release_variables` | empty | release into existing variables |
| | `diffusivity`, `rhie_chow_user_object` / `velocity` | none | transport |

Note **[code vs doc]**: the default `redox_class` differs between the two corrosion actions
(`purified_baseline` in `[CorrosionPlating]`, `oxidizing_fef2` in `[CorrosionPlatingFlow]`); always
set it explicitly.

### 5.8 Common pitfalls

1. **List every system in both places.** Each species is a separate linear system. Put it in
   `[Problem] linear_sys_names` *and* in `Executioner/system_names` (Transient fixed point) or
   `passive_scalar_systems` (SIMPLE), with one entry per system in
   `passive_scalar_equation_relaxation` and `passive_scalar_absolute_tolerance`. Use `verbose = true`
   to get the names and the parent-before-daughter order. When `track_removed_inventories = true`,
   remember the `_offgas` / `_plated` systems.
2. **Fixed-point sweeps are a model parameter.** `IterationCountConvergence` with
   `converge_at_max_iterations = true` always "converges" after `max_iterations`. Check the
   conservation postprocessors (section 2.7) and the time step (case 01: `dtmax`).
3. **Relaxation and the linear solver for closed loops.** For dissolved metals in a closed loop (case
   05) use a direct solve (`lu` + `mumps`), relaxation 1 on the passive scalars, a tight absolute
   tolerance, and enough `num_iterations`; start the metal near its expected level.
4. **Unused-parameter errors when overriding derived variables.** The inputs compute rates from
   top-level variables with `${fparse ...}` (case 04: `r_spray` from `eps_spray`). Overriding a
   *derived* variable on the command line (e.g. `r_spray=1e-3`) replaces its definition, so the
   variables it was built from are no longer used, and MOOSE stops with
   `unused parameter 'eps_spray'`. Override the *primary* variable (`eps_spray=0.2`) instead, or add
   `--allow-unused` (`-w`) if you really mean it. An unknown command-line name (a typo) gives
   `unused parameter` in `CLI_ARGS` too.
5. **One rate per (group, destination).** Combine rates of the same group and destination into one
   functor (case 04: `r_xe_offgas = r_spray + r_bubble`); a duplicate pair is an error.
6. **Wall plating needs a diffusivity.** `plating_boundaries` requires `diffusivity`; the boundary
   inventory is not in the `_plated` variables.
7. **`seed_fractions`.** Without it, every element is seeded with the full calibrated rate, which
   multiplies the total metal loss by the number of elements.
8. **The applied overpotential must be positive** (relative to $E_0$) to seed a dissolution rate; in
   the equilibrium mode it is not used by the wall law.
9. **Released variables.** With `release_variables`, the corrosion action adds only the wall BC; add
   the time derivative and the diffusion kernel yourself (the wall flux is assembled through the
   diffusion operator).
10. **Gold files and tests.** Each case has a `tests` file (`CSVDiff` against `gold/`, `abs_zero` 1e-9,
    `rel_err` 1e-5 for the light cases, 1e-4 for the MSFR) and is run by
    `./run_tests --test-root examples/training` (`--heavy` for case 05, 8 ranks). A second run of the
    same input (case 04 bubbles) uses `cli_args` with `Outputs/file_base`, and a profile check without
    a new run uses `should_execute = false`. Balance columns of nearly cancelling sums are compared
    with a loose `override_rel_err`. Regenerate gold files only after checking the physics, never to
    make a failing test pass.

## 6. Verification, calibration, validation

| label | definition | what it shows |
|---|---|---|
| **VERIFICATION** | the code solves its equations correctly: comparison with an analytic solution, an exact invariant, or an independent evaluation of the same model | numerical and implementation correctness, not physical truth |
| **CALIBRATION** | parameters were fitted to (or chosen to reproduce) the data the comparison uses; in-sample | the model *can* represent the data; no predictive claim for those values |
| **CODE-TO-CODE** (or code-to-model) | comparison with another code or another analyst's model of the same problem | consistency of models, not truth |
| **VALIDATION** | comparison with independent measured data that were not used to build or tune the model | predictive skill in that domain |
| prediction of a calibrated model | a quantity that was not fitted, from a model calibrated on others | weaker than validation when it comes from the same experiment |

How the five cases are classified (from the case headers and READMEs):

| case | checks | class |
|---|---|---|
| 01 fuel salt redox | U, F and redox invariants (1e-8, 1e-7, 1.9e-3); analytic buffer line ($t_\text{buffer}$ = 342 d, run 339 d); Henry equilibrium (1.001) | **VERIFICATION** |
| 02 xenon power history | Cs-137 burnup monitor 0.9987 (energy integral and decay); Xe worth and peaks are model results | **VERIFICATION** (+ model results; the Xe-137 caveat) |
| 03 corrosion models | (1) BV capsule: `faraday_ratio` = 1 is VERIFICATION; the comparison with NCL-16 is a **model-form check** (wrong model for a loop) | VERIFICATION + model-form |
| | (2) flow wall law: `inventory_balance_error` 7e-4 is VERIFICATION; $c_{\text{eq,Cr}}$, $\Delta H_\text{Cr}$ fitted to M-014 and M-012 (**CALIBRATION**); hot-leg loss and rate are **predictions of a calibrated model** (0.6x the measurements) | VERIFICATION + CALIBRATION + prediction |
| | (3) DRIDN: M-011 to M-015 are fit targets | **CALIBRATION** (in-sample) |
| | (4) MSTDB-TC: endpoint read from the calibration CSVs, not run | **CALIBRATION** (in-sample) |
| 04 MSRE loop | `xe_poison_ratio`, `Mo_salt_ratio`, Ru salt fraction, `rate_vs_correlation`, `faraday_ratio` | **VERIFICATION** |
| | a priori Xe poisoning 1.033 % vs ORNL's own 1.08 % | **CODE-TO-CODE** (code-to-model) |
| | bubble rate $r_b$ fitted to the measured 0.26-0.28 % | **CALIBRATION**; the approach to equilibrium after a power step is an unvalidated prediction |
| | noble metals vs ORNL-4865 Table 12.3 (salt, metal, graphite, off-gas fractions) | **VALIDATION** (the model fails it: salt too low, metal too high) |
| | salt Cr rise vs McCoy (TM-3063): M-003/M-004 are fit targets of the correlation | **CALIBRATION** check |
| 05 MSFR | flow and energy vs the VTB model | **CODE-TO-CODE** |
| | redox balance 5e-10, U total, Cr wall balance 5e-6, `Cr_vs_0d`, `Xe135_vs_0d` 1.023, `Mo99_vs_0d` 1.008, Mo-99 balance 3e-11 | **VERIFICATION** |
| | U(III)/U, Cr levels and fluxes, Xe stripping, Mo-99 partition | model results; no measured MSFR chemistry exists, so **no validation** |

The only independent measured data in the training are the MSRE xenon poisoning (used for
calibration of $r_b$ once the a priori model fails) and the MSRE noble-metal distribution (case 04).
All corrosion comparisons with NCL-16 and the MSRE are in-sample.

## 7. References

- E. L. Compere et al., "Fission product behavior in the Molten Salt Reactor Experiment",
  ORNL-4865 (1975). MSRE fission-product inventories, noble-metal deposition (Table 12.3).
- J. R. Engel, R. C. Steffy, "Xenon behavior in the Molten Salt Reactor Experiment", ORNL-TM-3464
  (1971). Spray stripping, salt-to-graphite mass transfer, measured and predicted Xe poisoning.
- H. E. McCoy, "MSRE Hastelloy N surveillance specimens - fourth group", ORNL-TM-3063 (1970/1971).
  Hastelloy N corrosion and salt Cr in the MSRE.
- "Compatibility of Hastelloy N and LiF-BeF2-UF4 in NCL-16", ORNL-TM-4188 (1972). The NCL-16 loop
  data (validation case ORNL-FL-04), https://www.osti.gov/servlets/purl/4597158.
- "Compatibility studies of potential molten-salt breeder reactor materials in molten fluoride
  salts", ORNL-TM-4286 (1972). Aggregate ORNL thermal-convection loop results (calibration source of
  the reduced correlation).
- INL Virtual Test Bed, MSFR model (`virtual_test_bed/msr/msfr`, `steady/run_ns.i`): EVOL geometry,
  LiF-ThF4-UF4 properties (Brovchenko 2019, as cited in `msfr_chemistry.i`), Griffin power shape.
- D. A. Brown et al., "ENDF/B-VIII.0: The 8th major release of the nuclear reaction data library",
  Nuclear Data Sheets 148, 1 (2018). Decay, fission-yield and capture data of
  `data/fission_product_database.json`.
- MSTDB-TC, Molten Salt Thermal Properties Database - Thermochemical, version 3.1 (ORNL; licensed,
  not redistributed). Files `MSTDB-TC_V3.1_Fluorides_No_Func.dat` and
  `MSTDB-TC_V3.1_Chlorides_No_Func.dat`, SHA-256 bound in `data/advanced_corrosion_models.json`.
- In the repository: `doc/content/theory/moltenSaltRadiolysis.md`, `fissionProductInventory.md`,
  `moltenSaltCorrosion.md`; `validation/corrosion/calibration/MODEL.md`, `ADVANCED_MODELS.md`,
  `PROVENANCE.md`; the five case READMEs.
