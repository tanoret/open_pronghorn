# Molten Salt Radiolysis

## Overview

Ionizing radiation deposits energy in a molten salt, generating reactive primary species (solvated
electrons and halogen radicals) that drive a network of homogeneous reactions. The products of
interest include the molecular halogens (Cl$_2$, F$_2$) that can evolve into a cover gas and the
redox states of dissolved metals (Zn, Cr, U). OpenPronghorn integrates this chemistry together with
species transport, treating each radiolysis species as a passive scalar.

## Governing equations

For a dissolved species $i$ with concentration $C_i$ (mol/m$^3$),

\begin{equation}
\frac{\partial C_i}{\partial t}
+ \nabla\!\cdot(\mathbf{u}\,C_i)
- \nabla\!\cdot(D_i^\text{eff} \nabla C_i)
= S_i + \sum_{\rho} \nu_{i\rho}\, r_\rho ,
\end{equation}

where $\mathbf{u}$ is the (optional) advecting velocity, $D_i^\text{eff}$ the effective diffusivity,
$S_i$ the radiolytic source, $r_\rho$ the rate of reaction $\rho$ and $\nu_{i\rho}$ the net
stoichiometric coefficient of $i$ in $\rho$. The advection and diffusion terms reuse the framework
`LinearFVScalarAdvection`/`LinearFVAdvection` and `LinearFVDiffusion` kernels; the source and
reaction terms are provided by [MSRReaction](MSRReaction.md) and `LinearFVSource`.

In a turbulent flow the effective diffusivity combines molecular and turbulent contributions,
$D_i^\text{eff} = D_i + \mu_t/(\rho\,\mathrm{Sc}_t)$, with $\mu_t$ the turbulent viscosity from the
k-$\epsilon$ model and $\mathrm{Sc}_t$ the turbulent Schmidt number; it is supplied as a functor
(e.g. a `ParsedFunctorMaterial`) to the action's `diffusivity` parameter. Because the velocity, the
temperature (energy conservation) and $\mu_t$ (turbulence) are all functors, the radiolysis transport
and chemistry couple directly to a full Navier-Stokes + energy + turbulence solve. The
`examples/msr_core` case demonstrates this complete coupling — linear momentum, mass and energy
conservation with k-$\epsilon$ turbulence driving the advection and molecular + turbulent diffusion of
a neutron-flux-shaped radiolytic source.

## Radiolytic source

The volumetric source of species $i$ follows from its radiolytic yield (G value, molecules per
100 eV) and the volumetric dose rate $\dot{D}$ (J/m$^3$/s):

\begin{equation}
S_i = \frac{G_i\, \dot{D}}{100\,\text{eV}\cdot N_A},
\end{equation}

with $N_A$ Avogadro's number. The 100 eV and Avogadro factors convert the yield to mol/m$^3$/s.

## Reaction kinetics

Each reaction is mass-action with an Arrhenius rate constant; see [MSRReaction](MSRReaction.md) for
the rate-constant forms and the implicit/explicit splitting used in the linear finite-volume
discretization. The chloride and fluoride reaction networks, metal redox templates, default G values
and Henry coefficients are stored in a JSON database (`data/msr_database.json`) that mirrors the
validated standalone model; the [MoltenSaltRadiolysis action](MoltenSaltRadiolysisAction.md) reads it
(or a user-supplied file) and builds the variables and kernels.

## Gas-liquid exchange

Dissolved diatomics exchange with a well-mixed headspace by a Henry-law mass-transfer model; see
[MSRGasExchange](MSRGasExchange.md).

## Numerical treatment

The species form segregated linear systems solved with a fixed-point (Gauss-Seidel) iteration that
converges the lagged mass-action coupling. Because the kinetics are integrated in time, a transient
executioner is used: a generic `Transient` with `multi_system_fixed_point` for chemistry-only or
prescribed-velocity problems, or `PIMPLE` when coupling to a segregated Navier-Stokes flow. A steady
dose-driven balance can be solved with `SIMPLE` by disabling the time derivatives
(`time_derivative = false`).

## Conservation and convergence

Each reaction's consumption is treated implicitly in the consumed species while its production is
lagged, so the rate "seen" by the reactants and the products is identical only when the segregated
fixed-point iteration is converged (lagged values equal the current iterate). At convergence the
scheme conserves the underlying elements exactly (for backward Euler); short of convergence there is
a conservation drift that decreases with the number of fixed-point sweeps. For a closed irradiated
cell the total-chromium conservation error scales roughly inversely with the sweep count:

| fixed-point sweeps | total-Cr conservation error |
| :----------------- | :-------------------------- |
| 1  | 1.5e-3 |
| 5  | 5.2e-4 |
| 20 | 1.5e-4 |
| 50 | 4.7e-5 |

(see `test/tests/msr/conservation`). In practice: iterate the `multi_system_fixed_point` coupling to
a real tolerance (or use enough sweeps) for the generic `Transient`, and tighten the
`passive_scalar_absolute_tolerance` / raise `num_iterations` for `SIMPLE`/`PIMPLE`. The radical
chemistry is stiff (rate constants up to ~1e10 m^3/mol/s, and the radical-forming reaction has a
pseudo-rate of order 1e11 1/s at bulk halide concentrations), so strongly-driven problems converge
slowly; for those, a monolithic Newton solve of the coupled chemistry is the more efficient and
exactly conservative choice.
