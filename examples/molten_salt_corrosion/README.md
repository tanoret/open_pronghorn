# Radiation-driven chromium corrosion in a flowing molten chloride salt

This example couples every part of the molten salt radiolysis capability in a single, realistic
multi-dimensional problem.

## Scenario

A section of a molten LiCl-KCl coolant channel carries dissolved Cr(II), a corrosion product shed by
structural alloys. The salt flows past a localized, intense gamma field (peaked near the reactor
core, represented by a Gaussian dose profile in the axial direction). In that region radiolysis
drives the chlorine radical chain

```
dose -> Cl_rad ;   Cl_rad + Cl_ion -> Cl2m_rad ;   Cl2m_rad + Cr(II) -> Cr(III) + 2 Cl-
```

so the salt oxidizes Cr(II) to Cr(III), while solvated electrons partially reduce the chromium back
(`e_sol + Cr(III) -> Cr(II)`, `e_sol + Cr(II) -> Cr(I)`). The steady solution is a spatial map of the
chromium oxidation state, i.e. where the salt becomes more oxidizing and corrosive.

## What it demonstrates

- laminar molten-salt flow solved with the linear finite-volume segregated **SIMPLE** algorithm and
  the Rhie-Chow mass flux;
- the radiolysis species created and transported by the `[MoltenSaltRadiolysis]` action as passive
  scalars (advection by the solved velocity field, plus molecular diffusion);
- a **spatially varying radiolytic dose** supplied as a functor (`dose_rate = dose_profile`);
- the full chloride + Cr mass-action chemistry;
- **steady-state** dose-driven operation (`time_derivative = false`, compatible with SIMPLE).

## Running and viewing

```
open_pronghorn-opt -i corrosion_channel.i
```

The Exodus output `corrosion_channel_out.e` contains the velocity, `Cr_II`, `Cr_III` and `Cl2m_rad`
fields; the CSV output reports outlet-averaged chromium oxidation state and a total-chromium
conservation check. With the near-core dose (1e7 J/m^3/s) the radiolytic radical fields are
localized in the irradiated region and oxidize a small fraction of the chromium, with total
chromium conserved to about 0.1 %.

## Convergence and conservation

The radical chemistry is stiff: at the bulk chloride concentration the radical-forming reaction has
a pseudo-rate of order 1e11 1/s, so the lagged segregated coupling converges slowly and both the
conservation error and the (small) net oxidation depend on how tightly the fixed point is iterated.
This example therefore uses a large iteration count and tight scalar absolute tolerances
(`num_iterations = 2000`, `passive_scalar_absolute_tolerance = 1e-11`) so that total chromium is
conserved to ~0.1 % — this is the "Tier-1" convergence control, validated quantitatively by the
closed-cell conservation study in `test/tests/msr/conservation`.

Because radiolytic oxidation (`Cl2m_rad + Cr(II)`) and reduction (`e_sol + Cr`) nearly balance, the
net steady oxidation at realistic dose is modest and continues to tighten with more iterations. For
quantitatively converged results in strongly-driven (high-dose) regimes, a monolithic Newton solve
of the chemistry (all species in one nonlinear system, exact per-reaction conservation) is the
appropriate tool; this example demonstrates the coupled multi-physics machinery and the spatial
structure of the radiation-driven chemistry.

## Notes

- To explore the transient build-up instead of the steady state, use the `PIMPLE` executioner with
  `time_derivative = true` (see `test/tests/msr/flow_pimple`).
