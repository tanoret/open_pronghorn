# Full MSR core: coupled flow, heat, turbulence and radiolysis

This is the fully-coupled, multiphysics MSR fuel-salt core example. A vertical channel of fluoride
fuel salt (FLiBe-relevant) flows upward through the active core, and everything below is solved
together on a linear finite-volume segregated (SIMPLE) basis:

- **Linear momentum + mass (pressure) conservation** — incompressible Navier-Stokes with the
  Rhie-Chow mass flux;
- **k-epsilon turbulence** — the turbulent viscosity `mu_t` (the repo's k-epsilon model with wall
  functions);
- **Energy conservation** — advected enthalpy, conduction plus turbulent conduction
  (`k_eff = k + mu_t cp / Pr_t`), and a nuclear heat source;
- **Radiolysis** — F2 produced by a **neutron-flux-shaped source**, transported by advection and by
  **molecular plus turbulent diffusion** (`D_eff = D + mu_t / (rho Sc_t)`), with the evolved gas
  phase rising buoyantly.

## Coupling

The nuclear heat and the radiolytic dose follow the same separable cosine flux shape
(`flux_shape`), peaked at the core mid-plane — a stand-in for the fundamental-mode neutron flux. The
energy equation supplies the temperature field `T_fluid`, which is passed to the radiolysis chemistry
(here it sets the gas-exchange equilibrium). The k-epsilon solution supplies `mu_t`, which enters both
the turbulent conductivity (energy) and the turbulent diffusivity (species) through
`ParsedFunctorMaterial` definitions. The species are advected by the same Rhie-Chow flow that the
momentum solve produces.

## Output

The CSV reports the outlet temperature, peak turbulent viscosity, and the F2 dissolved/gas state;
the Exodus output shows the velocity, `T_fluid`, `mu_t`, `F2_diss` and `F2_gas` fields. With the
default settings the salt heats by ~30 K across the core, the turbulent viscosity reaches
~700x molecular, and the radiolytic F2 is produced in the flux-peaked region, advected and turbulently
dispersed downstream, and partially carried off as gas.

## Notes

- The radical yields are lumped into the net F2 yield (`g_value_overrides`) so the chemistry is a
  source plus gas exchange; this keeps the already-large coupled solve inexpensive. Remove the
  `g_value_*` overrides to integrate the explicit fluoride radical network (a stiffer solve).
- The SIMPLE absolute tolerances are set just above the residual floors of each field (turbulence
  and energy residuals normalize to finite floors), so the coupled problem converges in ~100
  iterations rather than running to the iteration cap.
- For a transient power cycle (ramp-up / steady / ramp-down) with surface degassing, see
  `examples/operating_msr`; for steady radiation-driven corrosion, see
  `examples/molten_salt_corrosion`.
