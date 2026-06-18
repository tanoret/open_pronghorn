# CorrosionPlatingFlow

The `[CorrosionPlatingFlow]` block sets up salt-side molten salt corrosion of a wall on the linear
finite-volume segregated (SIMPLE/PIMPLE) basis, so the dissolved corrosion products are transported as
passive scalars by the same flow solve that carries the radiolysis chemistry and the energy equation.
It is the corrosion counterpart of [MoltenSaltRadiolysis](MoltenSaltRadiolysisAction.md) and the
natural way to couple corrosion into a flowing MSR (Navier-Stokes + k-epsilon + energy + radiolysis);
see `examples/flowing_msr_corrosion`.

For each tracked element the action creates a linear finite-volume cation concentration `c_<El>`, its
time derivative, advection by the solved flow (`rhie_chow_user_object`) or a prescribed `velocity`, and
molecular plus turbulent diffusion (`diffusivity`), and applies the Butler-Volmer electrode reaction at
the wall through [CorrosionLinearFVButlerVolmerBC](../linearfvbcs/CorrosionLinearFVButlerVolmerBC.md).
Optional inlet (Dirichlet) and outlet (outflow) boundaries are added when `inlet_boundary` /
`outlet_boundary` are supplied. Each generated solver system `c_<El>_sys` must be listed in
`[Problem] linear_sys_names` and in the executioner's `passive_scalar_systems`.

The exchange current is seeded from the validated effective corrosion correlation
([MoltenSaltCorrosionModel](../base/MoltenSaltCorrosionModel.md)) using the case feature selectors
(`material_class`, `salt_class`, `redox_class`, `flow_factor`, `delta_T_C`, `reference_temperature`),
so the wall reproduces the calibrated dissolution rate; with `temperature_dependent_kinetics = true`
the exchange current is Arrhenius in the local temperature field, so the wall corrodes faster where the
salt is hotter.

This action models the advective-diffusive transport of corrosion products and the wall kinetics. The
full electromigration / electric-potential (electrophoresis) solve and two-block (salt + solid)
interface kinetics are provided by the finite-element [CorrosionPlating](CorrosionPlatingAction.md)
action; both share the same calibrated kinetics and database.

## Example

!listing examples/flowing_msr_corrosion/flowing_msr_corrosion.i block=CorrosionPlatingFlow

!syntax parameters /CorrosionPlatingFlow

!syntax inputs /CorrosionPlatingFlow
