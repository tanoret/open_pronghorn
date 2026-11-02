# CorrosionPlatingFlow

The `[CorrosionPlatingFlow]` block sets up salt-side molten salt corrosion of a wall on the linear
finite-volume segregated (SIMPLE/PIMPLE) basis, so the dissolved corrosion products are transported as
passive scalars by the same flow solve that carries the radiolysis chemistry and the energy equation.
It is the corrosion counterpart of
[MoltenSaltRadiolysis](../../msr/actions/MoltenSaltRadiolysisAction.md) and the natural way to couple
corrosion into a flowing MSR (Navier-Stokes + k-epsilon + energy + radiolysis); see
`examples/flowing_msr_corrosion`.

For each tracked element the action creates a linear finite-volume cation concentration `c_<El>`, its
time derivative, advection by the solved flow (`rhie_chow_user_object`) or a prescribed `velocity`, and
molecular plus turbulent diffusion (`diffusivity`), and applies the Butler-Volmer electrode reaction at
the wall through [CorrosionLinearFVButlerVolmerBC](../linearfvbcs/CorrosionLinearFVButlerVolmerBC.md).
Optional inlet (Dirichlet) and outlet (outflow) boundaries are added when `inlet_boundary` /
`outlet_boundary` are supplied. Each generated solver system `c_<El>_sys` must be listed in
`[Problem] linear_sys_names` and in the executioner's `passive_scalar_systems`.

## Kinetics model

The default `kinetics_model = reduced_empirical` preserves the original action behavior. It seeds each
exchange current from the validated effective correlation
([MoltenSaltCorrosionModel](../base/MoltenSaltCorrosionModel.md)) using `material_class`, `salt_class`,
`redox_class`, `flow_factor`, `delta_T_C`, and `reference_temperature`. With
`temperature_dependent_kinetics = true`, the exchange current is then Arrhenius in the local
temperature field.

`kinetics_model = mstdb_tc_standard_state` evaluates one provenance-bound, static standard-state
endpoint described for
[MSTDBTCCorrosionUserObject](../userobjects/MSTDBTCCorrosionUserObject.md). It requires exactly
`elements = 'Cr Fe Ni'` and seeds each wall current from `front_rate_um_y * source_fraction[element]`,
so the total dissolution-front rate is not applied three times. This mode requires
`advanced_database`, external `fluoride_database` and `chloride_database` files,
`reference_cold_temperature`, `reference_exposure_time`, `area_to_salt_mass`, and
`inventory_coupling_factor`; `reference_temperature` is the hot-side temperature. The advanced JSON
must bind MSTDB-TC V3.1 and the SHA-256 hashes of both files. Hash and edition checks cannot be
disabled. The same artifact semantically binds the consumed base-model parameters, alloy densities,
and Cr/Fe/Ni transport and Butler-Volmer properties; a mismatch is rejected before Action objects are
generated. `allow_extrapolation` affects only out-of-interval Gibbs evaluation. If
`temperature_dependent_kinetics = true`, the static MSTDB seed is subsequently scaled in space with
the base corrosion database activation energy; this reduced closure does not reevaluate MSTDB at
each quadrature point. `inventory_coupling_factor` scales dissolved salt inventory only and does not
scale the static cold-capture mass-gain closure.

In MSTDB mode, `position_class`, `surface_class` (not exposed by this flow Action), and `delta_T_C`
do not govern the endpoint. The explicit `reference_temperature` and
`reference_cold_temperature` inputs define the hot- and cold-side temperatures, respectively.

The flow action's explicit `diffusivity` functor and `reference_concentration` are runtime transport
inputs that replace the corresponding database defaults. They do not disable semantic validation of
the source database or permit changes to bound electrochemical properties.

The MSTDB endpoint controls only the exchange-current magnitude and Cr/Fe/Ni source split. The wall
law's equilibrium/reversal behavior still uses `E0` and `c_ref` from the base `database`; MSTDB
affinities are not substituted as boundary equilibrium potentials.

The exchange-current seed uses the production Butler-Volmer bracket at
`applied_overpotential - E0`, including the boundary law's exponent clipping. Consequently custom
nonzero-`E0` databases reproduce the planned reference current at `c = c_ref`; the bundled zero-`E0`
database retains the original result. Numeric potential and inlet-concentration functors are emitted
with 17 significant digits so their runtime values match the values used while building the seed.

This action models the advective-diffusive transport of corrosion products and the wall kinetics. The
full electromigration / electric-potential (electrophoresis) solve and two-block (salt + solid)
interface kinetics are provided by the finite-element [CorrosionPlating](CorrosionPlatingAction.md)
action; both share the same kinetics options and databases.

## Example

!listing examples/flowing_msr_corrosion/flowing_msr_corrosion.i block=CorrosionPlatingFlow

!syntax parameters /CorrosionPlatingFlow

!syntax inputs /CorrosionPlatingFlow
