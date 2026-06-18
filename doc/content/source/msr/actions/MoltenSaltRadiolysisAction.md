# MoltenSaltRadiolysis

The `[MoltenSaltRadiolysis]` block assembles a complete molten salt radiolysis problem from the
built-in chemistry database. It is a C++ translation of the validated 0D model in
MoltenSaltRadiolysis (`msr_radiolysis`), extended so that the radiolysis species are passive scalars
that advect, diffuse and react on a finite-volume mesh.

Given a salt kernel (`chloride` or `fluoride`), an optional set of metals (`Zn`, `Cr`, `U`) and
optional gas-phase diatomics (`Cl2`, `F2`), the action creates one `MooseLinearVariableFVReal` per
species and the kernels that integrate the problem:

- `LinearFVTimeDerivative` for every species (transient mode only; see [#executioners] below),
- `LinearFVScalarAdvection` (when `rhie_chow_user_object` is supplied, coupling to a segregated
  Navier-Stokes flow) or `LinearFVAdvection` (when a constant `velocity` is supplied) for the
  dissolved species,
- `LinearFVDiffusion` for the dissolved species (when `diffusivity` is supplied),
- `LinearFVSource` for the radiolytic source $S_i = G_i\, \dot{D} / (100\,\text{eV}\cdot N_A)$,
- [MSRReaction](MSRReaction.md) for each species' contribution to every reaction it participates in,
- [MSRGasExchange](MSRGasExchange.md) for each gas-phase diatomic.

## Gas phase and off-gas

When `gas_species` lists diatomics (`Cl2`, `F2`), the action creates a gas-phase concentration for
each and couples it to the dissolved species through [MSRGasExchange](MSRGasExchange.md). The
gas-phase species are transported: they advect with the flow plus a buoyant rise velocity
(`gas_buoyancy_velocity`, a slip velocity directed against gravity) and disperse (`gas_dispersivity`),
so the evolved gas rises and can be vented at a free surface with a
`LinearFVAdvectionDiffusionOutflowBC`. Direct escape of the *dissolved* species at a free surface to
a swept cover gas is imposed with the [MSRDegassingBC](MSRDegassingBC.md) boundary condition. The two
together represent the volumetric (bubble) and surface (degassing) off-gas pathways; see
`examples/operating_msr` for a power-cycle off-gas demonstration.

## Chemistry database

The species lists, reaction networks, Arrhenius rate constants, default G values and Henry
coefficients are read from a JSON file at run time, selected with the `database` parameter. It
defaults to the bundled `data/msr_database.json` (a translation of the validated standalone model),
resolved through the application's registered data path, so no path is needed for the standard
networks. To change the chemistry without recompiling, point `database` at your own JSON file with
the same structure (`kernels`, `metals`, `G_values`, `henry`, `gas_pairs`).

## Species naming

Species are named with plain identifiers so that they are valid finite-volume variable names. The
mapping from the literature symbols is:

| Symbol | Variable | Symbol | Variable |
| :----- | :------- | :----- | :------- |
| $e_s^-$ | `e_sol` | $F^-$ | `F_ion` |
| $\text{Cl}^-$ | `Cl_ion` | $F^\bullet$ | `F_rad` |
| $\text{Cl}^\bullet$ | `Cl_rad` | $F_2$ (diss.) | `F2_diss` |
| $\text{Cl}_2^{\bullet -}$ | `Cl2m_rad` | $F_2$ (gas) | `F2_gas` |
| $\text{Cl}_3^-$ | `Cl3_ion` | $\text{Zn}^{2+}/\text{Zn}^+$ | `Zn_II`/`Zn_I` |
| $\text{Cl}_2$ (diss.) | `Cl2_diss` | $\text{Cr}^{2+}/\text{Cr}^{3+}/\text{Cr}^+$ | `Cr_II`/`Cr_III`/`Cr_I` |
| $\text{Cl}_2$ (gas) | `Cl2_gas` | $\text{U}^{4+}/\text{U}^{3+}$ | `U_IV`/`U_III` |

## Solver systems

Each species is placed in its own linear solver system named `<species>_sys` so that the segregated
fixed-point solve iterates the (lagged) mass-action coupling. These system names must be listed in
`[Problem] linear_sys_names` and in the executioner's `system_names` (or `passive_scalar_systems`
under a SIMPLE/PIMPLE executioner). Setting `verbose = true` prints the generated species and system
names.

## Executioners id=executioners

The action supports both steady and transient solves through the `time_derivative` parameter, which
controls whether a `LinearFVTimeDerivative` is added to each species:

| Use case | `time_derivative` | Executioner |
| :------- | :---------------- | :---------- |
| Pure chemistry / prescribed velocity, transient | `true` (default) | `Transient` with `multi_system_fixed_point = true` and a `[Convergence]` |
| Flow-coupled, transient | `true` (default) | `PIMPLE` (segregated Navier-Stokes); species listed in `passive_scalar_systems` |
| Flow-coupled, steady-state dose equilibrium | `false` | `SIMPLE` (segregated Navier-Stokes); species listed in `passive_scalar_systems` |

`SIMPLE` rejects time-derivative kernels, so steady flow-coupled cases must set
`time_derivative = false`. See `test/tests/msr/flow_simple` (SIMPLE) and `test/tests/msr/flow_pimple`
(PIMPLE) for complete production examples.

## Overrides

The default G values, initial concentrations and rate constants come from the database, selected by
`salt_type` and `metals`. The dose rate (`dose_rate`), temperature (`temperature`), initial
concentrations (`initial_condition_species`/`initial_condition_values`) and G values
(`g_value_species`/`g_value_overrides`) can be overridden in the input.

Both `dose_rate` and `temperature` are functors, so they may be constants or spatially/temporally
varying fields (functions, variables, etc.). A non-uniform `dose_rate` lets the radiolytic source
follow a realistic radiation field; see `examples/molten_salt_corrosion` for a flowing-salt example
in which a Gaussian near-core dose profile drives chromium oxidation.

## Example

!listing test/tests/msr/0d_disproportionation/disproportionation.i block=MoltenSaltRadiolysis

!syntax parameters /MoltenSaltRadiolysis

!syntax inputs /MoltenSaltRadiolysis
