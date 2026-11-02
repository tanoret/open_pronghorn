# CorrosionPlating

The `[CorrosionPlating]` block assembles a complete molten salt corrosion and plating problem from
the built-in corrosion database. It creates, for each tracked element, the salt-side cation and the
solid-side metal variables, the Nernst-Planck transport of the cations (diffusion plus
electromigration), the optional electric-potential (current-continuity) equation, and the
Butler-Volmer electrode kinetics. `kinetics_model` selects either the validated reduced-empirical
seed (the default) or a provenance-bound MSTDB-TC standard-state seed.

## Governing physics

For each dissolved cation $M^{z+}$ on the salt block the action solves the Nernst-Planck balance

!equation
\frac{\partial c_k}{\partial t} + \nabla\cdot\left(-D_k\nabla c_k - \frac{z_k F}{RT}D_k c_k\nabla\phi + \vec{u}\,c_k\right) = 0

through [CorrosionNernstPlanckFlux](../kernels/CorrosionNernstPlanckFlux.md) (diffusion + migration)
and [CorrosionAdvection](../kernels/CorrosionAdvection.md) (optional flow advection). When
`solve_potential = true` the salt electric potential is solved from current continuity,

!equation
\nabla\cdot\vec{i} = 0, \qquad \vec{i} = -F\sum_k z_k D_k \nabla c_k - \kappa\nabla\phi, \qquad \kappa = \frac{F^2}{RT}\sum_k z_k^2 D_k c_k,

through [CorrosionCurrentContinuity](../kernels/CorrosionCurrentContinuity.md), with the constant mode
pinned at `pin_potential_boundary`. Each tracked solid metal diffuses in the alloy through
[CorrosionSolidDiffusion](../kernels/CorrosionSolidDiffusion.md).

The metal/salt reaction $M(s)\rightleftharpoons M^{z+} + z\,e^-$ is imposed with the Butler-Volmer
kinetics

!equation
i_{BV} = i_0\left[\exp\!\left(\frac{\alpha_a z F\,\eta}{RT}\right) - \frac{c}{c_{ref}}\exp\!\left(-\frac{\alpha_c z F\,\eta}{RT}\right)\right], \qquad \eta = (\phi_{metal} - \phi_{salt}) - E_0,

as an interface kernel ([ButlerVolmerInterface](../interfacekernels/ButlerVolmerInterface.md)) in the
two-block topology or as a boundary condition ([ButlerVolmerBC](../bcs/ButlerVolmerBC.md)) in a
single-domain topology. The anodic branch carries the (constant) metal activity and the cathodic
branch the dissolved-ion activity $c/c_{ref}$, so fresh salt corrodes at a bounded calibrated rate
while the dissolved metal accumulates and the reaction approaches equilibrium.

## Topologies

| `topology` | Phases modeled | Butler-Volmer object |
| :--------- | :------------- | :------------------- |
| `two_block` | salt and solid as separate blocks | [ButlerVolmerInterface](../interfacekernels/ButlerVolmerInterface.md) on `interface_boundary` |
| `salt_only` | salt only (metal is external) | [ButlerVolmerBC](../bcs/ButlerVolmerBC.md) on `reaction_boundary` |
| `solid_only` | solid only (salt is external) | [ButlerVolmerBC](../bcs/ButlerVolmerBC.md) on `reaction_boundary` |

`applied_overpotential` always denotes the reference metal-minus-salt potential difference
$\phi_{metal}-\phi_{salt}$. In a single-domain topology, `applied_potential` is instead the absolute
potential of the absent counter phase. If that functor is omitted and the salt potential is not
solved, the action supplies `+applied_overpotential` as the external metal potential for `salt_only`,
or `-applied_overpotential` as the external salt potential for `solid_only`. When the salt potential
is solved and pinned to `pin_potential_value`, the generated metal absolute potential in `salt_only`
and `two_block` is instead `pin_potential_value + applied_overpotential`. The pin is therefore only
an absolute gauge offset: at the reference/pinned salt potential the generated interface, species
and charge boundaries, and diagnostics all retain the seeded metal-minus-salt difference. A solved
salt potential can subsequently change that difference locally away from the pin. Both scalar
inputs, and their generated absolute-potential sum, must be finite.

`solid_only` additionally requires `counter_concentration`; both the boundary flux and the
rate/recession diagnostics evaluate that same external-salt concentration functor. Because this
topology has no salt variable, `solve_potential = true` is rejected rather than silently creating no
potential equation. When an explicit `applied_potential` functor is supplied in `salt_only`, it is a
user-owned absolute metal potential and is not shifted by `pin_potential_value`; it can intentionally
differ from the scalar reference operating point used to seed `i0`, and the diagnostic follows that
functor exactly.

## Kinetics seeding and reproduction

With `kinetics_model = reduced_empirical`, the exchange current density is seeded so that, at the
imposed `applied_overpotential` and initial (reference) concentration, the Butler-Volmer current
Faradaically equals the calibrated dissolution rate computed by
[MoltenSaltCorrosionModel](../base/MoltenSaltCorrosionModel.md). This default retains the original
behavior: each selected element receives that total-rate seed. A 0D cell therefore reproduces the
reference rate exactly; see `validation/corrosion` for the reference-case reproduction.

For every kinetics model, the seed bracket is evaluated by the same production Butler-Volmer law as
the generated interface or boundary object. In particular, its activation overpotential is exactly
`applied_overpotential - E0`, and the same exponent clipping is applied. The bundled database has
the historical zero-$E_0$ values, for which this is identical to the original seed formula; a custom
nonzero-$E_0$ database no longer shifts the reproduced reference current.

With `kinetics_model = mstdb_tc_standard_state`, the action evaluates the static standard-state/Nernst
engineering model described for [MSTDBTCCorrosionUserObject](../userobjects/MSTDBTCCorrosionUserObject.md).
This mode requires exactly `elements = 'Cr Fe Ni'`. It multiplies the total dissolution-front rate by
each model source fraction before constructing the corresponding exchange current; the total rate is
therefore not duplicated across the three species.

This bridge uses the MSTDB-TC endpoint only for the exchange-current magnitude and Cr/Fe/Ni source
split. The subsequent reversible Butler-Volmer/Nernst boundary law continues to use each element's
`E0` and `c_ref` from the base `database`; the reported MSTDB affinities are diagnostics and are not
substituted for those equilibrium potentials. Consequently, this action is a static transport
coupling of the calibrated endpoint, not a local thermochemical-equilibrium solve.

The MSTDB mode also requires `advanced_database`, external `fluoride_database` and
`chloride_database` files, `reference_cold_temperature`, `reference_exposure_time`,
`area_to_salt_mass`, and `inventory_coupling_factor`. `reference_temperature` is the hot-side
temperature and must be at least the cold-side temperature. The advanced JSON binds calibration to
MSTDB-TC V3.1 and both raw-file SHA-256 hashes. The action always enforces that edition and those
hashes. It also verifies the calibrated semantic subset of the base corrosion database, including
the required correlation parameters, alloy densities, and Cr/Fe/Ni valence, molar mass,
diffusivity, $E_0$, transfer coefficients, and reference concentration. Changing any of those
properties requires a new provenance-bound calibration artifact. `allow_extrapolation` only changes
the Gibbs interval policy. In this static model, `inventory_coupling_factor` scales dissolved salt
inventory only, not cold-capture mass gain.

In MSTDB mode, `position_class`, `surface_class`, and `delta_T_C` do not govern the endpoint. The
explicit `reference_temperature` and `reference_cold_temperature` inputs define the hot- and
cold-side temperatures, respectively.

An explicitly supplied `salt_diffusivity` or `reference_concentration` replaces the corresponding
provenance-bound database default for that application. These documented runtime overrides do not
disable the semantic check of the source database or alter the bound valence, molar mass, $E_0$, or
transfer coefficients.

## Outputs

The action creates two boundary auxiliary fields through
[CorrosionRateAux](../auxkernels/CorrosionRateAux.md) or
[CorrosionTotalRateAux](../auxkernels/CorrosionTotalRateAux.md): the
instantaneous penetration rate `corrosion_rate_um_y` and the time-integrated recession/plating
thickness `recession_um`. The reduced-empirical path retains the historical
`recession_element` diagnostic. In MSTDB-TC mode, each Cr/Fe/Ni current is converted with its own
valence and molar mass before the three penetration-rate shares are summed; both fields therefore
report the total modeled front rate/recession rather than the default Cr share. In a single-domain
topology these diagnostics receive the same runtime concentration and counter-phase-potential
functors as the generated Butler-Volmer boundary.

## Database

The elements (valence, molar mass, diffusivity, $E_0$, transfer coefficients, reference
concentration), the engineering material tables (density, chromium fraction) and the calibrated
correlation parameters are read from a JSON file selected with `database`, defaulting to the bundled
`data/corrosion_database.json`. Point `database` at your own file to change the base chemistry without
recompiling. Raw MSTDB-TC files are external, permission-controlled runtime inputs and are not
distributed by OpenPronghorn.

## Example

!listing test/tests/corrosion/action_two_block/two_block.i block=CorrosionPlating

!syntax parameters /CorrosionPlating

!syntax inputs /CorrosionPlating
