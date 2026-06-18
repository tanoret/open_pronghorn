# CorrosionPlating

The `[CorrosionPlating]` block assembles a complete molten salt corrosion and plating problem from
the built-in corrosion database. It creates, for each tracked element, the salt-side cation and the
solid-side metal variables, the Nernst-Planck transport of the cations (diffusion plus
electromigration), the optional electric-potential (current-continuity) equation, and the
Butler-Volmer electrode kinetics, with the kinetics seeded from the validated effective-correlation
port [MoltenSaltCorrosionModel](MoltenSaltCorrosionModel.md).

## Governing physics

For each dissolved cation $M^{z+}$ on the salt block the action solves the Nernst-Planck balance

!equation
\frac{\partial c_k}{\partial t} + \nabla\cdot\left(-D_k\nabla c_k - \frac{z_k F}{RT}D_k c_k\nabla\phi + \vec{u}\,c_k\right) = 0

through [CorrosionNernstPlanckFlux](CorrosionNernstPlanckFlux.md) (diffusion + migration) and
[CorrosionAdvection](CorrosionAdvection.md) (optional flow advection). When `solve_potential = true`
the salt electric potential is solved from current continuity,

!equation
\nabla\cdot\vec{i} = 0, \qquad \vec{i} = -F\sum_k z_k D_k \nabla c_k - \kappa\nabla\phi, \qquad \kappa = \frac{F^2}{RT}\sum_k z_k^2 D_k c_k,

through [CorrosionCurrentContinuity](CorrosionCurrentContinuity.md), with the constant mode pinned at
`pin_potential_boundary`. Each tracked solid metal diffuses in the alloy through
[CorrosionSolidDiffusion](CorrosionSolidDiffusion.md).

The metal/salt reaction $M(s)\rightleftharpoons M^{z+} + z\,e^-$ is imposed with the Butler-Volmer
kinetics

!equation
i_{BV} = i_0\left[\exp\!\left(\frac{\alpha_a z F\,\eta}{RT}\right) - \frac{c}{c_{ref}}\exp\!\left(-\frac{\alpha_c z F\,\eta}{RT}\right)\right], \qquad \eta = (\phi_{metal} - \phi_{salt}) - E_0,

as an interface kernel ([ButlerVolmerInterface](ButlerVolmerInterface.md)) in the two-block topology
or as a boundary condition ([ButlerVolmerBC](ButlerVolmerBC.md)) in a single-domain topology. The
anodic branch carries the (constant) metal activity and the cathodic branch the dissolved-ion activity
$c/c_{ref}$, so fresh salt corrodes at a bounded calibrated rate while the dissolved metal accumulates
and the reaction approaches equilibrium.

## Topologies

| `topology` | Phases modeled | Butler-Volmer object |
| :--------- | :------------- | :------------------- |
| `two_block` | salt and solid as separate blocks | [ButlerVolmerInterface](ButlerVolmerInterface.md) on `interface_boundary` |
| `salt_only` | salt only (metal is external) | [ButlerVolmerBC](ButlerVolmerBC.md) on `reaction_boundary` |
| `solid_only` | solid only (salt is external) | [ButlerVolmerBC](ButlerVolmerBC.md) on `reaction_boundary` |

## Kinetics seeding and reproduction

The exchange current density is seeded so that, at the imposed `applied_overpotential` and the initial
(reference) concentration, the Butler-Volmer current Faradaically equals the calibrated dissolution
rate computed by [MoltenSaltCorrosionModel](MoltenSaltCorrosionModel.md) from the case feature
selectors (`material_class`, `salt_class`, `redox_class`, `position_class`, `flow_factor`,
`delta_T_C`, `reference_temperature`). As a result a 0D cell reproduces the reference rate exactly;
see `validation/corrosion` for the reproduction of all 76 reference cases.

## Outputs

The action creates two boundary auxiliary fields through [CorrosionRateAux](CorrosionRateAux.md): the
instantaneous penetration rate `corrosion_rate_um_y` and the time-integrated recession/plating
thickness `recession_um` on the recession-controlling element.

## Database

The elements (valence, molar mass, diffusivity, $E_0$, transfer coefficients, reference
concentration), the engineering material tables (density, chromium fraction) and the calibrated
correlation parameters are read from a JSON file selected with `database`, defaulting to the bundled
`data/corrosion_database.json`. Point `database` at your own file to change the chemistry without
recompiling.

## Example

!listing test/tests/corrosion/action_two_block/two_block.i block=CorrosionPlating

!syntax parameters /CorrosionPlating

!syntax inputs /CorrosionPlating
