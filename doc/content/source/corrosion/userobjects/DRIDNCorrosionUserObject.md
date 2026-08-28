# DRIDNCorrosionUserObject

`DRIDNCorrosionUserObject` advances the 19-state Dynamic Redox-Inventory-Depletion Network (DRIDN)
at `TIMESTEP_END`. Its state, elapsed physical time, and adaptive-integrator counters are restartable.
MOOSE time is interpreted in seconds and converted to years before the reduced network is advanced.
The execution schedule is fixed to exactly `TIMESTEP_END`; attempts to add or substitute execution
flags are rejected because they could advance the full MOOSE timestep more or less than once.

The complete `DRIDNModel::Context` is supplied explicitly: Cr/Fe/Ni alloy fractions, reaction
affinities, activity floor, cold-capture fractions, initial physical inventories, geometry,
inventory coupling, redox switches, mass-transfer limit, diffusion, and damage multipliers. No
governing term depends on a measurement ID, source citation, experiment family, or selected output.
The Nernst activity floor is distinct from actual dissolved inventory.

Arbitrary explicit physical contexts are accepted by design; DRIDN does not silently substitute the
calibration fixture. The frozen fit and validation metrics in `advanced_corrosion_models.json`,
however, apply only to the provenance-bound context table used for that calibration. Predictions for
other contexts are model applications or extrapolations, not additional validation evidence.

`charge_transfer_mode = legacy_irreversible` is the default because the stored parameter fit and
validation metrics were generated with the original positive-only charge-transfer law. That law
does not enforce detailed balance. `affinity_gated_dissolution` makes dissolution vanish at and below
the weighted equilibrium affinity, but is experimental and is **not** covered by the stored
calibration or validation metrics.

The adaptive RK4 step-doubling integrator rejects nonphysical intermediate states and restores the
salt-side species inventory identity after accepted steps. The reported `mass_balance_relative_error`
checks dissolved inventory, coupon deposits, bulk capture, and cumulative source for each modeled
element. It is not a proof of full alloy, charge, or redox-equivalent conservation.

Use [AdvancedCorrosionModelPostprocessor](../postprocessors/AdvancedCorrosionModelPostprocessor.md)
for endpoints, inventories, surface availability, integration statistics, and the mass-balance
diagnostic.

!syntax parameters /UserObjects/DRIDNCorrosionUserObject

!syntax inputs /UserObjects/DRIDNCorrosionUserObject

!syntax children /UserObjects/DRIDNCorrosionUserObject
