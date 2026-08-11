# Molten salt corrosion and plating validation

This directory reuses the complete validation set of the calibrated reference model
`msr_corrosion_plating_model` (the effective Butler-Volmer correlation fit to 76 molten salt
corrosion/plating cases and 43 detailed measurement targets) to validate the open_pronghorn
corrosion and plating framework.

## Reference dataset (`data/`)

Vendored verbatim from the reference model:

| File | Contents |
| :--- | :--- |
| `case_features.csv` | 76 cases (material/salt/redox/flow/temperature/time features) |
| `targets.csv` | 43 measurement targets with fit roles and bounds |
| `parameters.json` | the 59 fitted corrosion/plating parameters |
| `metrics.json` | the reference fit metrics |
| `validation_predictions.csv` | reference predictions for every target |
| `case_predictions_all_76_cases.csv` | reference predictions for every case |
| `ncl16_simplified_loop_simulation.csv` | the NCL-16 multi-segment loop simulation |

The frozen 59-parameter calibration has 28 direct/range calibration targets, with median factor
error 1.1166, 92.9 % (26/28) within a factor of 2, and 100 % (28/28) within a factor of 5.
The current production-model target matrix additionally applies the separate NCL-16 chromium
inventory correction to M-014, giving 27/28 direct/range targets within a factor of 2 and 28/28
within a factor of 5. M-005 is retained separately as a validation-only chromium-diffusivity target
and is within a factor of 2 of the independently supplied Hastelloy N diffusivity. Across the 29
quantitative targets reported by the frozen calibration metrics, 93.1 % are within a factor of 2
and 100 % are within a factor of 5.

## Layer 1 — term-for-term correlation reproduction (C++ unit tests)

`MoltenSaltCorrosionModel` is a faithful C++ port of the reference `MoltenSaltBVModel`. The unit
tests `unit/src/MoltenSaltCorrosionModelTest.C` and `unit/src/MoltenSaltCorrosionDataTest.C`
reproduce, to a relative tolerance of 1e-9:

- the corrosion rate, deposition rate, IGC depth, mass loss/gain, salt Cr ppm, Cr diffusivity,
  effective overpotential and every other `predict_response` branch across representative cases,
- the engineering tables (valences, molar masses, densities, chromium fractions), the current
  calibrated-parameter database, and the material-specific chromium solid-diffusivity properties,
- the NCL-16 loop simulation (final salt Cr ppm and cumulative dissolution/deposition).

Run with:

```
cd unit && make -j && ./open_pronghorn-unit-opt --gtest_filter='MoltenSaltCorrosion*'
```

## Layer 2 — end-to-end MOOSE reproduction (`run_corrosion_validation.py`)

`run_corrosion_validation.py` runs the mechanistic `CorrosionPlating` action (a 0D salt-only cell,
`case_0d.i`) for every one of the 76 cases and confirms the Butler-Volmer boundary current
Faradaically reproduces the reference dissolution rate. The exchange current is seeded from the
calibrated rate, so reproduction is essentially exact to floating-point precision (approximately
`1e-14` relative error in current validation runs). The script checks every case against a
configurable reproduction tolerance (default `1e-4` relative error) and reports the observed maximum
relative error.
```
python3 run_corrosion_validation.py            # all 76 cases
python3 run_corrosion_validation.py --limit 8  # quick subset
```

On local MPI/PETSc builds where the MOOSE process writes the `INITIAL` CSV output but lingers during
finalization, use `--case-timeout` to recover the completed CSV output case by case:

```
python3 run_corrosion_validation.py --case-timeout 3
```

This closes the loop from the validated 0D correlation to the spatial mechanistic framework: the same
calibrated kinetics that match the experimental targets drive the Nernst-Planck transport,
current-continuity potential and Butler-Volmer electrode reactions of the full open_pronghorn model.

## Experimental target matrix (`experimental_target_matrix.py`)

`run_corrosion_validation.py` is a code-reproduction check, not by itself an experimental validation.
The experimental validation check is the target matrix:

```
python3 experimental_target_matrix.py --output experimental_target_matrix.csv
```

It scores the 43 measurement rows in `validation_predictions.csv` against their experimental
relations. The 37 active calibration constraints exclude input-only rows, the auxiliary row, and
the validation-only M-005 chromium-diffusivity target. Direct targets pass when the model is within
a factor of 2, range targets pass when the prediction lies inside the extracted range, and
upper/lower targets pass against the corresponding bound with a small explicit digitization
tolerance. Validation-only quantitative targets are scored separately and do not enter the
calibration constraint count.

## Layer 3 — linear finite-volume Butler-Volmer wall (flow-coupled path)

The corrosion framework also provides a linear finite-volume path (`[CorrosionPlatingFlow]` +
`CorrosionLinearFVButlerVolmerBC`) so corrosion products are passive scalars in the same segregated
flow solve as the radiolysis and energy — the way corrosion couples into a flowing MSR
(`examples/flowing_msr_corrosion`). The linear-FV Butler-Volmer wall uses the identical kinetics and
exchange-current seeding as the finite-element path, so it reproduces the same calibrated rates:
`test/tests/corrosion/linearfv_wall` confirms the 304L hot-fluoride-loop case (16.821 um/y) in a
well-mixed cell, matching the finite-element `action_salt_only` result.
