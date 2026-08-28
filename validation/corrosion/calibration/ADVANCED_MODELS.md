# MSTDB-TC and DRIDN advanced corrosion calibration

This directory contains the auditable Python reference implementation and the
curated calibration snapshot used to migrate the MSTDB-TC standard-state/Nernst
and Dynamic Redox-Inventory-Depletion Network (DRIDN) models into OpenPronghorn.
The authoritative input data and shared effective-model coefficients are from
`pr_corrosion` commit `bbc8b69c87bc1c90af0d1cf1484c599e9553b625`.

The checked calibration fixes two material issues in the recovered development
package:

- `M-035`, `M-036`, and `M-037` use 973.15 K; `M-041` is a direct 19.7 um
  observation at 1079 h.
- For non-fuel FeX2/Fe buffering, the Fe reaction is the exact identity
  `FeX2 + Fe <=> Fe + FeX2`, so its affinity is zero. The former Python dict
  representation collapsed duplicate species keys and created a fictitious Fe
  affinity.

Direct, range, upper-bound, and lower-bound rows now use consistent
relation-aware residuals. A range produces zero residual inside its reported
interval and log error to the nearest bound outside it. Midpoint metrics are
retained only for historical comparison.

All inferred inventory/geometry choices are frozen in the measurement-keyed
`results/advanced/advanced_case_context.csv`. Advanced scripts merge this table
strictly and refuse incomplete context. The MSTDB/DRIDN advanced inventory and
dynamic paths do not branch on response names, source identifiers, or
experiment-family labels when constructing their physical context.

Context revision `explicit-geometry-v2` also separates the physical initial
Cr/Fe/Ni dissolved inventories from the Nernst activity floor and freezes the
inventory fallback scale plus transient-redox, stress, fluoride-impurity, and
chloride activation switches. These fields are required even where a selected
MSRE/loop closure makes the fallback scale inactive. The Python reference and
MOOSE parity fixture therefore consume the same physical context without
reconstructing switches from redox or salt labels.

The static MSTDB closure applies `inventory_coupling_factor` only to dissolved
salt inventory, not coupon mass gain. DRIDN applies it to the released source
before all downstream inventory sinks and deposits. These distinct frozen model
laws are exercised by coupling and conservation regressions.

`build_context` also fails closed if the normalized calibration row omits its
measurement ID, salt/redox/position classes, temperature, temperature span,
flow factor, or nonnegative exposure time. The checked `targets.csv` supplies
these fields for every one of the 32 supported advanced-model constraints;
there are no production fallbacks for missing values.

The reduced-mechanistic model retained in the tables is a legacy comparator,
not part of this context refactor. Its `salt_cr_ppm` path still chooses an
MSRE/loop scaling from `experiment_family`/`source_id` and does not consume the
explicit `area_to_salt_mass_cm2_g` value. Its scores must therefore not be used
as evidence that the legacy comparator is metadata-independent.

## MSTDB-TC input

MSTDB-TC is external and is intentionally not redistributed. Set
`MSTDB_TC_DIR` to an authorized directory containing the exact checked pair:

```text
MSTDB-TC_V3.1_Fluorides_No_Func.dat
MSTDB-TC_V3.1_Chlorides_No_Func.dat
```

The checked run used MSTDB-TC 3.1. File hashes and record counts are in
`advanced_provenance.json` and
`results/advanced/mstdb_thermochemistry_checks.json`. Calibration and replay
entry points reject another version, modified bytes, ambiguous duplicate V3.1
files, and above-interval extrapolation unless a caller explicitly opts in to
the latter with a documented basis.

`data/advanced_corrosion_models.json` also binds the semantic subset of the
base `data/corrosion_database.json` used by the advanced C++ paths: shared
transport/redox parameters, material densities, and Cr/Fe/Ni electrochemical
properties. The source base file SHA-256 is frozen in both artifacts.

## Reproduce the checked calibration

Run from `validation/corrosion/calibration`:

```bash
python -m pip install -r requirements.txt  # exact checked dependency versions
export MSTDB_TC_DIR=/authorized/path/to/MSTDB-TC
export PYTHONPATH=src
export MPLBACKEND=Agg

python scripts/run_mechanistic_validation.py \
  --targets data/targets.csv \
  --case-context results/advanced/advanced_case_context.csv \
  --parameters data/parameters.json \
  --outdir generated/advanced \
  --prior-strength 0.25 \
  --max-nfev 5000

python scripts/run_thermochemical_validation.py \
  --mstdb-dir "$MSTDB_TC_DIR" \
  --targets data/targets.csv \
  --case-context results/advanced/advanced_case_context.csv \
  --parameters data/parameters.json \
  --reduced-parameters generated/advanced/results/mechanistic_parameters.json \
  --outdir generated/advanced \
  --prior-strength 0.25 \
  --max-nfev 5000

python scripts/run_dynamic_network_validation.py \
  --mstdb-dir "$MSTDB_TC_DIR" \
  --targets data/targets.csv \
  --case-context results/advanced/advanced_case_context.csv \
  --parameters data/parameters.json \
  --reduced-parameters generated/advanced/results/mechanistic_parameters.json \
  --thermochemical-parameters generated/advanced/results/thermochemical_parameters.json \
  --dynamic-parameters results/advanced/dynamic_network_parameters.initial.json \
  --outdir generated/advanced \
  --prior-strength 0.35 \
  --max-nfev 35 \
  --fit-integration-steps 24 \
  --integration-steps 120 \
  --refit \
  --skip-plots
```

The initial DRIDN JSON is a warm start only. Every fitted parameter group is
re-optimized against the corrected targets; all six checked stages converged.
Generated products remain under the ignored `generated/` directory so a replay
cannot overwrite the curated snapshot.

Run tests without proprietary database files:

```bash
env -u MSTDB_TC_DIR PYTHONPATH=src MPLBACKEND=Agg \
  python -m unittest discover -s tests -v
```

Run the full database-backed integration suite:

```bash
MSTDB_TC_DIR="$MSTDB_TC_DIR" PYTHONPATH=src MPLBACKEND=Agg \
  python -m unittest discover -s tests -v
```

## Checked results

The benchmark contains 32 supported active constraints and 26 direct/range
observations.

| Model | Midpoint log RMSE | Relation-aware log RMSE | Median midpoint factor | Constraint pass |
|---|---:|---:|---:|---:|
| MSTDB-TC | 0.181078 | 0.109852 | 1.029108 | 93.75% |
| DRIDN | 0.200596 | 0.092065 | 1.052158 | 96.875% |

Both advanced models are within a factor of two on all 26 direct/range rows.
DRIDN's maximum endpoint elemental-balance error is `1.32e-15`; 120 versus 240
integration steps changes predictions by at most `6.27e-8` relative. No fitted
thermochemical or DRIDN coefficient is on a bound. Five non-fitted DRIDN
coefficients are intentionally fixed at their lower admissible values; this is
a model-definition/identifiability choice, not an optimizer bound hit.

The portable calibration consumed by C++ is
`data/advanced_corrosion_models.json`. Its lower-case model-law identifier is:

```text
mstdb-nst-v2-fe-identity_dridn-v2-explicit-geometry
```

Curated scalar parameters, metrics, constraint checks, convergence results,
calibration-stage records, endpoint state audits, and the compact per-case
prediction table are in `results/advanced/`. `cpp_parity_cases.csv` freezes full
explicit contexts and Python endpoint outputs for nine representative validation
cases plus one synthetic zero-element boundary case.
Every checked-in CSV or JSON artifact in `results/advanced/` is SHA-256 bound
under `curated_outputs` in `advanced_provenance.json`; the explicit case-context
table is intentionally bound there and again under `inputs` because it is both
a replay input and a published calibration artifact.

## Scope and limitations

- These are in-sample calibrated endpoint benchmarks, not independent
  validation datasets.
- The modern-family exclusion run is a conditional thermochemical-layer
  sensitivity calculation. The shared effective model was calibrated on all
  rows, so the report records `independent_holdout: false`; it must not be cited
  as out-of-sample validation of the complete stack.
- The Python migration reference uses the recovered always-positive irreversible
  charge-transfer branch. The MOOSE wrapper defaults to that calibrated legacy
  mode for parity; the standalone C++ core defaults to affinity-gated
  dissolution. The affinity-gated mode is not validated by this calibration and
  requires a separate benchmark before production claims.
- The internal MSTDB parser evaluates standard-state Gibbs polynomials, Nernst
  terms, and pure-phase saturation; it does not run a native SUBQ equilibrium.
  Fitted activity closures remain necessary because most legacy cases lack full
  initial compositions.
- New transient species and redox measurements are needed to identify the five
  fixed DRIDN feedback/capture coefficients.
- Extrapolation to new salt families, alloys, temperature windows, geometries,
  or flow regimes requires additional validation.
