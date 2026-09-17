# Corrosion and plating calibration

This directory is the complete, in-repository calibration package for the effective
molten-salt corrosion and plating model used by OpenPronghorn. It replaces the previous
dependency on the external `msr_corrosion_plating_model` working directory.

Everything needed by the calculation is local:

- `data/msr_corrosion_plating_validation_data.xlsx` is the source workbook.
- `data/workbook_exports/` contains readable CSV exports of the workbook's provenance,
  measurements, cases, variable definitions, completeness assessment, and digitization queue.
- `src/msr_corrosion_bv/` contains data ingestion, unit conversion, model, optimizer, and plotting
  code.
- reference/ is the frozen output of the current 59-parameter calibration. reference_legacy_61/ preserves the historical 61-parameter calibration outputs.
- `scripts/run_calibration.py` reruns the fit and checks it against the frozen reference.
- `tests/test_calibration.py` tests workbook ingestion, source traceability, numerical
  reproducibility, and the relationship between the 59 fitted parameters and the production
  OpenPronghorn database.

No network access or files outside this repository are required to rerun the calibration.

## Recommended local run

From the repository root:

```bash
conda run --no-capture-output -n full_anaconda \
  python validation/corrosion/calibration/scripts/run_calibration.py --no-plots
```

The runner resolves the workbook and reference files relative to its own location, so it does not
depend on the current working directory. By default it writes to
`validation/corrosion/calibration/generated/latest`, compares the new fit with
`reference/parameters.json`, and exits nonzero if the fitted vector drifts beyond the documented
tolerance.

To regenerate the diagnostic figures as well:

```bash
conda run --no-capture-output -n full_anaconda \
  python validation/corrosion/calibration/scripts/run_calibration.py
```

To use another Python environment:

```bash
python -m pip install -r validation/corrosion/calibration/requirements.txt
python validation/corrosion/calibration/scripts/run_calibration.py --no-plots
```

## Tests

```bash
conda run --no-capture-output -n full_anaconda \
  python -m unittest discover \
  -s validation/corrosion/calibration/tests \
  -p 'test_*.py' -v
```

The integration test performs a fresh nonlinear least-squares calibration. Small floating-point
differences across SciPy/BLAS versions are expected, so coefficient comparison uses an absolute
tolerance of `1e-3`. This is conservative relative to the observed cross-environment maximum
difference of approximately `3.2e-5`.

## What is calibrated

The workbook contains 76 source-mapped cases and 43 detailed measurement rows. Thirty-seven rows
are active calibration constraints:

- 25 direct numerical targets;
- 3 numerical ranges;
- 7 upper bounds;
- 1 lower bound; and
- 1 qualitative deposition ranking.

M-005, the historical chromium-diffusion target from ORNL-TM-3063, is retained as a quantitative
validation-only target and is not used to fit the corrosion/plating coefficients.

The ranking row produces two inequality residuals, so the 37 active rows produce 38 data residuals.
The objective also contains one regularizing prior residual for each of the 59 fitted parameters.

See [MODEL.md](MODEL.md) for the equations and objective and [PROVENANCE.md](PROVENANCE.md) for the
source inventory and the distinction between the original fit and later OpenPronghorn corrections.

## Frozen reference versus production database

`reference/parameters.json` contains exactly the current 59 fitted corrosion/plating coefficients.
The production file `data/corrosion_database.json` contains those same 59 values plus
`log_ncl16_cr_inventory_bonus`, a source-specific production correction, for 60 entries under
`calibrated_parameters`.

Chromium solid diffusivity is no longer part of that fitted parameter vector. Material-specific
diffusivities are stored separately under `solid_diffusivities` in the production database.
Hastelloy N chromium diffusivity uses the DeVan tracer-derived Arrhenius correlation, while
`generic_metal` retains the historical 61-parameter diffusion correlation as an explicitly labeled
legacy fallback.

`reference_legacy_61/` preserves the historical 61-parameter output snapshot. Reproducing that
historical optimizer exactly requires the corresponding historical code revision; the snapshot
alone is retained for traceability.

The calibration runner never overwrites the production database. Promoting a new fit must be an
explicit, reviewed action.