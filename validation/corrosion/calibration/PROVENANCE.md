# Calibration provenance

## Reproducible chain

The current calibration chain is preserved entirely in this directory:

```text
local workbook
  -> ingest.py: source rows, units, classes, target relations, and weights
  -> model.py: 59-parameter effective corrosion/plating correlation
              + material-specific solid diffusivities from data/corrosion_database.json
  -> calibrate.py: bounded robust least-squares fit with parameter priors
  -> reference/parameters.json: frozen 59-parameter fitted vector
  -> ../../../../data/corrosion_database.json:
       59 fitted corrosion/plating parameters
       + one NCL-16 production correction
       + independent solid-diffusivity properties
```

The tests enforce equality between the frozen local reference and the vendored
`validation/corrosion/data/parameters.json` copy. When a source audit corrected the MOD-F-07
experiment temperature from 650 C to 700 C, the workbook, frozen reference, and vendored parameters
were regenerated together.

`reference_legacy_61/` preserves the output snapshot from the historical 61-parameter formulation,
in which chromium solid diffusivity was included in the fitted parameter vector. Reproducing that
historical optimizer exactly also requires the corresponding historical code revision; the snapshot
is retained for traceability rather than as a standalone rerunnable definition of the old fit.

## Experimental sources used by active constraints

The 37 active calibration rows come from 14 source records:

| Source ID                 | Active rows | Experimental role                                                            |
| ------------------------- | ----------: | ---------------------------------------------------------------------------- |
| ORNL-TM-1997              |           1 | First MSRE Hastelloy N surveillance group                                    |
| ORNL-TM-3063              |           2 | Fourth MSRE surveillance group; depletion and salt Cr                        |
| ORNL-4865                 |           3 | MSRE metal transfer, noble-metal deposition ranking, and off-gas behavior    |
| ORNL-TM-4188              |           8 | NCL-16 mass loss/gain, corrosion, inventories, FeF2 acceleration, and attack |
| ORNL-TM-4189              |           3 | 9.2-year Hastelloy N loop and external oxidation                             |
| ORNL-TM-4286              |           3 | Fluoride-loop stainless-steel corrosion and void depth                       |
| ORNL-TM-5782              |           2 | Type 316 compatibility in Li2BeF4 before and after Be addition               |
| ORNL-TM-4271              |           1 | Hastelloy N mass transfer in NaBF4-NaF                                       |
| JNM-316H-FLINAK-2022      |           3 | Flowing 316H/FLiNaK mass loss, gain, and Cr depletion                        |
| ORNL-FLIBE-316H-2022      |           1 | Flowing 316H/FLiBe upper-bound attack severity                               |
| ENERGIES-FLINAK-2021      |           2 | Multi-alloy FLiNaK attack depth                                              |
| NPJ-HN-STRESS-FLINAK-2022 |           3 | Loaded/unloaded Hastelloy N and IGC in FLiNaK                                |
| NUC-SCI-TECH-316H-2024    |           1 | 316H IGC in impure chloride                                                  |
| ORNL-TM-6002              |           4 | Tellurium solubility, mass change, IGC, and redox threshold                  |

M-005, the ORNL-TM-3063 chromium-diffusion value at 650 C, remains in the target set as a
`validation_only` quantitative comparison. It is not an active calibration constraint. In the
current framework it independently checks the Hastelloy N chromium-diffusivity material property
against the historical surveillance-derived value; it is not itself a direct tracer measurement.

Full titles, URLs, extraction status, and source-quality classifications are retained in
`data/workbook_exports/source_index.csv`. The raw measurement-to-source mapping is in
`data/workbook_exports/detailed_measurements.csv` and the normalized model targets are in
`reference/targets.csv`.

## Current fit and source-audit changes

The frozen reference represents the current 59-parameter optimizer result, including both the
corrected 700 C MOD-F-07 experiment temperature and the source-audited M-041 direct target of
19.7 um at 1079 h.

Chromium solid diffusivity is no longer included in the fitted parameter vector. Hastelloy N
chromium diffusivity is supplied independently from the DeVan 51Cr radiotracer-derived Arrhenius
correlation stored under `solid_diffusivities` in `data/corrosion_database.json`.

Solution-treated SUS316 chromium volume/lattice diffusivity is supplied independently from the
Mizouchi et al. (2004) 51Cr radiotracer correlation,
`D = 1.13e-3 exp[-234000/(R T)] cm^2/s`, measured over 888-1173 K. The production database stores
this direct property under `stainless_316`. The `solid_diffusivity_fallbacks` table explicitly maps
`stainless_304`, `stainless_304l`, `stainless_316h`, and `stainless_316l` to that SUS316 correlation
when no exact grade-specific chromium lattice-diffusion property is installed. These mappings are
engineering fallbacks and do not imply that the source measurement was performed in those grades.

The `generic_metal` chromium entry preserves the historical diffusion correlation as an explicitly
labeled legacy fallback.

One OpenPronghorn production correction remains outside the 59-parameter optimizer:

1. `log_ncl16_cr_inventory_bonus` is a source-specific multiplier for the ORNL-TM-4188 NCL-16
   chromium inventory. It is the 60th entry under production `calibrated_parameters`; it is not a
   fitted material diffusivity.

The current target accounting distinguishes calibration from validation. The 37 active calibration
rows produce 38 data residuals because the qualitative deposition-ranking row encodes two
inequalities. M-005 is withheld from fitting and scored separately as validation-only. Quantitative
summary metrics may also be reported across the combined direct/range and validation-only targets,
but that combined score is descriptive and must not be presented as an independent-validation
metric.

The workbook retained here is the input to the current fit, so a clean recalibration reproduces the
current 59 coefficients and frozen reference metrics. The corrected M-041 measurement and the
validation-only M-005 classification are retained consistently in the workbook-derived target
mapping, frozen reference outputs, and vendored OpenPronghorn validation data.

Any future coefficient update should first update the workbook and its CSV exports, rerun the
optimizer, record the environment manifest, review the resulting frozen-reference changes, and
promote any production-database update explicitly.
