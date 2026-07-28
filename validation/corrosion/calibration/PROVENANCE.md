# Calibration provenance

## Reproducible chain

The original calibration chain is now preserved entirely in this directory:

```text
local workbook
  -> ingest.py: source rows, units, classes, target relations, and weights
  -> model.py: 61-parameter effective corrosion/plating correlation
  -> calibrate.py: bounded robust least-squares fit with parameter priors
  -> reference/parameters.json: frozen fitted vector
  -> ../../../../data/corrosion_database.json: production copy plus one later correction
```

The tests enforce equality between the frozen local reference and the vendored
`validation/corrosion/data/parameters.json` copy. When a source audit corrected the MOD-F-07
experiment temperature from 650 C to 700 C, the workbook, frozen reference, and vendored parameters
were regenerated together.

## Experimental sources used by active constraints

The 38 active rows come from 14 source records:

| Source ID | Active rows | Experimental role |
|---|---:|---|
| ORNL-TM-1997 | 1 | First MSRE Hastelloy N surveillance group |
| ORNL-TM-3063 | 3 | Fourth MSRE surveillance group; depletion, salt Cr, and Cr diffusion |
| ORNL-4865 | 3 | MSRE metal transfer, noble-metal deposition ranking, and off-gas behavior |
| ORNL-TM-4188 | 8 | NCL-16 mass loss/gain, corrosion, inventories, FeF2 acceleration, and attack |
| ORNL-TM-4189 | 3 | 9.2-year Hastelloy N loop and external oxidation |
| ORNL-TM-4286 | 3 | Fluoride-loop stainless-steel corrosion and void depth |
| ORNL-TM-5782 | 2 | Type 316 compatibility in Li2BeF4 before and after Be addition |
| ORNL-TM-4271 | 1 | Hastelloy N mass transfer in NaBF4-NaF |
| JNM-316H-FLINAK-2022 | 3 | Flowing 316H/FLiNaK mass loss, gain, and Cr depletion |
| ORNL-FLIBE-316H-2022 | 1 | Flowing 316H/FLiBe upper-bound attack severity |
| ENERGIES-FLINAK-2021 | 2 | Multi-alloy FLiNaK attack depth |
| NPJ-HN-STRESS-FLINAK-2022 | 3 | Loaded/unloaded Hastelloy N and IGC in FLiNaK |
| NUC-SCI-TECH-316H-2024 | 1 | 316H IGC in impure chloride |
| ORNL-TM-6002 | 4 | Tellurium solubility, mass change, IGC, and redox threshold |

Full titles, URLs, extraction status, and source-quality classifications are retained in
`data/workbook_exports/source_index.csv`. The raw measurement-to-source mapping is in
`data/workbook_exports/detailed_measurements.csv` and the normalized model targets are in
`reference/targets.csv`.

## Current fit and source-audit changes

The frozen reference represents the current 61-parameter optimizer result, including both the
corrected 700 C MOD-F-07 experiment temperature and the source-audited M-041 direct target of
19.7 um at 1079 h. One OpenPronghorn production change must not be mistaken for an output of that
optimization:

1. `log_ncl16_cr_inventory_bonus` was added to the production database as a source-specific
   multiplier for the ORNL-TM-4188 NCL-16 chromium inventory. It is the 62nd production entry.

The workbook retained here is the input to the current fit, so a clean recalibration reproduces the
current 61 coefficients and 94.7% constraint pass fraction. The corrected M-041 measurement is
retained consistently in the workbook, readable workbook exports, frozen reference, and vendored
OpenPronghorn target table.

Any future coefficient update should first update the workbook and its CSV exports, rerun the
optimizer, record the environment manifest, and review the resulting production-database change.
