# MoltenSaltCorrosionData

`MoltenSaltCorrosionDatabase` reads and queries the molten salt corrosion database
(`data/corrosion_database.json`) at run time. The database bundles three kinds of information:

- `elements`: the per-cation mechanistic properties (valence, molar mass, salt diffusivity, standard
  electrode potential, charge-transfer coefficients, reference concentration) consumed by the
  Nernst-Planck and Butler-Volmer objects,
- `densities_g_cm3` and `alloy_cr_wt_frac`: the engineering material tables from the reference model
  (`msr_corrosion_bv/chemistry.py`) used by the Faradaic conversions and the ported correlation,
- `calibrated_parameters`: the verbatim 61 fitted parameters of the reference effective Butler-Volmer
  correlation (`results/parameters.json`), keyed exactly as in that model so the C++ port reproduces
  it term for term.

Element and material lookups are case-insensitive and fall back to a `generic_metal` entry. Point the
`database` parameter of [CorrosionPlating](CorrosionPlatingAction.md) at your own JSON file with the
same structure to change the chemistry without recompiling.
