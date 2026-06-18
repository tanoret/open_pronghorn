# MoltenSaltCorrosionModel

`MoltenSaltCorrosionModel` is a faithful C++ port of the validated effective Butler-Volmer corrosion,
mass-transfer and plating correlation `MoltenSaltBVModel` from the reference model
`msr_corrosion_plating_model`. It reads the calibrated parameters and engineering tables from
[MoltenSaltCorrosionData](MoltenSaltCorrosionData.md) (the JSON corrosion database) and reproduces the
reference predictions term for term.

The model supplies:

- the corrosion rate (harmonic mean of a kinetic and a transport-limited branch), deposition rate,
  IGC depth, mass loss/gain, salt chromium ppm, chromium diffusivity, effective overpotential and the
  full `predict_response` dispatch,
- the Faradaic conversions between current density and penetration rate (shared with the
  Butler-Volmer objects through `CorrosionChemistry`),
- the seeds for the mechanistic exchange current density and solid-state chromium diffusivity used by
  the [CorrosionPlating](CorrosionPlatingAction.md) action,
- the simplified multi-segment loop simulation (`simulateLoop`).

It is exercised directly by the unit tests `unit/src/MoltenSaltCorrosionModelTest.C` and
`unit/src/MoltenSaltCorrosionDataTest.C`, which reproduce the reference case and target predictions
and the NCL-16 loop to a relative tolerance of 1e-9. See `validation/corrosion` for the end-to-end
MOOSE reproduction of all 76 reference cases.
