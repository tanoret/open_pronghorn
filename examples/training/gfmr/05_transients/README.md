# Stage 5: transients — NOT DONE

This stage is planned but **not implemented or validated**.

## Plan

Use the homogenized full-core mesh (`../01_mesh/gfmr_homog_3d.e`) with:

- time-dependent 8-group diffusion (`ArrayTimeDerivative` with 1/v from `xs/gfmr_xs.i`);
- 6 delayed-neutron precursor groups advected with the salt (β_eff 681 pcm, decay constants in
  `xs/gfmr_xs.i`), including the loss of precursors to the external loop;
- porous-medium salt and graphite energy, with feedback from the OpenMC coefficients
  (fuel −2.7 pcm/K, graphite −4.1 pcm/K);
- a lumped loop with a transit time (not in the specification; must be assumed).

Scenarios: rod withdrawal (use the full-insertion worth: the 36 cm insertion is worth
< 50 pcm), loss of flow, overcooling, loss of heat sink.

## Inputs already available

- Kinetics: β_eff and generation time from OpenMC IFP (`../02_neutronics/xs/openmc_summary.json`).
- Rod worth: 4844 pcm total, with diffusion calibrated to it (`--rod-factor 0.602`).
- Steady-state temperatures: stage 3.
