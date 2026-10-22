# GFMR training: mesh, neutronics and core thermal-hydraulics of a FUNaK graphite MSR

This training builds a multiphysics model of the **Generic FUNaK-Fuelled Thermal Molten Salt
Reactor (GFMR)**, a public reference design published by Saltfoss Energy ApS
([specification](https://saltfossenergy.github.io/GFMR/), Apache 2.0,
DOI 10.5281/zenodo.23348200). The GFMR is a 150 MWth graphite-moderated reactor. Its fuel salt is
26.4 UF4 – 24.7 KF – 48.9 NaF mol% (2.12 % enriched) and flows upward through 379 channels in a
hexagonal graphite lattice. Nine hafnium control rods sit in Hastelloy N guide tubes inside fuel
channels.

The specification is a **reference problem without published results**. Numbers produced here
are code results, cross-checked between independent methods (Monte Carlo against diffusion, CFD
against correlations). Nothing here is validated against measurements.

## Stages

| Stage | Directory | What it does | Status |
|---|---|---|---|
| 1 | [01_mesh](01_mesh/README.md) | Reactor-module mesh generators: 2-D cross section, 3-D 60° CFD sector, homogenized full core | complete |
| 2 | [02_neutronics](02_neutronics/README.md) | OpenMC continuous-energy reference, 8-group cross sections, MOOSE multigroup diffusion on the CFD mesh | complete |
| 3 | [03_cfd](03_cfd/README.md) | 3-D conjugate k-epsilon CFD of the sector with the diffusion power map | complete |
| 4 | [04_chemistry](04_chemistry/README.md) | redox, guide-tube corrosion, I-135 / Xe-135 on the CFD flow, with a 0-D loop balance | complete (redox and corrosion from the 0-D loop: the 3-D fields are below the numerical floor) |
| 5 | [05_transients](05_transients/README.md) | 0.3 $ rod step: point kinetics with precursor drift and two-temperature feedback | one scenario (loss of flow, overcooling, loss of heat sink need a loop model) |

[THEORY.md](THEORY.md) is the theory manual for all five stages.

## Headline results

| Quantity | Value | Stage |
|---|---|---|
| k-eff, nominal (OpenMC) | 0.99510 ± 0.00014 (≈ 490 pcm subcritical as specified) | 2 |
| Total rod worth | 4844 ± 21 pcm (7.1 $) | 2 |
| Temperature coefficients, fuel / graphite | −2.7 / −4.1 pcm/K | 2 |
| β_eff / prompt generation time | 681 pcm / 0.61 ms | 2 |
| Diffusion vs OpenMC k | +115 pcm | 2 |
| Channel Δp vs SPEC | within 1–4 % | 3 |
| Graphite moderator average / maximum | 721 / 830 °C (reflector 859 °C) | 3 |
| Reactivity lost to precursor drift | 295 pcm (0.43 $) | 5 |
| 0.3 $ rod step, unprotected | peak 2.78 × nominal at 15 s, settles at 1.38 × | 5 |

## Shared data

[common/gfmr.py](common/gfmr.py) holds every specification value and its derived quantities. The
mesh, the OpenMC model and the MOOSE inputs all take their numbers from this one file. Running
`python common/gfmr.py` prints the consistency checks of the specification itself:

| Check | Result |
|---|---|
| Channel count from the decoded lattice | 82 inner + 288 outer + 9 rod = 379 (SPEC 379) |
| Mass flow from the zone velocities | 2294 kg/s (SPEC 2300) |
| Smooth-pipe friction vs SPEC pressure drops | 14.3 / 6.0 / 33.4 kPa vs 14.3 / 6.0 / 34 kPa |
| Core temperature rise at 2300 kg/s | 72.5 K |

The lattice map is not tabulated in the specification. It was decoded by sampling the colour of
the published top view at every lattice position, and it reproduces the stated counts exactly.

## Assumptions (not given in the specification)

| Quantity | Value used | Basis |
|---|---|---|
| Salt thermal conductivity | 0.8 W/m/K | typical NaF-KF-UF4 / NaF-ZrF4-UF4 fuel-salt data (0.6–1.0) |
| Meaning of the 650 °C fuel temperature | core-average salt temperature | gives T_in 613.8 °C, T_out 686.2 °C |
| Graphite heat capacity | 1900 J/kg/K | nuclear graphite near 800 °C |
| Axial reflectors / plena | none (vacuum top and bottom) | not in the specification |

## Requirements

- OpenPronghorn built with this branch (the C++ object `LinearFVTurbulentWallCHTBC` is new).
- For stage 2's OpenMC step only: OpenMC 0.16 with MPI and the ENDF/B-VIII.0 HDF5 library. The
  generated cross sections are committed in `02_neutronics/xs/`, so the MOOSE stages run without
  OpenMC.

## Computing note

On the INL workstation a login session is capped at 28 cores
(`systemctl show user-$UID.slice -p CPUQuotaPerSecUSec`), whatever the machine size. Run times
quoted in the READMEs are for 24 MPI ranks under that cap.
