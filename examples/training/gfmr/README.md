# GFMR training: a graphite-moderated MSR from mesh to transient

This hands-on training builds a full multiphysics model of one reactor with OpenPronghorn. The
reactor is the **Generic FUNaK-Fuelled Thermal Molten Salt Reactor (GFMR)**, a public reference
design by Saltfoss Energy ApS ([specification](https://saltfossenergy.github.io/GFMR/), Apache 2.0,
DOI 10.5281/zenodo.23348200).

The GFMR in one paragraph: **150 MWth**. The fuel salt is FUNaK (26.4 UF₄ – 24.7 KF – 48.9 NaF
mol %, 2.12 % enriched). It flows upward at 2300 kg/s through **379 channels** of 4 cm diameter in
a 14 cm hexagonal graphite lattice, 3.6 m tall, inside a graphite reflector of 3.5 m diameter. The
inner channels run faster (1.70 m/s) than the outer ones (1.03 m/s) to flatten the outlet
temperature. Nine hafnium rods in Hastelloy N guide tubes sit inside fuel channels, inserted 36 cm.

The five stages each feed the next:

```
 01_mesh ──► 02_neutronics ──► 03_cfd ──► 04_chemistry
  meshes     OpenMC + 8-group   3-D k-ε      redox, corrosion,
             diffusion: power   + graphite   I-135 / Xe-135
                 │
                 └──────── kinetics, feedback ──► 05_transients
                                                  0.3 $ rod step
```

- [THEORY.md](THEORY.md) is the theory manual: the equations, the methods and their limits for
  every stage.
- Each stage has its own `README.md` with the model, the expected results, the checks, and
  exercises.

## The stages

| Stage | What you build | New tools | Run time (24 ranks) | Time in class |
|---|---|---|---|---|
| [01 mesh](01_mesh/README.md) | 2-D cross section, exact 60° 3-D sector with every channel | Reactor-module mesh generators, `HexagonMeshTrimmer`, extrusion with block swaps | < 1 min each | 30 min |
| [02 neutronics](02_neutronics/README.md) | OpenMC reference; 8-group diffusion on the CFD mesh | OpenMC MGXS / IFP; `ArrayDiffusion`, `Eigenvalue` | 6 min (diffusion); OpenMC optional | 60 min |
| [03 core CFD](03_cfd/README.md) | 3-D conjugate CFD of the sector, driven by the diffusion power | linear FV `SIMPLE`, k-ε, `LinearFVTurbulentWallCHTBC` (new) | 19 min | 60 min |
| [04 chemistry](04_chemistry/README.md) | fuel-salt chemistry on the frozen flow, plus a 0-D loop | `[MoltenSaltRadiolysis]`, `[CorrosionPlatingFlow]`, `[FissionProductInventory]` | 20 min | 45 min |
| [05 transient](05_transients/README.md) | unprotected 0.3 $ rod step | scalar ODEs, `Transient` (BDF2) | seconds, 1 rank | 30 min |

A half-day session runs stages 1–3 live and stages 4–5 from their committed results. A full day
runs everything and adds the matching [THEORY.md](THEORY.md) section before each stage.

## Setup

```bash
ml use.moose moose-dev-openmpi
moose-dev-shell
cd ~/projects/open_pronghorn
make -j 24                                   # builds LinearFVTurbulentWallCHTBC as well
export EXE=$PWD/open_pronghorn-opt
cd examples/training/gfmr
python common/gfmr.py                        # the specification checks (below)
```

OpenMC is **not** needed: the cross sections it produced are committed in
`02_neutronics/xs/`. To rerun OpenMC itself you need OpenMC 0.16 with MPI and the ENDF/B-VIII.0
HDF5 library (see stage 2).

## Running the whole chain

Each stage reads a file the previous stage writes. The three hand-off files are kept in the
working copy so that you can start at any stage. They are generated, not tracked:
if you cloned fresh, run the stages in order.

| Step | Command | Writes | Read by |
|---|---|---|---|
| 1 | `cd 01_mesh && $EXE -i gfmr_mesh_3d.i --mesh-only gfmr_sector_3d_cfd.e` | `gfmr_sector_3d_cfd.e` | 2, 3 |
| 2 | `cd ../02_neutronics && mpirun -n 24 $EXE -i gfmr_diffusion_sector.i -w Outputs/file_base=sector_power` | `sector_power.e` | 3 |
| 2 (optional) | build the rod meshes (stage-1 README), then `./run_rod_cases.sh 24` | rod worths | — |
| 3 | `cd ../03_cfd && mpirun -n 24 $EXE -i gfmr_core_cfd.i` | `gfmr_core_cfd_out.e` | 4 |
| 4 | `cd ../04_chemistry && python loop_balance.py && mpirun -n 24 $EXE -i gfmr_chemistry.i` | results | — |
| 5 | `cd ../05_transients && $EXE -i gfmr_rod_step.i` | results | — |

The `-w` in step 2 is needed because `xs/gfmr_xs.i` defines constants for every region and each
input uses only some of them. Committed result tables are in `03_cfd/results/`,
`04_chemistry/results/` and `05_transients/results/`, so the reported numbers can be checked
without rerunning.

## Headline results

| Quantity | Value | Stage |
|---|---|---|
| k-eff, nominal (OpenMC) | 0.99510 ± 0.00014 | 2 |
| Total rod worth | 4844 ± 21 pcm (7.1 $) | 2 |
| Temperature coefficients, fuel / graphite | −2.7 / −4.1 pcm/K | 2 |
| β_eff / prompt generation time | 681 pcm / 0.61 ms | 2 |
| Heat deposited directly in the graphite | 7.9 % | 2 |
| Channel Δp vs specification | within 1–4 % | 3 |
| Graphite moderator average / maximum | 721 / 830 °C (reflector 859 °C) | 3 |
| Mixed outlet of the inner / outer channels | 685.4 / 686.8 °C | 3 |
| Reductant needed to hold U(III)/U at 1 % | 0.54 mol eq/day | 4 |
| Reactivity lost to precursor drift | 295 pcm (0.43 $) | 5 |
| 0.3 $ rod step, unprotected | peak 2.78 × nominal at 15 s, settles at 1.38 × | 5 |

### The thread to pull on: is the specified core critical?

Three results from different stages answer it together, and they are the main discussion of
the day:

1. The static OpenMC core is **490 pcm subcritical** at the specified temperatures (stage 2).
2. Circulating fuel loses **295 pcm** of delayed neutrons to the loop (stage 5).
3. The CFD graphite runs **80 K cooler** than the 800 °C the specification assumes, which is
   worth about **+330 pcm** with the graphite coefficient (stage 3).

So as specified the core is about 450 pcm short of critical at operating conditions. The
specification gives no excess-reactivity target, and the enrichment was not adjusted. Each number
can be traced to its stage README.

## What each check really is

The GFMR has no measured data. Every check is therefore labelled **V** (verification: analytic,
or exact conservation), **C** (calibration: reproduces what it was tuned to) or **B**
(code-to-code benchmark). None of them is a validation.

| Stage | Checks | Key result |
|---|---|---|
| spec | V: channel count, mass flow, Δp from friction | 379 channels; 2294 vs 2300 kg/s; Δp within 2 % |
| 1 | V: block areas against the analytic geometry | exact to the mesh tolerance |
| 2 | B: diffusion vs OpenMC (k, power split, axial shape); **C**: rod worth | +115 pcm; split within 0.5 points; axial peak 1.557 both; Hf factor 0.602 fitted |
| 3 | V: energy balance, interface fluxes; B: Δp vs the specification | 4.5×10⁻⁵; fluxes to 10⁻⁵; Δp within 1–4 % |
| C++ | V: salt pipe vs Gnielinski (correlation, not measurement) | Nu ratio 0.93 |
| 4 | V: 3-D I-135 / Xe-135 core rise vs the 0-D balance | within 4–5 % |
| 5 | V: null transient; final state vs analytic feedback balance | steady to 8 digits; 1.3789 vs 1.3797 |

Two results are **not** checks:

- **The Hf rod worth in diffusion is calibrated.** Diffusion overestimates a black absorber,
  so its cross sections are scaled by 0.602 to match OpenMC.
- **The 3-D redox and corrosion fields are below the numerical floor.** The advection error
  (10⁻⁵) is larger than their real change in one 3-second core pass, so stage 4 takes those
  results from the 0-D loop instead.

## Shared data

[common/gfmr.py](common/gfmr.py) holds every specification value and its derived quantities. The
mesh, the OpenMC model and the MOOSE inputs all take their numbers from this one file.
`python common/gfmr.py` checks the specification against itself:

| Check | Result |
|---|---|
| Channel count from the decoded lattice | 82 inner + 288 outer + 9 rod = 379 (SPEC 379) |
| Mass flow from the zone velocities | 2294 kg/s (SPEC 2300) |
| Smooth-pipe friction vs SPEC pressure drops | 14.3 / 6.0 / 33.4 kPa vs 14.3 / 6.0 / 34 kPa |
| Core temperature rise at 2300 kg/s | 72.5 K |

The specification does not tabulate the lattice map. It was decoded by sampling the colour of the
published top view at every lattice position, and it reproduces the stated channel counts exactly.

## Assumptions (not given in the specification)

| Quantity | Value used | Stage | Basis |
|---|---|---|---|
| Salt thermal conductivity | 0.8 W/m/K | 3 | typical fluoride fuel salts (0.6–1.0) |
| Meaning of the 650 °C fuel temperature | core-average salt temperature | 3 | gives T_in 613.8 °C, T_out 686.2 °C |
| Axial reflectors, plena | none (vacuum top and bottom) | 2 | not in the specification |
| Graphite heat capacity | 1900 J/kg/K | 5 | nuclear graphite near 800 °C |
| External loop | ex-core volume 2 × core salt (9.1 s loop) | 4, 5 | MSRE-like; the specification has no loop |
| Redox control, Xe stripping | U(III)/U 1 %, 1 day; 10 % Xe per pass | 4 | MSRE practice |
| Corrosion data | `fluoride_fuel` / `hastelloy_n` classes | 4 | no FUNaK-specific data exist |

## Known limitations

- **Neutronics.** The diffusion rod worth is calibrated (above). The power map is frozen: there is
  no temperature feedback in the steady CFD (one-way coupling).
- **CFD.** Standard k-ε with wall functions, at 0.93 × Gnielinski on the salt pipe. The graphite
  heat follows the E > 5.53 keV flux rather than a gamma transport solution. The reflector surface
  is adiabatic, which makes its 859 °C maximum a sensitivity rather than a prediction.
- **Chemistry.** The loop is assumed. Redox and corrosion come from the 0-D loop.
- **Transient.** Point kinetics with lumped temperatures and a fixed inlet temperature. Loss of
  flow, overcooling and loss of heat sink are **not modelled**: they need a loop thermal model,
  which the specification does not define. Stage 5 shows how the same input extends to them.

## Regression test

The new C++ object has a heavy CSVDiff test:

```bash
cd ~/projects/open_pronghorn
./run_tests -p 4 --heavy --re=salt_pipe_heat_transfer
```

## Computing note

On the INL workstation a login session is capped at 28 cores
(`systemctl show user-$UID.slice -p CPUQuotaPerSecUSec`), whatever the machine size. The run times
above are for 24 MPI ranks under that cap.
