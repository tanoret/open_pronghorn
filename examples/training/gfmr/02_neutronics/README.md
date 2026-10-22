# Stage 2: GFMR neutronics — OpenMC reference and multigroup diffusion on the CFD mesh

The neutronics is built in two steps:

1. A **continuous-energy Monte Carlo** model (OpenMC 0.16, ENDF/B-VIII.0) of the full 3-D core.
   It is the reference for k-eff, rod worth, reactivity coefficients, kinetics parameters and the
   power distribution, and it produces the 8-group cross sections.
2. A **MOOSE 8-group diffusion** model with those cross sections, solved directly on the 3-D CFD
   sector mesh. Every salt channel, guide tube, absorber and graphite region has its own
   constants, so the power density lands on the CFD cells without any mapping.

## Files

| File | Purpose |
|---|---|
| `openmc/build_model.py` | OpenMC geometry: hexagonal lattice of unit-cell universes, SPEC materials |
| `openmc/run_reference.py` | the reference cases (below), cell-domain MGXS, IFP kinetics |
| `openmc/postprocess.py` | summary, channel / axial power, multigroup constants → `xs/gfmr_8g.json` |
| `make_diffusion_input.py` | `xs/gfmr_8g.json` → `xs/gfmr_xs.i` (HIT include, SI units) |
| `gfmr_diffusion_rz.i` | 2-D R-Z diffusion, homogenized lattice (quick check) |
| `gfmr_diffusion_sector.i` | 3-D heterogeneous diffusion on the CFD sector mesh |
| `run_rod_cases.sh` | sector diffusion at nominal, rods-out and rods-in |
| `xs/` | committed cross sections and OpenMC results, so the MOOSE runs need no OpenMC |

## OpenMC reference

Each case uses 2×10⁵ particles × 200 active batches (4×10⁷ histories). The source and statepoint
files are in `/data/tanome/gfmr_runs/openmc_reference` and are not committed.

| Case | k-eff |
|---|---|
| nominal (rods inserted 36 cm, fuel 650 °C, graphite 800 °C) | 0.99510 ± 0.00014 |
| rods withdrawn | 0.99525 ± 0.00015 |
| rods fully inserted (3.6 m) | 0.94948 ± 0.00013 |
| fuel +100 K | 0.99241 ± 0.00014 |
| graphite +100 K | 0.99103 ± 0.00015 |

| Derived quantity | Value |
|---|---|
| Total rod worth (9 rods) | 4844 ± 21 pcm (7.1 $) |
| Worth of the nominal 36 cm insertion | 16 ± 21 pcm (not significant) |
| Fuel temperature coefficient | −2.7 ± 0.2 pcm/K |
| Graphite temperature coefficient | −4.1 ± 0.2 pcm/K |
| β_eff (IFP) | 681 ± 5 pcm |
| Prompt generation time (IFP) | 0.605 ms |
| Heat deposited in the graphite (coupled n-γ) | 7.9 % |
| Power split inner / outer / rod channels | 31.2 / 66.1 / 2.7 % |
| Channel / axial power peaking | 1.57 / 1.56 |

### Findings to discuss

- **The core as specified is about 490 pcm subcritical at nominal conditions.** The specification
  gives no excess-reactivity target, and the enrichment has not been adjusted. Reasons to
  discuss: no axial reflector or plena (vacuum at the top and bottom of the graphite), and 1 appm
  boron in the graphite.
- **The rods at 36 cm are worth almost nothing.** The axial flux at the top of a bare core is
  low. The rods' worth comes from deeper insertion.
- Both temperature coefficients are negative. The graphite coefficient dominates.

### Making it run fast

Three changes were needed to make the campaign affordable:

- **Lattice geometry.** Explicit hexagonal prisms with the reflector defined as "outside all
  397 prisms" cost 70 s per batch. A `HexLattice` of universes runs 12× faster.
- **MPI ranks instead of threads.** Many MGXS tallies make threads contend for the tally
  arrays. With 7 ranks × 4 threads the run is 3× faster than with 28 threads.
- **Source sampling.** Only ~7 % of the volume is salt, so the uniform initial source needs
  `source_rejection_fraction = 0.01`.

## Diffusion on the sector mesh

`gfmr_diffusion_sector.i` solves the 8-group eigenvalue problem with CFEM (`ArrayDiffusion`,
`ArrayReaction`, `ArrayVacuumBC`). It has 3.6×10⁶ unknowns and takes 6 min on 24 ranks.
Each mesh block takes the constants of the matching OpenMC cell region (`regions` in
`gfmr_8g.json`). The vacuum condition is Marshak (α = 1/2) on the top, bottom and reflector side;
the symmetry planes are reflective. The eigenvector is normalized to 25 MW (150 MW / 6).

### The black-absorber calibration

Diffusion theory overestimates the absorption of a nearly black absorber. With the raw OpenMC
constants the total rod worth comes out 21 % high (5868 vs 4844 pcm). The Hf removal matrix is
therefore scaled by **0.602** (`make_diffusion_input.py --rod-factor`), which reproduces the
OpenMC total worth. **This is a calibration, not a prediction.**

### Comparison with OpenMC

On the final mesh (with the 3 mm wall ring of the CFD); the coarser mesh without the ring gives
the same k within 1 pcm.

| Quantity | Diffusion | OpenMC |
|---|---|---|
| k-eff nominal | 0.99624 (+115 pcm) | 0.99510 |
| k-eff rods out / rods in | 0.99670 / 0.95078 (+146 / +144 pcm) | 0.99525 / 0.94948 |
| Total rod worth | 4845 pcm (calibrated) | 4844 ± 21 pcm |
| Power split inner / outer / rod | 30.8 / 66.6 / 2.6 % | 31.2 / 66.1 / 2.7 % |
| Axial peaking | 1.557 | 1.557 |
| Worth of the 36 cm insertion | 46 pcm | 16 ± 21 pcm |
| Peak salt power density | 212 MW/m³ | — |

The k bias of about +130 pcm is the same with the rods in and out, so it comes from the
homogenization of the lattice and not from the rods.

## How to run

```bash
# OpenMC (optional: the constants are committed)
cd openmc
python run_reference.py --out RUNDIR --ranks 7 --threads 4
python postprocess.py --ref RUNDIR --out ../xs
cd ..
python make_diffusion_input.py --rod-factor 0.602
# diffusion
mpirun -n 24 open_pronghorn-opt -i gfmr_diffusion_rz.i -w
./run_rod_cases.sh 24
```

The `-w` flag is needed because `xs/gfmr_xs.i` defines constants for every region, and each input
uses only some of them.

## Exercises

1. Compare the R-Z and sector k-eff. What does the R-Z homogenization miss?
2. Rerun the rods-in case with `--rod-factor 1.0` and explain the difference in worth.
3. From the two temperature coefficients and β_eff, estimate the prompt temperature rise that
   would compensate a 100 pcm reactivity insertion.
