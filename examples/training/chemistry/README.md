# OpenPronghorn chemistry and radiochemistry training

A one-day, hands-on training on the molten salt reactor (MSR) chemistry capabilities of
OpenPronghorn: radiolysis, fission-product inventory, and wall corrosion and plating. The day builds
from a one-cell input with an analytic answer to a 2D, 24-rank MSFR primary loop that combines all
three. Each case is a validated input with a regression test, a plot script and a README with
exercises and solutions.

Slides: `doc/content/training/chemistry/index.md` (MooseDocs reveal.js; see *Building the slides*
below).

## Agenda

| time | module | cases | topics |
|---|---|---|---|
| 08:30–09:00 | **0. Setup** | | environment, running, reading CSV output, `plot.py` |
| 09:00–10:15 | **1. Radiolysis** | [01](01_f2_accumulation) [02](02_chloride_radiolysis) [03](03_cr_redox) [04](04_fuel_salt_redox_offgas) | G values, radical chains, metal redox and invariants, redox buffer, Henry-law off-gas |
| 10:15–10:30 | *break* | | |
| 10:30–12:00 | **2. Fission-product inventory** | [05](05_fp_buildup) [06](06_xenon_transient) [07](07_xenon_stripping) [08](08_noble_metal_partition) [09](09_power_history) | Bateman chains, Xe-135 and the xenon pit, sparging, noble-metal plating, power histories |
| 12:00–13:00 | *lunch* | | |
| 13:00–14:15 | **3. Corrosion and plating** | [10](10_corrosion_coupon_0d) [11](11_cr_depletion_two_block) [12](12_alloy_salt_redox_scan) [13](13_dridn_loop) | calibrated Butler-Volmer, Nernst limit, Cr-depletion layer, alloy/salt/redox, long-term loops |
| 14:15–14:30 | *break* | | |
| 14:30–16:30 | **4. Coupled, reactor-grade** | [14](14_irradiated_pool_corrosion) [15](15_core_channel_plating) [16](16_msfr_loop_capstone) | radiolysis-driven corrosion, turbulent heated channel, MSFR loop capstone |
| 16:30–17:00 | **Wrap-up** | | validation status, limits of the models, Q&A |

## The three actions

| action | creates | used in |
|---|---|---|
| `[MoltenSaltRadiolysis]` | radiolysis species, reaction network and G values from `data/msr_database.json` | 01–04, 14–16 |
| `[FissionProductInventory]` | nuclide inventories, yields, decay/capture links and removal pathways from `data/fission_product_database.json` | 05–09, 15, 16 |
| `[CorrosionPlating]` / `[CorrosionPlatingFlow]` | Butler-Volmer wall kinetics seeded from the calibrated corrosion correlation (`data/corrosion_database.json`) | 10–12, 14–16 |

All three build **linear finite-volume** systems (one per species). They are combined in one input
by listing every system in `[Problem] linear_sys_names` and in the executioner, and coupled through
shared variables (for example, corrosion releasing Cr directly into the radiolysis Cr(II)).

## Setup

```bash
ml use.moose moose-dev-openmpi
moose-dev-shell                       # interactive shell with the MOOSE environment
cd ~/projects/open_pronghorn
make -j 24                            # if open_pronghorn-opt is not built yet
export OPEN_PRONGHORN_EXE=$PWD/open_pronghorn-opt
cd examples/training/chemistry
```

Run a case and plot it:

```bash
cd 01_f2_accumulation
$OPEN_PRONGHORN_EXE -i f2_accumulation.i           # cases 01-14: serial, seconds
python plot.py                                      # writes ../figures/case01_*.png
cd ../16_msfr_loop_capstone
mpirun -n 24 $OPEN_PRONGHORN_EXE -i msfr_loop_capstone.i   # cases 15-16: parallel
```

Any input variable can be overridden on the command line (`dose_rate=2e6`,
`Executioner/end_time=3600`). The exercises use this throughout. `python plot_all.py` regenerates
every figure; it runs an input first if its output is missing.

## Case summary and validation

| case | runtime | validation check | result |
|---|---|---|---|
| 01 F2 accumulation | 2 s | analytic $S\,t$ | 3e-15 |
| 02 chloride radical chain | 6 s | dose exponents 1/2 and 1 | 0.50, 1.00 |
| 03 Cr redox | 14 s | Cr conservation, redox balance | 5e-7, 2e-6 |
| 04 fuel-salt redox + off-gas | 5 min | U conservation, redox balance | 2e-6, 2e-4 |
| 05 FP build-up | 3 s | analytic saturation / growth | 1e-13, 2e-7 |
| 06 xenon transient | 9 s | analytic equilibrium and post-scram | 1e-4, 2e-4 |
| 07 xenon stripping | 8 s | analytic with stripping; MSRE poisoning | 2.1% -> 0.19% |
| 08 noble-metal partition | 18 s | analytic salt inventory; MSRE salt fractions (ORNL-4865) | exact; < 1% in salt |
| 09 power history | 27 s | Cs-137 burnup monitor | 0.9987 |
| 10 corrosion coupon | 3 s | reference rate, Nernst limit, Faraday | 1.000, 0.9997, 1.000 |
| 11 Cr depletion | 3 s | Cr conservation, $\sqrt{Dt}$ growth | 1e-13, 1.06 |
| 12 alloy/salt/redox scan | 5 s | 7 measured loop/capsule rates | all within 2x (GM 1.18) |
| 13 DRIDN loop | < 1 s | NCL-16 salt Cr (~500 ppm), mass balance | 533 ppm, 4e-16 |
| 14 irradiated pool | 51 s | no corrosion without dose | 0 |
| 15 core channel (24 ranks) | 12 s | energy, Mo-99 balance; mesh study | 3e-4, 4e-4; < 3% |
| 16 MSFR loop (24 ranks) | 1.5 min | energy balance; Xe-135, Mo-99 vs 0D; mesh study | 6e-5; 1.00002, 1.0004 |

Every case is a CSVDiff regression test:

```bash
cd ~/projects/open_pronghorn
./run_tests -j 24 --test-root examples/training             # cases 01-03, 05-15
./run_tests -j 24 --test-root examples/training --heavy     # cases 04 and 16
```

## Building the slides

```bash
cd ~/projects/open_pronghorn/doc
./moosedocs.py build --config training_chemistry.yml --serve
```

Then open http://127.0.0.1:8000/training/chemistry/index.html.

## Known limitations (discussed in the wrap-up)

- Radiolysis G values and rate constants are literature values for idealized melts. Real salts
  carry impurities (moisture, oxides) that dominate the redox chemistry.
- The corrosion correlation is calibrated to the loop and capsule database. Outside its range (new
  alloys, chlorides, high dose) treat the rates as estimates. Coupling radiolysis to corrosion (case
  14) is a model hypothesis, not yet validated against irradiated-loop data.
- The fission-product database tracks selected chains with one-group cross sections. Noble-metal
  plating uses a mass-transfer coefficient; the surface chemistry (sticking, re-entrainment,
  graphite) is not modeled.
- Cases 15 and 16 are steady. The inventory equilibrium is exact when the loop time is much
  shorter than the half-lives. Power transients in a resolved loop need a transient executioner.
