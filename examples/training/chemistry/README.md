# OpenPronghorn chemistry and radiochemistry training

A hands-on training on the molten salt reactor (MSR) chemistry capabilities of OpenPronghorn:
fuel-salt redox and radiolysis, the fission-product inventory, and wall corrosion and plating. It
has **five cases of increasing complexity**: a single-cell redox problem, then a month-long xenon
transient, then corrosion in a real loop seen by four models, then the **MSRE** compared with its
measured chemistry, and finally the **MSFR** primary circuit on the Virtual Test Bed geometry.

- [THEORY.md](THEORY.md) is the theory and usage manual: the equations, parameters and pitfalls of
  every model.
- Each case has its own `README.md` with the physics, a step-by-step walk through the input, the
  expected results, the checks, and exercises with solutions.
- Slides: `doc/content/training/chemistry/index.md` (MooseDocs reveal.js; see *Building the slides*).

## The five cases

| case | geometry | new in this case | physics | time |
|---|---|---|---|---|
| [01 fuel-salt redox](01_fuel_salt_redox) | 0D, serial, 1.5 min | `[MoltenSaltRadiolysis]`, hand-added kernels | fission oxidation, the U(III) buffer, F2 to the cover gas | 45 min |
| [02 xenon and power history](02_xenon_power_history) | 0D, serial, 30 s | `[FissionProductInventory]`, functors of time | Bateman chains, Xe-135 worth, stripping, burnup | 40 min |
| [03 one loop, four corrosion models](03_corrosion_models) | 0D/1D, serial, 20 s | `[CorrosionPlating]`, `[CorrosionPlatingFlow]`, DRIDN | Butler-Volmer, the Nernst limit, hot/cold transfer, model comparison | 50 min |
| [04 MSRE fuel loop](04_msre_loop) | multi-region, serial, 16 s | all three physics in one input | Xe with spray, graphite and bubbles; noble metals; corrosion; against the MSRE record | 50 min |
| [05 MSFR primary circuit](05_msfr_reactor) | 2D-RZ, 8 ranks, ~7 min | CFD (SIMPLE, k-epsilon, energy), flow restart | everything at reactor scale | 60 min |

A suggested half-day plan: setup (20 min); cases 01–02 (1 h 30); break; cases 03–05 (2 h 30);
wrap-up (20 min). For a full day, add the theory sections of [THEORY.md](THEORY.md) before each case.

## Setup

```bash
ml use.moose moose-dev-openmpi
moose-dev-shell                       # interactive shell with the MOOSE environment
cd ~/projects/open_pronghorn
make -j 24                            # if open_pronghorn-opt is not built yet
export OPEN_PRONGHORN_EXE=$PWD/open_pronghorn-opt
cd examples/training/chemistry
```

Run a case and plot it (matplotlib is available in the MOOSE environment):

```bash
cd 01_fuel_salt_redox
$OPEN_PRONGHORN_EXE -i fuel_salt_redox.i
python plot.py                                       # writes ../figures/case01_*.png
cd ../05_msfr_reactor
mpirun -n 8 $OPEN_PRONGHORN_EXE -i msfr_flow.i       # step 1: flow and energy (~1 min)
mpirun -n 8 $OPEN_PRONGHORN_EXE -i msfr_chemistry.i  # step 2: chemistry (~6 min)
```

You can override any input variable on the command line, for example `P_fission=4e7` or
`Executioner/end_time=3600`. The exercises use this throughout. `python plot_all.py`
regenerates every figure, and runs an input first if its output is missing.

## What each check really is

Every check is labelled **V** (verification: analytic or exact conservation), **C** (calibration:
reproduces data the model was fitted to), **B** (code-to-code benchmark) or **VAL** (validation:
independent measured data).

| case | checks | key result |
|---|---|---|
| 01 | V: U, F, redox balances; analytic buffer line | buffer spent at day 339 (analytic 342) |
| 02 | V: analytic xenon equilibrium, Cs-137 burnup | -1066 pcm; 0.9987 (0.967 with Xe-137) |
| 03 | V: balances; C: all four loop models on NCL-16 | flow wall law: loss 0.63x measured; DRIDN Fe wrong sign |
| 04 | V; B: ORNL design model; **VAL**: Xe and noble metals; C: bubbles, Cr | Xe 1.03% (ORNL 1.08%), **4x above the measured 0.27%**; salt noble metals below range |
| 05 | B: flow vs VTB; V: all balances, 2D = 0D | mdot 19,270 vs 19,220 kg/s; T_max 210 K low |

The only true validations, in case 04, **fail**, in the same way ORNL's own analyses did. That is
the main lesson: a model that reproduces its calibration, or another code, is not validated until
it is compared with independent measurements.

Every case is a CSVDiff regression test:

```bash
cd ~/projects/open_pronghorn
./run_tests -j 24 --test-root examples/training             # cases 01-04
./run_tests -j 24 --test-root examples/training --heavy     # case 05 (8 ranks)
```

## Building the slides

```bash
cd ~/projects/open_pronghorn/doc
./moosedocs.py build --config training_chemistry.yml --serve
```

Then open http://127.0.0.1:8000/training/chemistry/index.html.

## Known limitations

The [THEORY.md](THEORY.md) sections list the limitations of each model. The most important are:

- **Radiolysis.** G values are assumed gamma yields, while fission fragments dominate the dose in
  fuel salt.
- **Corrosion.**
  - The element-resolved equilibria are user assumptions.
  - Wall comproportionation and film mass transfer are missing.
  - MSTDB-TC needs licensed data.
- **Fission products.**
  - Noble metals are a perfect mass-transfer sink.
  - Untracked precursors such as Xe-137 are lumped.
  - Wall-plated inventories (`plating_boundaries`) are not stored in the `_plated` variables.
- **MSFR.** Its peak temperature depends on the turbulence model (±200 K), and U-235 fast yields
  stand in for U-233.
