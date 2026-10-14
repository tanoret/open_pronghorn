# Module 0: Setup id=module0

!---

## Agenda

| time | module | cases |
| :- | :- | :- |
| 08:30 | 0. Setup | |
| 09:00 | 1. Radiolysis | 01–04 |
| 10:30 | 2. Fission-product inventory | 05–09 |
| 13:00 | 3. Corrosion and plating | 10–13 |
| 14:30 | 4. Coupled, reactor-grade | 14–16 |
| 16:30 | Wrap-up | |

All material is in `examples/training/chemistry`: one directory per case with the input, a
`README.md` (exercises and solutions), a `plot.py` and a regression test.

!---

## Environment

```bash
ml use.moose moose-dev-openmpi
moose-dev-shell
cd ~/projects/open_pronghorn/examples/training/chemistry
export OPEN_PRONGHORN_EXE=~/projects/open_pronghorn/open_pronghorn-opt

cd 01_f2_accumulation
$OPEN_PRONGHORN_EXE -i f2_accumulation.i
python plot.py
```

- Cases 01–14 run serially in seconds; 15–16 use `mpirun -n 24`.
- Override any input variable from the command line: `dose_rate=2e6`.

!---

## Three actions, one solver basis

| action | what it builds | data |
| :- | :- | :- |
| `[MoltenSaltRadiolysis]` | radiolysis species, reactions, G values | `msr_database.json` |
| `[FissionProductInventory]` | nuclides, yields, decay/capture, removal | `fission_product_database.json` |
| `[CorrosionPlating(Flow)]` | Butler-Volmer walls seeded from the corrosion correlation | `corrosion_database.json` |

- Every species is a **linear finite-volume** system.
- List the systems in `[Problem] linear_sys_names` and in the executioner.
- The actions share variables: corrosion can release Cr straight into the radiolysis Cr(II).
- The same input runs 0D, 1D scans, or 2D with the SIMPLE flow solution.
