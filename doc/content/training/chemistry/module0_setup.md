# Setup id=setup

!---

## Five cases, increasing complexity

| case | model | new in this case | time |
| :- | :- | :- | :- |
| 01 fuel-salt redox | 0D | `[MoltenSaltRadiolysis]`, hand-added kernels | 45 min |
| 02 xenon and power history | 0D, transient | `[FissionProductInventory]`, functors of time | 40 min |
| 03 one loop, four corrosion models | 0D/1D walls | `[CorrosionPlating(Flow)]`, DRIDN, model comparison | 50 min |
| 04 MSRE fuel loop | multi-region | all actions together, compared with reactor data | 50 min |
| 05 MSFR primary circuit | 2D-RZ CFD, 8 ranks | flow restart, all actions at reactor scale | 60 min |

All material is in `examples/training/chemistry`: one directory per case (input, `README.md` with
exercises and solutions, `plot.py`, regression test). The theory and usage manual of every model is
`THEORY.md`.

!---

## Environment

```bash
ml use.moose moose-dev-openmpi
moose-dev-shell
cd ~/projects/open_pronghorn/examples/training/chemistry
export OPEN_PRONGHORN_EXE=~/projects/open_pronghorn/open_pronghorn-opt

cd 01_fuel_salt_redox
$OPEN_PRONGHORN_EXE -i fuel_salt_redox.i
python plot.py
```

- Cases 01–04 run serially in seconds to a minute and a half. Case 05 runs on 8 ranks in two
  steps (`mpirun -n 8`, ~7 min).
- Any input variable can be overridden from the command line, for example `P_fission=4e7`.

!---

## Three actions, one solver basis

| action | what it builds | data |
| :- | :- | :- |
| `[MoltenSaltRadiolysis]` | radiolysis species, reactions, G values, gas exchange | `msr_database.json` |
| `[FissionProductInventory]` | nuclides, yields, decay/capture, removal pathways | `fission_product_database.json` |
| `[CorrosionPlating(Flow)]` | Butler-Volmer walls seeded from the corrosion correlation | `corrosion_database.json` |

- Every species is a **linear finite-volume** system.
- List the systems in `[Problem] linear_sys_names` and in the executioner.
- The actions share variables, and you can add your own kernels to them (case 01).
- The same input runs in 0D, in 1D, or in 2D with the SIMPLE flow solution.
