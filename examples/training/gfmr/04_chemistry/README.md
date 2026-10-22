# Stage 4: fuel-salt chemistry on the 3-D core flow

`gfmr_chemistry.i` solves the fuel-salt chemistry as passive scalars on the converged stage-3
flow. Turbulence and temperature are frozen. Momentum and pressure are re-solved from the
converged restart, which only rebuilds the face fluxes that advect the species.

| Action | What it models here |
|---|---|
| `[MoltenSaltRadiolysis]` (`salt_type = fluoride`, `metals = U`) | the U(IV)/U(III) redox buffer of FUNaK (26.4 mol % UF4) |
| `[CorrosionPlatingFlow]` (`hastelloy_n`, `fluoride_fuel`) | Cr / Fe / Ni corrosion and plating of the 9 Hastelloy N guide tubes, the only metal in the core |
| `[FissionProductInventory]` (I-135, Xe-135) | fission yields, decay, Xe-135 burnout in the thermal flux, and uptake by the graphite channel walls |

## Why a loop balance is needed

The GFMR core is **once through**: the salt crosses it in 3.0 s and leaves. Everything that
builds up over minutes to days happens in the whole circulating loop. The 3-D core therefore
needs inlet values, and those come from a loop balance. The specification defines **no loop**,
so `loop_balance.py` assumes one. Every number below marked ASSUMED is a choice, not data:

| Loop quantity | Value | Status |
|---|---|---|
| Core salt volume, core transit | 1.70 m³, 3.02 s | from the stage-3 CFD |
| Ex-core volume | 2 × core salt | ASSUMED (MSRE: core ≈ 1/3 of the fuel salt) |
| Loop time | 9.06 s | follows |
| U(III)/U held by redox control | 1 % | ASSUMED (MSRE practice) |
| Redox-control time constant | 1 day | ASSUMED |
| Xe stripping (He sparging) | 10 % per loop pass | ASSUMED (MSRE spray ring: 10–15 %) |
| Fission oxidation | 0.8 equivalents per fission | ORNL MSRE experience |
| Radiolytic dose | 0.87 × fission power density | fragments + β + part of γ |

`python loop_balance.py` prints the loop quantities and writes `inlet_values.i`, which
`gfmr_chemistry.i` includes.

## Results

Run on 24 ranks, ≈ 20 min for 60 SIMPLE iterations. The species residuals reach 10⁻⁹–10⁻¹⁰.

### Verification: the 3-D core against the 0-D balance

| Quantity | 3-D (flow-weighted, last minus first layer) | 0-D core pass | Ratio |
|---|---|---|---|
| I-135 rise across the core | 5.28×10⁻⁷ mol/m³ | 5.09×10⁻⁷ | 1.04 |
| Xe-135 rise across the core | 3.05×10⁻⁷ mol/m³ | 3.20×10⁻⁷ | 0.95 |

### Numerical floor (read before using the redox and corrosion fields)

A source-free species must leave the core at its inlet value. `F_ion` has no net source here,
and it comes out with a **mixed relative error of 1.1×10⁻⁵**, and **up to 0.6 % in single
cells** (corner cells of the channel inlets).

That floor is **larger than the physical per-pass changes of the redox state and of the
corrosion products**:

| Quantity | Physical change per core pass (0-D) | Relative to the value |
|---|---|---|
| U(III) oxidized by fission | 1.0×10⁻⁵ mol/m³ | 1.1×10⁻⁷ |
| Cr exchanged at the 9 guide tubes | < 10⁻⁶ mol/m³ | < 10⁻⁶ |

**The 3-D U(III) and Cr fields are therefore numerical noise at the 10⁻⁵ level and carry no
physical information.** The redox and corrosion results below come from the loop balance. The
3-D model is meaningful for the fission products, whose per-pass rise (1.6×10⁻⁴ relative for
I-135) is above the floor.

### Physical results

| Result | Value | Source |
|---|---|---|
| U(III)/U | 1.00 % | held by the assumed redox control |
| Fission oxidation, whole loop | 0.54 mol eq/day | 0-D: 0.8 eq × 4.68×10¹⁸ fissions/s |
| Reductant needed to hold U(III)/U | 0.54 mol eq/day: 43 g/day of U metal (U + 3 UF₄ → 4 UF₃) or 2.4 g/day of Be | 0-D |
| F₂ in the salt | 0 | charge-consistent radiolysis (net G = 0) |
| I-135 in the salt (loop average) | 3.26×10⁻³ mol/m³ | 0-D, verified in 3-D (above) |
| Xe-135 in the salt (loop average) | 8.47×10⁻⁶ mol/m³ | 0-D, verified in 3-D (above) |
| Xe-135 removal rates (1/s) | decay 2.1×10⁻⁵, burnout 8.1×10⁻⁵, stripping 1.2×10⁻² | 0-D |
| Xe-135 poisoning, Σ_a,Xe / Σ_a,salt (thermal) | 2.9×10⁻⁴ (order 30 pcm) | 0-D, rough one-group |
| Guide-tube corrosion | below the numerical floor: the salt is held at the Cr equilibrium of the loop | see above |

**Stripping dominates the xenon balance.** At 10 % per pass the removal rate is 150 times the
burnout rate, so the xenon poison is small (order 30 pcm). The result depends almost entirely on
the assumed stripping efficiency (exercise 2).

**The graphite uptake of xenon is included as a channel-wall sink.** It uses a turbulent mass
transfer coefficient and a perfect-sink graphite, which is an upper bound. In graphite MSRs this
uptake is the main xenon-poisoning concern; the MSRE sealed its graphite to limit it.

### Limits

- **No GFMR-specific chemistry data exist.** FUNaK has no U(III)/U(IV) data comparable to FLiBe;
  the `fluoride_fuel` corrosion class and the FLiBe-like equilibria are stand-ins.
- The loop is assumed. Redox control, stripping and the ex-core volume set the answers.
- The 3-D redox and corrosion fields are below the advection floor (above).

## How to run

```bash
python loop_balance.py                               # 0-D loop, writes inlet_values.i
mpirun -n 24 open_pronghorn-opt -i gfmr_chemistry.i  # needs ../03_cfd/gfmr_core_cfd_out.e
```

## Exercises

1. Remove the redox control and integrate the 0-D U(III) balance: how many days until U(III)/U
   falls below 0.1 %?
2. Change the stripping efficiency from 10 % to 1 % per pass in `loop_balance.py`: how do the
   Xe-135 concentration and its reactivity change?
3. Explain why a source-free species can show a 10⁻⁵ error at the outlet, and why that hides
   the U(III) change in 3-D but not the I-135 change.
