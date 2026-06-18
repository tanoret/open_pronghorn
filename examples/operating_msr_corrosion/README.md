# Coupled radiolysis + corrosion over an MSR power cycle

`operating_msr_corrosion.i` couples radiolysis and corrosion through the **dissolved chromium redox
couple** in a molten chloride coolant pool, over a reactor power cycle (ramp up, hold at full power,
ramp down). It is the chemistry-coupling counterpart of the spatially-resolved
`examples/flowing_msr_corrosion`.

## The coupling

Two physics act on the same dissolved-chromium species:

- **Corrosion** — the structural-alloy (316 stainless) walls dissolve through a Butler-Volmer
  reaction ([CorrosionLinearFVButlerVolmerBC](../../doc/content/source/corrosion/linearfvbcs/CorrosionLinearFVButlerVolmerBC.md)),
  releasing chromium into the salt **as Cr(II)**. The release rate is seeded from the validated
  corrosion correlation. The action targets the radiolysis-tracked `Cr_II` variable directly
  (`release_variables = 'Cr_II'`), so the corrosion product *is* the species the radiolysis network
  acts on.

- **Radiolysis** — irradiation produces the chlorine radical chain (here the net oxidant `Cl2m_rad`,
  Cl2.-, is produced proportional to the reactor power). The database reactions
  `Cl2m_rad + Cr(II) -> Cr(III)` and `e_sol + Cr -> ...` then oxidize and reduce the dissolved
  chromium.

The coupling is two-way:

1. **corrosion -> radiolysis:** the wall reaction feeds Cr(II) into the salt, which is both the
   substrate that the radiolytic oxidant transmutes and a sink (getter) for that oxidant;
2. **radiolysis -> corrosion:** radiolysis sets the salt's Cr(III)/Cr(II) oxidation state — its
   oxidizing power — which is what drives corrosion (represented here by the oxidizing redox class of
   the correlation).

## What the cycle shows

| Quantity | Behavior |
| :------- | :------- |
| `reactor_power_fraction` | ramps 0 -> 1, holds, ramps down |
| `Cl2m_rad_max` | radiolytic oxidant, tracks the power |
| `redox_ratio_CrIII_CrII` | the coupling observable: rises from 0 to ~1.1 as the radiolytic oxidant converts the corrosion-product Cr(II) to Cr(III) |
| `Cr_total` | grows steadily from the wall corrosion source |

So the salt oxidation state responds to the reactor power (radiolytic oxidation of the corrosion
product), while corrosion continually replenishes the reduced Cr(II) — exactly the radiolysis-driven
redox dynamics that govern chromium corrosion in an irradiated molten salt.

```
open_pronghorn-opt -i operating_msr_corrosion.i
```

The radical chain is lumped into a net oxidant yield (Cl2.- produced directly, as the other examples
lump radicals into the net product) to remove the stiff solvated-electron kinetics; remove the
`g_value_overrides` to integrate the explicit radical network at a smaller time step. For the fully
flow-resolved coupling (Navier-Stokes + turbulence + energy + radiolysis + corrosion in one solve)
see `examples/flowing_msr_corrosion`.
