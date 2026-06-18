# Molten salt corrosion and plating examples

These examples exercise the `[CorrosionPlating]` action, which sets up a mechanistic molten salt
corrosion/plating problem (per-element salt and solid species, Nernst-Planck transport with
electromigration, an optional electric-potential solve, and Butler-Volmer electrode kinetics) and
seeds the electrode kinetics from the validated effective-correlation port
(`MoltenSaltCorrosionModel`).

## corrosion_coupon.i — two-block coupon (interface kinetics)

A Hastelloy N coupon (solid block) immersed in initially clean fluoride salt (salt block). The solid
starts at the **alloy composition** — constant molar concentrations of Cr, Fe and Ni (≈ 7, 4 and
71 wt %) — and the salt starts essentially free of dissolved metal. For each element the action
creates the salt-side cation `c_<El>` and the solid-side metal `cs_<El>`, and a `ButlerVolmerInterface`
kernel dissolves each metal across the shared interface at the calibrated rate. As corrosion proceeds:

- the salt concentrations `c_Cr/c_Fe/c_Ni` rise from near zero and **diffuse away from the wall into
  the bulk salt** (`salt_Cr_ppm` reaches the tens-of-ppm scale over two years, matching the reference
  inventory predictions),
- a **chromium-depletion layer** grows in the solid near the interface (`solid_Cr_min` drops from
  11968 toward ~6000 mol/m³ as Cr leaves faster than solid-state diffusion replenishes it),
- `recession_um` records the Faradaic-equivalent penetration, and `cr_total` stays flat (chromium is
  conserved: every mole leaving the solid appears in the salt).

The 0.5 m salt reservoir backing the 1 mm coupon sets a realistic reservoir-to-area ratio — a small
closed bath would concentrate the corrosion product unrealistically fast.

```
open_pronghorn-opt -i corrosion_coupon.i
```

## flow_coupled_corrosion.i — single-domain wall (boundary kinetics + flow)

Clean salt enters a heated channel and flows past a corroding alloy wall (the bottom boundary), which
sheds chromium through a `ButlerVolmerBC` at the calibrated rate (≈ 50 µm/y for this 316 stainless
oxidizing case, reported by `corrosion_rate_um_y`; the wall recedes accordingly through
`recession_um`). The shed chromium is advected downstream by a prescribed developed channel velocity
profile (read by the corrosion advection kernel) and spread by dispersion, giving a steady downstream
plume that the flow keeps dilute (the outlet rises only to ≈ 1 ppm). The plume reaches steady state
within a few flow-through times while the wall keeps receding at the slow corrosion pace.

This is the recommended **one-way** flow coupling: the flow drives the dissolved species, which do not
feed back on the flow. For a fully resolved CFD velocity — e.g. the segregated Navier-Stokes solution
in `examples/molten_salt_corrosion/corrosion_channel.i` (which produces `vel_x`/`vel_y` through a
`RhieChowMassFlux` and a `SIMPLE` executioner) — transfer those velocity variables into this system
with a `MultiApp`/functor transfer and point `velocity_x`/`velocity_y` at the transferred fields. The
corrosion transport is a monolithic Newton solve and the segregated flow is a separate linear solve,
so they are coupled through lagged velocity functors rather than a single monolithic system.

```
open_pronghorn-opt -i flow_coupled_corrosion.i
```

## A note on the electrode kinetics

The Butler-Volmer reaction is written in concentration-explicit form: the **anodic** (dissolution)
branch carries the solid metal activity (constant), so the corrosion rate is set by the overpotential
and stays bounded for fresh salt, while the **cathodic** (plating) branch carries the dissolved-ion
activity `c/c_ref`, so plating accelerates as the salt loads and the reaction reaches equilibrium when
`c = c_ref·exp[(α_a+α_c) z F η / R T]`. This is why fresh salt corrodes at the calibrated rate (rather
than runaway) and the dissolved metal accumulates gradually.

## Validation

The calibrated correlation that seeds the kinetics is reproduced term for term by the C++ port and
its unit tests (`unit/src/MoltenSaltCorrosion*Test.C`) and by the MOOSE-level validation in
`validation/corrosion/`, which regenerates the reference effective-Butler-Volmer predictions for all
76 cases / 43 targets and the NCL-16 loop simulation.
