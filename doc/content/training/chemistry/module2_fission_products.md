# Module 2: Fission-product inventory id=module2

!---

## The `[FissionProductInventory]` action

!equation
\frac{dN_i}{dt} = \frac{Y_i F}{N_A} + \sum_j b_{j\to i}\lambda_j N_j + \sigma_j\phi N_j
 - \left(\lambda_i + \sigma_i\phi + \sum_k r_{ik}\right)N_i

- ENDF/B-VIII.0 yields (U-235 thermal/fast), NUBASE2020 half-lives, one-group cross sections.
- Decay/capture links between tracked nuclides are found automatically. Untracked precursors are
  lumped in (the `effective` yield).
- Removal by chemical group: `noble_gas`, `noble_metal`, `halogen`, `salt_seeking`. The
  destinations `offgas` and `plated` become inventories that keep decaying.
- Systems are ordered parent-before-daughter, so one sweep per step is exact.

!---

## Case 05: saturation versus growth

!row!
!col! width=40%
- 20 MW/m$^3$ for one year.
- I-131 (8 d) saturates; its activity $= Y F$ measures the power.
- Cs-137 (30 y) grows linearly and measures the burnup.
- Both match the analytic solution to $10^{-7}$ or better.
!col-end!

!col! width=60%
!media case05_fp_buildup.png style=width:100%; alt=I-131 saturating and Cs-137 growing linearly
!col-end!
!row-end!

!---

## Case 06: the xenon pit

!row!
!col! width=40%
- I-135 -> Xe-135 with burnout ($\sigma = 2.65\times10^6$ b).
- 4 days at power, then a scram.
- Peak **8.3 h after the scram, 1.89x** the equilibrium.
- Analytic agreement: $10^{-4}$ (equilibrium), $2\times10^{-4}$ (post-scram).
!col-end!

!col! width=60%
!media case06_xenon_transient.png style=width:100%; alt=I-135 and Xe-135 inventories and xenon worth through a scram
!col-end!
!row-end!

!---

## Case 07: stripping the poison

!row!
!col! width=40%
- `removal_groups = 'noble_gas'` with `destinations = 'offgas'`.
- 5 cells scan the stripping rate $r$.
- Poisoning goes from 2.1% to **0.19%** at $r = 10^{-3}$ 1/s (MSRE: a few tenths of a percent).
- Stripping pays off once $r > \lambda_{Xe} + \sigma\phi$.
!col-end!

!col! width=60%
!media case07_xenon_stripping.png style=width:100%; alt=Xe-135 poisoning versus time and versus stripping rate
!col-end!
!row-end!

!---

## Case 08: where the noble metals go

!row!
!col! width=40%
- Two pathways: plated ($r_p$) and off-gas ($r_o$).
- **< 1% stays in the salt** (MSRE, ORNL-4865).
- Plated Te-132 produces **iodine on the walls**: removed inventories must keep decaying.
!col-end!

!col! width=60%
!media case08_noble_metal_partition.png style=width:100%; alt=Mo-99 in salt, plated and off-gas, and partition bars for three nuclides
!col-end!
!row-end!

!---

## Case 09: a month of operation

!row!
!col! width=40%
- Start-up, load follow, shutdown.
- Xe worth is -1066 pcm at power, and **-1325 pcm** at the post-shutdown peak with stripping
  (-4030 pcm without).
- Cs-137 tracks the energy produced to 0.9987.
!col-end!

!col! width=60%
!media case09_power_history.png style=width:90%; alt=Power history, xenon worth and normalized I-131 and Cs-137
!col-end!
!row-end!
