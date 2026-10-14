# Module 4: Coupled, reactor-grade id=module4

!---

## Case 14: radiation switches corrosion on

!row!
!col! width=40%
- 316 walls, chloride pool, **pre-equilibrated** salt: no corrosion without dose.
- At power, Cl$_2^-$ oxidizes Cr(II) to Cr(III) (84%), and the walls corrode at **35 µm/y**.
- After shutdown the salt remembers: 33 µm/y.
- `release_variables = 'Cr_II'` couples the two actions.
!col-end!

!col! width=60%
!media case14_irradiated_pool.png style=width:85%; alt=Power, Cr(II)/Cr(III) and wall corrosion rate over 11 days
!col-end!
!row-end!

!---

## Case 15: a heated turbulent channel

!row!
!col! width=35%
- SIMPLE + k-epsilon + energy, with all three actions on the flow.
- $k_m$ from Chilton-Colburn: $1.3\times10^{-5}$ m/s.
- Energy balance $3\times10^{-4}$; Mo-99 balance $4\times10^{-4}$.
- Only $10^{-4}$ plates per pass: plating is a many-pass process.
!col-end!

!col! width=65%
!media case15_core_channel.png style=width:100%; alt=Wall temperature, corrosion rate and mid-height Cr and Mo-99 profiles
!col-end!
!row-end!

!---

## Case 16: the MSFR primary loop

- Pump-driven closed loop, conjugate heat transfer, chloride fuel salt, 24 ranks, 1.5 min.
- Radiolysis + corrosion + I/Xe with **top-leg sparging** + Mo/Ru/Te **plating on the HX**.
- Removal is placed with spatial functors (`strip_rate`, `plate_coefficient`).

!media case16_msfr_loop.png style=width:75%;margin-left:auto;margin-right:auto;display:block; alt=Top-leg temperature and loop inventories relative to the 0D prediction

!---

## Case 16: validation against the 0D cases

| quantity | 60x60 | 90x90 |
| :- | -: | -: |
| energy balance error | 6e-5 | 2e-4 |
| Xe-135 / 0D | 1.00002 | 1.00002 |
| Mo-99 / 0D | 1.0004 | 1.0005 |
| Xe-135 stripped | 83% | 83% |
| Mo-99 plated | 94% | 94% |
| core-leg velocity (m/s) | 1.56 | 1.68 |

- The loop time ($\sim$5 s) is much shorter than the half-lives, so the inventory is uniform to
  $7\times10^{-4}$ and equals the 0D prediction.
- Closed-loop corrosion settles at the Nernst equilibrium (Cr(III)/Cr(II) = 5.1).
