# Module 3: Corrosion and plating id=module3

!---

## Calibrated Butler-Volmer walls

!equation
i = i_0\left[a_s\,e^{\alpha z f\eta} - \frac{c}{c_{ref}}e^{-\alpha z f\eta}\right], \qquad
J = \frac{i}{zF}

- $i_0$ is **seeded** so the initial rate reproduces the calibrated correlation for the alloy, salt,
  redox, position, flow and $\Delta T$ classes.
- $a_s = c_s/c_{s,0}$, the alloy-surface activity (two-block, `solid_activity = true`).
- Topologies: `salt_only`, `solid_only`, `two_block` (`[CorrosionPlating]`). For flow problems,
  the linear-FV `[CorrosionPlatingFlow]`.

!---

## Case 10: a static capsule is self-limiting

!row!
!col! width=40%
- 304L in fuel salt (ORNL-FL-01): **16.82 µm/y** initially.
- The salt saturates at $c_{eq} = c_{ref}e^{zf\eta}$ = 11.2 mol/m$^3$.
- The recession stops at **0.66 µm** (the Faraday limit), not 5.5 µm.
- Sustained attack needs a sink.
!col-end!

!col! width=60%
!media case10_corrosion_coupon.png style=width:100%; alt=Corrosion rate decaying and recession approaching the Faraday limit
!col-end!
!row-end!

!---

## Case 11: the Cr-depletion layer

!row!
!col! width=40%
- The alloy is resolved as a second block (1 µm cells).
- The rate falls from 23.4 to **3.2 µm/y** in a year.
- Surface Cr is at 48% of the bulk; the zone below 95% is 125 µm deep ($\sim\sqrt{4Dt}$).
- Cr balance: $10^{-13}$.
!col-end!

!col! width=60%
!media case11_cr_depletion.png style=width:100%; alt=Corrosion rate versus time and Cr profile in the alloy
!col-end!
!row-end!

!---

## Case 12: alloy, salt and redox against loop data

!row!
!col! width=45%
- 7 measured loop/capsule rates from 1 to 217 µm/y.
- One input, CLI overrides, a Python driver.
- **All within 2x**, geometric-mean error 1.18.
- Stainless corrodes ~10x faster than Hastelloy N; moisture or stress add another ~10x.
!col-end!

!col! width=55%
!media case12_parity.png style=width:85%; alt=Parity plot of model versus measured corrosion rates
!col-end!
!row-end!

!---

## Case 13: a three-year thermal-convection loop

!row!
!col! width=40%
- The DRIDN reduced-order model: kinetics, inventory, cold-leg capture, depletion.
- ORNL NCL-16: **533 ppm** Cr against ~500 measured.
- The hot leg loses, the cold leg gains; mass balance $4\times10^{-16}$.
!col-end!

!col! width=60%
!media case13_dridn_loop.png style=width:100%; alt=Dissolved Cr, Ni, Fe and hot/cold-leg mass change over three years
!col-end!
!row-end!
