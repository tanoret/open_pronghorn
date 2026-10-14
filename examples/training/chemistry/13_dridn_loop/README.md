# Case 13: three years of a Hastelloy N thermal-convection loop (DRIDN)

**Module 3 (corrosion) | ~15 min | 0D, serial, < 1 s**

## Objective

Use a reduced-order long-term model when resolving the physics is out of reach: a multi-year loop
with hot-leg dissolution, salt inventory, cold-leg plating and alloy depletion.

## Physics

Case 10 showed that a closed capsule stops corroding at the Nernst limit. A **loop** never gets
there: the temperature dependence of the equilibrium makes the hot leg dissolve and the cold leg
plate, so mass is carried continuously from hot to cold.

The Dynamic Redox-Inventory-Depletion Network (DRIDN) integrates, over years:

- the coupled Cr/Fe/Ni dissolution kinetics and their thermodynamic driving forces,
- the dissolved inventory of the salt (the area/salt-mass ratio and the capacity of the loop),
- cold-leg capture (plating) of the dissolved metals,
- depletion of the alloy surface and the attack morphology (front and intergranular depth).

The physical case is ORNL thermal-convection loop NCL-16 (ORNL-TM-4188, validation target M-014):
purified fuel salt, 704 C hot leg, 166 C loop temperature span, 3.37 years. The measured salt
chromium increase was **~500 ppm**.

## The input, step by step

[dridn_loop.i](dridn_loop.i):

- `[Problem] solve = false`: no PDE is solved. The `DRIDNCorrosionUserObject` integrates its ODE
  network internally once per time step (40 steps over 3.37 years).
- The parameters are grouped in the input: alloy composition, calibrated kinetics, loop mass
  transfer and cold capture, salt inventory, and alloy-side depletion.
- `AdvancedCorrosionModelPostprocessor` exposes the model state: dissolved ppm, front rate and
  depth, mass loss and gain, and the mass-balance error.

## Run

```bash
open_pronghorn-opt -i dridn_loop.i
python plot.py
```

## Expected results

![](../figures/case13_dridn_loop.png)

- Salt Cr after 3.37 y: **532.7 ppm**, against ~500 ppm measured.
- The hot leg loses 3.13 mg/cm^2 and the cold leg gains 1.61 mg/cm^2. The difference is the
  dissolved inventory.
- Front rate 1.54 µm/y (compare with case 12's purified Hastelloy N, 1.0–1.4 µm/y), with 37 µm
  of intergranular attack.
- `mass_balance_error` stays below $10^{-15}$.

## Exercises

1. Turn the cold-leg capture off (`"UserObjects/dridn/cold_capture_fraction=0 0 0"`). What
   changes, and what does not?
2. Why does the dissolved Cr level off instead of growing linearly?
3. Which parameters of this model are calibrated and which are loop design data? (Read the input
   comments.)

<details><summary>Solutions</summary>

1. The hot-leg loss is unchanged (3.13 mg/cm^2, front rate 1.54 µm/y): in this calibration the
   attack is set by the hot-leg kinetics and mass transfer, not by the cold leg. What changes is
   where the metal ends up. Nothing plates (`mass_gain` = 0), so the salt holds more Fe
   (106 -> 133 ppm) and Ni (378 -> 488 ppm). Cr is nearly unchanged (529 ppm) because it is the
   least captured element.
2. The salt Cr rises to 169 ppm within a month and 535 ppm within the first year, then stays flat
   (533 ppm at 2.4 y, 533 ppm at 3.4 y). From then on the hot-leg release is balanced by the
   cold-leg capture, and the salt is a conveyor carrying the metal from hot to cold.
3. Calibrated: `log_exchange_offsets`, `affinity_baseline`, `selectivity_scale`,
   `mass_transfer_rate`, `cold_capture_fraction` and `deposit_area_factor`. Design data:
   `mass_fractions`, `density`, `area_to_salt_mass` and the exposure time.

</details>
