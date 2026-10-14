# Wrap-up id=wrapup

!---

## What we validated

| check | cases | result |
| :- | :- | :- |
| analytic solutions | 01, 02, 05, 06, 07, 08 | $10^{-4}$ – $10^{-15}$ |
| exact invariants | 03, 04, 11, 13 | $10^{-6}$ – $10^{-16}$ |
| measured data | 07 (MSRE), 08 (MSRE), 12 (7 loops), 13 (NCL-16) | within 2x |
| balances | 15, 16 | < $5\times10^{-4}$ |
| 2D = 0D equilibrium | 16 | $4\times10^{-4}$ |

All 16 cases are regression tests: `./run_tests --test-root examples/training [--heavy]`.

!---

## Limits of the models

- **Radiolysis**: literature G values and rate constants for idealized melts. Impurities (moisture,
  oxides) dominate real salt redox.
- **Corrosion**: calibrated to the loop/capsule database. New alloys, chlorides or high dose are
  extrapolations. The radiolysis-corrosion coupling (case 14) is a hypothesis, not yet validated.
- **Fission products**: selected chains, one-group cross sections. Plating is a mass-transfer
  coefficient (no sticking, re-entrainment or graphite).
- **Coupled loops**: steady. Power transients in a resolved loop need a transient executioner or a
  frozen-flow restart.

!---

## Where to go next

- Theory: `doc/content/theory/fissionProductInventory.md` and the action pages under
  `doc/content/source/`.
- Validation: `validation/corrosion` and `validation/msr`.
- Larger examples: `examples/operating_msr_corrosion`, `examples/msfr_loop`,
  `examples/flowing_msr_corrosion`.

!style halign=center fontsize=150%
Questions?
