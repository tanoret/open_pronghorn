# Wrap-up id=wrapup

!---

## Verification, calibration, validation

| kind | meaning | in this training |
| :- | :- | :- |
| **Verification** | analytic solutions and exact conservation | all five cases |
| **Calibration** | reproduces data the model was fitted to | 03 (NCL-16), 04 (bubble rate, Cr) |
| **Code-to-code** | another code or model for the same problem | 04 (ORNL design model), 05 (VTB) |
| **Validation** | independent measured data | 04: MSRE Xe and noble metals (**both fail**) |

All five cases are regression tests: `./run_tests --test-root examples/training [--heavy]`.

!---

## Limits of the models

- **Radiolysis and redox:** G values are assumed gamma yields, while fission fragments dominate the
  dose. Fission oxidation is ~0.8 eq/fission. Impurities are not modelled.
- **Fission products:** selected chains with one-group cross sections; lumped precursors (Xe-137).
  Noble metals are a perfect mass-transfer sink, which fails against the MSRE.
- **Corrosion:** the correlation, MSTDB-TC and DRIDN are calibrated in-sample. The equilibria of the
  element-resolved wall law are assumptions. Wall comproportionation and film mass transfer are
  missing.
- **Reactor models:** the MSRE is lumped (loop time much shorter than the chemistry times). The MSFR
  is steady 2D-RZ with one lumped loop, and its T$_{max}$ is uncertain by ±200 K.

!---

## Where to go next

- `examples/training/chemistry/THEORY.md`: the equations, usage and pitfalls of every model.
- Validation data: `validation/corrosion` (76 cases) and `validation/msr`.
- The VTB MSFR and MSRE models: `virtual_test_bed/msr`.

!style halign=center fontsize=150%
Questions?
