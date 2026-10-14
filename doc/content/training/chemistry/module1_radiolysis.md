# Module 1: Radiolysis id=module1

!---

## Case 01: the G value and its units

!row!
!col! width=45%
- One cell, constant dose rate $D$, direct F$_2$ yield only:

!equation
\frac{dC}{dt} = \frac{G\,D}{100\,\mathrm{eV}\,N_A}

- $D = 10^6$ W/m$^3$: gamma heating in an MSR core.
- Result: matches $S\,t$ to $3\times10^{-15}$.
!col-end!

!col! width=55%
!media case01_f2_accumulation.png style=width:100%; alt=F2 concentration growing linearly, model on top of the analytic line
!col-end!
!row-end!

!---

## Case 01: the input

!listing examples/training/chemistry/01_f2_accumulation/f2_accumulation.i block=MoltenSaltRadiolysis

- `g_value_species`/`g_value_overrides` switch the radical yields off.
- One linear system per species: `linear_sys_names` in `[Problem]` and `system_names` in the
  executioner.

!---

## Case 02: the chloride radical chain

!row!
!col! width=40%
- LiCl-KCl, 400 C, G = 0.3 for e$^-$ and Cl.
- 4 independent cells form a **dose scan**, $10^2$–$10^5$ W/m$^3$.
- Cl$_2^-$ $\to \sqrt{S/2k}$ (slope 0.50).
- Cl$_3^-$ $\to S/2$ (slope 1.00).
- Fixed-point sweeps couple the segregated systems.
!col-end!

!col! width=60%
!media case02_chloride_radiolysis.png style=width:100%; alt=Radical concentrations versus time and dose scaling of Cl2- and Cl3-
!col-end!
!row-end!

!---

## Case 03: chromium redox under gamma

!row!
!col! width=40%
- `metals = 'Cr'` adds Cr(I)/(II)/(III).
- Invariants:

  - Cr conservation: $5\times10^{-7}$
  - redox balance: $2\times10^{-6}$

- The salt **stays reduced**: Cr(III)/Cr $= 4\times10^{-6}$.
- Oxidizing equivalents accumulate as Cl$_2$.
!col-end!

!col! width=60%
!media case03_cr_redox.png style=width:100%; alt=Chromium oxidation states and conservation errors versus time
!col-end!
!row-end!

!---

## Case 04: the uranium redox buffer and the cover gas

!row!
!col! width=40%
- FLiBe-UF$_4$, U(III)/U = 1% (MSRE).
- The direct F$_2$ yield slowly oxidizes the salt, and the **buffer is spent in 1.6 h**.
- Then F$_2$ appears in the salt and the headspace (Henry's law, `gas_species = 'F2'`).
- Exercise: `U_III_0=0`.
!col-end!

!col! width=60%
!media case04_fuel_salt_redox.png style=width:100%; alt=U(III) fraction falling and F2 rising in salt and cover gas
!col-end!
!row-end!
