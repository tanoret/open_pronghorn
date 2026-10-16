# The five cases id=cases

!---

## Case 01: fission oxidation and the U(III) buffer

!row!
!col! width=40%
- MSRE fuel salt (F$^-$ 8.0e4, U 490 mol/m$^3$, U(III)/U = 1%).
- Radiolysis recombines; the net oxidant is **fission**: ~0.8 eq/fission, i.e. 1.66e-7 eq/(m$^3$ s)
  at 4 MW/m$^3$.
- The buffer is spent at **day 339** (analytic 342), on the time scale of MSRE's Be additions.
  F$_2$ then passes to the cover gas (Henry's law).
- New: `[MoltenSaltRadiolysis]` plus your own `LinearFVSource` kernels.
- **Verification:** U, F and redox balances.
!col-end!

!col! width=60%
!media case01_fuel_salt_redox.png style=width:100%; alt=U(III) fraction falling over a year, F2 appearing in salt and cover gas, and the conservation errors
!col-end!
!row-end!

!---

## Case 02: xenon over a month of operation

!row!
!col! width=40%
- New: `[FissionProductInventory]` with I-135/Xe-135, I-131 and Cs-137; power and flux as functions
  of time; off-gas stripping.
- Xe worth: -1066 pcm at power and **-1325 pcm** at the peak after shutdown (-4030 pcm without
  stripping).
- Cs-137 as a burnup monitor: 0.9987. But untracked Xe-137 hides a ~3% stripping loss (0.967).
- **Verification:** analytic equilibria.
!col-end!

!col! width=60%
!media case02_xenon_power_history.png style=width:90%; alt=Power history, xenon worth and normalized I-131 and Cs-137 over a month
!col-end!
!row-end!

!---

## Case 03: one loop, four corrosion models

!row!
!col! width=40%
- ORNL loop NCL-16 (Hastelloy N, 704/538 C, 3.37 y) seen by four models:
  - a capsule Butler-Volmer wall stalls at the Nernst limit;
  - the **element-resolved flow wall law** (hot walls dissolve, cold walls plate): loss 1.84 vs 2.9
    mg/cm$^2$, Fe −127 vs −100 ppm;
  - **DRIDN** fits all targets, but its Fe has the wrong sign;
  - the **MSTDB-TC** endpoint fits them too, but needs licensed data.
- **Calibration**, not validation: NCL-16 tuned every loop model.
!col-end!

!col! width=60%
!media case03_corrosion_models.png style=width:100%; alt=NCL-16 measured values against four corrosion models, salt Cr versus time and wall mass changes
!col-end!
!row-end!

!---

## Case 04: the MSRE against its record, xenon

!row!
!col! width=40%
- A lumped MSRE loop: core, spray stripping, graphite, Hastelloy N; all actions together.
- **A priori Xe-135: 1.03%**, matching ORNL's design model (1.08%, code-to-model).
- **Measured: 0.26–0.28%.** The validation fails by 4x: circulating helium bubbles were missing.
- One fitted bubble-stripping rate reproduces 0.270% (**calibration**).
!col-end!

!col! width=60%
!media case04_msre_xenon.png style=width:100%; alt=MSRE xenon poisoning versus time, a priori and with bubble stripping, against the measured band
!col-end!
!row-end!

!---

## Case 04: noble metals and corrosion in the MSRE

!row!
!col! width=40%
- Mass-transfer-limited deposition puts **0.06–0.5%** of the noble metals in the salt, against
  **1–20% measured** (ORNL-4865). The validation fails, as ORNL's own analysis did.
- Corrosion: 3.0 µm/y, and the salt Cr rises +50 ppm in 6 weeks (McCoy: +60 ppm). This is a
  calibration target of the correlation.
- Lesson: a model is validated only against data it was not tuned to.
!col-end!

!col! width=60%
!media case04_msre_noble_metals.png style=width:100%; alt=Model noble-metal shares in salt, metal, graphite and off-gas against the MSRE ranges
!col-end!
!row-end!

!---

## Case 05: the MSFR primary circuit (VTB geometry)

- EVOL/VTB MSFR: 3 GW, LiF-ThF$_4$-UF$_4$, 2D-RZ toroidal core, pump and HX, Griffin power.
  New: CFD (SIMPLE, k-epsilon, energy), a flow restart, 8 ranks.
- **Code-to-code against the VTB:** mass flow 19,270 vs 19,220 kg/s; loop 997 to 902 K
  (VTB dT 97.6 K). T$_{max}$ is 1095 vs 1304 K (k-epsilon against a mixing length).

!media case05_msfr_flow.png style=width:90%;margin-left:auto;margin-right:auto;display:block; alt=MSFR core-cavity temperature and velocity profiles and the external-loop temperature

!---

## Case 05: chemistry at reactor scale

- **Redox:** fission oxidation is balanced by a processing-unit control, giving U(III)/U = 0.96%.
- **Corrosion:** walls above ~965 K dissolve and the cold HX ducts plate. Cr settles at 49 wppm, and
  0.52 kg/y of Cr moves from hot to cold.
- **Fission products:** 99.94% of the Xe-135 is stripped by bubbling; Mo-99 plates 99.5% in the HX.
- **Verification** only: all balances are below 1e-5 and the 2D result equals 0D to 2%. No MSFR
  chemistry has been measured.

!media case05_msfr_chemistry.png style=width:90%;margin-left:auto;margin-right:auto;display:block; alt=MSFR fission products around the loop, wall temperatures and the Cr wall flux sign along the walls
