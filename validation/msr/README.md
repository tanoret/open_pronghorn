# Molten Salt Radiolysis Validation Suite

This directory ports the literature validation suite of the standalone `MoltenSaltRadiolysis`
model. Each case keeps the original `manifest.yaml` (citation, experimental conditions, digitization
metadata) and any digitized data under `data/`, mirroring the source repository.

## Status

Three cases ship with digitized quantitative data and are implemented as full validation cases
(MOOSE input + VnV script comparing against the measurement): the Cr(II)/Cr(III) and Zn(II) kinetics
in chloride and the F2 yield in fluoride. The remaining cases are scoping/stub manifests in the
source model itself (`traces: []` and no digitized figures): their citations, conditions and any
auxiliary data are preserved so they can be promoted to quantitative validation cases once their
figures are digitized.

| Case | Kernel | Status | Title |
| :--- | :----- | :----- | :---- |
| `cr_licl_kcl/iwamatsu_2026_pccp` | chloride | **quantitative (9 traces)** | Cr(II)/Cr(III) radiation-induced redox kinetics in LiCl-KCl |
| `zn_licl_kcl/horne_2022_pccp` | chloride | **quantitative (rate)** | Zn2+ + e_s- pseudo-first-order rate in LiCl-KCl (400 C) |
| `eflibe_f2_yield/davis_2022_nse` | fluoride | **quantitative (G value)** | Radiolytic F2 yield for FLiBe-UF4 |
| `licl_kcl_baseline/hagiwara_1987_rpc` | chloride | stub/scoping | Baseline pulse radiolysis of LiCl-KCl at 400 C |
| `ealkali_chlorides/horne_2018_jpcc` | chloride | stub/scoping | Foundational e_s- characterization in molten alkali chlorides |
| `cl3_chlorobasicity/horne_2025_pccp` | chloride | stub/scoping | Cl2 vs Cl3- branching by chlorobasicity (MD) |
| `i_licl_kcl/horne_2023_pccp` | chloride | stub/scoping | Iodide-impurity impact on radiolytic transients in LiCl-KCl |
| `nd_licl_kcl/horne_2026_ic` | chloride | stub/scoping | Nd(II)/Nd(III) influence on radiolytic transients in LiCl-KCl |
| `lanthanide_epsilon_licl_kcl/moon_chidambaram_2022_pnse` | chloride | stub/scoping | Molar absorption coefficients of trivalent lanthanides in LiCl-KCl |
| `zncl2_neat/gibson_2023_jpcb` | chloride | stub/scoping | Excess-electron states in molten ZnCl2 (MD) |
| `kcl_mgcl2_solid/ramos_ballesteros_2022_jpcc` | chloride | stub/scoping | Long-lived transients in irradiated solid KCl-MgCl2 |
| `mcfr_uc13_irradiation/karlsson_2022_inl` | chloride | stub/scoping | Gamma irradiation of NaCl-UCl3 for the MCFR |
| `oxidants_halide_melts/makarov_1982_rcb` | chloride | stub/scoping | Oxidizing agents from radiolysis of alkali halide melts |
| `alkali_halide_review/pikaev_1982_rpc` | chloride | stub/scoping | Solvated electron in irradiated alkali halide melts (review) |
| `libr_kbr_pulse_rad/sawamura_1990_rpc` | (bromide) | stub/scoping | Pulse radiolysis of LiBr-KBr melts |
| `flibe_msre/haubenreich_msr_69_46` | fluoride | stub/scoping | Steady-state F2 buildup in stored MSRE fuel salt |
| `flibe_msre_f2/heron_1990_redds` | fluoride | stub/scoping | F2 generation by gamma radiolysis of MSRE-composition salt |
| `flinak_pulse_rad/akiyama_1994_jnst` | fluoride | stub/scoping | Short-lived species in pulse-irradiated FLiNaK |

## Quantitative cases

- **`cr_licl_kcl/iwamatsu_2026_pccp`** reproduces the transient-absorption decay of the solvated
  electron captured by Cr(II) and Cr(III) in molten LiCl-KCl. `cr_pulse.i` seeds a radiolytic pulse
  and integrates the chloride + Cr chemistry; `cr_vnv.py` rescales the modeled `e_sol` decay to the
  measured maximum (the absorbance is proportional to `[e_sol]` up to an unknown scale factor) and
  bounds the scale-free relative RMSE against each digitized trace. All nine traces (four Cr(II) and
  five Cr(III) concentrations) are validated, with scale-free relative RMSE of 0.05-0.14.

- **`zn_licl_kcl/horne_2022_pccp`** reproduces the pseudo-first-order rate of e_sol capture by Zn(II).
  `zn_vnv.py` extracts `k_obs` from the modeled decay and compares it to the digitized Fig. 4B
  observation at 400 C and 9.41 mM Zn(II) (agreement within ~4 %).

- **`eflibe_f2_yield/davis_2022_nse`** recovers the radiolytic F2 G value from the modeled production
  rate and compares it to the Davis et al. Table III FLiBe-UF4 measurement (G = 0.005).

## Running

```
./run_tests --spec-file validation -i validation/msr
```

A complete multi-dimensional application of the capability — radiation-driven Cr corrosion in a
flowing salt with a spatially varying dose — is provided as an example in
`examples/molten_salt_corrosion`.
