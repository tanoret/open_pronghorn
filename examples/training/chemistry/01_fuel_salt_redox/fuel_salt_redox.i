# =====================================================================================================
# Training case 01: fission oxidation of an MSRE fuel salt, the U(IV)/U(III) redox buffer and the
# cover gas.
#
# Salt: MSRE fuel salt 65LiF-29.1BeF2-5ZrF4-0.9UF4 at 600 C (2.27 g/cm^3): F- ~ 8.0e4 mol/m^3,
# U ~ 490 mol/m^3. MSRE kept U(III)/U ~ 1 % (4.9 mol/m^3 of U(III)) to protect the Hastelloy N.
#
# Where the oxidation comes from. Radiolysis only moves charge around: every F atom it makes has a
# matching electron (G(e) = G(F)), and in the molten salt they recombine (net G(F2) ~ 0; MSRE saw F2
# only from frozen salt). The NET oxidant in a fuel salt is FISSION: a fissioned UF4 releases four
# fluorides, but the fission products bind fewer of them (~3.2), so each fission leaves ~0.8
# oxidizing equivalents in the salt. With F = P / (200 MeV N_A) fissions [mol/m^3/s],
#   R_ox = 0.8 F = 0.8 * 4e6 / (200 MeV * N_A) = 1.66e-7 eq/m^3/s     (MSRE loop average 4 MW/m^3)
# The excess fluorine is modelled explicitly and charge-consistently as nascent F2, i.e. a fluoride
# sink -R_ox on F- and an F2 source +R_ox/2 (fluorine and charge are both conserved), which is then
# scavenged by the existing U(III) chemistry of the database:
#   F2 + 2 U(III) -> 2 U(IV) + 2 F-                (k = 5e4 m^6/mol^2/s)
# Net effect: U(III) -> U(IV) at R_ox, fluoride unchanged. The 4.9 mol/m^3 buffer lasts
#   t_buffer = 4.9 / 1.66e-7 s = 342 days,
# the months-to-year time scale on which MSRE added beryllium to restore U(III):
#   Be + 2 U(IV) -> Be(2+) + 2 U(III)              ('Be_rate', 0 by default; exercise 3)
# Only after the buffer is spent does F2 survive and pass to the cover gas by Henry's law:
#   d C_liq / dt = -kLa (C_liq - kH R T C_gas),  d C_gas / dt = (V_liq / V_gas) kLa (C_liq - kH R T C_gas)
# kH(F2) = 1e-6 mol/(m^3 Pa) is an ASSUMED value (no measurement; F2 is reactive), kH R T = 7.26e-3.
# kLa = 1e-3 1/s (loop-averaged, MSRE pump-bowl stripping) and V_liq/V_gas = 20 (2 m^3 of salt under
# a ~0.1 m^3 gas space) are MSRE-like; the headspace here is closed (no He purge).
#
# Radiolysis: the action still builds e_sol, F_rad and their reactions (G = 0.2 each, assumed, gamma
# values). The base run uses dose_rate = 0 because the radiolytic e/F cycle is closed (net zero) and
# turns over ~4e5 times faster than fission oxidation, so it would force second-scale steps over a
# year. Exercise 1 switches it on for one hour (D ~ 0.875 x 4 MW/m^3 = 3.5e6 W/m^3) to show that it
# does not consume the buffer. G(F2) is set to 0 (charge-balanced molten-salt radiolysis).
#
# Results (gold, 450 days): U(III)/U falls linearly from 1 % (0.708 % at day 100, matching the
# analytic line to 1e-12) and drops below 0.01 % at day 339 (analytic 342); then F2 accumulates
# (0.098 mol/m^3 in the salt, 13.4 mol/m^3 in the gas at day 450) with henry_ratio = 1.001.
# Max errors: uranium 1.1e-8, fluorine 1.2e-7, redox 1.9e-3 (after the buffer is spent).
#
# VERIFICATION invariants:
#   * uranium:   U(III) + U(IV) = U_0 (the model does not burn uranium; at 4 MW/m^3 burnup would be
#                ~1.6 % of U in 450 days)                                         ('uranium_error')
#   * fluorine:  F- + F_rad + 2 (F2_liq + F2_gas / (V_l/V_g)) = F_0              ('fluorine_error')
#   * redox:     (U(III) + e_sol) - (F_rad + 2 F2_liq + 2 F2_gas / (V_l/V_g))
#                    = U(III)_0 + 2 Be_rate t - R_ox t                            ('redox_error')
# =====================================================================================================

T = 873.15 # K
F0 = 8.0e4 # mol/m^3 fluoride in MSRE fuel salt
U_tot = 490.0 # mol/m^3 uranium (0.9 mol% UF4)
U_III_0 = 4.9 # mol/m^3 (U(III)/U = 1 %)
U_IV_0 = ${fparse U_tot - U_III_0}
P_fission = 4.0e6 # W/m^3, MSRE loop-averaged fission power density (8 MW / ~2 m^3)
eq_per_fission = 0.8 # oxidizing equivalents left in the salt per fission
R_ox = ${fparse eq_per_fission * P_fission / (200e6 * 1.602176634e-19) / 6.02214076e23} # eq/m^3/s
Be_rate = 0 # mol Be/m^3/s added (each Be reduces 2 U(IV))
dose_rate = 0 # W/m^3 radiolysis (exercise 1: 3.5e6, i.e. ~0.875 x P_fission)
G_F = 0.2 # F_rad yield of the database (molecules / 100 eV, assumed gamma value)
kLa = 1.0e-3 # 1/s, loop-averaged gas-liquid mass transfer (MSRE-like)
volume_ratio = 20 # V_liq / V_gas (MSRE-like)
kH = 1.0e-6 # mol/(m^3 Pa), F2 Henry coefficient of the database (assumed)

[Mesh]
  [cell]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys U_IV_sys U_III_sys F2_gas_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = fluoride
  metals = 'U'
  gas_species = 'F2'
  temperature = ${T}
  dose_rate = ${dose_rate}
  # Molten-salt radiolysis is charge balanced: no direct F2 yield without its electrons
  g_value_species = 'F2_diss'
  g_value_overrides = '0'
  kLa = ${kLa}
  volume_ratio = ${volume_ratio}
  initial_condition_species = 'F_ion U_IV U_III'
  initial_condition_values = '${F0} ${U_IV_0} ${U_III_0}'
[]

[LinearFVKernels]
  # Fission oxidation: the excess fluoride of each fissioned UF4 becomes nascent F2 (F- -> 1/2 F2 + e,
  # the electron being taken up by the fission products). Fluorine and charge are conserved.
  [fission_F_ion_sink]
    type = LinearFVSource
    variable = F_ion
    source_density = ${fparse -R_ox}
  []
  [fission_F2_source]
    type = LinearFVSource
    variable = F2_diss
    source_density = ${fparse R_ox / 2}
  []
  # The database radiolysis source creates F_rad (G = 0.2) without removing the F- it comes from;
  # this sink closes the fluorine balance when the dose is on (zero in the base run).
  [radiolysis_F_ion_sink]
    type = LinearFVSource
    variable = F_ion
    source_density = ${fparse -G_F * dose_rate / (100 * 1.602176634e-19) / 6.02214076e23}
  []
  # Beryllium addition: Be + 2 U(IV) -> Be(2+) + 2 U(III)
  [Be_U_III]
    type = LinearFVSource
    variable = U_III
    source_density = ${fparse 2 * Be_rate}
  []
  [Be_U_IV]
    type = LinearFVSource
    variable = U_IV
    source_density = ${fparse -2 * Be_rate}
  []
[]

[Postprocessors]
  [c_e_sol]
    type = ElementAverageValue
    variable = e_sol
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_F_rad]
    type = ElementAverageValue
    variable = F_rad
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_F_ion]
    type = ElementAverageValue
    variable = F_ion
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_F2_diss]
    type = ElementAverageValue
    variable = F2_diss
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_F2_gas]
    type = ElementAverageValue
    variable = F2_gas
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_U_III]
    type = ElementAverageValue
    variable = U_III
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_U_IV]
    type = ElementAverageValue
    variable = U_IV
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [U_III_fraction]
    type = ParsedPostprocessor
    expression = 'c_U_III / (c_U_III + c_U_IV)'
    pp_names = 'c_U_III c_U_IV'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Analytic buffer for comparison: U(III)_0 + (2 Be_rate - R_ox) t, floored at zero
  [U_III_analytic]
    type = ParsedPostprocessor
    expression = 'max(${U_III_0} + (${fparse 2 * Be_rate - R_ox}) * t, 0)'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Gas-liquid equilibrium indicator: 1 when C_liq = kH R T C_gas
  [henry_ratio]
    type = ParsedPostprocessor
    expression = 'if(c_F2_gas > 0, c_F2_diss / (${kH} * 8.314462618 * ${T} * c_F2_gas), 0)'
    pp_names = 'c_F2_diss c_F2_gas'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [uranium_error]
    type = ParsedPostprocessor
    expression = 'abs(c_U_III + c_U_IV - ${U_tot}) / ${U_tot}'
    pp_names = 'c_U_III c_U_IV'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [fluorine_error]
    type = ParsedPostprocessor
    expression = 'abs(c_F_ion + c_F_rad + 2 * (c_F2_diss + c_F2_gas / ${volume_ratio}) - ${F0}) / ${F0}'
    pp_names = 'c_F_ion c_F_rad c_F2_diss c_F2_gas'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Redox balance, relative to the initial buffer
  [redox_error]
    type = ParsedPostprocessor
    expression = 'abs(c_U_III + c_e_sol - c_F_rad - 2 * (c_F2_diss + c_F2_gas / ${volume_ratio})'
                 ' - ${U_III_0} - (${fparse 2 * Be_rate - R_ox}) * t) / ${U_III_0}'
    pp_names = 'c_U_III c_e_sol c_F_rad c_F2_diss c_F2_gas'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 30
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys U_IV_sys U_III_sys F2_gas_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  end_time = 3.888e7 # 450 days
  # Variable time step: grows from 1 s to dtmax
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.4
  []
  # Time-step study (max redox_error): 1.9e-3 at dtmax = 12 h, 2e-5 at 6 h (2x runtime), 2.3e-2 at 1 d.
  # The error comes from the lagged F2 / U(III) coupling when the buffer runs out.
  dtmax = 43200 # 12 h
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
