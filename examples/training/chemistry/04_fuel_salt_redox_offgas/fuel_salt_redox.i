# =====================================================================================================
# Training case 04: radiolytic F2, the U(IV)/U(III) redox buffer and the cover-gas headspace.
#
# In a fluoride fuel salt (FLiBe-UF4, 600 C) gamma radiolysis creates electrons and fluorine atoms.
# The salt chemistry is controlled by the uranium couple, which MSRE kept at U(III)/U_total ~ 1 %
# to protect the Hastelloy N:
#   e_sol + U(IV)       -> U(III)                (reduces the salt)
#   F2    + 2 U(III)    -> 2 U(IV) + 2 F-        (U(III) scavenges radiolytic fluorine)
#   2 F_rad -> F2,  e_sol + F_rad -> F-,  plus a direct F2 yield G = 0.005
# and dissolved F2 exchanges with the cover gas through Henry's law
#   d C_liq / dt = -kLa (C_liq - kH R T C_gas),  d C_gas / dt = (V_liq / V_gas) kLa (C_liq - kH R T C_gas)
#
# Run the base case, then the exercise with 'U_III_0=0' (no redox buffer): without U(III) the
# fluorine is not scavenged and accumulates in the salt and in the headspace until gas-liquid
# equilibrium C_liq = kH R T C_gas (kH R T = 7.26e-3 at 600 C) is reached.
#
# Invariants checked in the postprocessors:
#   * uranium conservation: U(III) + U(IV) = U_0
#   * redox balance. Radiolysis makes one electron per F atom (G = 0.2 each) plus a direct F2 yield
#     (G = 0.005, i.e. 0.01 oxidizing equivalents per 100 eV with no matching electron). Counting
#     oxidizing equivalents (F_rad, 2 F2 in salt and gas) minus reducing ones (e_sol, U(III) gained),
#       ox - red = 2 G_F2 D t / (100 eV N_A),
#     so the salt is slowly oxidized at ~1e-3 mol/m^3/s: the 1 % U(III) buffer is used up in ~1.6 h
#     and only then does F2 build up in the salt and the cover gas. (F- is a 1e5 mol/m^3 reservoir
#     that the source does not deplete, so it is not a conserved quantity here.)
# =====================================================================================================

dose_rate = 1.0e6 # W/m^3, gamma heating in the core
T = 873.15 # K
F0 = 1.0e5 # mol/m^3 fluoride in FLiBe
U_IV_0 = 600.0 # mol/m^3 (~ 1 mol% UF4)
U_III_0 = 6.0 # mol/m^3 (U(III)/U = 1 %)
kLa = 0.1 # 1/s, gas-liquid mass transfer (helium bubbling)
volume_ratio = 0.2 # V_liq / V_gas
G_F2 = 0.005 # direct F2 yield (molecules / 100 eV)

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
  kLa = ${kLa}
  volume_ratio = ${volume_ratio}
  initial_condition_species = 'F_ion U_IV U_III'
  initial_condition_values = '${F0} ${U_IV_0} ${U_III_0}'
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
  # Gas-liquid equilibrium indicator: 1 when C_liq = kH R T C_gas
  [henry_ratio]
    type = ParsedPostprocessor
    expression = 'if(c_F2_gas > 0, c_F2_diss / (1e-6 * 8.314462618 * ${T} * c_F2_gas), 0)'
    pp_names = 'c_F2_diss c_F2_gas'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [uranium_error]
    type = ParsedPostprocessor
    expression = 'abs(c_U_III + c_U_IV - ${fparse U_IV_0 + U_III_0}) / ${fparse U_IV_0 + U_III_0}'
    pp_names = 'c_U_III c_U_IV'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [net_oxidation]
    type = ParsedPostprocessor
    expression = '2 * ${G_F2} * ${dose_rate} / (100 * 1.602176634e-19) / 6.02214076e23 * t'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [redox_error]
    type = ParsedPostprocessor
    expression = 'if(net_oxidation > 0, abs(c_F_rad + 2 * (c_F2_diss + c_F2_gas / ${volume_ratio})'
                 ' - c_e_sol - (c_U_III - ${U_III_0}) - net_oxidation) / net_oxidation, 0)'
    pp_names = 'c_F_rad c_F2_diss c_F2_gas c_e_sol c_U_III net_oxidation'
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
  end_time = 2.16e4 # 6 hours
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-6
    growth_factor = 1.4
  []
  # The F2 / U(III) scavenging couples through lagged partners in the fixed point: once the buffer
  # is gone the redox balance error is ~ 2e-4 at dtmax = 10 s (1.4e-2 at 60 s, 2e-6 at 2 s)
  dtmax = 10
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
