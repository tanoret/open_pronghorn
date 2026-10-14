# =====================================================================================================
# Training case 05: fission-product build-up in a circulating fuel salt -- saturation vs growth.
#
# A well-mixed volume of fuel salt is irradiated for one year at a loop-averaged power density of
# 20 MW/m^3 (the order of a thermal MSR: 8 MW of fission power in ~0.4 m^3 of salt). For a nuclide
# with no tracked precursor, no capture (no 'neutron_flux' given) and no removal the inventory obeys
#
#   dN/dt = Y F / N_A - lambda N   ->   N(t) = Y F / (N_A lambda) (1 - exp(-lambda t))
#
# * I-131 (8 d) saturates within a few weeks: its activity = Y F (the decay rate equals the
#   production rate) is a measure of the power, not of the burnup;
# * Cs-137 and Sr-90 (30 y and 29 y) grow linearly over the year: they measure the burnup.
#
# F = P / E_f with E_f = 200 MeV; the action converts the power density itself. Yields: ENDF/B-VIII.0
# U-235 thermal cumulative (the precursors are not tracked, so they are lumped in).
# =====================================================================================================

power_density = 2.0e7 # W/m^3
day = 86400

[Mesh]
  [cell]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  linear_sys_names = 'I131_sys Cs137_sys Sr90_sys Kr85_sys'
[]

[FissionProductInventory]
  nuclides = 'I131 Cs137 Sr90 Kr85'
  power_density = ${power_density}
  # Print the yields, decay constants and the system ordering the action generated
  verbose = true
[]

[Postprocessors]
  [time_d]
    type = ParsedPostprocessor
    expression = 't / ${day}'
    use_t = true
  []
  # Inventories [mol/m^3]
  [N_I131]
    type = ElementAverageValue
    variable = I131
  []
  [N_Cs137]
    type = ElementAverageValue
    variable = Cs137
  []
  [N_Sr90]
    type = ElementAverageValue
    variable = Sr90
  []
  [N_Kr85]
    type = ElementAverageValue
    variable = Kr85
  []
  # Activities A = lambda N N_A [Bq/m^3]
  [A_I131]
    type = ParsedPostprocessor
    expression = 'log(2) / 693377.3 * N_I131 * 6.02214076e23'
    pp_names = 'N_I131'
  []
  [A_Cs137]
    type = ParsedPostprocessor
    expression = 'log(2) / 949252600.0 * N_Cs137 * 6.02214076e23'
    pp_names = 'N_Cs137'
  []
  # Production rate of I-131 [Bq/m^3 at saturation]: Y F
  [I131_saturation_activity]
    type = ParsedPostprocessor
    expression = '0.0289069 * ${power_density} / (200e6 * 1.602176634e-19)'
  []
  # Analytic Cs-137 inventory for comparison (continuous solution)
  [N_Cs137_analytic]
    type = ParsedPostprocessor
    expression = 'lam := log(2) / 949252600.0; 0.0618832 * ${power_density} / (200e6 * 1.602176634e-19) '
                 '/ 6.02214076e23 / lam * (1 - exp(-lam * t))'
    use_t = true
  []
  [N_I131_analytic]
    type = ParsedPostprocessor
    expression = 'lam := log(2) / 693377.3; 0.0289069 * ${power_density} / (200e6 * 1.602176634e-19) '
                 '/ 6.02214076e23 / lam * (1 - exp(-lam * t))'
    use_t = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'I131_sys Cs137_sys Sr90_sys Kr85_sys'
  # Second-order time integration: with one-day steps the 8-day I-131 build-up is within 0.5 % of
  # the analytic curve after a week and within 0.01 % once saturated (exercise: try dt = 6 h)
  scheme = bdf2
  dt = ${day}
  end_time = ${fparse 365 * day}
  # The inventories are small numbers in mol/m^3: converge the linear solves tightly
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
