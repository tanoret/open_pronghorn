# Noble-metal (A = 99) production in a 0D salt volume with plating on metal surfaces and
# entrainment into the off-gas. The plated and off-gas inventories keep decaying, so the chain
# Mo-99 -> Tc-99m -> Tc-99 continues on the walls and in the off-gas system. Mass balance:
#   d/dt sum_{salt, plated, offgas} (Mo99 + Tc99m + Tc99) = (Y_Mo + Y_Tcm + Y_Tc) F / N_A
#                                                          - lambda_Tc sum Tc99
# The volumetric plating rate is k_m A / V = 1e-5 m/s * 10 1/m = 1e-4 1/s.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  linear_sys_names = 'Mo99_sys Tc99m_sys Tc99_sys Mo99_offgas_sys Tc99m_offgas_sys Tc99_offgas_sys
                      Mo99_plated_sys Tc99m_plated_sys Tc99_plated_sys'
[]

[FissionProductInventory]
  nuclides = 'Mo99 Tc99m Tc99'
  power_density = 2e7
  removal_groups = 'noble_metal noble_metal'
  removal_rates = '1e-4 2e-5'
  removal_destinations = 'plated offgas'
[]

[Postprocessors]
  [salt]
    type = ParsedPostprocessor
    expression = 'a + b + c'
    pp_names = 'a b c'
  []
  [plated]
    type = ParsedPostprocessor
    expression = 'ap + bp + cp'
    pp_names = 'ap bp cp'
  []
  [offgas]
    type = ParsedPostprocessor
    expression = 'ao + bo + co'
    pp_names = 'ao bo co'
  []
  [a]
    type = ElementAverageValue
    variable = Mo99
    outputs = none
  []
  [b]
    type = ElementAverageValue
    variable = Tc99m
    outputs = none
  []
  [c]
    type = ElementAverageValue
    variable = Tc99
    outputs = none
  []
  [ap]
    type = ElementAverageValue
    variable = Mo99_plated
    outputs = none
  []
  [bp]
    type = ElementAverageValue
    variable = Tc99m_plated
    outputs = none
  []
  [cp]
    type = ElementAverageValue
    variable = Tc99_plated
    outputs = none
  []
  [ao]
    type = ElementAverageValue
    variable = Mo99_offgas
    outputs = none
  []
  [bo]
    type = ElementAverageValue
    variable = Tc99m_offgas
    outputs = none
  []
  [co]
    type = ElementAverageValue
    variable = Tc99_offgas
    outputs = none
  []
  [Tc99_total]
    type = ParsedPostprocessor
    expression = 'c + cp + co'
    pp_names = 'c cp co'
  []
[]

[Executioner]
  type = Transient
  system_names = 'Mo99_sys Tc99m_sys Tc99_sys Mo99_offgas_sys Tc99m_offgas_sys Tc99_offgas_sys
                  Mo99_plated_sys Tc99m_plated_sys Tc99_plated_sys'
  dt = 7200
  num_steps = 120
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
