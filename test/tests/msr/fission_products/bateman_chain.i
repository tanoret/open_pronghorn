# Decay of an initial Mo-99 inventory through Tc-99m to Tc-99 (no fission source, no flux).
# The systems are solved parent before daughter, so each implicit-Euler step is exact for the
# lower-triangular Bateman system in a single sweep:
#   Mo_{n+1}   = Mo_n / (1 + lambda_Mo dt)
#   Tcm_{n+1}  = (Tcm_n + dt b_m lambda_Mo Mo_{n+1}) / (1 + lambda_m dt)
#   Tc_{n+1}   = (Tc_n + dt (b_g lambda_Mo Mo_{n+1} + b lambda_m Tcm_{n+1})) / (1 + lambda_Tc dt)
# and the total inventory is conserved to round-off (Tc-99 decay to Ru-99 is ~1e-13 1/s).

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  linear_sys_names = 'Mo99_sys Tc99m_sys Tc99_sys'
[]

[FissionProductInventory]
  nuclides = 'Mo99 Tc99m Tc99'
  fission_rate_density = 0
  initial_condition_nuclides = 'Mo99'
  initial_condition_values = '1'
[]

[Postprocessors]
  [N_Mo99]
    type = ElementAverageValue
    variable = Mo99
  []
  [N_Tc99m]
    type = ElementAverageValue
    variable = Tc99m
  []
  [N_Tc99]
    type = ElementAverageValue
    variable = Tc99
  []
  [total]
    type = ParsedPostprocessor
    expression = 'N_Mo99 + N_Tc99m + N_Tc99'
    pp_names = 'N_Mo99 N_Tc99m N_Tc99'
  []
[]

[Executioner]
  type = Transient
  system_names = 'Mo99_sys Tc99m_sys Tc99_sys'
  dt = 3600
  num_steps = 100
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
