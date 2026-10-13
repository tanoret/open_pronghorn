# Saturation of I-131 and Cs-137 from a constant fission rate (no flux, no removal):
#   N(t) -> Y F / (N_A lambda), implicit Euler: N_{n+1} = (N_n + dt Y F / N_A) / (1 + lambda dt)
# The cumulative yields are used since the precursors (Te-131, Xe-137) are not tracked.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  linear_sys_names = 'I131_sys Cs137_sys'
[]

[FissionProductInventory]
  nuclides = 'I131 Cs137'
  fission_rate_density = 6.2e17
[]

[Postprocessors]
  [N_I131]
    type = ElementAverageValue
    variable = I131
  []
  [N_Cs137]
    type = ElementAverageValue
    variable = Cs137
  []
[]

[Executioner]
  type = Transient
  system_names = 'I131_sys Cs137_sys'
  dt = 86400
  num_steps = 60
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
