# Ru-103 produced in a 1 m slab and plated on both walls with J = k_m N. With a large diffusivity
# the slab is well mixed, so the steady inventory approaches Y F / (N_A (lambda + 2 k_m / L)).

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 20
  []
[]
[Problem]
  linear_sys_names = 'Ru103_sys'
[]
[FissionProductInventory]
  nuclides = 'Ru103'
  fission_rate_density = 6.2e17
  diffusivity = 1
  plating_boundaries = 'left right'
  plating_mass_transfer_coefficient = 1e-4
[]
[Postprocessors]
  [N]
    type = ElementAverageValue
    variable = Ru103
  []
[]
[Executioner]
  type = Transient
  system_names = 'Ru103_sys'
  dt = 1e5
  num_steps = 200
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]
[Outputs]
  csv = true
[]
