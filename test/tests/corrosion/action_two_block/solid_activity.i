# Two-block coupon with the alloy-surface activity (solid_activity = true): the anodic
# Butler-Volmer branch is scaled by c_s / c_s,0 at the interface. With a slow solid diffusivity the
# surface depletes, so the dissolution slows and the surface concentration stays positive, while the
# total chromium is still conserved. Compare with two_block.i (unit metal activity).

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 40
    xmax = 2.0
  []
  [solid_block]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    bottom_left = '1.0 0 0'
    top_right = '2.0 0 0'
  []
  [names]
    type = RenameBlockGenerator
    input = solid_block
    old_block = '0 1'
    new_block = 'salt solid'
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = names
    primary_block = salt
    paired_block = solid
    new_boundary = interface
  []
[]

[CorrosionPlating]
  topology = two_block
  salt_block = 'salt'
  solid_block = 'solid'
  interface_boundary = 'interface'
  elements = 'Cr'
  temperature = 923.15
  reference_temperature = 923.15
  material_class = hastelloy_n
  salt_class = fluoride_fuel
  redox_class = oxidizing_fef2
  applied_overpotential = 0.1
  default_salt_concentration = 1.0
  default_solid_concentration = 1.0
  solid_diffusivity = 1.0e-9
  solid_activity = true
  transient = true
[]

[Postprocessors]
  [cr_salt]
    type = ElementIntegralVariablePostprocessor
    variable = c_Cr
    block = salt
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_solid]
    type = ElementIntegralVariablePostprocessor
    variable = cs_Cr
    block = solid
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_total]
    type = LinearCombinationPostprocessor
    pp_names = 'cr_salt cr_solid'
    pp_coefs = '1.0 1.0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Lowest solid concentration (the surface cell): must stay positive
  [surface_cs]
    type = ElementExtremeValue
    variable = cs_Cr
    block = solid
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [recession_um]
    type = SideAverageValue
    variable = recession_um
    boundary = interface
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 5.0e4
  num_steps = 4
  automatic_scaling = true
  nl_abs_tol = 1.0e-12
[]

[Outputs]
  csv = true
[]
