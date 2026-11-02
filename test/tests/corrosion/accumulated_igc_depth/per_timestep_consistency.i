[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Problem]
  solve = false
[]

[AuxVariables]
  [corrosion_rate_um_s]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 1
  []
  [corrosion_depth_um]
    family = MONOMIAL
    order = CONSTANT
  []
  [igc_depth_um]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [corrosion_depth]
    type = VariableTimeIntegrationAux
    variable = corrosion_depth_um
    variable_to_integrate = corrosion_rate_um_s
    order = 1
    execute_on = TIMESTEP_END
  []
  [igc_depth]
    type = MoltenSaltIGCDepthAux
    variable = igc_depth_um
    corrosion_depth = corrosion_depth_um
    salt_class = fluoride_fuel
    redox_class = purified_baseline
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [corrosion_depth]
    type = ElementAverageValue
    variable = corrosion_depth_um
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [igc_depth]
    type = ElementAverageValue
    variable = igc_depth_um
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  dt = 1
  num_steps = 3
[]

[Outputs]
  csv = true
  execute_on = 'INITIAL TIMESTEP_END'
[]
