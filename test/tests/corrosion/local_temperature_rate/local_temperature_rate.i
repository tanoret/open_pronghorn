[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 2
  []
[]

# This is an initialization-only diagnostic test; no nonlinear solve is needed.
[Problem]
  solve = false
[]

[Variables]
  [temperature]
    type = MooseVariableFVReal
  []
[]

[AuxVariables]
  [corrosion_rate_um_y]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  [local_temperature]
    type = ParsedFunction
    expression = 'if(x < 0.5, 850, 1050)'
  []
[]

[FVICs]
  [temperature]
    type = FVFunctionIC
    variable = temperature
    function = local_temperature
  []
[]

[AuxKernels]
  [corrosion_rate]
    type = MoltenSaltCorrosionRateAux
    variable = corrosion_rate_um_y
    temperature = temperature
    material_class = hastelloy_n
    salt_class = fluoride_fuel
    redox_class = purified_baseline
    flow_factor = 0.75
    delta_T_C = 50
    execute_on = INITIAL
  []
[]

[Postprocessors]
  [temperature_cold]
    type = PointValue
    variable = temperature
    point = '0.25 0 0'
    execute_on = INITIAL
  []
  [temperature_hot]
    type = PointValue
    variable = temperature
    point = '0.75 0 0'
    execute_on = INITIAL
  []
  [rate_cold]
    type = PointValue
    variable = corrosion_rate_um_y
    point = '0.25 0 0'
    execute_on = INITIAL
  []
  [rate_hot]
    type = PointValue
    variable = corrosion_rate_um_y
    point = '0.75 0 0'
    execute_on = INITIAL
  []
[]
[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
  execute_on = INITIAL
[]
