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
  [corrosion_depth_um]
    type = MooseVariableFVReal
  []
[]

[AuxVariables]
  [igc_depth_um]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  [accumulated_depth]
    type = ParsedFunction
    expression = 'if(x < 0.5, 4, 9)'
  []
[]

[FVICs]
  [corrosion_depth]
    type = FVFunctionIC
    variable = corrosion_depth_um
    function = accumulated_depth
  []
[]

[AuxKernels]
  [igc_depth]
    type = MoltenSaltIGCDepthAux
    variable = igc_depth_um
    corrosion_depth = corrosion_depth_um
    salt_class = fluoride_fuel
    redox_class = purified_baseline
    execute_on = INITIAL
  []
[]

[Postprocessors]
  [corrosion_depth_low]
    type = PointValue
    variable = corrosion_depth_um
    point = '0.25 0 0'
    execute_on = INITIAL
  []
  [corrosion_depth_high]
    type = PointValue
    variable = corrosion_depth_um
    point = '0.75 0 0'
    execute_on = INITIAL
  []
  [igc_depth_low]
    type = PointValue
    variable = igc_depth_um
    point = '0.25 0 0'
    execute_on = INITIAL
  []
  [igc_depth_high]
    type = PointValue
    variable = igc_depth_um
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
