# Minimal valid direct CorrosionRateAux input used by the exclusivity error regressions.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Variables]
  [c]
  []
[]

[AuxVariables]
  [rate]
  []
[]

[AuxKernels]
  [rate_aux]
    type = CorrosionRateAux
    variable = rate
    boundary = left
    concentration = c
    valence = 2
    molar_mass = 51.9961
    density = 8
    exchange_current_density = 1
    temperature = 923.15
    execute_on = INITIAL
  []
[]

[Executioner]
  type = Steady
[]
