# Gauge-invariance regression for the Action-generated default external metal potential. This uses
# the reduced model so the same shared potential wiring is covered independently of MSTDB fixtures.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  solve = false
[]

[CorrosionPlating]
  topology = salt_only
  reaction_boundary = left
  elements = 'Cr'
  temperature = 923.15
  reference_temperature = 923.15
  applied_overpotential = 0.1
  solve_potential = true
  pin_potential_boundary = right
  pin_potential_value = 0.0375
  transient = true
  verbose = true
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
[]
