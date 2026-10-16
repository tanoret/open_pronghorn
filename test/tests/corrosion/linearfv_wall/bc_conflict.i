# Error check: the linear FV Butler-Volmer wall takes its overpotential either from an equilibrium
# concentration or from the metal/salt potentials, not both.
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 2
  []
[]
[Problem]
  linear_sys_names = 'c_sys'
[]
[Variables]
  [c]
    type = MooseLinearVariableFVReal
    solver_sys = c_sys
  []
[]
[LinearFVKernels]
  [diff]
    type = LinearFVDiffusion
    variable = c
    diffusion_coeff = 1
  []
[]
[LinearFVBCs]
  [wall]
    type = CorrosionLinearFVButlerVolmerBC
    variable = c
    boundary = left
    valence = 2
    exchange_current_density = 1e-3
    temperature = 900
    metal_potential = 0.1
    equilibrium_concentration = 10
  []
[]
[Executioner]
  type = Steady
  system_names = 'c_sys'
[]
