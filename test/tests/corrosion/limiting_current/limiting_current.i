# Mechanistic verification of the diffusion-limited (limiting) current.
#
# A 1D salt cell of length L has a fast cathodic Butler-Volmer electrode at x = 0 (consuming the
# cation) and a bulk reservoir held at c_bulk at x = L. When the electrode kinetics are fast, the
# surface concentration is driven toward zero and the steady current saturates at the analytic
# diffusion-limited value
#       i_lim = z F D c_bulk / L .
# The test reports the surface concentration and the relative error of the computed current against
# i_lim; with fast kinetics the surface concentration is near zero and the error is below 0.1%.

L = 1.0e-3        # cell length [m]
D = 1.0e-9        # cation diffusivity [m^2/s]
c_bulk = 1.0      # bulk concentration [mol/m^3]
z = 2             # cation charge number

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 100
    xmax = ${L}
  []
[]

[Variables]
  [c_Cr]
    initial_condition = ${c_bulk}
  []
[]

[Kernels]
  [diffusion]
    type = CorrosionNernstPlanckFlux
    variable = c_Cr
    diffusivity = ${D}
    valence = ${z}
    temperature = 923.15
  []
[]

[BCs]
  [bulk]
    type = ADDirichletBC
    variable = c_Cr
    boundary = right
    value = ${c_bulk}
  []
  [electrode]
    type = ButlerVolmerBC
    variable = c_Cr
    boundary = left
    flux_type = species
    metal_domain = false
    valence = ${z}
    exchange_current_density = 1.0
    alpha_a = 0.5
    alpha_c = 1.0
    temperature = 923.15
    c_ref = 1.0
    # The metal is held strongly cathodic so the electrode consumes the cation.
    applied_potential = -0.3
  []
[]

[Postprocessors]
  [c_surface]
    type = SideAverageValue
    variable = c_Cr
    boundary = left
  []
  [i_limiting]
    type = ParsedPostprocessor
    expression = '${z} * 96485.33212 * ${D} * ${c_bulk} / ${L}'
  []
  [i_actual]
    type = ParsedPostprocessor
    pp_names = 'c_surface'
    expression = '${z} * 96485.33212 * ${D} * (${c_bulk} - c_surface) / ${L}'
  []
  [relative_error]
    type = ParsedPostprocessor
    pp_names = 'c_surface'
    expression = 'c_surface / ${c_bulk}'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = basic
  nl_abs_tol = 1.0e-12
[]

[Outputs]
  csv = true
[]
