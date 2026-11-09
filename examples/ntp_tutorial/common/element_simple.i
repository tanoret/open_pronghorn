# SIMPLE (steady) executioner for the segregated element model; include after
# element_physics_segregated.i and its dimension file (which name the linear systems), instead of
# element_pimple.i.
#
# A SIMPLE solve converges the flow, the gas energy and the particles with the lagged state held
# fixed; the fixed-point iterations update the lagged state (and an orifice design) between solves.
# Parameters the including input must set:
#   simple_iterations     largest number of SIMPLE iterations per solve, e.g. 300
#   simple_fixed_point    number of fixed-point iterations (solves), e.g. 40

[Executioner]
  type = SIMPLE
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  energy_equation_relaxation = 1.0
  num_iterations = ${simple_iterations}
  continue_on_max_its = true

  pressure_absolute_tolerance = 1e-8
  momentum_absolute_tolerance = 1e-8
  energy_absolute_tolerance = 1e-8
  solid_energy_absolute_tolerance = 1e-8
  momentum_l_abs_tol = 1e-13
  pressure_l_abs_tol = 1e-13
  energy_l_abs_tol = 1e-13
  solid_energy_l_abs_tol = 1e-13
  momentum_l_tol = 0
  pressure_l_tol = 0
  energy_l_tol = 0
  solid_energy_l_tol = 0
  momentum_petsc_options_iname = '-pc_type'
  momentum_petsc_options_value = 'lu'
  pressure_petsc_options_iname = '-pc_type'
  pressure_petsc_options_value = 'lu'
  energy_petsc_options_iname = '-pc_type'
  energy_petsc_options_value = 'lu'
  solid_energy_petsc_options_iname = '-pc_type'
  solid_energy_petsc_options_value = 'lu'

  # the lagged state (and an orifice design) is updated after every solve
  fixed_point_max_its = ${simple_fixed_point}
  accept_on_max_fixed_point_iteration = true
  disable_fixed_point_residual_norm_check = true
[]
