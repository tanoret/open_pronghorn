# PIMPLE (pseudo-transient) executioner for the segregated element model; include after
# element_physics_segregated.i and its dimension file (which name the linear systems), instead of
# element_simple.i.
#
# The flow Physics add the time derivatives rho du/dt and rho dh/dt when the executioner is
# transient; the particle temperature gets one here. The segregated pressure equation has no
# compressibility term, so the march is a pseudo-transient towards the steady state, not a physical
# transient: each step runs a few PIMPLE iterations, the lagged density state (p_rho, T_rho) is
# updated at the end of the step, and the step grows while the iterations converge. The run stops
# when the solution no longer changes (steady_state_detection).
#
# Parameters the including input must set:
#   rho_cp_solid   pseudo heat capacity of the particles [J/(m^3 K)], e.g. 5e4: shapes the path only
#   pimple_dt      first time step [s], e.g. 1e-3
#   pimple_dtmax   largest time step [s]
#   pimple_end     end time [s]
#   pimple_steady_state_detection   stop once the solution no longer changes (true or false)

[LinearFVKernels]
  [particle_storage]
    type = LinearFVTimeDerivative
    variable = T_solid
    factor = ${rho_cp_solid}
  []
[]

[Executioner]
  type = PIMPLE

  momentum_equation_relaxation = 0.8
  pressure_variable_relaxation = 0.4
  energy_equation_relaxation = 1.0
  num_iterations = 30
  num_piso_iterations = 0
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

  end_time = ${pimple_end}
  dtmax = ${pimple_dtmax}
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = ${pimple_dt}
    optimal_iterations = 15
    iteration_window = 5
    growth_factor = 1.5
    cutback_factor = 0.5
  []
  steady_state_detection = ${pimple_steady_state_detection}
  steady_state_tolerance = 1e-8
[]
