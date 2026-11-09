# Axisymmetric (r, z) part of the segregated element model; include after
# element_physics_segregated.i. The including input sets, before including the model:
#
#   segregated_momentum_systems = 'u_system v_system'
#   segregated_velocities = 'vel_x vel_y'
#   segregated_inlet_velocity = '0 u_inlet'
#   segregated_inlet_fraction = 1
#   segregated_walls = 'axis outer_wall closed_ends'
#   segregated_wall_types = 'symmetry noslip noslip'
#   segregated_wall_energy_types = 'heatflux heatflux heatflux'
#   segregated_wall_energy_functors = '0 0 0'
#   segregated_initial_velocity = '0 0 0'

[LinearFVKernels]
  # The segregated flow Physics does not add the axisymmetric hoop term of the radial momentum
  # equation (mu u_r / r^2); without it the radial momentum next to the axis is ill posed.
  [radial_hoop_stress]
    type = LinearFVRZViscousSource
    variable = vel_x
    mu = mu
    momentum_component = x
  []
[]

[Executioner]
  rhie_chow_user_object = ins_rhie_chow_interpolator
  momentum_systems = 'u_system v_system'
  pressure_system = pressure_system
  energy_system = energy_system
  solid_energy_system = solid_energy_system
[]
