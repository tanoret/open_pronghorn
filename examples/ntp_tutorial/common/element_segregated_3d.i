# Three-dimensional part of the segregated element model (05_element_3d/element_mesh_3d.i, the
# element axis along y); include after element_physics_segregated.i. The including input sets,
# before including the model:
#
#   segregated_momentum_systems = 'u_system v_system w_system'
#   segregated_velocities = 'vel_x vel_y vel_z'
#   segregated_inlet_velocity = '0 u_inlet 0'
#   segregated_inlet_fraction = 0.5          (the mesh is half of the element)
#   segregated_walls = 'axis outer_wall closed_ends symmetry_0 symmetry_pi'
#   segregated_wall_types = 'symmetry noslip noslip symmetry symmetry'
#   segregated_wall_energy_types = 'heatflux heatflux heatflux heatflux heatflux'
#   segregated_wall_energy_functors = '0 0 0 0 0'
#   segregated_initial_velocity = '0 segregated_plug_velocity 0'   (SIMPLE; PIMPLE may start at rest)

[UserObjects]
  [ins_rhie_chow_interpolator]
    w = vel_z
  []
[]

[Executioner]
  rhie_chow_user_object = ins_rhie_chow_interpolator
  momentum_systems = 'u_system v_system w_system'
  pressure_system = pressure_system
  energy_system = energy_system
  solid_energy_system = solid_energy_system
[]
