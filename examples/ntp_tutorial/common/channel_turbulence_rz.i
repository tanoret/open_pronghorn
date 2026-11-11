# RZ part of channel_turbulence.i: include after it. The velocities of the turbulence objects, the
# radial hoop stress with its compressible (2/3 div u) part and the turbulent viscosity, and the
# solid walls of the channels (the axis is a symmetry line).
channel_walls = 'outer_wall closed_ends'
segregated_plug_start = '0 segregated_plug_velocity'

[AuxKernels]
  [mu_t]
    u = vel_x_channel
    v = vel_y_channel
  []
  [yplus]
    u = vel_x
    v = vel_y
  []
[]

[LinearFVKernels]
  [TKE_source_sink]
    u = vel_x_channel
    v = vel_y_channel
  []
  [TKED_source_sink]
    u = vel_x_channel
    v = vel_y_channel
  []
  [radial_hoop_stress]
    use_deviatoric_terms = true
    mu = mu_eff
    u = vel_x
    v = vel_y
  []
[]

[LinearFVBCs]
  [mu_t_walls]
    u = vel_x
    v = vel_y
  []
[]
