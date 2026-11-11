# 3-D part of channel_turbulence.i: include after it. The third velocity component in the
# turbulence objects and the momentum terms, and the solid walls of the channels (the axis
# cylinder is a symmetry surface, as are the cut planes).
channel_walls = 'outer_wall closed_ends'
segregated_plug_start = '0 segregated_plug_velocity 0'

[AuxKernels]
  [mu_t]
    u = vel_x_channel
    v = vel_y_channel
    w = vel_z_channel
  []
  [yplus]
    u = vel_x
    v = vel_y
    w = vel_z
  []
[]

[LinearFVKernels]
  [TKE_source_sink]
    u = vel_x_channel
    v = vel_y_channel
    w = vel_z_channel
  []
  [TKED_source_sink]
    u = vel_x_channel
    v = vel_y_channel
    w = vel_z_channel
  []
  [channel_stress_x]
    w = vel_z
    gradient_functors := 'vel_x_channel vel_y_channel vel_z_channel'
  []
  [channel_stress_y]
    w = vel_z
    gradient_functors := 'vel_x_channel vel_y_channel vel_z_channel'
  []
  [channel_stress_z]
    type = NTPLinearFVViscousStressCorrection
    variable = vel_z
    momentum_component = z
    u = vel_x
    v = vel_y
    w = vel_z
    mu = mu_eff
    gradient_functors = 'vel_x_channel vel_y_channel vel_z_channel'
    block = 'plenum chamber'
  []
  [channel_mu_t_diffusion_z]
    type = LinearFVDiffusion
    variable = vel_z
    diffusion_coeff = mu_t_all
    coeff_interp_method = channel_harmonic
    use_nonorthogonal_correction = false
  []
[]

[LinearFVBCs]
  [mu_t_walls]
    u = vel_x
    v = vel_y
    w = vel_z
  []
[]
