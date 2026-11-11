# Turbulence and the compressible viscous stress of the free channels (plenum, chamber) of the
# segregated element model, resolved radially: include after element_physics_segregated.i and its
# dimension file, with several cells across the channels (n_pl, n_ch, e.g. 8), and the dimension
# part channel_turbulence_rz.i or channel_turbulence_3d.i after this file.
#
# The plenum (Re ~ 1e5 - 1e6) and the chamber (Re ~ 1e5 - 5e5, Mach up to 0.2, the gas at 3000 K
# and a hundred times lighter than at the inlet) are turbulent. In them:
#  - k-epsilon of OpenPronghorn (kEpsilonTKESourceSink, kEpsilonTKEDSourceSink, kEpsilonViscosity),
#    variant 'channel_turbulence_variant', with the compressible production
#      G_k = mu_t 2 S:S - 2/3 (rho k div u + mu_t (div u)^2)
#    and the dilatational dissipation gamma_M = rho C_M k eps / c^2 (use_compressibility). The
#    density and viscosity are those of the momentum equation (the lagged state), the speed of
#    sound the hydrogen one at that state; the strain rate is the one of the channel flow
#    (NTPChannelVelocity: no gradient through the porous interfaces);
#  - wall functions on the solid walls of the channels, 'channel_walls' (the plenum outer wall and
#    the closed ends); the frit interfaces are permeable: k and epsilon leave the plenum with the
#    gas into the cold frit and enter the chamber from the hot frit;
#  - the momentum stress mu_eff (grad u + grad u^T - 2/3 div u I), mu_eff = mu + mu_t: the Physics
#    assembles mu grad u, channel_mu_t_diffusion adds mu_t grad u (with the wall-function value of
#    mu_t on the walls), NTPLinearFVViscousStressCorrection the transpose and dilatation part with
#    mu_eff; -2/3 rho k I is absorbed in the pressure;
#  - the turbulent conduction k_t = mu_t cp / Pr_t in the gas energy equation;
#  - no body-force duct friction in the channels: the walls carry the shear.
# The packed beds carry no RANS turbulence (mu_t = 0 there): their momentum and heat exchange are
# the packed-bed correlations.

channel_turbulence_variant = Standard
channel_wall_treatment = neq
channel_inflow_intensity = 0.05
channel_Pr_t = 0.9
channel_C_M = 1.0
# initial state: a cold-plenum guess (rho 72 kg/m^3, u ~ 20 m/s, 5 % intensity, l ~ 1 mm)
channel_k_init = 1.5
channel_eps_init = 1.2e3

# the walls carry the shear of the channels
f_ann := 0
channel_wall_friction := false
# a flowing start: the resolved channels diverge from rest
segregated_initial_velocity := '${segregated_plug_start}'

[Problem]
  linear_sys_names := '${segregated_momentum_systems} pressure_system energy_system
                       solid_energy_system TKE_system TKED_system'
[]

[UserObjects]
  # velocity with a cell gradient that ignores the porous neighbours (vel_x_channel, ...)
  [channel_velocity]
    type = NTPChannelVelocity
    velocity = '${segregated_velocities}'
    channel_blocks = 'plenum chamber'
  []
[]

[Variables]
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    block = 'plenum chamber'
    initial_condition = ${channel_k_init}
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    block = 'plenum chamber'
    initial_condition = ${channel_eps_init}
  []
[]

[AuxVariables]
  [mu_t]
    type = MooseLinearVariableFVReal
    block = 'plenum chamber'
    initial_condition = ${fparse 0.09 * 72 * channel_k_init^2 / channel_eps_init}
  []
  [wall_distance]
    type = MooseVariableFVReal
    block = 'plenum chamber'
    initial_condition = 1
  []
  [yplus]
    type = MooseVariableFVReal
    block = 'plenum chamber'
  []
[]

[AuxKernels]
  [wall_distance]
    type = WallDistanceAux
    variable = wall_distance
    walls = '${channel_walls}'
    block = 'plenum chamber'
    execute_on = 'INITIAL'
  []
  [mu_t]
    type = kEpsilonViscosity
    variable = mu_t
    block = 'plenum chamber'
    k_epsilon_variant = ${channel_turbulence_variant}
    tke = TKE
    epsilon = TKED
    rho = rho
    mu = mu
    walls = '${channel_walls}'
    wall_treatment = ${channel_wall_treatment}
    bulk_wall_treatment = false
    mu_t_ratio_max = 1e5
    execute_on = 'NONLINEAR'
  []
  [yplus]
    type = RANSYPlusAux
    variable = yplus
    block = 'plenum chamber'
    tke = TKE
    rho = rho
    mu = mu
    walls = '${channel_walls}'
    wall_treatment = ${channel_wall_treatment}
    execute_on = 'TIMESTEP_END'
  []
[]

[FunctorMaterials]
  # mu_t on every block (0 in the packed beds)
  [mu_t_porous]
    type = GenericFunctorMaterial
    prop_names = 'mu_t_all'
    prop_values = '0'
    block = 'cold_frit fuel_bed hot_frit'
  []
  [mu_t_channels]
    type = ParsedFunctorMaterial
    property_name = mu_t_all
    functor_names = mu_t
    expression = 'mu_t'
    block = 'plenum chamber'
  []
  [mu_eff]
    type = ParsedFunctorMaterial
    property_name = mu_eff
    functor_names = 'mu mu_t_all'
    expression = 'mu + mu_t_all'
  []
  # turbulent conduction of the enthalpy equation, k_t / cp = mu_t / Pr_t
  [k_t_by_cp]
    type = ParsedFunctorMaterial
    property_name = k_t_by_cp
    functor_names = 'mu_t_all'
    expression = 'mu_t_all / ${channel_Pr_t}'
  []
  # inflow turbulence: at the inlet, intensity channel_inflow_intensity of the inlet speed and a
  # length scale of a tenth of the plenum gap; from the hot frit, the same intensity of the speed
  # of the gas leaving the frit and the hot-frit particle diameter. k = 3/2 (I U)^2,
  # eps = C_mu^(3/4) k^(3/2) / l.
  [k_inlet]
    type = ParsedFunctorMaterial
    property_name = k_inlet
    functor_names = u_inlet
    expression = '1.5 * (${channel_inflow_intensity} * u_inlet)^2'
  []
  [eps_inlet]
    type = ParsedFunctorMaterial
    property_name = eps_inlet
    functor_names = k_inlet
    expression = '0.1643 * k_inlet^1.5 / ${fparse 0.1 * (r5_in - r4_in)}'
  []
  [k_frit]
    type = ParsedFunctorMaterial
    property_name = k_frit
    functor_names = 'G_x G_y G_z rho'
    expression = '1.5 * (${channel_inflow_intensity} * sqrt(G_x * G_x + G_y * G_y + G_z * G_z) / rho)^2
                  + 1e-8'
  []
  [eps_frit]
    type = ParsedFunctorMaterial
    property_name = eps_frit
    functor_names = k_frit
    expression = '0.1643 * k_frit^1.5 / ${d_hf}'
  []
[]

[FVInterpolationMethods]
  [channel_upwind]
    type = FVAdvectedUpwind
  []
  # The turbulent diffusion of momentum and enthalpy acts on every block, with mu_t = 0 in the
  # packed beds: harmonic faces carry no turbulent flux through a porous interface, and every face
  # flux enters its two cells with opposite signs. Restricted to the channel blocks, a diffusion
  # kernel of a variable defined on all blocks would add the flux of an interface face to the
  # channel cell only, a source.
  [channel_harmonic]
    type = FVHarmonicAverage
  []
[]

[LinearFVKernels]
  # ---- momentum: mu_t grad u (implicit, all blocks), and the transpose / dilatation stress with
  # mu_eff (channel faces only)
  [channel_mu_t_diffusion_x]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = mu_t_all
    coeff_interp_method = channel_harmonic
    use_nonorthogonal_correction = false
  []
  [channel_mu_t_diffusion_y]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = mu_t_all
    coeff_interp_method = channel_harmonic
    use_nonorthogonal_correction = false
  []
  [channel_stress_x]
    type = NTPLinearFVViscousStressCorrection
    variable = vel_x
    momentum_component = x
    u = vel_x
    v = vel_y
    mu = mu_eff
    gradient_functors = 'vel_x_channel vel_y_channel'
    block = 'plenum chamber'
  []
  [channel_stress_y]
    type = NTPLinearFVViscousStressCorrection
    variable = vel_y
    momentum_component = y
    u = vel_x
    v = vel_y
    mu = mu_eff
    gradient_functors = 'vel_x_channel vel_y_channel'
    block = 'plenum chamber'
  []
  # ---- gas energy: turbulent conduction
  [channel_turbulent_conduction]
    type = LinearFVDiffusion
    variable = h
    diffusion_coeff = k_t_by_cp
    coeff_interp_method = channel_harmonic
    use_nonorthogonal_correction = false
  []
  # ---- turbulent kinetic energy
  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
    rhie_chow_user_object = ins_rhie_chow_interpolator
    advected_interp_method_name = channel_upwind
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = mu
    use_nonorthogonal_correction = false
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = mu_t
    scaling_coeff = 1.0
    use_nonorthogonal_correction = false
  []
  [TKE_source_sink]
    type = kEpsilonTKESourceSink
    variable = TKE
    k_epsilon_variant = ${channel_turbulence_variant}
    epsilon = TKED
    rho = rho
    mu = mu
    mu_t = mu_t
    walls = '${channel_walls}'
    wall_treatment = ${channel_wall_treatment}
    use_compressibility = true
    speed_of_sound = speed_of_sound
    C_M = ${channel_C_M}
    C_pl = 10
  []
  # ---- dissipation rate
  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    rhie_chow_user_object = ins_rhie_chow_interpolator
    advected_interp_method_name = channel_upwind
    walls = '${channel_walls}'
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = mu
    use_nonorthogonal_correction = false
    walls = '${channel_walls}'
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = mu_t
    scaling_coeff = 1.3
    use_nonorthogonal_correction = false
    walls = '${channel_walls}'
  []
  [TKED_source_sink]
    type = kEpsilonTKEDSourceSink
    variable = TKED
    k_epsilon_variant = ${channel_turbulence_variant}
    tke = TKE
    rho = rho
    mu = mu
    mu_t = mu_t
    walls = '${channel_walls}'
    wall_treatment = ${channel_wall_treatment}
    use_compressibility = true
    speed_of_sound = speed_of_sound
    C_M = ${channel_C_M}
    C_pl = 10
  []
[]

[LinearFVBCs]
  # wall-function value of mu_t on the solid walls: the wall shear of the momentum equations
  [mu_t_walls]
    type = LinearFVTurbulentViscosityWallFunctionBC
    variable = mu_t
    boundary = '${channel_walls}'
    rho = rho
    mu = mu
    tke = TKE
    wall_treatment = ${channel_wall_treatment}
  []
  # k and epsilon leave the plenum with the gas, and enter the chamber from the hot frit
  [TKE_to_frit]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = TKE
    boundary = plenum_frit
    use_two_term_expansion = false
  []
  [TKED_to_frit]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = TKED
    boundary = plenum_frit
    use_two_term_expansion = false
  []
  [TKE_from_frit]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = TKE
    boundary = chamber_frit
    functor = k_frit
  []
  [TKED_from_frit]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = TKED
    boundary = chamber_frit
    functor = eps_frit
  []
  [TKE_inlet]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = TKE
    boundary = inlet
    functor = k_inlet
  []
  [TKED_inlet]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = TKED
    boundary = inlet
    functor = eps_inlet
  []
  [TKE_outlet]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = TKE
    boundary = outlet
    use_two_term_expansion = false
  []
  [TKED_outlet]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = TKED
    boundary = outlet
    use_two_term_expansion = false
  []
[]

[Executioner]
  # k and epsilon are solved after the energy in every SIMPLE / PIMPLE iteration
  turbulence_systems = 'TKE_system TKED_system'
  turbulence_equation_relaxation = '0.3 0.3'
  turbulence_field_relaxation = '0.3 0.3'
  # a stagnant corner (the closed end of the plenum) converges k and epsilon slowly: 1e-6 is a
  # converged flow (the inlet pressure no longer changes)
  turbulence_absolute_tolerance = '1e-6 1e-6'
  turbulence_l_abs_tol = 1e-13
  turbulence_l_tol = 0
  turbulence_petsc_options_iname = '-pc_type'
  turbulence_petsc_options_value = 'lu'
[]
