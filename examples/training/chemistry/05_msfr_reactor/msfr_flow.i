# =====================================================================================================
# Training case 05 (capstone), step 1 of 2: the MSFR primary circuit -- flow, k-epsilon and energy.
#
# The EVOL / Virtual Test Bed MSFR (virtual_test_bed/msr/msfr, steady/run_ns.i): 2D-RZ toroidal core
# cavity (H ~ 2.6 m, R ~ 2.3 m) with one lumped external loop (hot leg -> pump -> HX -> cold leg) in
# revolution around the axis, LiF-ThF4-UF4 fuel salt, 3 GWth. The mesh (8272 cells) and the Griffin
# fission power density come from the VTB coupled steady state, copied here as msfr_vtb_mesh_power.e
# (= virtual_test_bed/msr/msfr/steady/restart/run_ns_coupled_restart.e).
# Differences from the VTB (nonlinear INSFV, capped mixing length):
#   * linear FV + SIMPLE so that the OpenPronghorn chemistry actions can be coupled (step 2),
#   * k-epsilon with wall functions instead of the capped mixing length,
#   * Forchheimer HX friction F |u| u linearised as a Darcy term with a lagged speed,
#   * pump force -24000 N/m^3 (VTB -20000) to recover the VTB mass flow rate with k-epsilon.
# Results (8 ranks, 1000 SIMPLE iterations, ~1 min; code-to-code BENCHMARK against the VTB, not a
# validation): mdot 19270 kg/s (VTB 19220), bulk hot / cold leg (HX inlet / outlet, mass-flow
# weighted) 997 / 902 K (VTB dT 97.6 K), Tmax 1095 K (VTB 1304 K, see README), |u|max 2.07 m/s
# (VTB 2.03), energy balance 2.2e-3. Uniform refinement (33088 cells): mdot +6.7 %, dT 89 K, Tmax 1102 K.
# Writes msfr_flow_restart.e, the converged fields read by msfr_chemistry.i.
# =====================================================================================================

# --- fuel salt LiF-ThF4-UF4 (VTB / Brovchenko 2019, @ ~1000 K) ---
rho = 4284            # kg/m^3
cp = 1594             # J/kg/K
drho_dT = 0.882       # kg/m^3/K
mu = 0.0166           # Pa s
k_salt = 1.7          # W/m/K
alpha_b = ${fparse drho_dT / rho}
Pr_t = 0.9

# --- VTB operating point ---
T_HX = 873.15         # K heat-exchanger secondary (ambient) temperature
h_hx = ${fparse 600 * 20e3}   # W/m^3/K volumetric HX conductance (VTB)
F_hx = 4.0e3          # kg/m^4 Forchheimer coefficient of the HX (VTB, non-standard formulation)
pump_force = -24000.  # N/m^3 pump body force (downward); VTB uses -20000 with the mixing length, tuned here
                      # to the VTB mass flow rate (~19200 kg/s) with k-epsilon wall functions

# --- k-epsilon ---
C_mu = 0.09
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
u_scale = 1.0
intensity = 0.05
k_init = '${fparse 1.5*(intensity*u_scale)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / 0.1}'
walls = 'shield_wall reflector_wall'
wall_treatment = 'eq_newton'
fluid = 'fuel pump hx'

[Mesh]
  [fmg]
    type = FileMeshGenerator
    # VTB mesh carrying the Griffin-Pronghorn coupled steady state (power density)
    file = 'msfr_vtb_mesh_power.e'
    use_for_exodus_restart = true
  []
  [hx_top]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y > 0'
    included_subdomains = 'hx'
    included_neighbors = 'fuel'
    fixed_normal = true
    normal = '0 1 0'
    new_sideset_name = 'hx_top'
    input = fmg
  []
  [hx_bot]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y < -0.6'
    included_subdomains = 'hx'
    included_neighbors = 'fuel'
    fixed_normal = true
    normal = '0 -1 0'
    new_sideset_name = 'hx_bot'
    input = hx_top
  []
  coord_type = RZ
  rz_coord_axis = Y
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system TKE_system TKED_system energy_system'
  previous_nl_solution_required = true
  allow_initial_conditions_with_restart = true  # only power_density is read from the file
[]

[FVInterpolationMethods]
  [upwind]
    type = FVAdvectedUpwind
  []
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method = 'upwind'
  advected_interp_method_name = 'upwind'
[]

[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
    block = ${fluid}
  []
[]

[Functions]
  [pump_profile]
    type = ConstantFunction
    value = ${pump_force}
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_condition = 1e-6
    block = ${fluid}
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_condition = 1e-6
    block = ${fluid}
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 0
    block = ${fluid}
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    initial_condition = ${k_init}
    block = ${fluid}
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    initial_condition = ${eps_init}
    block = ${fluid}
  []
  [T_fluid]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
    initial_condition = ${T_HX}
    block = ${fluid}
  []
[]

[AuxVariables]
  [wall_distance]
    type = MooseVariableFVReal
    initial_condition = 0.1
    block = ${fluid}
  []
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse rho * C_mu * k_init^2 / eps_init}'
    block = ${fluid}
  []
  [power_density]
    type = MooseVariableFVReal
    initial_from_file_var = power_density
    block = ${fluid}
  []
[]

[AuxKernels]
  [compute_wall_distance]
    type = WallDistanceAux
    variable = wall_distance
    walls = ${walls}
    execute_on = 'INITIAL'
    block = ${fluid}
  []
  [compute_mu_t]
    type = kEpsilonViscosityAux
    variable = mu_t
    C_mu = ${C_mu}
    tke = TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    bulk_wall_treatment = false
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    mu_t_ratio_max = 1e5
    execute_on = 'NONLINEAR'
    block = ${fluid}
  []
[]

[FunctorMaterials]
  [k_eff_mat]
    type = ParsedFunctorMaterial
    property_name = k_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${k_salt} + mut*${cp}/${Pr_t}'
    block = ${fluid}
  []
  [mu_eff_mat]
    type = ParsedFunctorMaterial
    property_name = mu_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${mu} + mut'
    block = ${fluid}
  []
  [alpha_mat]
    type = GenericFunctorMaterial
    prop_names = 'alpha_b'
    prop_values = '${alpha_b}'
    block = ${fluid}
  []
  # Lagged speed for the linearised Forchheimer friction F |u| u of the heat exchanger
  [speed_mat]
    type = ParsedFunctorMaterial
    property_name = speed_lag
    functor_names = 'vel_x vel_y'
    functor_symbols = 'u v'
    expression = 'sqrt(u*u + v*v)'
    block = 'hx'
  []
  [forch_mat]
    type = GenericVectorFunctorMaterial
    prop_names = 'forch'
    prop_values = '${F_hx} ${F_hx} ${F_hx}'
    block = 'hx'
  []
  [hx_removal_mat]
    type = ParsedFunctorMaterial
    property_name = hx_removal
    functor_names = 'T_fluid'
    functor_symbols = 'T'
    expression = '${h_hx} * (T - ${T_HX})'
    block = 'hx'
  []
[]

[LinearFVKernels]
  # --- momentum ---
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'x'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
    block = ${fluid}
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [u_rz]
    type = LinearFVRZViscousSource
    variable = vel_x
    mu = 'mu_eff'
    momentum_component = 'x'
    block = ${fluid}
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
    block = ${fluid}
  []
  [u_friction]
    type = LinearFVMomentumFriction
    variable = vel_x
    mu = 'speed_lag'
    Darcy_name = 'forch'
    momentum_component = 'x'
    block = 'hx'
  []
  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'y'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
    block = ${fluid}
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
    block = ${fluid}
  []
  [v_friction]
    type = LinearFVMomentumFriction
    variable = vel_y
    mu = 'speed_lag'
    Darcy_name = 'forch'
    momentum_component = 'y'
    block = 'hx'
  []
  [v_buoyancy]
    type = LinearFVMomentumBoussinesq
    variable = vel_y
    T_fluid = T_fluid
    gravity = '0 -9.81 0'
    rho = ${rho}
    alpha_name = 'alpha_b'
    ref_temperature = ${T_HX}
    momentum_component = 'y'
    block = ${fluid}
  []
  [pump]
    type = LinearFVSource
    variable = vel_y
    source_density = pump_profile
    block = 'pump'
  []
  [p_diffusion]
    type = LinearFVPressureCorrectionDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [HbyA_divergence]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = true
    block = ${fluid}
  []

  # --- k-epsilon ---
  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
    block = ${fluid}
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [TKE_source_sink]
    type = LinearFVTKESourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    block = ${fluid}
  []
  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = ${walls}
    block = ${fluid}
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    walls = ${walls}
    block = ${fluid}
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    walls = ${walls}
    block = ${fluid}
  []
  [TKED_source_sink]
    type = LinearFVTKEDSourceSink
    variable = TKED
    u = vel_x
    v = vel_y
    tke = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    block = ${fluid}
  []

  # --- energy ---
  [energy_advection]
    type = LinearFVEnergyAdvection
    variable = T_fluid
    advected_quantity = temperature
    cp = ${cp}
    block = ${fluid}
  []
  [energy_conduction]
    type = LinearFVDiffusion
    variable = T_fluid
    diffusion_coeff = 'k_eff'
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [fission_heat]
    type = LinearFVSource
    variable = T_fluid
    source_density = power_density
    block = ${fluid}
  []
  [hx_sink]
    type = LinearFVReaction
    variable = T_fluid
    coeff = ${h_hx}
    block = 'hx'
  []
  [hx_secondary]
    type = LinearFVSource
    variable = T_fluid
    source_density = ${fparse h_hx * T_HX}
    block = 'hx'
  []
[]

[LinearFVBCs]
  [walls_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_x
    functor = 0.0
  []
  [walls_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_y
    functor = 0.0
  []
  [axis_u]
    type = LinearFVVelocitySymmetryBC
    boundary = fluid_symmetry
    variable = vel_x
    u = vel_x
    v = vel_y
    momentum_component = x
  []
  [axis_v]
    type = LinearFVVelocitySymmetryBC
    boundary = fluid_symmetry
    variable = vel_y
    u = vel_x
    v = vel_y
    momentum_component = y
  []
  [axis_p]
    type = LinearFVPressureSymmetryBC
    boundary = fluid_symmetry
    variable = pressure
    HbyA_flux = 'HbyA'
  []
  [walls_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = ${walls}
    variable = mu_t
    u = vel_x
    v = vel_y
    rho = ${rho}
    mu = ${mu}
    tke = TKE
    wall_treatment = ${wall_treatment}
  []
[]

[Postprocessors]
  [max_v]
    type = ElementExtremeValue
    variable = vel_x
    value_type = max
    block = ${fluid}
  []
  [max_vy]
    type = ElementExtremeValue
    variable = vel_y
    value_type = max
    block = ${fluid}
  []
  [min_vy]
    type = ElementExtremeValue
    variable = vel_y
    value_type = min
    block = ${fluid}
  []
  [T_max]
    type = ElementExtremeValue
    variable = T_fluid
    block = ${fluid}
  []
  [T_min]
    type = ElementExtremeValue
    variable = T_fluid
    value_type = min
    block = ${fluid}
  []
  [T_avg]
    type = ElementAverageValue
    variable = T_fluid
    block = ${fluid}
  []
  [p_max]
    type = ElementExtremeValue
    variable = pressure
    block = ${fluid}
  []
  [p_min]
    type = ElementExtremeValue
    variable = pressure
    value_type = min
    block = ${fluid}
  []
  [mu_t_max]
    type = ElementExtremeValue
    variable = mu_t
    block = ${fluid}
  []
  [flow_hx_top]
    type = VolumetricFlowRate
    boundary = 'hx_top'
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = 1
    rhie_chow_user_object = rc
    execute_on = 'TIMESTEP_END'
  []
  [mdot]
    type = ParsedPostprocessor
    expression = '-${rho} * flow_hx_top'
    pp_names = 'flow_hx_top'
  []
  [Tflow_hx_top]
    type = VolumetricFlowRate
    boundary = 'hx_top'
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = T_fluid
    rhie_chow_user_object = rc
    execute_on = 'TIMESTEP_END'
    outputs = none
  []
  [flow_hx_bot]
    type = VolumetricFlowRate
    boundary = 'hx_bot'
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = 1
    rhie_chow_user_object = rc
    execute_on = 'TIMESTEP_END'
    outputs = none
  []
  [Tflow_hx_bot]
    type = VolumetricFlowRate
    boundary = 'hx_bot'
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = T_fluid
    rhie_chow_user_object = rc
    execute_on = 'TIMESTEP_END'
    outputs = none
  []
  [T_hx_in]
    type = ParsedPostprocessor
    expression = 'Tflow_hx_top / flow_hx_top'
    pp_names = 'Tflow_hx_top flow_hx_top'
  []
  [T_hx_out]
    type = ParsedPostprocessor
    expression = 'Tflow_hx_bot / flow_hx_bot'
    pp_names = 'Tflow_hx_bot flow_hx_bot'
  []
  [fission_power]
    type = ElementIntegralVariablePostprocessor
    variable = power_density
    block = ${fluid}
  []
  [hx_removal_power]
    type = ElementIntegralFunctorPostprocessor
    functor = hx_removal
    block = 'hx'
  []
  [energy_balance_error]
    type = ParsedPostprocessor
    pp_names = 'fission_power hx_removal_power'
    expression = 'abs(fission_power - hx_removal_power) / fission_power'
  []
  [salt_volume]
    type = VolumePostprocessor
    block = ${fluid}
  []
[]

[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  energy_system = 'energy_system'
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.3 0.3'
  turbulence_field_relaxation = '0.3 0.3'
  energy_equation_relaxation = 0.9
  num_iterations = 1000
  pressure_absolute_tolerance = 1e-6
  momentum_absolute_tolerance = 1e-6
  turbulence_absolute_tolerance = '1e-6 1e-6'
  energy_absolute_tolerance = 1e-6
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'
  turbulence_l_max_its = 30
  pin_pressure = true
  pressure_pin_value = 1e5
  pressure_pin_point = '0.5 0.0 0.0'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'
  print_fields = false
  continue_on_max_its = true
[]

[VectorPostprocessors]
  # Core mid-plane (y = 0) from the axis to the shield wall
  [midplane]
    type = LineValueSampler
    warn_discontinuous_face_values = false  # FV cell values; a point may sit on a face
    variable = 'T_fluid vel_x vel_y mu_t'
    start_point = '0.01 0 0'
    end_point = '1.05 0 0'
    num_points = 53
    sort_by = x
    execute_on = TIMESTEP_END
  []
  # Core axis (x = 0.05) from the bottom to the top of the cavity
  [axis]
    type = LineValueSampler
    warn_discontinuous_face_values = false  # FV cell values; a point may sit on a face
    variable = 'T_fluid vel_x vel_y mu_t'
    start_point = '0.05 -0.79 0'
    end_point = '0.05 0.79 0'
    num_points = 80
    sort_by = y
    execute_on = TIMESTEP_END
  []
  # External loop (x = 2.13), top to bottom: hot leg, pump, HX, cold leg
  [loop]
    type = LineValueSampler
    warn_discontinuous_face_values = false  # FV cell values; a point may sit on a face
    variable = 'T_fluid vel_x vel_y mu_t'
    start_point = '2.13 1.0 0'
    end_point = '2.13 -1.0 0'
    num_points = 101
    sort_by = y
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  csv = true
  [restart]
    # converged flow / k-epsilon / energy fields (and the Griffin power) read by msfr_chemistry.i
    type = Exodus
    file_base = msfr_flow_restart
    execute_on = FINAL
  []
[]
