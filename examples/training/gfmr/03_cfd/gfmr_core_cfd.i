################################################################################
# GFMR training -- Stage 3: 3-D conjugate-heat-transfer CFD of the core (60-degree sector)
#
# Linear finite volume, SIMPLE, on the CFD sector mesh (../01_mesh/gfmr_sector_3d_cfd.e: every
# fuel channel, the control-rod annuli, the graphite moderator / reflector and the guide tubes).
#
#   salt      : incompressible turbulent flow, standard k-epsilon with equilibrium wall functions,
#               energy with turbulent conduction (Pr_t = 0.85), FUNaK properties at 650 C (SPEC)
#   graphite  : conduction, k = 25 W/m/K (SPEC); the guide tube uses the same k (simplification:
#               Hastelloy N ~21 W/m/K; the tube carries no net heat, see README)
#   walls     : LinearFVTurbulentWallCHTBC on both sides (Jayatilleke thermal wall function for
#               Pr ~ 17; verified in test/tests/kEpsilon/salt_pipe_heat_transfer)
#   power     : the 8-group diffusion power map (../02_neutronics, calibrated against OpenMC),
#               150 MW total: 92.1 % deposited in the salt, 7.9 % in the graphite (OpenMC
#               coupled neutron-photon heating), the graphite share distributed like the fast
#               flux
#   inlets    : the SPEC zone velocities 1.70 / 1.03 / 1.74 m/s at T_in = 613.8 C (the 650 C
#               SPEC fuel temperature taken as the core-average salt temperature, 2300 kg/s)
#
# The channels are hydraulically separate (no cross flow in the graphite), so each zone keeps
# its specified velocity; the computed pressure drops are compared with the SPEC values.
#
# Run: mpirun -n 24 open_pronghorn-opt -i gfmr_core_cfd.i
################################################################################

# ---- FUNaK at 650 C (SPEC correlations), k ASSUMED (not specified)
rho = 4073.5
mu = 0.01494
cp = 900.0
k_salt = 0.8
Pr_t = 0.85
k_graphite = 25.0

T_in = ${fparse 613.8 + 273.15}
u_inner = 1.70
u_outer = 1.03
u_cr = 1.74

# ---- heating split (OpenMC coupled n-gamma): fraction deposited in the graphite
f_graphite = 0.0787
power_sector = ${fparse 150e6 / 6}

# ---- k-epsilon
C_mu = 0.09
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
intensity = 0.05
k_init = '${fparse 1.5 * (intensity * u_inner)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / (0.07 * 0.04)}'
wall_treatment = 'eq_newton'

salt = 'salt_inner_tri salt_inner salt_outer_tri salt_outer salt_cr'
solid = 'graphite graphite_plug_tri graphite_plug reflector guide_tube'
walls = 'wall_inner wall_outer wall_cr wall_guide'

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = ../01_mesh/gfmr_sector_3d_cfd.e
  []
  # the absorber and helium are not part of the thermal model (stagnant, inside the guide tube)
  [drop_rod]
    type = BlockDeletionGenerator
    input = fmg
    block = 'absorber_tri absorber helium helium_tri'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system w_system pressure_system TKE_system TKED_system energy_system solid_energy_system'
  previous_nl_solution_required = true
  material_coverage_check = false
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method_name = 'upwind'
[]

[FVInterpolationMethods]
  [upwind]
    type = FVAdvectedUpwind
  []
[]

[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    w = vel_z
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
    block = ${salt}
  []
  # power map of the diffusion solution (salt power density and graphite fast flux)
  [power_map]
    type = SolutionUserObject
    mesh = ../02_neutronics/sector_power.e
    system_variables = 'power_density fast_flux'
    timestep = LATEST
  []
[]

[Functions]
  [q_salt_fn]
    type = SolutionFunction
    solution = power_map
    from_variable = power_density
  []
  [fast_fn]
    type = SolutionFunction
    solution = power_map
    from_variable = fast_flux
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_condition = 0
    block = ${salt}
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_condition = 0
    block = ${salt}
  []
  [vel_z]
    type = MooseLinearVariableFVReal
    solver_sys = w_system
    initial_condition = ${u_outer}
    block = ${salt}
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 0
    block = ${salt}
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    initial_condition = ${k_init}
    block = ${salt}
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    initial_condition = ${eps_init}
    block = ${salt}
  []
  [T_fluid]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
    initial_condition = ${T_in}
    block = ${salt}
  []
  [T_solid]
    type = MooseLinearVariableFVReal
    solver_sys = solid_energy_system
    initial_condition = ${fparse T_in + 60}
    block = ${solid}
  []
[]

[AuxVariables]
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse rho * C_mu * k_init^2 / eps_init}'
    block = ${salt}
  []
  [wall_distance]
    type = MooseVariableFVReal
    block = ${salt}
  []
  [q_salt]
    type = MooseVariableFVReal
    block = ${salt}
  []
  [fast_flux]
    type = MooseVariableFVReal
    block = 'graphite graphite_plug_tri graphite_plug reflector'
  []
[]

[AuxKernels]
  [wall_distance]
    type = WallDistanceAux
    variable = wall_distance
    walls = ${walls}
    execute_on = INITIAL
    block = ${salt}
  []
  [mu_t]
    type = kEpsilonViscosityAux
    variable = mu_t
    C_mu = ${C_mu}
    tke = TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    w = vel_z
    bulk_wall_treatment = false
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
    block = ${salt}
  []
  [q_salt]
    type = FunctionAux
    variable = q_salt
    function = q_salt_fn
    execute_on = INITIAL
    block = ${salt}
  []
  [fast_flux]
    type = FunctionAux
    variable = fast_flux
    function = fast_fn
    execute_on = INITIAL
    block = 'graphite graphite_plug_tri graphite_plug reflector'
  []
[]

[Postprocessors]
  # ---- power normalization: salt 92.1 % of the sector power, graphite 7.9 % shaped like the
  # fast flux (the gamma and neutron heating of the moderator follow the flux closely)
  [q_salt_integral]
    type = ElementIntegralVariablePostprocessor
    variable = q_salt
    block = ${salt}
    execute_on = INITIAL
  []
  [fast_integral]
    type = ElementIntegralVariablePostprocessor
    variable = fast_flux
    block = 'graphite graphite_plug_tri graphite_plug reflector'
    execute_on = INITIAL
  []
  [salt_scale]
    type = ParsedPostprocessor
    expression = '${power_sector} * (1 - ${f_graphite}) / q_salt_integral'
    pp_names = 'q_salt_integral'
    execute_on = INITIAL
  []
  [graphite_scale]
    type = ParsedPostprocessor
    expression = '${power_sector} * ${f_graphite} / fast_integral'
    pp_names = 'fast_integral'
    execute_on = INITIAL
  []
[]

[FunctorMaterials]
  [k_eff]
    type = ParsedFunctorMaterial
    property_name = k_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${k_salt} + mut * ${cp} / ${Pr_t}'
    block = ${salt}
  []
  [q_salt_scaled]
    type = ParsedFunctorMaterial
    property_name = q_salt_W
    functor_names = 'q_salt salt_scale'
    functor_symbols = 'q s'
    expression = 'q * s'
    block = ${salt}
  []
  [q_graphite]
    type = ParsedFunctorMaterial
    property_name = q_graphite_W
    functor_names = 'fast_flux graphite_scale'
    functor_symbols = 'f s'
    expression = 'f * s'
    block = 'graphite graphite_plug_tri graphite_plug reflector'
  []
  [rhocpT]
    type = ParsedFunctorMaterial
    property_name = rhocpT
    functor_names = 'T_fluid'
    expression = '${rho} * ${cp} * T_fluid'
    block = ${salt}
  []
[]

[LinearFVKernels]
  # ---------------------------------------------------------------- momentum
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = mu_t
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = 'x'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
    block = ${salt}
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = ${salt}
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
    block = ${salt}
  []
  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    mu = mu_t
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = 'y'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
    block = ${salt}
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = ${salt}
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
    block = ${salt}
  []
  [w_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_z
    mu = mu_t
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = 'z'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
    block = ${salt}
  []
  [w_diffusion]
    type = LinearFVDiffusion
    variable = vel_z
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = ${salt}
  []
  [w_pressure]
    type = LinearFVMomentumPressure
    variable = vel_z
    pressure = pressure
    momentum_component = 'z'
    block = ${salt}
  []
  [p_diffusion]
    type = LinearFVPressureCorrectionDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
    block = ${salt}
  []
  [HbyA_divergence]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = true
    block = ${salt}
  []

  # ---------------------------------------------------------------- k-epsilon
  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
    block = ${salt}
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = ${salt}
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = mu_t
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = false
    block = ${salt}
  []
  [TKE_source_sink]
    type = LinearFVTKESourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    w = vel_z
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = mu_t
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    block = ${salt}
  []
  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = ${walls}
    block = ${salt}
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    walls = ${walls}
    block = ${salt}
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = mu_t
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    walls = ${walls}
    block = ${salt}
  []
  [TKED_source_sink]
    type = LinearFVTKEDSourceSink
    variable = TKED
    u = vel_x
    v = vel_y
    w = vel_z
    tke = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = mu_t
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    block = ${salt}
  []

  # ---------------------------------------------------------------- energy
  [energy_advection]
    type = LinearFVEnergyAdvection
    variable = T_fluid
    advected_quantity = temperature
    cp = ${cp}
    block = ${salt}
  []
  [energy_conduction]
    type = LinearFVDiffusion
    variable = T_fluid
    diffusion_coeff = k_eff
    use_nonorthogonal_correction = false
    block = ${salt}
  []
  [fission_heat]
    type = LinearFVSource
    variable = T_fluid
    source_density = q_salt_W
    block = ${salt}
  []
  [solid_conduction]
    type = LinearFVDiffusion
    variable = T_solid
    diffusion_coeff = ${k_graphite}
    use_nonorthogonal_correction = false
    block = ${solid}
  []
  [graphite_heat]
    type = LinearFVSource
    variable = T_solid
    source_density = q_graphite_W
    block = 'graphite graphite_plug_tri graphite_plug reflector'
  []
[]

[LinearFVBCs]
  # ---------------------------------------------------------------- inlets (SPEC zone velocities)
  [inlet_vel_x]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_inner inlet_outer inlet_cr'
    variable = vel_x
    functor = 0
  []
  [inlet_vel_y]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_inner inlet_outer inlet_cr'
    variable = vel_y
    functor = 0
  []
  [inlet_w_inner]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_inner'
    variable = vel_z
    functor = ${u_inner}
  []
  [inlet_w_outer]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_outer'
    variable = vel_z
    functor = ${u_outer}
  []
  [inlet_w_cr]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_cr'
    variable = vel_z
    functor = ${u_cr}
  []
  [inlet_TKE]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_inner inlet_outer inlet_cr'
    variable = TKE
    functor = ${k_init}
  []
  [inlet_TKED]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_inner inlet_outer inlet_cr'
    variable = TKED
    functor = ${eps_init}
  []
  [inlet_T]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_inner inlet_outer inlet_cr'
    variable = T_fluid
    functor = ${T_in}
  []
  # ---------------------------------------------------------------- outlet
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'outlet'
    variable = pressure
    functor = 0
  []
  [outlet_vel_x]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = vel_x
    use_two_term_expansion = false
  []
  [outlet_vel_y]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = vel_y
    use_two_term_expansion = false
  []
  [outlet_vel_z]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = vel_z
    use_two_term_expansion = false
  []
  [outlet_TKE]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = TKE
    use_two_term_expansion = false
  []
  [outlet_TKED]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = TKED
    use_two_term_expansion = false
  []
  [outlet_T_fluid]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = T_fluid
    use_two_term_expansion = false
  []
  # ---------------------------------------------------------------- walls
  [wall_vel_x]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls}'
    variable = vel_x
    functor = 0
  []
  [wall_vel_y]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls}'
    variable = vel_y
    functor = 0
  []
  [wall_vel_z]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls}'
    variable = vel_z
    functor = 0
  []
  [wall_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = '${walls}'
    variable = mu_t
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    tke = TKE
    wall_treatment = ${wall_treatment}
  []
  # conjugate heat transfer through the Jayatilleke wall function, one pair per interface
  [fluid_wall_inner]
    type = LinearFVTurbulentWallCHTBC
    boundary = 'wall_inner'
    variable = T_fluid
    side = fluid
    thermal_conductivity = k_eff
    other_temperature = interface_temperature_solid_wall_inner
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    cp = ${cp}
    k_fluid = ${k_salt}
    Pr_t = ${Pr_t}
    wall_treatment = ${wall_treatment}
  []
  [solid_wall_inner]
    type = LinearFVTurbulentWallCHTBC
    boundary = 'wall_inner'
    variable = T_solid
    side = solid
    thermal_conductivity = ${k_graphite}
    other_temperature = interface_temperature_fluid_wall_inner
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    cp = ${cp}
    k_fluid = ${k_salt}
    Pr_t = ${Pr_t}
    wall_treatment = ${wall_treatment}
  []
  [fluid_wall_outer]
    type = LinearFVTurbulentWallCHTBC
    boundary = 'wall_outer'
    variable = T_fluid
    side = fluid
    thermal_conductivity = k_eff
    other_temperature = interface_temperature_solid_wall_outer
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    cp = ${cp}
    k_fluid = ${k_salt}
    Pr_t = ${Pr_t}
    wall_treatment = ${wall_treatment}
  []
  [solid_wall_outer]
    type = LinearFVTurbulentWallCHTBC
    boundary = 'wall_outer'
    variable = T_solid
    side = solid
    thermal_conductivity = ${k_graphite}
    other_temperature = interface_temperature_fluid_wall_outer
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    cp = ${cp}
    k_fluid = ${k_salt}
    Pr_t = ${Pr_t}
    wall_treatment = ${wall_treatment}
  []
  [fluid_wall_cr]
    type = LinearFVTurbulentWallCHTBC
    boundary = 'wall_cr'
    variable = T_fluid
    side = fluid
    thermal_conductivity = k_eff
    other_temperature = interface_temperature_solid_wall_cr
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    cp = ${cp}
    k_fluid = ${k_salt}
    Pr_t = ${Pr_t}
    wall_treatment = ${wall_treatment}
  []
  [solid_wall_cr]
    type = LinearFVTurbulentWallCHTBC
    boundary = 'wall_cr'
    variable = T_solid
    side = solid
    thermal_conductivity = ${k_graphite}
    other_temperature = interface_temperature_fluid_wall_cr
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    cp = ${cp}
    k_fluid = ${k_salt}
    Pr_t = ${Pr_t}
    wall_treatment = ${wall_treatment}
  []
  [fluid_wall_guide]
    type = LinearFVTurbulentWallCHTBC
    boundary = 'wall_guide'
    variable = T_fluid
    side = fluid
    thermal_conductivity = k_eff
    other_temperature = interface_temperature_solid_wall_guide
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    cp = ${cp}
    k_fluid = ${k_salt}
    Pr_t = ${Pr_t}
    wall_treatment = ${wall_treatment}
  []
  [solid_wall_guide]
    type = LinearFVTurbulentWallCHTBC
    boundary = 'wall_guide'
    variable = T_solid
    side = solid
    thermal_conductivity = ${k_graphite}
    other_temperature = interface_temperature_fluid_wall_guide
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    cp = ${cp}
    k_fluid = ${k_salt}
    Pr_t = ${Pr_t}
    wall_treatment = ${wall_treatment}
  []
  # ---------------------------------------------------------------- outer surfaces
  [solid_adiabatic]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    boundary = 'top bottom reflector_outer symmetry'
    variable = T_solid
    functor = 0
    diffusion_coeff = ${k_graphite}
  []
  # symmetry planes cut no salt channel (the mirror lines pass between and through the
  # channel centres; a cut channel is a half channel with a slip / adiabatic plane)
  [sym_vel_x]
    type = LinearFVVelocitySymmetryBC
    boundary = 'symmetry'
    variable = vel_x
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = x
  []
  [sym_vel_y]
    type = LinearFVVelocitySymmetryBC
    boundary = 'symmetry'
    variable = vel_y
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = y
  []
  [sym_vel_z]
    type = LinearFVVelocitySymmetryBC
    boundary = 'symmetry'
    variable = vel_z
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = z
  []
  [sym_p]
    type = LinearFVPressureSymmetryBC
    boundary = 'symmetry'
    variable = pressure
    HbyA_flux = HbyA
  []
  [sym_TKE]
    type = LinearFVAdvectionDiffusionScalarSymmetryBC
    boundary = 'symmetry'
    variable = TKE
  []
  [sym_TKED]
    type = LinearFVAdvectionDiffusionScalarSymmetryBC
    boundary = 'symmetry'
    variable = TKED
  []
  [sym_T_fluid]
    type = LinearFVAdvectionDiffusionScalarSymmetryBC
    boundary = 'symmetry'
    variable = T_fluid
  []
[]

[Postprocessors]
  [heat_salt]
    type = ElementIntegralFunctorPostprocessor
    functor = q_salt_W
    block = ${salt}
  []
  [heat_graphite]
    type = ElementIntegralFunctorPostprocessor
    functor = q_graphite_W
    block = 'graphite graphite_plug_tri graphite_plug reflector'
  []
  [enthalpy_in]
    type = VolumetricFlowRate
    boundary = 'inlet_inner inlet_outer inlet_cr'
    vel_x = vel_x
    vel_y = vel_y
    vel_z = vel_z
    advected_quantity = rhocpT
    subtract_mesh_velocity = false
  []
  [enthalpy_out]
    type = VolumetricFlowRate
    boundary = outlet
    vel_x = vel_x
    vel_y = vel_y
    vel_z = vel_z
    advected_quantity = rhocpT
    subtract_mesh_velocity = false
  []
  [mdot_sector]
    type = VolumetricFlowRate
    boundary = outlet
    vel_x = vel_x
    vel_y = vel_y
    vel_z = vel_z
    advected_quantity = ${rho}
    subtract_mesh_velocity = false
  []
  [mdot_core]
    type = ParsedPostprocessor
    expression = '6 * mdot_sector'
    pp_names = 'mdot_sector'
  []
  [energy_balance]
    type = ParsedPostprocessor
    expression = '(enthalpy_out + enthalpy_in - heat_salt - heat_graphite) / (heat_salt + heat_graphite)'
    pp_names = 'enthalpy_in enthalpy_out heat_salt heat_graphite'
  []
  [T_outlet_mixed]
    type = ParsedPostprocessor
    expression = 'enthalpy_out / (${cp} * mdot_sector) - 273.15'
    pp_names = 'enthalpy_out mdot_sector'
  []
  [T_salt_max]
    type = ElementExtremeValue
    variable = T_fluid
    block = ${salt}
  []
  [T_graphite_max]
    type = ElementExtremeValue
    variable = T_solid
    block = 'graphite graphite_plug_tri graphite_plug reflector'
  []
  [T_graphite_avg]
    type = ElementAverageValue
    variable = T_solid
    block = 'graphite graphite_plug_tri graphite_plug'
  []
  [T_guide_max]
    type = ElementExtremeValue
    variable = T_solid
    block = guide_tube
  []
  [p_inlet_inner]
    type = SideAverageValue
    variable = pressure
    boundary = inlet_inner
  []
  [p_inlet_outer]
    type = SideAverageValue
    variable = pressure
    boundary = inlet_outer
  []
  [p_inlet_cr]
    type = SideAverageValue
    variable = pressure
    boundary = inlet_cr
  []
[]

[Executioner]
  type = SIMPLE
  momentum_systems = 'u_system v_system w_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  energy_system = 'energy_system'
  solid_energy_system = 'solid_energy_system'
  momentum_l_abs_tol = 1e-10
  pressure_l_abs_tol = 1e-10
  turbulence_l_abs_tol = 1e-10
  energy_l_abs_tol = 1e-10
  solid_energy_l_abs_tol = 1e-10
  momentum_l_tol = 1e-4
  pressure_l_tol = 1e-4
  turbulence_l_tol = 1e-4
  energy_l_tol = 1e-6
  solid_energy_l_tol = 1e-6
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.5 0.5'
  turbulence_field_relaxation = '0.5 0.5'
  energy_equation_relaxation = 0.95
  num_iterations = 600
  pressure_absolute_tolerance = 1e-7
  momentum_absolute_tolerance = 1e-7
  turbulence_absolute_tolerance = '1e-7 1e-7'
  energy_absolute_tolerance = 1e-7
  solid_energy_absolute_tolerance = 1e-7
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'
  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'
  solid_energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  solid_energy_petsc_options_value = 'hypre boomeramg'
  cht_interfaces = 'wall_inner wall_outer wall_cr wall_guide'
  max_cht_fpi = 1
  cht_solid_flux_relaxation = '1 1 1 1'
  cht_fluid_flux_relaxation = '1 1 1 1'
  cht_solid_temperature_relaxation = '1 1 1 1'
  cht_fluid_temperature_relaxation = '1 1 1 1'
  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = FINAL
  [console]
    type = Console
    execute_postprocessors_on = 'INITIAL FINAL'
  []
[]
