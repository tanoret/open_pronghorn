################################################################################
# Turbulent molten-salt pipe flow with conjugate heat transfer to a graphite wall
# (verification of LinearFVTurbulentWallCHTBC)
#
# 2-D R-Z, a 4 cm FUNaK channel (GFMR inner zone: 1.7 m/s, Re ~ 1.9e4, Pr ~ 17) inside a
# graphite sleeve heated uniformly. The fluid-side wall condition is the Jayatilleke thermal
# wall function, applied on both sides of the interface (LinearFVTurbulentWallCHTBC).
#
# Reported: the local Nusselt number Nu = q'' D / (k (T_w - T_b)) at z = 1.61 m (L/D = 40), compared in the
# test with the Gnielinski correlation (Nu = 194 at Re 1.85e4, Pr 16.8): Nu_CFD = 180, ratio 0.93.
################################################################################

D = 0.04
R = ${fparse D / 2}
R_out = 0.035
L = 2.0
u_in = 1.70
T_in = 900.0

rho = 4073.5
mu = 0.01494
cp = 900.0
k = 0.8
k_s = 25.0
Pr_t = 0.85
# uniform heating of the graphite (W/m^3) giving ~ 120 kW/m^2 at the wall
q_s = ${fparse 1.2e5 * 2 * R / (R_out^2 - R^2)}

C_mu = 0.09
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
intensity = 0.05
k_init = '${fparse 1.5*(intensity * u_in)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / (0.07 * D)}'
wall_treatment = 'eq_newton'

nr = 8
nr_s = 4
nz = 100

[Mesh]
  [fluid]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${R}
    ymin = 0
    ymax = ${L}
    nx = ${nr}
    ny = ${nz}
    # uniform: first cell centre at y = 1.25 mm, y+ ~ 33 (log layer, as wall functions need)
    subdomain_ids = 0
    subdomain_name = fluid
    boundary_name_prefix = fluid
  []
  [solid]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = ${R}
    xmax = ${R_out}
    ymin = 0
    ymax = ${L}
    nx = ${nr_s}
    ny = ${nz}
    subdomain_ids = 1
    subdomain_name = solid
    boundary_name_prefix = solid
    boundary_id_offset = 10
  []
  [stitch]
    type = StitchMeshGenerator
    inputs = 'fluid solid'
    stitch_boundaries_pairs = 'fluid_right solid_left'
    clear_stitched_boundary_ids = true
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = stitch
    primary_block = fluid
    paired_block = solid
    new_boundary = interface
  []
  coord_type = RZ
  rz_coord_axis = Y
  allow_renumbering = false
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system TKE_system TKED_system energy_system solid_energy_system'
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
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
    block = fluid
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_condition = 0
    block = fluid
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_condition = ${u_in}
    block = fluid
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 0
    block = fluid
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    initial_condition = ${k_init}
    block = fluid
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    initial_condition = ${eps_init}
    block = fluid
  []
  [T_fluid]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
    initial_condition = ${T_in}
    block = fluid
  []
  [T_solid]
    type = MooseLinearVariableFVReal
    solver_sys = solid_energy_system
    initial_condition = ${T_in}
    block = solid
  []
[]

[AuxVariables]
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse rho * C_mu * k_init^2 / eps_init}'
    block = fluid
  []
  [wall_distance]
    type = MooseVariableFVReal
    block = fluid
  []
[]

[AuxKernels]
  [wall_distance]
    type = WallDistanceAux
    variable = wall_distance
    walls = 'interface'
    execute_on = INITIAL
    block = fluid
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
    bulk_wall_treatment = false
    walls = 'interface'
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
    block = fluid
  []
[]

[FunctorMaterials]
  [k_eff]
    type = ParsedFunctorMaterial
    property_name = k_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${k} + mut * ${cp} / ${Pr_t}'
    block = fluid
  []
  [rhocpT]
    type = ParsedFunctorMaterial
    property_name = rhocpT
    functor_names = 'T_fluid'
    expression = '${rho} * ${cp} * T_fluid'
    block = fluid
  []
[]

[LinearFVKernels]
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = mu_t
    u = vel_x
    v = vel_y
    momentum_component = 'x'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
    block = fluid
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = fluid
  []
  [u_rz]
    type = LinearFVRZViscousSource
    variable = vel_x
    mu = ${fparse mu}
    momentum_component = 'x'
    block = fluid
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
    block = fluid
  []
  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    mu = mu_t
    u = vel_x
    v = vel_y
    momentum_component = 'y'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
    block = fluid
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = fluid
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
    block = fluid
  []
  [p_diffusion]
    type = LinearFVPressureCorrectionDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
    block = fluid
  []
  [HbyA_divergence]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = true
    block = fluid
  []

  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
    block = fluid
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = fluid
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = mu_t
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = false
    block = fluid
  []
  [TKE_source_sink]
    type = LinearFVTKESourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = mu_t
    walls = 'interface'
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    block = fluid
  []
  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = 'interface'
    block = fluid
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    walls = 'interface'
    block = fluid
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = mu_t
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    walls = 'interface'
    block = fluid
  []
  [TKED_source_sink]
    type = LinearFVTKEDSourceSink
    variable = TKED
    u = vel_x
    v = vel_y
    tke = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = mu_t
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = 'interface'
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    block = fluid
  []

  [energy_advection]
    type = LinearFVEnergyAdvection
    variable = T_fluid
    advected_quantity = temperature
    cp = ${cp}
    block = fluid
  []
  [energy_conduction]
    type = LinearFVDiffusion
    variable = T_fluid
    diffusion_coeff = k_eff
    use_nonorthogonal_correction = false
    block = fluid
  []
  [solid_conduction]
    type = LinearFVDiffusion
    variable = T_solid
    diffusion_coeff = ${k_s}
    use_nonorthogonal_correction = false
    block = solid
  []
  [solid_heating]
    type = LinearFVSource
    variable = T_solid
    source_density = ${q_s}
    block = solid
  []
[]

[LinearFVBCs]
  [inlet_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = fluid_bottom
    variable = vel_x
    functor = 0
  []
  [inlet_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = fluid_bottom
    variable = vel_y
    functor = ${u_in}
  []
  [inlet_TKE]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = fluid_bottom
    variable = TKE
    functor = ${k_init}
  []
  [inlet_TKED]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = fluid_bottom
    variable = TKED
    functor = ${eps_init}
  []
  [inlet_T]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = fluid_bottom
    variable = T_fluid
    functor = ${T_in}
  []
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = fluid_top
    variable = pressure
    functor = 0
  []
  [outlet_vel_x]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = fluid_top
    variable = vel_x
    use_two_term_expansion = false
  []
  [outlet_vel_y]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = fluid_top
    variable = vel_y
    use_two_term_expansion = false
  []
  [outlet_TKE]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = fluid_top
    variable = TKE
    use_two_term_expansion = false
  []
  [outlet_TKED]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = fluid_top
    variable = TKED
    use_two_term_expansion = false
  []
  [outlet_T_fluid]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = fluid_top
    variable = T_fluid
    use_two_term_expansion = false
  []
  [axis_u]
    type = LinearFVVelocitySymmetryBC
    boundary = fluid_left
    variable = vel_x
    u = vel_x
    v = vel_y
    momentum_component = x
  []
  [axis_v]
    type = LinearFVVelocitySymmetryBC
    boundary = fluid_left
    variable = vel_y
    u = vel_x
    v = vel_y
    momentum_component = y
  []
  [axis_p]
    type = LinearFVPressureSymmetryBC
    boundary = fluid_left
    variable = pressure
    HbyA_flux = 'HbyA'
  []
  [axis_TKE]
    type = LinearFVAdvectionDiffusionScalarSymmetryBC
    boundary = fluid_left
    variable = TKE
  []
  [axis_TKED]
    type = LinearFVAdvectionDiffusionScalarSymmetryBC
    boundary = fluid_left
    variable = TKED
  []
  [axis_T_fluid]
    type = LinearFVAdvectionDiffusionScalarSymmetryBC
    boundary = fluid_left
    variable = T_fluid
  []
  [wall_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = interface
    variable = vel_x
    functor = 0
  []
  [wall_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = interface
    variable = vel_y
    functor = 0
  []
  [wall_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = interface
    variable = mu_t
    u = vel_x
    v = vel_y
    rho = ${rho}
    mu = ${mu}
    tke = TKE
    wall_treatment = ${wall_treatment}
  []
  # Wall-function CHT on both sides: the fluid sees the solid surface temperature, the solid
  # sees the near-wall fluid temperature through the wall function in series with its half cell
  [fluid_wall_T]
    type = LinearFVTurbulentWallCHTBC
    side = fluid
    boundary = interface
    variable = T_fluid
    u = vel_x
    v = vel_y
    rho = ${rho}
    mu = ${mu}
    cp = ${cp}
    k_fluid = ${k}
    Pr_t = ${Pr_t}
    wall_treatment = ${wall_treatment}
    thermal_conductivity = k_eff
    other_temperature = interface_temperature_solid_interface
  []
  [solid_wall_T]
    type = LinearFVTurbulentWallCHTBC
    side = solid
    boundary = interface
    variable = T_solid
    u = vel_x
    v = vel_y
    rho = ${rho}
    mu = ${mu}
    cp = ${cp}
    k_fluid = ${k}
    Pr_t = ${Pr_t}
    wall_treatment = ${wall_treatment}
    thermal_conductivity = ${k_s}
    other_temperature = interface_temperature_fluid_interface
  []
  [solid_insulated]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    boundary = 'solid_right solid_top solid_bottom'
    diffusion_coeff = ${k_s}
    variable = T_solid
    functor = 0
  []
[]

[Postprocessors]
  [T_bulk_probe]
    # mixed-mean temperature at L/D = 40 (line integral of rho cp T v r dr / of rho cp v r dr)
    type = Receiver
  []
  [T_wall]
    type = PointValue
    variable = T_solid
    point = '${fparse R + 1e-4} 1.61 0'
  []
  [T_center]
    type = PointValue
    variable = T_fluid
    point = '1e-4 1.61 0'
  []
  [heat_in]
    type = ElementIntegralFunctorPostprocessor
    functor = ${q_s}
    block = solid
  []
  [enthalpy_in]
    type = VolumetricFlowRate
    boundary = fluid_bottom
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = rhocpT
    subtract_mesh_velocity = false
  []
  [enthalpy_out]
    type = VolumetricFlowRate
    boundary = fluid_top
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = rhocpT
    subtract_mesh_velocity = false
  []
  [energy_balance]
    type = ParsedPostprocessor
    # the inflow rate is negative (outward-normal convention)
    expression = '(enthalpy_out + enthalpy_in - heat_in) / heat_in'
    pp_names = 'enthalpy_in enthalpy_out heat_in'
  []
[]

[VectorPostprocessors]
  [radial_solid]
    type = LineValueSampler
    variable = 'T_solid'
    # cell centres of the uniform solid mesh
    start_point = '${fparse R + (R_out - R) / (2 * nr_s)} 1.61 0'
    end_point = '${fparse R_out - (R_out - R) / (2 * nr_s)} 1.61 0'
    num_points = ${nr_s}
    sort_by = x
    execute_on = FINAL
  []
  [radial]
    type = LineValueSampler
    variable = 'T_fluid vel_y'
    # cell centres of the uniform fluid mesh
    start_point = '${fparse R / (2 * nr)} 1.61 0'
    end_point = '${fparse R - R / (2 * nr)} 1.61 0'
    num_points = ${nr}
    sort_by = x
    execute_on = FINAL
  []
[]

[Executioner]
  type = SIMPLE
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  energy_system = 'energy_system'
  solid_energy_system = 'solid_energy_system'
  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  turbulence_l_abs_tol = 1e-12
  energy_l_abs_tol = 1e-12
  solid_energy_l_abs_tol = 1e-12
  momentum_l_tol = 0
  pressure_l_tol = 0
  turbulence_l_tol = 0
  energy_l_tol = 0
  solid_energy_l_tol = 0
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.5 0.5'
  turbulence_field_relaxation = '0.5 0.5'
  energy_equation_relaxation = 0.95
  num_iterations = 2000
  pressure_absolute_tolerance = 1e-8
  momentum_absolute_tolerance = 1e-8
  turbulence_absolute_tolerance = '1e-8 1e-8'
  energy_absolute_tolerance = 1e-8
  solid_energy_absolute_tolerance = 1e-8
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
  cht_interfaces = 'interface'
  max_cht_fpi = 1
  cht_solid_flux_relaxation = 1.0
  cht_fluid_flux_relaxation = 1.0
  cht_solid_temperature_relaxation = 1.0
  cht_fluid_temperature_relaxation = 1.0
  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  csv = true
  execute_on = FINAL
[]
