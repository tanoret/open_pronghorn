# =====================================================================================================
# Molten salt fast reactor (MSFR) PRIMARY LOOP -- a pump-driven closed fuel-salt circuit around a
# central solid reflector, enclosed by a solid reflector/vessel, with conjugate heat transfer.
#
# Unlike a channel pushed by inlet/outlet boundary conditions, the fuel salt CIRCULATES around a
# closed loop driven by a PUMP (a momentum source) and assisted by buoyancy. The loop is a rectangular
# salt annulus around a central solid reflector block, enclosed by an outer solid reflector/vessel --
# the surrounding reactor structure is modeled as real solid SUBDOMAINS, not boundary conditions.
#
#   * core leg (left)   : the active core -- volumetric fission heat (and later radiolysis),
#   * top leg           : hot leg,
#   * heat-exchanger leg (right) : heat removed to the secondary circuit (Newton cooling toward the
#                                  cold-leg temperature, which also sets the loop temperature level),
#   * bottom leg        : the pump (a momentum body force closing the clockwise circulation),
#   * solids            : the central reflector and the outer reflector/vessel conduct heat
#                         (conjugate heat transfer with the salt across the fluid-solid interfaces).
#
# The model integrates, in one segregated (SIMPLE) linear finite-volume solve: turbulent flow
# (k-epsilon) driven by the pump and buoyancy; conjugate heat transfer (fluid + solid temperatures
# coupled at the interfaces); chloride radiolysis (the oxidant Cl2.- produced in the core leg); and
# corrosion of the solid walls (a Butler-Volmer reaction at the true fluid-solid interfaces releasing
# chromium as Cr(II), which the radiolytic oxidant then oxidizes to Cr(III)).
#
# Run:  open_pronghorn-opt -i msfr_loop.i
# View: the Exodus output -- the recirculating velocity, T_fluid/T_solid (conjugate heat transfer),
#       mu_t, Cl2m_rad, and the dissolved chromium Cr(II)/Cr(III) released from the walls.
# =====================================================================================================

L = 3.0               # m, outer size of the reactor block
salt_in = 0.5         # m, outer solid frame thickness
core_w = 0.5          # m, salt leg width

# --- NaCl-UCl3 fuel salt (~650-750 C) ---
rho_0 = 3300          # kg/m^3 reference density
mu = 4.0e-3           # Pa s
cp = 1000             # J/kg/K
k_salt = 0.4          # W/m/K
alpha_b = 2.0e-4      # 1/K thermal expansion coefficient
T_cold = 923.15       # K cold-leg (heat-exchanger) temperature, 650 C
T_ref = 973.15        # K reference temperature for buoyancy

# --- Ni-based reflector / vessel structure ---
k_solid = 23.0        # W/m/K (nickel alloy)

# --- operating point ---
q_core = 8.0e7        # W/m^3 volumetric fission heat in the core leg
h_hx = 1.0e6          # W/m^3/K volumetric heat-transfer coefficient of the heat exchanger
pump_force = 1.5e3    # N/m^3 pump body force in the bottom leg

# --- k-epsilon turbulence closure ---
C_mu = 0.09
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
intensity = 0.05
u_scale = 1.5         # m/s nominal loop speed (for the turbulence initial conditions)
k_init = '${fparse 1.5*(intensity*u_scale)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / (0.07*core_w)}'
walls = 'inner_wall outer_wall'
wall_treatment = 'eq_newton'

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${L}
    ymin = 0
    ymax = ${L}
    nx = 60
    ny = 60
  []
  [reflector]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    bottom_left = '1.0 1.0 0'
    top_right = '2.0 2.0 0'
  []
  [vessel]
    type = ParsedSubdomainMeshGenerator
    input = reflector
    combinatorial_geometry = 'x < ${salt_in} | x > ${fparse L - salt_in} | y < ${salt_in} | y > ${fparse L - salt_in}'
    block_id = 2
  []
  [names]
    type = RenameBlockGenerator
    input = vessel
    old_block = '0 1 2'
    new_block = 'salt reflector vessel'
  []
  [inner_wall]
    type = SideSetsBetweenSubdomainsGenerator
    input = names
    primary_block = salt
    paired_block = reflector
    new_boundary = inner_wall
  []
  [outer_wall]
    type = SideSetsBetweenSubdomainsGenerator
    input = inner_wall
    primary_block = salt
    paired_block = vessel
    new_boundary = outer_wall
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system TKE_system TKED_system energy_system '
                     'solid_energy_system '
                     'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys '
                     'Cr_II_sys Cr_III_sys Cr_I_sys c_Fe_sys c_Ni_sys'
  previous_nl_solution_required = true
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method = 'upwind'
[]

[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho_0}
    p_diffusion_kernel = p_diffusion
    block = 'salt'
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_condition = 1e-6
    block = 'salt'
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_condition = 1e-6
    block = 'salt'
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 0
    block = 'salt'
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    initial_condition = ${k_init}
    block = 'salt'
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    initial_condition = ${eps_init}
    block = 'salt'
  []
  [T_fluid]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
    initial_condition = ${T_cold}
    block = 'salt'
  []
  [T_solid]
    type = MooseLinearVariableFVReal
    solver_sys = solid_energy_system
    initial_condition = ${T_cold}
    block = 'reflector vessel'
  []
[]

[FunctorMaterials]
  # Boussinesq density for buoyancy.
  [rho_mat]
    type = ParsedFunctorMaterial
    property_name = rho
    functor_names = 'T_fluid'
    functor_symbols = 'T'
    expression = '${rho_0} * (1 - ${alpha_b} * (T - ${T_ref}))'
    block = 'salt'
  []
  # Placeholder material on the solids (MOOSE requires every block to carry a material once any does).
  [solid_mat]
    type = GenericFunctorMaterial
    prop_names = 'rho'
    prop_values = '10000'
    block = 'reflector vessel'
  []
  # Effective (molecular + turbulent) thermal conductivity of the salt.
  [k_eff_mat]
    type = ParsedFunctorMaterial
    property_name = k_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${k_salt} + mut*${cp}/0.9'
    block = 'salt'
  []
  # Effective (molecular + turbulent) diffusivity for the radiolysis and corrosion species.
  [D_eff_mat]
    type = ParsedFunctorMaterial
    property_name = D_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '2.0e-9 + mut/(${rho_0}*0.7)'
    block = 'salt'
  []
  # Fission heat and heat-exchanger removal as functors, for the energy-balance verification.
  [fission_functor_mat]
    type = ParsedFunctorMaterial
    property_name = fission_functor
    expression = 'if(x < ${fparse salt_in + core_w}, ${q_core} * cos(pi*(y-0.5*${L})/(1.4*${L})), 0.0)'
    block = 'salt'
  []
  [hx_removal_mat]
    type = ParsedFunctorMaterial
    property_name = hx_removal
    functor_names = 'T_fluid'
    functor_symbols = 'T'
    expression = 'if(x > ${fparse L - salt_in - core_w}, ${h_hx} * (T - ${T_cold}), 0.0)'
    block = 'salt'
  []
[]

[AuxVariables]
  [wall_distance]
    type = MooseVariableFVReal
    initial_condition = 0.1
    block = 'salt'
  []
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse rho_0 * C_mu * k_init^2 / eps_init}'
    block = 'salt'
  []
[]

[AuxKernels]
  [compute_wall_distance]
    type = WallDistanceAux
    variable = wall_distance
    walls = ${walls}
    execute_on = 'INITIAL'
    block = 'salt'
  []
  [compute_mu_t]
    type = kEpsilonViscosity
    variable = mu_t
    C_mu = ${C_mu}
    tke = TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho_0}
    u = vel_x
    v = vel_y
    bulk_wall_treatment = false
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    mu_t_ratio_max = 1e5
    execute_on = 'NONLINEAR'
    wall_distance = wall_distance
    block = 'salt'
  []
[]

[Functions]
  # Pump body force: bottom leg, in -x, closing the clockwise circulation.
  [pump_profile]
    type = ParsedFunction
    expression = 'if(y < ${fparse salt_in + core_w}, -${pump_force}, 0.0)'
  []
  # Fission heat: the core leg (left), a cosine axial shape peaked at the core mid-height.
  [fission_shape]
    type = ParsedFunction
    expression = 'if(x < ${fparse salt_in + core_w}, ${q_core} * cos(pi*(y-0.5*${L})/(1.4*${L})), 0.0)'
  []
  # Heat-exchanger coverage: the HX leg (right).
  [hx_mask]
    type = ParsedFunction
    expression = 'if(x > ${fparse L - salt_in - core_w}, ${h_hx}, 0.0)'
  []
  [hx_source]
    type = ParsedFunction
    expression = 'if(x > ${fparse L - salt_in - core_w}, ${fparse h_hx * T_cold}, 0.0)'
  []
  # Radiolytic dose: peaked in the core leg, following the same neutron-flux shape as the fission heat.
  [dose_field]
    type = ParsedFunction
    expression = 'if(x < ${fparse salt_in + core_w}, ${fparse 5.0e5} * cos(pi*(y-0.5*${L})/(1.4*${L})), 0.0)'
  []
[]

[MoltenSaltRadiolysis]
  block = 'salt'
  salt_type = chloride
  metals = 'Cr'                       # Cr(II)/Cr(III)/Cr(I) and their radiolytic reactions
  temperature = T_fluid
  dose_rate = dose_field              # radiolytic oxidant produced in the core leg
  g_value_species = 'e_sol Cl_rad Cl2m_rad'
  g_value_overrides = '0.0 0.0 0.30'  # lump the fast kinetics into the net oxidant Cl2.-
  rhie_chow_user_object = 'rc'        # advected around the loop by the solved flow
  diffusivity = D_eff
  time_derivative = false             # steady operating point (SIMPLE)
  initial_condition_species = 'Cl_ion Cr_II'
  initial_condition_values = '2.0e4 3.0'
[]

[CorrosionPlatingFlow]
  block = 'salt'
  elements = 'Cr Fe Ni'
  # Chromium is released into the radiolysis-tracked Cr(II) (coupling corrosion to radiolysis); iron
  # and nickel have no chloride radiolysis chemistry, so they get their own corrosion-product scalars
  # c_Fe and c_Ni (whose systems c_Fe_sys / c_Ni_sys are listed in [Problem] and the executioner).
  release_variables = 'Cr_II none none'
  reaction_boundary = 'inner_wall outer_wall'   # the true fluid-solid interfaces
  temperature = T_fluid
  reference_temperature = 973.15
  material_class = stainless_316
  salt_class = chloride
  redox_class = oxidizing_fef2
  applied_overpotential = 0.1
  reference_concentration = 3.0
  temperature_dependent_kinetics = true
  rhie_chow_user_object = 'rc'
  diffusivity = D_eff
  time_derivative = false             # steady operating point (SIMPLE); no time integrator
[]

[LinearFVKernels]
  # --- momentum: viscous stress + pressure + buoyancy + pump ---
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'x'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
    block = 'salt'
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = ${mu}
    block = 'salt'
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
    block = 'salt'
  []
  [pump]
    type = LinearFVSource
    variable = vel_x
    source_density = pump_profile
    block = 'salt'
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
    block = 'salt'
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = ${mu}
    block = 'salt'
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
    block = 'salt'
  []
  [buoyancy]
    type = LinearFVMomentumBuoyancy
    variable = vel_y
    rho = 'rho'
    reference_rho = ${rho_0}
    gravity = '0 -9.81 0'
    momentum_component = 'y'
    block = 'salt'
  []
  [p_diffusion]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
    block = 'salt'
  []
  [HbyA_divergence]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = true
    block = 'salt'
  []

  # --- k-epsilon turbulence ---
  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
    block = 'salt'
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = 'salt'
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = false
    block = 'salt'
  []
  [TKE_source_sink]
    type = kEpsilonTKESourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    epsilon = TKED
    rho = ${rho_0}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    block = 'salt'
  []
  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = ${walls}
    block = 'salt'
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    walls = ${walls}
    block = 'salt'
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    walls = ${walls}
    block = 'salt'
  []
  [TKED_source_sink]
    type = kEpsilonTKEDSourceSink
    variable = TKED
    u = vel_x
    v = vel_y
    tke = TKE
    rho = ${rho_0}
    mu = ${mu}
    mu_t = 'mu_t'
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    wall_distance = wall_distance
    block = 'salt'
  []

  # --- fluid energy: advection + conduction + fission heat + heat-exchanger Newton cooling ---
  [energy_advection]
    type = LinearFVEnergyAdvection
    variable = T_fluid
    advected_quantity = temperature
    cp = ${cp}
    block = 'salt'
  []
  [energy_conduction]
    type = LinearFVDiffusion
    variable = T_fluid
    diffusion_coeff = 'k_eff'
    use_nonorthogonal_correction = false
    block = 'salt'
  []
  [fission_heat]
    type = LinearFVSource
    variable = T_fluid
    source_density = fission_shape
    block = 'salt'
  []
  [hx_sink]
    type = LinearFVReaction
    variable = T_fluid
    coeff = hx_mask
    block = 'salt'
  []
  [hx_secondary]
    type = LinearFVSource
    variable = T_fluid
    source_density = hx_source
    block = 'salt'
  []

  # --- solid energy: conduction in the reflector and vessel ---
  [solid_conduction]
    type = LinearFVDiffusion
    variable = T_solid
    diffusion_coeff = ${k_solid}
    use_nonorthogonal_correction = false
    block = 'reflector vessel'
  []
[]

[LinearFVBCs]
  # No-slip at the fluid-solid interfaces (the salt is enclosed; no inlet or outlet).
  [walls_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inner_wall outer_wall'
    variable = vel_x
    functor = 0.0
  []
  [walls_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inner_wall outer_wall'
    variable = vel_y
    functor = 0.0
  []
  [walls_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = 'inner_wall outer_wall'
    variable = mu_t
    u = vel_x
    v = vel_y
    rho = ${rho_0}
    mu = ${mu}
    tke = TKE
    wall_treatment = ${wall_treatment}
  []
  # Conjugate heat transfer across the fluid-solid interfaces.
  [fluid_to_solid]
    type = LinearFVConvectiveHeatTransferBC
    variable = T_fluid
    T_fluid = T_fluid
    T_solid = T_solid
    boundary = 'inner_wall outer_wall'
    h = 1.0e4
  []
  [solid_to_fluid_inner]
    type = LinearFVConvectiveHeatTransferBC
    variable = T_solid
    T_fluid = T_fluid
    T_solid = T_solid
    boundary = 'inner_wall outer_wall'
    h = 1.0e4
  []
[]

[Postprocessors]
  [vel_y_core_leg_max]
    type = ElementExtremeValue
    variable = vel_y
    value_type = max
    block = 'salt'
  []
  [vel_y_hx_leg_min]
    type = ElementExtremeValue
    variable = vel_y
    value_type = min
    block = 'salt'
  []
  [T_fluid_max]
    type = ElementExtremeValue
    variable = T_fluid
    value_type = max
    block = 'salt'
  []
  [T_fluid_min]
    type = ElementExtremeValue
    variable = T_fluid
    value_type = min
    block = 'salt'
  []
  [T_solid_max]
    type = ElementExtremeValue
    variable = T_solid
    value_type = max
    block = 'reflector vessel'
  []
  # --- radiolysis + corrosion coupling ---
  [Cl2m_rad_max]
    type = ElementExtremeValue
    variable = Cl2m_rad
    value_type = max
    block = 'salt'
  []
  [Cr_II_avg]
    type = ElementAverageValue
    variable = Cr_II
    block = 'salt'
  []
  [Cr_III_avg]
    type = ElementAverageValue
    variable = Cr_III
    block = 'salt'
  []
  [redox_ratio_CrIII_CrII]
    type = ParsedPostprocessor
    pp_names = 'Cr_II_avg Cr_III_avg'
    expression = 'Cr_III_avg / (Cr_II_avg + 1e-30)'
  []
  [Fe_avg]
    type = ElementAverageValue
    variable = c_Fe
    block = 'salt'
  []
  [Ni_avg]
    type = ElementAverageValue
    variable = c_Ni
    block = 'salt'
  []
  # --- energy balance verification: at steady state the fission power equals the heat-exchanger
  #     removal (the loop is closed and the walls are adiabatic on the outside). ---
  [fission_power]
    type = ElementIntegralFunctorPostprocessor
    functor = fission_functor
    block = 'salt'
  []
  [hx_removal_power]
    type = ElementIntegralFunctorPostprocessor
    functor = hx_removal
    block = 'salt'
  []
  [energy_balance_error]
    type = ParsedPostprocessor
    pp_names = 'fission_power hx_removal_power'
    expression = 'abs(fission_power - hx_removal_power) / fission_power'
  []
[]

[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  energy_system = 'energy_system'
  solid_energy_system = 'solid_energy_system'
  passive_scalar_systems = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys '
                           'Cr_II_sys Cr_III_sys Cr_I_sys c_Fe_sys c_Ni_sys'
  momentum_equation_relaxation = 0.5
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.3 0.3'
  turbulence_field_relaxation = '0.3 0.3'
  energy_equation_relaxation = 0.9
  passive_scalar_equation_relaxation = '0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9'
  num_iterations = 800
  pressure_absolute_tolerance = 1e-7
  momentum_absolute_tolerance = 1e-7
  turbulence_absolute_tolerance = '1e-6 1e-6'
  energy_absolute_tolerance = 1e-5
  solid_energy_absolute_tolerance = 1e-5
  passive_scalar_absolute_tolerance = '1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9'
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'
  passive_scalar_petsc_options_iname = '-pc_type -pc_hypre_type'
  passive_scalar_petsc_options_value = 'hypre boomeramg'
  turbulence_l_max_its = 30
  pin_pressure = true
  pressure_pin_value = 0.0
  pressure_pin_point = '0.75 0.75 0.0'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'
  solid_energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  solid_energy_petsc_options_value = 'hypre boomeramg'
  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  exodus = true
  csv = true
[]
