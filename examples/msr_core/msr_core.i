# =====================================================================================================
# Full MSR core example: coupled flow, heat, turbulence and radiolysis in a fuel-salt channel.
#
# A vertical fuel-salt channel (fluoride / FLiBe-relevant) carries salt upward through the active
# core. The model solves, fully coupled, on a linear finite-volume segregated (SIMPLE) basis:
#
#   * linear momentum + mass (pressure) conservation  -- Navier-Stokes with Rhie-Chow,
#   * k-epsilon turbulence                            -- turbulent viscosity mu_t,
#   * energy conservation                             -- nuclear heating, conduction + turbulent
#                                                        conduction, advected enthalpy,
#   * radiolysis                                      -- F2 produced with a neutron-flux-shaped
#                                                        source, transported by advection and by
#                                                        molecular + turbulent diffusion, with the
#                                                        evolved gas phase rising buoyantly.
#
# The radiolytic source and the nuclear heat both follow the same separable cosine flux shape
# (peaked at the core mid-plane). The energy equation supplies the temperature field that drives the
# chemistry (here the gas-exchange equilibrium). The species turbulent diffusivity is mu_t/(rho Sc_t)
# from the k-epsilon solution, added to the molecular diffusivity.
#
# To keep the (otherwise stiff) chemistry inexpensive the radical yields are lumped into the net F2
# yield (G = 0.005); remove the g_value overrides to integrate the explicit radical network.
#
# Run:  open_pronghorn-opt -i msr_core.i      View: the Exodus output (vel, T_fluid, mu_t, F2_diss, F2_gas).
# =====================================================================================================

# --- geometry ---
W = 0.3               # m, channel half-pitch (radial)
H = 2.0               # m, active core height (axial)

# --- molten fluoride salt properties (FLiBe-relevant, ~700 C) ---
rho = 2000            # kg/m^3
mu = 1.0e-3           # Pa s (molecular)
cp = 2000             # J/kg/K
k_mol = 1.0           # W/m/K (molecular conductivity)
D_mol = 1.0e-9        # m^2/s (molecular species diffusivity)
Pr_t = 0.9            # turbulent Prandtl number
Sc_t = 0.7            # turbulent Schmidt number

# --- operating point ---
u_in = 0.5            # m/s inlet velocity (upward)
T_in = 900            # K inlet temperature
q_peak = 5.0e7        # W/m^3 peak nuclear heat density
dose_peak = 1.0e5     # J/m^3/s peak radiolytic dose

# --- k-epsilon closure ---
C_mu = 0.09
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
intensity = 0.05
k_init = '${fparse 1.5*(intensity*u_in)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / (0.07*W)}'
walls = 'left right'
wall_treatment = 'eq_newton'
advected_interp_method = 'upwind'

[Mesh]
  [core]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${W}
    ymin = 0
    ymax = ${H}
    nx = 6
    ny = 24
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system TKE_system TKED_system energy_system '
                     'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys F2_gas_sys'
  previous_nl_solution_required = true
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method = ${advected_interp_method}
[]

[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_condition = 0
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_condition = ${u_in}
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 1e-8
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    initial_condition = ${k_init}
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    initial_condition = ${eps_init}
  []
  [T_fluid]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
    initial_condition = ${T_in}
  []
[]

[FunctorMaterials]
  [cp_mat]
    type = GenericFunctorMaterial
    prop_names = 'cp'
    prop_values = '${cp}'
  []
  # Effective (molecular + turbulent) conductivity for the energy equation.
  [k_eff_mat]
    type = ParsedFunctorMaterial
    property_name = k_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${k_mol} + mut*${cp}/${Pr_t}'
  []
  # Effective (molecular + turbulent) diffusivity for the radiolysis species.
  [D_eff_mat]
    type = ParsedFunctorMaterial
    property_name = D_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${D_mol} + mut/(${rho}*${Sc_t})'
  []
[]

[Functions]
  # Separable cosine flux shape, peaked at the core mid-plane (the neutron-flux profile that both
  # the nuclear heat and the radiolytic source follow).
  [flux_shape]
    type = ParsedFunction
    expression = 'cos(pi*(x-0.5*${W})/(1.4*${W})) * cos(pi*(y-0.5*${H})/(1.4*${H}))'
  []
  [dose_field]
    type = ParsedFunction
    expression = '${dose_peak} * f'
    symbol_names = 'f'
    symbol_values = 'flux_shape'
  []
[]

[MoltenSaltRadiolysis]
  salt_type = fluoride
  gas_species = 'F2'
  temperature = T_fluid              # chemistry driven by the solved temperature field
  dose_rate = dose_field             # neutron-flux-shaped radiolytic source
  g_value_species = 'e_sol F_rad'    # lump radicals into the net F2 yield
  g_value_overrides = '0.0 0.0'
  rhie_chow_user_object = 'rc'       # advection by the solved turbulent flow
  diffusivity = D_eff                # molecular + turbulent diffusion
  time_derivative = false            # steady-state core (SIMPLE executioner)
  kLa = 0.02
  gas_buoyancy_velocity = '0 0.1 0'  # gas rises relative to the salt
  gas_dispersivity = D_eff
[]

[LinearFVKernels]
  # --- momentum (turbulent stress mu_t + molecular mu) + pressure ---
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'x'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = '${mu}'
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
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
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = '${mu}'
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
  []
  [p_diffusion]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
  []
  [HbyA_divergence]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = true
  []

  # --- k-epsilon turbulence ---
  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = false
  []
  [TKE_source_sink]
    type = kEpsilonTKESourceSink
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
  []
  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = ${walls}
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    walls = ${walls}
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    walls = ${walls}
  []
  [TKED_source_sink]
    type = kEpsilonTKEDSourceSink
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
    wall_distance = wall_distance
  []

  # --- energy: advected enthalpy + (molecular + turbulent) conduction + nuclear heat ---
  [energy_advection]
    type = LinearFVEnergyAdvection
    variable = T_fluid
    advected_quantity = temperature
    cp = ${cp}
  []
  [energy_conduction]
    type = LinearFVDiffusion
    variable = T_fluid
    diffusion_coeff = 'k_eff'
    use_nonorthogonal_correction = false
  []
  [nuclear_heat]
    type = LinearFVSource
    variable = T_fluid
    source_density = flux_shape
    scaling_factor = ${q_peak}
  []
[]

[LinearFVBCs]
  # --- flow ---
  [inlet_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    variable = vel_x
    functor = 0.0
  []
  [inlet_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    variable = vel_y
    functor = ${u_in}
  []
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
  [outlet_u]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = top
    variable = vel_x
    use_two_term_expansion = false
  []
  [outlet_v]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = top
    variable = vel_y
    use_two_term_expansion = false
  []
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = top
    variable = pressure
    functor = 0.0
  []
  # --- turbulence ---
  [inlet_TKE]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    variable = TKE
    functor = ${k_init}
  []
  [outlet_TKE]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = top
    variable = TKE
    use_two_term_expansion = false
  []
  [inlet_TKED]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    variable = TKED
    functor = ${eps_init}
  []
  [outlet_TKED]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = top
    variable = TKED
    use_two_term_expansion = false
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
  # --- energy ---
  [inlet_T]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    variable = T_fluid
    functor = ${T_in}
  []
  [outlet_T]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = top
    variable = T_fluid
    use_two_term_expansion = false
  []
  # --- radiolysis species: clean salt enters at the bottom; everything leaves at the top outlet ---
  [inlet_e_sol]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    variable = e_sol
    functor = 0
  []
  [outlet_e_sol]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = top
    variable = e_sol
    use_two_term_expansion = false
  []
  [inlet_F_ion]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    variable = F_ion
    functor = 0
  []
  [outlet_F_ion]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = top
    variable = F_ion
    use_two_term_expansion = false
  []
  [inlet_F_rad]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    variable = F_rad
    functor = 0
  []
  [outlet_F_rad]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = top
    variable = F_rad
    use_two_term_expansion = false
  []
  [inlet_F2_diss]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    variable = F2_diss
    functor = 0
  []
  [outlet_F2_diss]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = top
    variable = F2_diss
    use_two_term_expansion = false
  []
  [inlet_F2_gas]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    variable = F2_gas
    functor = 0
  []
  [outlet_F2_gas]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = top
    variable = F2_gas
    use_two_term_expansion = false
  []
[]

[AuxVariables]
  [wall_distance]
    type = MooseVariableFVReal
    initial_condition = 1.0
  []
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse rho * C_mu * k_init^2 / eps_init}'
  []
[]

[AuxKernels]
  [compute_wall_distance]
    type = WallDistanceAux
    variable = wall_distance
    walls = ${walls}
    execute_on = 'INITIAL'
  []
  [compute_mu_t]
    type = kEpsilonViscosity
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
    mu_t_ratio_max = 1e20
    execute_on = 'NONLINEAR'
    wall_distance = wall_distance
  []
[]

[Postprocessors]
  [T_outlet]
    type = SideAverageValue
    variable = T_fluid
    boundary = top
  []
  [mu_t_max]
    type = ElementExtremeValue
    variable = mu_t
    value_type = max
  []
  [F2_diss_max]
    type = ElementExtremeValue
    variable = F2_diss
    value_type = max
  []
  [F2_diss_outlet]
    type = SideAverageValue
    variable = F2_diss
    boundary = top
  []
  [F2_gas_outlet]
    type = SideAverageValue
    variable = F2_gas
    boundary = top
  []
[]

[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  energy_system = 'energy_system'
  passive_scalar_systems = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys F2_gas_sys'
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.3 0.3'
  turbulence_field_relaxation = '0.3 0.3'
  energy_equation_relaxation = 0.9
  passive_scalar_equation_relaxation = '0.9 0.9 0.9 0.9 0.9'
  num_iterations = 500
  momentum_absolute_tolerance = 1e-5
  pressure_absolute_tolerance = 1e-5
  turbulence_absolute_tolerance = '1e-4 1e-4'
  energy_absolute_tolerance = 1e-3
  passive_scalar_absolute_tolerance = '1e-9 1e-9 1e-9 1e-9 1e-9'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'
  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'
  passive_scalar_petsc_options_iname = '-pc_type -pc_hypre_type'
  passive_scalar_petsc_options_value = 'hypre boomeramg'
  momentum_l_max_its = 300
  pressure_l_max_its = 300
  turbulence_l_max_its = 30
  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  exodus = true
  csv = true
[]
