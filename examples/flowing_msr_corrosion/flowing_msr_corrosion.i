# =====================================================================================================
# Complete flowing-MSR example: coupled flow, nuclear power, heat, turbulence, radiolysis AND
# corrosion in a fuel-salt channel, all solved together on one linear finite-volume segregated
# (SIMPLE) basis.
#
# A vertical fuel-salt channel (fluoride / FLiBe-relevant) carries salt upward through the active
# core between two structural-alloy (Hastelloy N) walls. The model solves, fully coupled:
#
#   * linear momentum + mass (pressure) conservation  -- Navier-Stokes with Rhie-Chow,
#   * k-epsilon turbulence                            -- turbulent viscosity mu_t,
#   * energy conservation                             -- nuclear heating, conduction + turbulent
#                                                        conduction, advected enthalpy,
#   * radiolysis                                      -- F2 produced with a neutron-flux-shaped
#                                                        source, transported by advection and by
#                                                        molecular + turbulent diffusion, with the
#                                                        evolved gas phase rising buoyantly,
#   * corrosion                                       -- the radiolytically oxidizing salt corrodes
#                                                        the alloy walls; chromium, iron and nickel
#                                                        shed by a Butler-Volmer wall reaction enter
#                                                        the salt and are carried up the core by the
#                                                        same turbulent flow.
#
# The nuclear heat and the radiolytic source follow the same separable cosine flux shape (peaked at
# the core mid-plane). The energy equation supplies the temperature field, which drives both the
# chemistry and the Arrhenius corrosion kinetics -- so the wall corrodes fastest where the salt is
# hottest. The corrosion products share the same turbulent diffusivity mu_t/(rho Sc_t) as the
# radiolysis species. The radiolytically oxidizing environment is represented by the oxidizing redox
# class of the corrosion correlation; a two-way coupling (the dissolved Cr feeding the radiolysis
# redox chemistry) is a natural extension.
#
# Run:  open_pronghorn-opt -i flowing_msr_corrosion.i
# View: the Exodus output (vel, T_fluid, mu_t, F2_diss, F2_gas, c_Cr/c_Fe/c_Ni corrosion plumes).
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
                     'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys F2_gas_sys '
                     'c_Cr_sys c_Fe_sys c_Ni_sys'
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
  # Effective (molecular + turbulent) diffusivity for the radiolysis and corrosion species.
  [D_eff_mat]
    type = ParsedFunctorMaterial
    property_name = D_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${D_mol} + mut/(${rho}*${Sc_t})'
  []
  # Wall corrosion rate [um/y] evaluated at the local temperature: the calibrated reference rate
  # (7.435 um/y at 950 K for this oxidizing Hastelloy N case) scaled by the Arrhenius factor of the
  # exchange current (Ea = 38.92 kJ/mol). Valid where the dissolved metal is dilute (c << c_ref), as
  # in this single-pass core; it is the rate the Butler-Volmer wall boundary imposes.
  [corrosion_rate_mat]
    type = ParsedFunctorMaterial
    property_name = corrosion_rate_functor
    functor_names = 'T_fluid'
    functor_symbols = 'T'
    expression = '7.43499 * exp(38921.89 / 8.31446261815324 * (1.0/950.0 - 1.0/T))'
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

[CorrosionPlatingFlow]
  elements = 'Cr Fe Ni'
  temperature = T_fluid               # Arrhenius corrosion driven by the solved temperature field
  reference_temperature = 950         # anchor temperature for the calibrated rate [K]
  material_class = hastelloy_n        # the structural alloy of the core walls
  salt_class = fluoride_fuel
  redox_class = oxidizing_fef2        # the radiolytically oxidizing core salt
  applied_overpotential = 0.12
  temperature_dependent_kinetics = true
  reaction_boundary = ${walls}        # the alloy walls corrode
  rhie_chow_user_object = 'rc'        # corrosion products advected by the solved turbulent flow
  advected_interp_method = ${advected_interp_method}
  diffusivity = D_eff                 # molecular + turbulent diffusion of the corrosion products
  inlet_boundary = 'bottom'           # fresh salt enters clean
  inlet_concentration = 0.0
  outlet_boundary = 'top'
  reference_concentration = 1.0
  time_derivative = false             # steady-state core (SIMPLE executioner)
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
  [corrosion_rate_um_y]
    type = MooseLinearVariableFVReal
    initial_condition = 0
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
  [compute_corrosion_rate]
    type = FunctorAux
    variable = corrosion_rate_um_y
    functor = corrosion_rate_functor
    execute_on = 'TIMESTEP_END'
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
  # --- corrosion product transport ---
  [Cr_outlet]
    type = SideAverageValue
    variable = c_Cr
    boundary = top
  []
  [Cr_max]
    type = ElementExtremeValue
    variable = c_Cr
    value_type = max
  []
  [Fe_outlet]
    type = SideAverageValue
    variable = c_Fe
    boundary = top
  []
  [Ni_outlet]
    type = SideAverageValue
    variable = c_Ni
    boundary = top
  []
  # Wall corrosion rate at the hottest point of the core (Arrhenius peak).
  [corrosion_rate_max]
    type = ElementExtremeValue
    variable = corrosion_rate_um_y
    value_type = max
  []
[]

[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  energy_system = 'energy_system'
  passive_scalar_systems = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys F2_gas_sys '
                           'c_Cr_sys c_Fe_sys c_Ni_sys'
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.3 0.3'
  turbulence_field_relaxation = '0.3 0.3'
  energy_equation_relaxation = 0.9
  passive_scalar_equation_relaxation = '0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9'
  num_iterations = 500
  momentum_absolute_tolerance = 1e-5
  pressure_absolute_tolerance = 1e-5
  turbulence_absolute_tolerance = '1e-4 1e-4'
  energy_absolute_tolerance = 1e-3
  passive_scalar_absolute_tolerance = '1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9'
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
