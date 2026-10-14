# =====================================================================================================
# Training case 16 (capstone): the MSFR primary loop -- flow, heat, radiolysis, corrosion AND the
# fission-product inventory with xenon stripping and heat-exchanger plating, on 24 ranks.
#
# A closed, pump-driven NaCl-UCl3 fuel-salt loop around a central reflector inside a reflector/vessel
# (examples/msfr_loop): k-epsilon flow driven by the pump and buoyancy, conjugate heat transfer with
# the solids, fission heat in the core leg (left) and heat removal in the heat-exchanger leg (right).
# On top of it everything from the day, in one steady segregated SIMPLE solve:
#
#   * chloride radiolysis with the Cr(II)/Cr(III) redox (module 1, case 14),
#   * Butler-Volmer corrosion of the 316 walls releasing Cr, Fe, Ni (module 3),
#   * the fission-product inventory (module 2) with the fast yields:
#       - I-135 -> Xe-135, with the xenon stripped by sparging in the top (hot) leg only,
#       - the noble metals Mo-99 (-> Tc-99m), Ru-103 and Te-132 plating on the heat-exchanger walls
#         with a Chilton-Colburn mass-transfer coefficient (case 15) and nowhere else.
#
# Why steady state is enough for the inventory: the loop time (~ 5 s) is much shorter than the
# half-lives (Xe-135 9.1 h, Mo-99 66 h), so the inventories are nearly uniform around the loop and
# their equilibrium is what the 0D cases 07 and 08 predicted. That is the validation:
#
#   'Xe135_vs_0d'  : loop-average Xe-135 / (P_I + P_Xe) / (lambda_Xe + R/V),  R = int r_strip dV
#   'Mo99_vs_0d'   : loop-average Mo-99  / P_Mo / (lambda_Mo + G/V),          G = int k_m dA
#   'Xe135_stripped_fraction', 'Mo99_plated_fraction': where the inventory goes,
#   'energy_balance_error' : fission power = heat-exchanger removal              (< 1 %)
#
# The Xe-135 burnout (sigma phi) is left out: the 2.6e6 b capture cross section of the database is a
# thermal value and does not apply to a fast spectrum (neutron_flux = 0).
# Corrosion in a closed loop: the dissolved metals have no outlet, so their steady state is the wall
# (Nernst) equilibrium, c ~ c_ref exp(z f eta) ~ 30 mol/m^3 per metal, with the radiolytic oxidant
# keeping Cr(III)/Cr(II) at ~ 5. Hot-wall dissolution then balances cold-wall plating -- the
# temperature-gradient mass transfer of case 13. The salt starts from that state (see the ICs).
# Results (60x60, 24 ranks, 1 min 35 s): T 971 -> 1036 K, energy balance 6e-5, Cr(III)/Cr(II) 5.1;
# Xe135_vs_0d = 1.00002 and Mo99_vs_0d = 1.0004 with the inventories uniform to 7e-4 around the loop;
# 83 % of the Xe-135 is stripped (the rest decays) and 94 % of the Mo-99 plates on the HX walls.
# Run:              mpirun -n 24 open_pronghorn-opt -i msfr_loop_capstone.i
# Mesh sensitivity: ... Mesh/gen/nx=90 Mesh/gen/ny=90 (see README)
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

# --- fission-product removal ---
r_sparge = 3.0e-4     # 1/s xenon stripping rate in the sparged top leg
D_mol = 2.0e-9        # m^2/s molecular diffusivity of the dissolved species
u_hx = 1.5            # m/s heat-exchanger leg velocity, D_h = core_w for the Sherwood number
# Chilton-Colburn: Sh = 0.023 Re^0.8 Sc^(1/3); Re ~ 6e5, Sc ~ 600 -> k_m ~ 2e-5 m/s
k_plate = ${fparse 0.023 * (rho_0 * u_hx * core_w / mu)^0.8 * (mu / (rho_0 * D_mol))^(1/3) * D_mol / core_w}

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
                     'Cr_II_sys Cr_III_sys Cr_I_sys c_Fe_sys c_Ni_sys '
                     'I135_sys Xe135_sys Mo99_sys Tc99m_sys Ru103_sys Te132_sys'
  previous_nl_solution_required = true
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
  [strip_rate_mat]
    type = ParsedFunctorMaterial
    property_name = strip_loss
    functor_names = 'strip_rate Xe135'
    functor_symbols = 'r N'
    expression = 'r * N'
    block = 'salt'
  []
  [plating_flux_mat]
    type = ParsedFunctorMaterial
    property_name = Mo99_plating_flux
    functor_names = 'plate_coefficient Mo99'
    functor_symbols = 'k N'
    expression = 'k * N'
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
  # Xenon stripping by sparging: the top (hot) leg only
  [strip_rate]
    type = ParsedFunction
    expression = 'if(y > ${fparse L - salt_in - core_w}, ${r_sparge}, 0.0)'
  []
  # Noble-metal plating: the heat-exchanger leg walls only
  [plate_coefficient]
    type = ParsedFunction
    expression = 'if(x > ${fparse L - salt_in - core_w} - 1e-6, ${k_plate}, 0.0)'
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
  # The loop is closed: the dissolved metals have no outlet and settle where the hot-wall dissolution
  # balances the cold-wall plating (Nernst back-reaction). The segregated iterations approach that
  # state only slowly (the wall exchange is the only coupling), so the salt starts from it: the
  # values are the converged loop averages of a previous run.
  initial_condition_species = 'Cl_ion Cr_II Cr_III'
  initial_condition_values = '2.0e4 30.38 156.2'
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
  initial_concentration = 30.4        # Fe and Ni near their wall equilibrium (see above)
[]

[FissionProductInventory]
  block = 'salt'
  nuclides = 'I135 Xe135 Mo99 Tc99m Ru103 Te132'
  yield_set = U235_fast
  power_density = fission_shape        # the same fission heat as the energy equation
  removal_groups = 'noble_gas'
  removal_rates = 'strip_rate'
  removal_destinations = 'offgas'
  track_removed_inventories = false    # the stripped and plated amounts follow from the rates
  plating_boundaries = 'inner_wall outer_wall'
  plating_mass_transfer_coefficient = plate_coefficient
  rhie_chow_user_object = 'rc'
  diffusivity = D_eff
  time_derivative = false              # equilibrium inventory (loop time << half-lives)
  verbose = true
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
    type = LinearFVPressureCorrectionDiffusion
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
  # --- fission-product inventory: equilibrium against the 0D well-mixed loop (cases 07, 08) ---
  [salt_volume]
    type = VolumePostprocessor
    block = 'salt'
    outputs = none
  []
  [strip_integral]
    type = FunctionElementIntegral
    function = strip_rate
    block = 'salt'
    outputs = none
  []
  [plating_conductance]
    type = SideIntegralFunctorPostprocessor
    functor = plate_coefficient
    boundary = 'inner_wall outer_wall'
    outputs = none
  []
  [Xe135_avg]
    type = ElementAverageValue
    variable = Xe135
    block = 'salt'
  []
  [Mo99_avg]
    type = ElementAverageValue
    variable = Mo99
    block = 'salt'
  []
  [Ru103_avg]
    type = ElementAverageValue
    variable = Ru103
    block = 'salt'
  []
  [Te132_avg]
    type = ElementAverageValue
    variable = Te132
    block = 'salt'
  []
  # Production rates [mol/s per m depth]: cumulative mass-135 yield (I-135 + direct Xe-135) and Mo-99
  [P135]
    type = ParsedPostprocessor
    expression = 'fission_power * 0.0660108 * 5.182134828e-14'
    pp_names = 'fission_power'
    outputs = none
  []
  [P_Mo99]
    type = ParsedPostprocessor
    expression = 'fission_power * 0.0594313 * 5.182134828e-14'
    pp_names = 'fission_power'
    outputs = none
  []
  [Xe135_vs_0d]
    type = ParsedPostprocessor
    expression = 'Xe135_avg / (P135 / salt_volume / (2.106574218e-05 + strip_integral / salt_volume))'
    pp_names = 'Xe135_avg P135 salt_volume strip_integral'
  []
  [Mo99_vs_0d]
    type = ParsedPostprocessor
    expression = 'Mo99_avg / (P_Mo99 / salt_volume / (2.918347331e-06 + plating_conductance / salt_volume))'
    pp_names = 'Mo99_avg P_Mo99 salt_volume plating_conductance'
  []
  [Xe135_stripped]
    type = ElementIntegralFunctorPostprocessor
    functor = strip_loss
    block = 'salt'
    outputs = none
  []
  [Xe135_stripped_fraction]
    type = ParsedPostprocessor
    expression = 'Xe135_stripped / P135'
    pp_names = 'Xe135_stripped P135'
  []
  [Mo99_plated]
    type = SideIntegralFunctorPostprocessor
    functor = Mo99_plating_flux
    boundary = 'inner_wall outer_wall'
    outputs = none
  []
  [Mo99_plated_fraction]
    type = ParsedPostprocessor
    expression = 'Mo99_plated / P_Mo99'
    pp_names = 'Mo99_plated P_Mo99'
  []
  # Loop non-uniformity of the inventory (max/min - 1): small because the loop time << half-life
  [Xe135_max]
    type = ElementExtremeValue
    variable = Xe135
    value_type = max
    block = 'salt'
    outputs = none
  []
  [Xe135_min]
    type = ElementExtremeValue
    variable = Xe135
    value_type = min
    block = 'salt'
    outputs = none
  []
  [Xe135_nonuniformity]
    type = ParsedPostprocessor
    expression = 'Xe135_max / Xe135_min - 1'
    pp_names = 'Xe135_max Xe135_min'
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
                           'Cr_II_sys Cr_III_sys Cr_I_sys c_Fe_sys c_Ni_sys '
                           'I135_sys Xe135_sys Mo99_sys Tc99m_sys Ru103_sys Te132_sys'
  momentum_equation_relaxation = 0.5
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.3 0.3'
  turbulence_field_relaxation = '0.3 0.3'
  energy_equation_relaxation = 0.9
  # No relaxation for the metals and fission products: relaxation 0.9 stalls the closed-loop metal
  # inventories (no outlet), and the linear fission-product equations need none
  passive_scalar_equation_relaxation = '0.9 0.9 0.9 0.9 0.9 0.9 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 '
                                       '1.0 1.0 1.0'
  num_iterations = 800
  pressure_absolute_tolerance = 1e-7
  momentum_absolute_tolerance = 1e-7
  turbulence_absolute_tolerance = '1e-6 1e-6'
  energy_absolute_tolerance = 1e-5
  solid_energy_absolute_tolerance = 1e-5
  passive_scalar_absolute_tolerance = '1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9'
                                      ' 1e-14 1e-14 1e-14 1e-14 1e-14 1e-14'
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

[VectorPostprocessors]
  # Centreline of the sparged top (hot) leg, core side to heat-exchanger side
  [loop_line]
    type = LineValueSampler
    variable = 'T_fluid Xe135 Mo99 Cr_III'
    start_point = '0.775 2.275 0'
    end_point = '2.225 2.275 0'
    num_points = 30
    sort_by = x
    execute_on = FINAL
  []
[]

[Outputs]
  exodus = true
  [out]
    type = CSV
    execute_vector_postprocessors_on = FINAL
  []
[]
