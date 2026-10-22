################################################################################
# GFMR training -- Stage 4: fuel-salt chemistry on the 3-D core flow
#
# Restarts the converged CFD of stage 3 (../03_cfd/gfmr_core_cfd_out.e: velocity, pressure,
# k-epsilon, salt and graphite temperatures) and solves ONLY the chemistry, as passive scalars
# on the frozen flow (should_solve_* = false for the flow and energy):
#
#   [MoltenSaltRadiolysis]   U(IV)/U(III) redox buffer of the FUNaK salt (26.4 mol% UF4). The
#                            radiolytic dose is 0.87 x the salt power density (fragments + beta +
#                            part of the gamma). Charge-consistent radiolysis: the e_sol / F pairs
#                            recombine in the hot salt, so the net radiolytic yields are zero
#                            (G = 0), and the salt is oxidized by FISSION (0.8 eq per fission).
#   [CorrosionPlatingFlow]   element-resolved Cr / Fe / Ni corrosion and plating of the
#                            Hastelloy N guide tubes, the only metal in the core (9 rods).
#   [FissionProductInventory] I-135 / Xe-135 with neutron burnout and graphite permeation.
#
# The once-through core needs inlet values. The GFMR specification has no loop, so the inlets
# come from the loop balance of loop_balance.py (ASSUMED ex-core volume 2 x core salt, 9.1 s loop
# time, 1 day redox control, 10 % Xe stripping per pass); see README.
#
# Run: mpirun -n 24 open_pronghorn-opt -i gfmr_chemistry.i
################################################################################

!include inlet_values.i

# ---- FUNaK at 650 C (SPEC), as in stage 3
rho = 4073.5
mu = 0.01494

# ---- transport
D_mol = 2.0e-9        # m^2/s molecular diffusivity of dissolved species in a fluoride salt
Sc_t = 0.85           # turbulent Schmidt number (= Pr_t of the energy equation)
dose_fraction = 0.87  # radiolytic dose / fission power density in a fuel salt

# ---- corrosion of the Hastelloy N guide tube (same database class and equilibria as the
# MSFR capstone of examples/training/chemistry; GFMR-specific FUNaK data do not exist)
T_ref = ${fparse 650 + 273.15}
c_eq_Cr = 4.0         # mol/m^3 Cr equilibrium at U(III)/U ~ 1 % (fluoride fuel class)
c_eq_Fe = 0.05
c_eq_Ni = 1.0e-4
f_Cr = 0.90           # selective Cr leaching of the Ni alloy
f_Fe = 0.08
f_Ni = 0.02
dH_diss = 50.0        # kJ/mol dissolution enthalpy of the equilibria (assumed, all elements)

# ---- graphite uptake of Xe-135: k_m to the channel walls (Sh = 0.023 Re^0.8 Sc^(1/3)) times a
# graphite acceptance of 1 (the graphite is a perfect xenon sink, an upper bound of the
# permeation; MSRE graphite was sealed to limit it). Inner channel: Re 1.85e4, D = 4 cm
Sc_mol = ${fparse mu / (rho * D_mol)}
k_m_graphite = ${fparse 0.023 * 18500^0.8 * Sc_mol^(1/3) * D_mol / 0.04}

salt = 'salt_inner_tri salt_inner salt_outer_tri salt_outer salt_cr'
channel_walls = 'wall_inner wall_outer wall_cr'
walls = 'wall_inner wall_outer wall_cr wall_guide'
inlets = 'inlet_inner inlet_outer inlet_cr'

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = ../03_cfd/gfmr_core_cfd_out.e
    use_for_exodus_restart = true
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system w_system pressure_system TKE_system TKED_system energy_system '
                     'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys U_IV_sys U_III_sys '
                     'c_Cr_sys c_Fe_sys c_Ni_sys I135_sys Xe135_sys'
  previous_nl_solution_required = true
  allow_initial_conditions_with_restart = true
  material_coverage_check = false
  kernel_coverage_check = false
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
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_from_file_var = vel_x
    block = ${salt}
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_from_file_var = vel_y
    block = ${salt}
  []
  [vel_z]
    type = MooseLinearVariableFVReal
    solver_sys = w_system
    initial_from_file_var = vel_z
    block = ${salt}
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_from_file_var = pressure
    block = ${salt}
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    initial_from_file_var = TKE
    block = ${salt}
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    initial_from_file_var = TKED
    block = ${salt}
  []
  [T_fluid]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
    initial_from_file_var = T_fluid
    block = ${salt}
  []
[]

[AuxVariables]
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_from_file_var = mu_t
    block = ${salt}
  []
  [q_salt]
    type = MooseVariableFVReal
    initial_from_file_var = q_salt
    block = ${salt}
  []
[]

[FunctorMaterials]
  [D_eff]
    type = ParsedFunctorMaterial
    property_name = D_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${D_mol} + mut / (${rho} * ${Sc_t})'
    block = ${salt}
  []
  # q_salt of the CFD is the diffusion shape; scaled to 92.1 % of 25 MW in stage 3
  [power_density]
    type = ParsedFunctorMaterial
    property_name = power_density
    functor_names = 'q_salt'
    functor_symbols = 'q'
    expression = 'q * 0.9213'
    block = ${salt}
  []
  [dose_rate]
    type = ParsedFunctorMaterial
    property_name = dose_rate
    functor_names = 'power_density'
    functor_symbols = 'q'
    expression = '${dose_fraction} * q'
    block = ${salt}
  []
  # thermal neutron flux in the salt, shaped like the power (Xe-135 burnout)
  [thermal_flux]
    type = ParsedFunctorMaterial
    property_name = thermal_flux
    functor_names = 'power_density'
    functor_symbols = 'q'
    expression = '${phi_thermal} * q / ${fparse 150e6 * 0.9213 / 1.70462}'
    block = ${salt}
  []
  [no_removal]
    type = GenericFunctorMaterial
    prop_names = 'zero_rate'
    prop_values = '0'
    block = ${salt}
  []
[]

[MoltenSaltRadiolysis]
  block = ${salt}
  salt_type = fluoride
  metals = 'U'
  temperature = T_fluid
  dose_rate = dose_rate
  # charge-consistent radiolysis: zero net radiolytic oxidation (see header and README)
  g_value_species = 'e_sol F_rad F2_diss'
  g_value_overrides = '0 0 0'
  diffusivity = D_eff
  time_derivative = false
  initial_condition_species = 'F_ion U_IV U_III'
  initial_condition_values = '${F_total} ${U_IV_in} ${U_III_in}'
  verbose = true
[]

[CorrosionPlatingFlow]
  block = ${salt}
  elements = 'Cr Fe Ni'
  reaction_boundary = 'wall_guide'
  # inlets: per-element loop equilibria, set below in [LinearFVBCs] (the action applies one
  # inlet value to every element, so the inlets are not created here)
  outlet_boundary = 'outlet'
  temperature = T_fluid
  reference_temperature = ${T_ref}
  material_class = hastelloy_n
  salt_class = fluoride_fuel
  redox_class = msre_or_fuel_baseline
  reference_concentration = 1.0
  temperature_dependent_kinetics = true
  equilibrium_concentrations = '${c_eq_Cr} ${c_eq_Fe} ${c_eq_Ni}'
  dissolution_enthalpies = '${dH_diss}'
  seed_fractions = '${f_Cr} ${f_Fe} ${f_Ni}'
  diffusivity = D_eff
  time_derivative = false
  initial_concentration = ${c_eq_Cr}
  verbose = true
[]

[FissionProductInventory]
  block = ${salt}
  nuclides = 'I135 Xe135'
  yield_set = U235_thermal
  power_density = power_density
  neutron_flux = thermal_flux
  removal_groups = 'noble_gas'
  removal_rates = 'zero_rate'
  removal_destinations = 'offgas'
  track_removed_inventories = false
  plating_boundaries = ${channel_walls}
  plating_mass_transfer_coefficient = ${k_m_graphite}
  diffusivity = D_eff
  time_derivative = false
  verbose = true
[]

[LinearFVKernels]
  # ---- fission oxidation: 0.8 equivalents per fission convert U(III) to U(IV)
  [fission_oxidation_UIII]
    type = LinearFVSource
    variable = U_III
    source_density = power_density
    scaling_factor = ${fparse -fission_eq_per_J}
    block = ${salt}
  []
  [fission_oxidation_UIV]
    type = LinearFVSource
    variable = U_IV
    source_density = power_density
    scaling_factor = ${fission_eq_per_J}
    block = ${salt}
  []
  # ---- momentum and pressure (re-solved from the converged restart; see the executioner)
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
[]

[LinearFVBCs]
  # ---- corrosion products entering the core: each element at its loop equilibrium at the
  # core-average salt temperature (the closed loop holds the salt there; see README)
  [in_c_Cr]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${inlets}
    variable = c_Cr
    functor = ${c_eq_Cr}
  []
  [in_c_Fe]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${inlets}
    variable = c_Fe
    functor = ${c_eq_Fe}
  []
  [in_c_Ni]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${inlets}
    variable = c_Ni
    functor = ${c_eq_Ni}
  []
  # ---- species entering the core: the loop-balance inlet values (loop_balance.py)
  [in_e_sol]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${inlets}'
    variable = e_sol
    functor = 0
  []
  [out_e_sol]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = e_sol
    use_two_term_expansion = false
  []
  [in_F_ion]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${inlets}'
    variable = F_ion
    functor = ${F_total}
  []
  [out_F_ion]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = F_ion
    use_two_term_expansion = false
  []
  [in_F_rad]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${inlets}'
    variable = F_rad
    functor = 0
  []
  [out_F_rad]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = F_rad
    use_two_term_expansion = false
  []
  [in_F2_diss]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${inlets}'
    variable = F2_diss
    functor = 0
  []
  [out_F2_diss]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = F2_diss
    use_two_term_expansion = false
  []
  [in_U_IV]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${inlets}'
    variable = U_IV
    functor = ${U_IV_in}
  []
  [out_U_IV]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = U_IV
    use_two_term_expansion = false
  []
  [in_U_III]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${inlets}'
    variable = U_III
    functor = ${U_III_in}
  []
  [out_U_III]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = U_III
    use_two_term_expansion = false
  []
  [in_I135]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${inlets}'
    variable = I135
    functor = ${I135_in}
  []
  [out_I135]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = I135
    use_two_term_expansion = false
  []
  [in_Xe135]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${inlets}'
    variable = Xe135
    functor = ${Xe135_in}
  []
  [out_Xe135]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = Xe135
    use_two_term_expansion = false
  []
  # ---- momentum and pressure (as in stage 3)
  [in_vel_x]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${inlets}'
    variable = vel_x
    functor = 0
  []
  [in_vel_y]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${inlets}'
    variable = vel_y
    functor = 0
  []
  [in_w_inner]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_inner'
    variable = vel_z
    functor = 1.7
  []
  [in_w_outer]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_outer'
    variable = vel_z
    functor = 1.03
  []
  [in_w_cr]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_cr'
    variable = vel_z
    functor = 1.74
  []
  [out_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'outlet'
    variable = pressure
    functor = 0
  []
  [out_vel_x]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = vel_x
    use_two_term_expansion = false
  []
  [wall_vel_x]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls}'
    variable = vel_x
    functor = 0
  []
  [out_vel_y]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = vel_y
    use_two_term_expansion = false
  []
  [wall_vel_y]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls}'
    variable = vel_y
    functor = 0
  []
  [out_vel_z]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = vel_z
    use_two_term_expansion = false
  []
  [wall_vel_z]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls}'
    variable = vel_z
    functor = 0
  []
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
[]

[Postprocessors]
  # ---- redox
  [U_III_out]
    type = SideAverageValue
    variable = U_III
    boundary = outlet
  []
  [U_IV_out]
    type = SideAverageValue
    variable = U_IV
    boundary = outlet
  []
  [U_III_fraction_out]
    type = ParsedPostprocessor
    expression = 'U_III_out / ${U_total}'
    pp_names = 'U_III_out'
  []
  # U(III) oxidized in one core pass (area-averaged outlet minus inlet)
  [dU_III_core]
    type = ParsedPostprocessor
    expression = '${U_III_in} - U_III_out'
    pp_names = 'U_III_out'
  []
  [U_total_out]
    type = ParsedPostprocessor
    expression = 'U_III_out + U_IV_out'
    pp_names = 'U_III_out U_IV_out'
  []
  [F2_max]
    type = ElementExtremeValue
    variable = F2_diss
    block = ${salt}
  []
  # ---- corrosion of the guide tubes
  [Cr_out]
    type = SideAverageValue
    variable = c_Cr
    boundary = outlet
  []
  [Fe_out]
    type = SideAverageValue
    variable = c_Fe
    boundary = outlet
  []
  [Ni_out]
    type = SideAverageValue
    variable = c_Ni
    boundary = outlet
  []
  # net metal exchange of the guide tubes (positive = dissolution into the salt), from the
  # concentration rise across the core: mdot_sector / rho * (c_out - c_in) [mol/s per sector]
  [Cr_release_sector]
    type = ParsedPostprocessor
    expression = '(Cr_out - ${c_eq_Cr}) * 2294.13 / 6 / ${rho}'
    pp_names = 'Cr_out'
  []
  # ---- xenon
  [Xe135_out]
    type = SideAverageValue
    variable = Xe135
    boundary = outlet
  []
  [Xe135_core_avg]
    type = ElementAverageValue
    variable = Xe135
    block = ${salt}
  []
  [I135_out]
    type = SideAverageValue
    variable = I135
    boundary = outlet
  []
  # ---- 0D loop-balance targets (loop_balance.py) for the verification
  [Xe135_loop_0d]
    type = Receiver
    default = ${Xe135_loop}
  []
  [I135_loop_0d]
    type = Receiver
    default = ${I135_loop}
  []
  [core_transit_0d]
    type = Receiver
    default = ${t_core}
  []
  [loop_time_0d]
    type = Receiver
    default = ${t_loop}
  []
[]

[Executioner]
  type = SIMPLE
  momentum_systems = 'u_system v_system w_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  energy_system = 'energy_system'
  passive_scalar_systems = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys U_IV_sys U_III_sys '
                           'c_Cr_sys c_Fe_sys c_Ni_sys I135_sys Xe135_sys'
  # Turbulence and temperature are the converged fields of stage 3 (frozen). Momentum and
  # pressure are re-solved from the converged restart: they stay converged (residual ~1e-6 at the
  # first iteration) and this builds the Rhie-Chow face fluxes that advect the species (with
  # should_solve_momentum = false the face flux is never initialized).
  should_solve_turbulence = false
  should_solve_energy = false
  passive_scalar_equation_relaxation = '1 1 1 1 1 1 1 1 1 1 1'
  passive_scalar_absolute_tolerance = '1e-10 1e-10 1e-10 1e-10 1e-10 1e-10 1e-10 1e-10 1e-10 1e-10 1e-10'
  passive_scalar_l_tol = 1e-10
  passive_scalar_l_abs_tol = 1e-14
  passive_scalar_l_max_its = 500
  passive_scalar_petsc_options_iname = '-pc_type -pc_hypre_type'
  passive_scalar_petsc_options_value = 'hypre boomeramg'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'
  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  energy_equation_relaxation = 1.0
  turbulence_equation_relaxation = '1 1'
  num_iterations = 60
  pressure_absolute_tolerance = 1e-30
  momentum_absolute_tolerance = 1e-30
  energy_absolute_tolerance = 1e-30
  turbulence_absolute_tolerance = '1e-30 1e-30'
  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = FINAL
[]
