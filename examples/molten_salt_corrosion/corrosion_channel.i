# =====================================================================================================
# Realistic example: radiation-driven chromium corrosion in a flowing molten chloride salt.
#
# A section of a molten LiCl-KCl coolant channel carries dissolved Cr(II) (a corrosion product shed by
# structural alloys). The salt flows past a localized, intense gamma field (peaked near the reactor
# core, modeled by a Gaussian dose profile). Radiolysis there generates the chlorine radical chain
# (Cl_rad -> Cl2m_rad) whose oxidant Cl2m_rad oxidizes Cr(II) -> Cr(III), while solvated electrons
# partially reduce the chromium back. The steady solution is the spatial map of the chromium
# oxidation state, i.e. where the salt becomes more oxidizing/corrosive.
#
# Physics demonstrated:
#   * laminar molten-salt flow (linear FV segregated SIMPLE + Rhie-Chow),
#   * radiolysis species transported as passive scalars (advection by the solved flow + diffusion),
#   * a spatially varying radiolytic dose source,
#   * the full chloride + Cr mass-action chemistry,
#   * steady-state (dose-driven) operation (time_derivative = false, compatible with SIMPLE).
#
# Run:   open_pronghorn-opt -i corrosion_channel.i
# View:  the Exodus output shows vel_x, Cr_II, Cr_III and Cl2m_rad fields.
# =====================================================================================================

# --- molten LiCl-KCl properties (approximate, ~500 C) ---
rho = 1800.0          # kg/m^3
mu = 0.0025           # Pa s
T = 773.15            # K
D_species = 1.0e-9    # m^2/s molecular diffusivity
inlet_velocity = 0.02 # m/s

# --- chemistry inlet state ---
Cl_bulk = 2.0e4       # mol/m^3 chloride (eutectic)
Cr_inlet = 1.0        # mol/m^3 dissolved Cr(II) corrosion product

# --- radiation field ---
dose_peak = 1.0e7     # J/m^3/s peak volumetric dose rate near the core
dose_center = 0.1     # m, axial location of the dose peak
dose_width = 0.025    # m, Gaussian half-width of the dose region

advected_interp_method = 'upwind'

[Mesh]
  [channel]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0.0
    xmax = 0.2
    ymin = 0.0
    ymax = 0.02
    nx = 40
    ny = 8
  []
[]

[Functions]
  [dose_profile]
    type = ParsedFunction
    expression = 'dose_peak * exp(-((x - dose_center)^2) / (2.0 * dose_width^2))'
    symbol_names = 'dose_peak dose_center dose_width'
    symbol_values = '${dose_peak} ${dose_center} ${dose_width}'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system e_sol_sys Cl_ion_sys Cl_rad_sys '
                     'Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys Cr_II_sys Cr_III_sys Cr_I_sys'
  previous_nl_solution_required = true
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
    initial_condition = ${inlet_velocity}
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    initial_condition = 0.0
    solver_sys = v_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    initial_condition = 0.0
    solver_sys = pressure_system
  []
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  metals = 'Cr'
  temperature = ${T}
  dose_rate = dose_profile          # spatially varying near-core dose field
  rhie_chow_user_object = 'rc'      # advection by the solved flow
  advected_interp_method = ${advected_interp_method}
  diffusivity = ${D_species}
  time_derivative = false           # steady-state corrosion map (SIMPLE)
  initial_condition_species = 'Cl_ion Cr_II'
  initial_condition_values = '${Cl_bulk} ${Cr_inlet}'
[]

[LinearFVKernels]
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    mu = ${mu}
    u = vel_x
    v = vel_y
    momentum_component = 'x'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
  []
  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    mu = ${mu}
    u = vel_x
    v = vel_y
    momentum_component = 'y'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
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
[]

[LinearFVBCs]
  # --- flow ---
  [inlet_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = left
    variable = vel_x
    functor = ${inlet_velocity}
  []
  [inlet_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = left
    variable = vel_y
    functor = 0.0
  []
  [walls_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom'
    variable = vel_x
    functor = 0.0
  []
  [walls_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom'
    variable = vel_y
    functor = 0.0
  []
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = right
    variable = pressure
    functor = 0.0
  []
  [outlet_u]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x
    boundary = right
    use_two_term_expansion = false
  []
  [outlet_v]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y
    boundary = right
    use_two_term_expansion = false
  []

  # --- species: bulk salt with dissolved Cr(II) enters; everything leaves through the outlet ---
  [inlet_e_sol]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = e_sol
    boundary = left
    functor = 0
  []
  [inlet_Cl_ion]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl_ion
    boundary = left
    functor = ${Cl_bulk}
  []
  [inlet_Cl_rad]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl_rad
    boundary = left
    functor = 0
  []
  [inlet_Cl2m_rad]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl2m_rad
    boundary = left
    functor = 0
  []
  [inlet_Cl3_ion]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl3_ion
    boundary = left
    functor = 0
  []
  [inlet_Cl2_diss]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl2_diss
    boundary = left
    functor = 0
  []
  [inlet_Cr_II]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cr_II
    boundary = left
    functor = ${Cr_inlet}
  []
  [inlet_Cr_III]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cr_III
    boundary = left
    functor = 0
  []
  [inlet_Cr_I]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cr_I
    boundary = left
    functor = 0
  []
  [outlet_e_sol]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = e_sol
    boundary = right
    use_two_term_expansion = false
  []
  [outlet_Cl_ion]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl_ion
    boundary = right
    use_two_term_expansion = false
  []
  [outlet_Cl_rad]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl_rad
    boundary = right
    use_two_term_expansion = false
  []
  [outlet_Cl2m_rad]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl2m_rad
    boundary = right
    use_two_term_expansion = false
  []
  [outlet_Cl3_ion]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl3_ion
    boundary = right
    use_two_term_expansion = false
  []
  [outlet_Cl2_diss]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl2_diss
    boundary = right
    use_two_term_expansion = false
  []
  [outlet_Cr_II]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cr_II
    boundary = right
    use_two_term_expansion = false
  []
  [outlet_Cr_III]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cr_III
    boundary = right
    use_two_term_expansion = false
  []
  [outlet_Cr_I]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cr_I
    boundary = right
    use_two_term_expansion = false
  []
[]

[Postprocessors]
  [Cr_III_outlet]
    type = SideAverageValue
    variable = Cr_III
    boundary = right
  []
  [Cr_II_outlet]
    type = SideAverageValue
    variable = Cr_II
    boundary = right
  []
  [Cr_I_outlet]
    type = SideAverageValue
    variable = Cr_I
    boundary = right
  []
  [Cl2_diss_max]
    type = ElementExtremeValue
    variable = Cl2_diss
    value_type = max
  []
  # Total chromium leaving the channel; should match the 1 mol/m^3 entering (conservation check).
  [Cr_total_outlet]
    type = ParsedPostprocessor
    expression = 'Cr_II_outlet + Cr_III_outlet + Cr_I_outlet'
    pp_names = 'Cr_II_outlet Cr_III_outlet Cr_I_outlet'
  []
  # Oxidized fraction Cr(III) / Cr_total at the outlet: the corrosion-state indicator.
  [oxidized_fraction_outlet]
    type = ParsedPostprocessor
    expression = 'Cr_III_outlet / (Cr_total_outlet + 1e-30)'
    pp_names = 'Cr_III_outlet Cr_total_outlet'
  []
[]

[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  passive_scalar_systems = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys '
                           'Cr_II_sys Cr_III_sys Cr_I_sys'
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  passive_scalar_equation_relaxation = '0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9'
  # The stiff radical chemistry makes the segregated coupling converge slowly, so the species
  # systems are iterated to a tight absolute tolerance (Tier-1 convergence control). This keeps
  # total chromium conserved to ~0.1 %; see test/tests/msr/conservation for the conservation study.
  num_iterations = 2000
  pressure_absolute_tolerance = 1e-9
  momentum_absolute_tolerance = 1e-9
  passive_scalar_absolute_tolerance = '1e-11 1e-11 1e-11 1e-11 1e-11 1e-11 1e-11 1e-11 1e-11'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  passive_scalar_petsc_options_iname = '-pc_type -pc_hypre_type'
  passive_scalar_petsc_options_value = 'hypre boomeramg'
  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  exodus = true
  csv = true
[]
