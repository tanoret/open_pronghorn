# Production transient flow-coupled radiolysis. A laminar channel flow is advanced in time with the
# linear finite-volume segregated PIMPLE algorithm (momentum + pressure + Rhie-Chow mass flux); the
# chloride radiolysis species are transported as passive scalars (advection by the solved velocity
# via Rhie-Chow, plus molecular diffusion) with a uniform radiolytic dose source and the full
# mass-action chemistry. The radiolysis action is run in transient mode (the default), so each
# species carries a time derivative consistent with the PIMPLE time integration.

mu = 2.6
rho = 1.0
advected_interp_method = 'upwind'

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '0.5 0.5'
    dy = '0.2'
    ix = '5 5'
    iy = '4'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system e_sol_sys Cl_ion_sys Cl_rad_sys '
                     'Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys'
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
    initial_condition = 0.5
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_condition = 0.0
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 0.2
  []
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  temperature = 673.15
  dose_rate = 1.0e6
  rhie_chow_user_object = 'rc'
  advected_interp_method = ${advected_interp_method}
  diffusivity = 1.0e-3
[]

[LinearFVKernels]
  [u_time]
    type = LinearFVTimeDerivative
    variable = vel_x
    factor = ${rho}
  []
  [v_time]
    type = LinearFVTimeDerivative
    variable = vel_y
    factor = ${rho}
  []
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
  [inlet_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = vel_x
    functor = '1.0'
  []
  [inlet_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = vel_y
    functor = '0.0'
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
    boundary = 'right'
    variable = pressure
    functor = 0.0
  []
  [outlet_u]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x
    use_two_term_expansion = false
    boundary = right
  []
  [outlet_v]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y
    use_two_term_expansion = false
    boundary = right
  []

  # Clean salt enters at the inlet; species leave through the outlet.
  [inlet_e_sol]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = e_sol
    boundary = left
    functor = 0
  []
  [outlet_e_sol]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = e_sol
    boundary = right
    use_two_term_expansion = false
  []
  [inlet_Cl_ion]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl_ion
    boundary = left
    functor = 0
  []
  [outlet_Cl_ion]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl_ion
    boundary = right
    use_two_term_expansion = false
  []
  [inlet_Cl_rad]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl_rad
    boundary = left
    functor = 0
  []
  [outlet_Cl_rad]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl_rad
    boundary = right
    use_two_term_expansion = false
  []
  [inlet_Cl2m_rad]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl2m_rad
    boundary = left
    functor = 0
  []
  [outlet_Cl2m_rad]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl2m_rad
    boundary = right
    use_two_term_expansion = false
  []
  [inlet_Cl3_ion]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl3_ion
    boundary = left
    functor = 0
  []
  [outlet_Cl3_ion]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl3_ion
    boundary = right
    use_two_term_expansion = false
  []
  [inlet_Cl2_diss]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl2_diss
    boundary = left
    functor = 0
  []
  [outlet_Cl2_diss]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl2_diss
    boundary = right
    use_two_term_expansion = false
  []
[]

[Postprocessors]
  [Cl2_diss_avg]
    type = ElementAverageValue
    variable = Cl2_diss
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cl2_diss_outlet]
    type = SideAverageValue
    variable = Cl2_diss
    boundary = right
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = PIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  passive_scalar_systems = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys'
  momentum_equation_relaxation = 0.8
  pressure_variable_relaxation = 0.3
  passive_scalar_equation_relaxation = '0.9 0.9 0.9 0.9 0.9 0.9'
  num_iterations = 50
  pressure_absolute_tolerance = 1e-9
  momentum_absolute_tolerance = 1e-9
  passive_scalar_absolute_tolerance = '1e-9 1e-9 1e-9 1e-9 1e-9 1e-9'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  passive_scalar_petsc_options_iname = '-pc_type -pc_hypre_type'
  passive_scalar_petsc_options_value = 'hypre boomeramg'
  print_fields = false
  continue_on_max_its = true
  dt = 0.5
  num_steps = 6
  num_piso_iterations = 0
[]

[Outputs]
  csv = true
[]
