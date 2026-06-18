# Flow-coupled demonstration: chloride radiolysis species transported through a 1D channel by a
# prescribed constant velocity, with molecular diffusion, a uniform radiolytic dose source and the
# full mass-action chemistry. Clean salt enters at the left (all species zero); the dose generates
# radicals that react and accumulate (notably dissolved Cl2) as the fluid moves downstream. This
# exercises the species advection + diffusion + source + reaction kernels created by the action.

velocity = '0.1 0 0'

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 20
    xmax = 1.0
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  temperature = 673.15
  dose_rate = 1.0e6
  velocity = ${velocity}
  advected_interp_method = upwind
  diffusivity = 1.0e-3
[]

[LinearFVBCs]
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
    type = PointValue
    variable = Cl2_diss
    point = '0.975 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 5
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = 20.0
  dt = 2.0
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
