# Tier-1 conservation verification. A closed, well-mixed cell (no flow, no in/out-flux) of LiCl-KCl
# with dissolved Cr is irradiated, so radiolysis continuously cycles the chromium oxidation state
# (Cl2m_rad oxidizes Cr(II)->Cr(III); e_sol reduces it). No reaction creates or destroys chromium,
# so the total chromium Cr_II + Cr_III + Cr_I must stay equal to its initial value. With the
# segregated reaction coupling iterated to convergence each step (multi_system_fixed_point with a
# residual/solution-change Convergence), the lagged-linearized scheme conserves total chromium to a
# tight tolerance. The 'total_Cr' postprocessor is the conservation check.

n_fp = 50   # fixed-point sweeps per step (overridden in the convergence study)

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
    xmax = 1.0
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys '
                     'Cr_II_sys Cr_III_sys Cr_I_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  metals = 'Cr'
  temperature = 673.15
  dose_rate = 1.0e6
  initial_condition_species = 'Cl_ion Cr_II'
  initial_condition_values = '2.0e4 1.0'
[]

[Postprocessors]
  [c_Cr_II]
    type = ElementAverageValue
    variable = Cr_II
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Cr_III]
    type = ElementAverageValue
    variable = Cr_III
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Cr_I]
    type = ElementAverageValue
    variable = Cr_I
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_Cr]
    type = ParsedPostprocessor
    expression = 'c_Cr_II + c_Cr_III + c_Cr_I'
    pp_names = 'c_Cr_II c_Cr_III c_Cr_I'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Absolute conservation error |total_Cr - 1| (initial total chromium = 1 mol/m^3).
  [conservation_error]
    type = ParsedPostprocessor
    expression = 'abs(total_Cr - 1.0)'
    pp_names = 'total_Cr'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = ${n_fp}
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys '
                 'Cr_II_sys Cr_III_sys Cr_I_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = 1.0
  dt = 0.1
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
