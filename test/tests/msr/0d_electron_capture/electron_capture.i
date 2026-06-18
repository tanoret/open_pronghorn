# 0D verification of the metal electron-capture chemistry against
# test_pseudo_first_order_electron_capture. With e_sol consumed only by e_sol + Zn_II -> Zn_I
# (all chloride radical channels start at zero and Cl_rad is never produced), and Zn_II held nearly
# constant at 100 mol/m^3, the solvated electron decays as e(t) = e0 exp(-k [Zn_II] t) with
# k(673.15 K) = 2.4e10 exp(-3.56e4/(R T)) = 4.1476e7. At t = 2e-10 s, e ~ 0.436 (pseudo-first order);
# the fully coupled result is slightly higher because Zn_II is depleted a little.

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
                     'Zn_II_sys Zn_I_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  metals = 'Zn'
  temperature = 673.15
  dose_rate = 0.0
  initial_condition_species = 'e_sol Zn_II Cl_ion'
  initial_condition_values = '1.0 100.0 5000.0'
[]

[Postprocessors]
  [c_e_sol]
    type = ElementAverageValue
    variable = e_sol
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Zn_I]
    type = ElementAverageValue
    variable = Zn_I
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 20
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys '
                 'Zn_II_sys Zn_I_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = 2.0e-10
  dt = 5.0e-12
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
