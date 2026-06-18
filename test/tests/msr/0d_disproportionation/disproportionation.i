# 0D verification of the MoltenSaltRadiolysis action against the analytic second-order
# disproportionation 2 Cl2m_rad -> Cl3_ion + Cl_ion. With no other species present, Cl2m_rad is
# consumed only by this reaction (nothing produces Cl_rad, so the reverse-forming reaction stays
# inactive), giving dC/dt = -2 k C^2 and C(t) = C0 / (1 + 2 k C0 t). The database rate at the
# reference temperature 673.15 K is k = 2.2e6, so C reaches C0/2 at t_half = 1/(2 k C0) = 2.2727e-7 s.
# This mirrors MoltenSaltRadiolysis tests.py::test_second_order_disproportionation.

k = 2.2e6
C0 = 1.0
t_half = '${fparse 1.0 / (2.0 * k * C0)}'

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
    xmax = 1.0
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  temperature = 673.15
  dose_rate = 0.0
  initial_condition_species = 'Cl2m_rad'
  initial_condition_values = '${C0}'
[]

[Postprocessors]
  [c_Cl2m_rad]
    type = ElementAverageValue
    variable = Cl2m_rad
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
  system_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = ${t_half}
  dt = '${fparse t_half / 200.0}'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
