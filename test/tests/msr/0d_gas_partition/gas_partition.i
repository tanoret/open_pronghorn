# 0D verification of MSRGasExchange against test_gas_partition_equilibrium. The dissolved Cl2 has no
# chemical source (all chloride radical channels start at zero), so it only exchanges with the
# headspace. Carrying the gas inventory as a concentration C_gas = n_gas / V_gas, with all of it
# initially in the gas (C_gas = 1, C_liq = 0) and V_liq = V_gas, total-mole conservation
# C_liq + C_gas = 1 and equilibrium C_liq = kH R T C_gas give
#   C_liq_eq = kH R T / (1 + kH R T) = 0.0529994 mol/m^3   (kH = 1e-5, T = 673.15).

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
                     'Cl2_gas_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  gas_species = 'Cl2'
  temperature = 673.15
  dose_rate = 0.0
  kLa = 10.0
  volume_ratio = 1.0
  initial_condition_species = 'Cl2_gas'
  initial_condition_values = '1.0'
[]

[Postprocessors]
  [c_Cl2_diss]
    type = ElementAverageValue
    variable = Cl2_diss
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Cl2_gas]
    type = ElementAverageValue
    variable = Cl2_gas
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
                 'Cl2_gas_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = 10.0
  dt = 0.25
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
