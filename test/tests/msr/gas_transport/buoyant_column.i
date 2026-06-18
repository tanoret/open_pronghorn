# Verification of gas-phase buoyant transport and venting. A vertical column initially holds a
# uniform dissolved-gas-phase F2 inventory. With the gas-exchange turned off (kLa = 0) the gas
# phase only transports: it rises by buoyancy (an upward slip/advection velocity) and disperses,
# leaving through the vent at the top (outflow BC). The gas is therefore swept out of the domain -
# the bottom of the column clears first as the gas rises, and the total gas inventory decays toward
# zero on the residence-time scale L/v_buoy = 1.0/0.1 = 10 s.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 20
    xmax = 0.1
    ymax = 1.0
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys F2_gas_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = fluoride
  gas_species = 'F2'
  temperature = 673.15
  dose_rate = 0.0
  kLa = 0.0                          # gas exchange off: isolate gas-phase transport
  gas_buoyancy_velocity = '0 0.1 0'  # buoyant rise (upward)
  gas_dispersivity = 1.0e-4
  initial_condition_species = 'F2_gas'
  initial_condition_values = '1.0'
[]

[LinearFVBCs]
  [vent_F2_gas]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = F2_gas
    boundary = top
    use_two_term_expansion = false
  []
[]

[Postprocessors]
  [F2_gas_total]
    type = ElementAverageValue
    variable = F2_gas
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [F2_gas_bottom]
    type = PointValue
    variable = F2_gas
    point = '0.05 0.05 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [F2_gas_top]
    type = PointValue
    variable = F2_gas
    point = '0.05 0.95 0'
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
  system_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys F2_gas_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = 15.0
  dt = 0.5
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
