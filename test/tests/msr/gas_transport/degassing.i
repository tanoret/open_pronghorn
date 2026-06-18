# Verification of the MSRDegassingBC surface outgassing boundary condition. A single well-mixed cell
# of fluoride salt holds dissolved F2 with no radiolysis (dose = 0, so no radicals form and no
# reactions fire). The F2 escapes through the top free surface at the convective mass-transfer rate
# J = k_surf * (C - kH * p_cover) with a clean cover gas (p_cover = 0), so a mass balance over the
# cell (area A_top = 1, volume V = 1) gives the analytic decay
#   C(t) = C0 * exp(-(k_surf * A_top / V) t) = exp(-0.1 t)  ->  C(10) = exp(-1) = 0.367879.

k_surf = 0.1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 1
    xmax = 1.0
    ymax = 1.0
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = fluoride
  temperature = 673.15
  dose_rate = 0.0
  diffusivity = 1.0e-6
  initial_condition_species = 'F2_diss'
  initial_condition_values = '1.0'
[]

[LinearFVBCs]
  [degas_F2]
    type = MSRDegassingBC
    variable = F2_diss
    boundary = top
    mass_transfer_coefficient = ${k_surf}
    henry_coefficient = 0.0
    cover_gas_pressure = 0.0
  []
[]

[Postprocessors]
  [c_F2_diss]
    type = ElementAverageValue
    variable = F2_diss
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 10
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = 10.0
  dt = 0.05
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
