# Pulse-radiolysis validation: solvated-electron decay by capture on Cr(II) and Cr(III) in molten
# LiCl-KCl. Reproduces the Iwamatsu et al. PCCP 2026 transient-absorption experiment (case
# cr_licl_kcl_iwamatsu_2026). A radiolytic pulse seeds e_sol and Cl_rad; the solvated electron is
# then consumed by e_sol + Cr_II -> Cr_I and/or e_sol + Cr_III -> Cr_II, so the (scale-free) e_sol
# concentration decays on the tens-of-nanoseconds timescale. The trace metal concentrations are
# overridable from the command line so the validation harness can sweep the measured traces.

Cr_II_initial = 0.0    # mol/m^3, set per Cr(II) trace
Cr_III_initial = 0.0   # mol/m^3, set per Cr(III) trace

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
  dose_rate = 0.0
  initial_condition_species = 'e_sol Cl_rad Cl_ion Cr_II Cr_III'
  initial_condition_values = '1.0e-3 1.0e-3 2.0e4 ${Cr_II_initial} ${Cr_III_initial}'
[]

[Postprocessors]
  [c_e_sol]
    type = ElementAverageValue
    variable = e_sol
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
                 'Cr_II_sys Cr_III_sys Cr_I_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = 2.1e-8
  dt = 1.0e-10
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
