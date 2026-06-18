# Pulse-radiolysis validation: solvated-electron capture by Zn(II) in molten LiCl-KCl, reproducing
# the Iwamatsu/Horne et al. PCCP 2022 measurement (case zn_licl_kcl_horne_2022). With [Zn_II] held
# nearly constant (>> e_sol), the solvated electron decays pseudo-first-order, e(t) = e0 exp(-k_obs t)
# with k_obs = k(T) [Zn_II]. The validation extracts k_obs from the modeled decay and compares it to
# the digitized Fig. 4B observation (data/vision_fig4B_Zn2_pseudo1st_order.csv) at 400 C and
# 9.41 mM Zn(II): k_obs = 4.05e8 1/s.

Zn_II_initial = 9.41   # mol/m^3 (9.41 mM ZnCl2 in LiCl-KCl)

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
  initial_condition_species = 'e_sol Cl_rad Cl_ion Zn_II'
  initial_condition_values = '1.0e-3 1.0e-3 2.0e4 ${Zn_II_initial}'
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
                 'Zn_II_sys Zn_I_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = 6.0e-9
  dt = 2.0e-11
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
