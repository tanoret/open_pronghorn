# =====================================================================================================
# Training case 03: chromium redox in an irradiated chloride salt (LiCl-KCl, 400 C).
#
# Dissolved chromium is the main corrosion product of Ni-Cr alloys. Under gamma irradiation the
# radiolysis products cycle its oxidation state:
#   e_sol    + Cr(III) -> Cr(II)              e_sol    + Cr(II) -> Cr(I)
#   Cl2m_rad + Cr(II)  -> Cr(III) + 2 Cl-     Cl2m_rad + Cr(III) -> Cr(II) + Cl2
#   Cr(III)  + Cr(I)   -> 2 Cr(II)
# The Cr(III)/Cr(II) ratio sets the salt redox potential, i.e. how corrosive the salt is.
#
# Two exact invariants are monitored:
#   * chromium is neither created nor destroyed: Cr(I) + Cr(II) + Cr(III) = Cr_0     ('cr_error')
#   * radiolysis makes electrons and holes in equal numbers, and every reaction conserves charge, so
#     the reducing equivalents (e_sol, Cr(I) = Cr(II) + e) equal the oxidizing ones (Cl_rad, Cl2m_rad,
#     2 Cl3-, 2 Cl2, Cr(III) = Cr(II) - e). 'redox_balance' is their difference relative to the
#     electrons produced, G_e D t / (100 eV N_A); most of them recombine with Cl atoms.
#
# Result to discuss: the salt stays reduced (Cr(III)/Cr ~ 4e-6) because e_sol reduces Cr(III) and
# Cr(III) + Cr(I) comproportionate as fast as Cl2- oxidizes Cr(II), while Cl2 slowly accumulates.
# =====================================================================================================

dose_rate = 1.0e4 # W/m^3
Cr0 = 1.0 # mol/m^3 initially all Cr(II) (~ 50 wppm Cr in LiCl-KCl)
n_fp = 30 # fixed-point sweeps per time step

[Mesh]
  [cell]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
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
  dose_rate = ${dose_rate}
  initial_condition_species = 'Cl_ion Cr_II'
  initial_condition_values = '2.0e4 ${Cr0}'
[]

[Postprocessors]
  [c_e_sol]
    type = ElementAverageValue
    variable = e_sol
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Cl_rad]
    type = ElementAverageValue
    variable = Cl_rad
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Cl2m_rad]
    type = ElementAverageValue
    variable = Cl2m_rad
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Cl3_ion]
    type = ElementAverageValue
    variable = Cl3_ion
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Cl2_diss]
    type = ElementAverageValue
    variable = Cl2_diss
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Cr_I]
    type = ElementAverageValue
    variable = Cr_I
    execute_on = 'INITIAL TIMESTEP_END'
  []
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
  [Cr_III_fraction]
    type = ParsedPostprocessor
    expression = 'c_Cr_III / ${Cr0}'
    pp_names = 'c_Cr_III'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Invariant 1: chromium conservation
  [cr_error]
    type = ParsedPostprocessor
    expression = 'abs(c_Cr_I + c_Cr_II + c_Cr_III - ${Cr0}) / ${Cr0}'
    pp_names = 'c_Cr_I c_Cr_II c_Cr_III'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Invariant 2: electrons produced = holes produced. Electrons are carried by e_sol and by
  # Cr(I) (Cr(II) + e); holes by Cl_rad, Cl2m_rad (one each), Cl3- and Cl2 (two each) and Cr(III).
  [electrons_produced]
    type = ParsedPostprocessor
    expression = '0.3 * ${dose_rate} / (100 * 1.602176634e-19) / 6.02214076e23 * t'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [redox_balance]
    type = ParsedPostprocessor
    expression = 'if(electrons_produced > 0, '
                 'abs((c_e_sol + c_Cr_I) - (c_Cl_rad + c_Cl2m_rad + 2 * c_Cl3_ion + 2 * c_Cl2_diss '
                 '+ c_Cr_III)) / electrons_produced, 0)'
    pp_names = 'c_e_sol c_Cr_I c_Cl_rad c_Cl2m_rad c_Cl3_ion c_Cl2_diss c_Cr_III '
               'electrons_produced'
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
  end_time = 1e3
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-6
    growth_factor = 1.3
  []
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
