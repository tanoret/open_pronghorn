# Fluoride F2-yield validation against Davis et al. NSE 2022, Table III (case eflibe_f2_yield_davis_2022).
# Davis reports the radiolytic F2 production yield for FLiBe-UF4 as G(F2) = 0.005 +/- 0.001
# molecules per 100 eV. Here the direct radiolytic F2 source is isolated (the e_sol and F_rad yields
# are overridden to zero, so dissolved F2 grows only from its own G value) and accumulates linearly;
# the validation recovers the effective G value from the modeled production rate and checks it
# against the measured FLiBe-UF4 yield. This ties the bundled fluoride F2 G value to the literature
# measurement and verifies the source-term implementation.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
    xmax = 1.0
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = fluoride
  temperature = 333.15        # Davis irradiation temperature ~60 C
  dose_rate = 1.0e6           # J/m^3/s (magnitude cancels in the recovered G value)
  g_value_species = 'e_sol F_rad'
  g_value_overrides = '0.0 0.0'
[]

[Postprocessors]
  [c_F2_diss]
    type = ElementAverageValue
    variable = F2_diss
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys'
  scheme = implicit-euler
  start_time = 0.0
  end_time = 10.0
  dt = 1.0
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
