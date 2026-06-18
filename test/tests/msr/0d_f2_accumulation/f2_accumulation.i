# 0D verification of the radiolytic source against test_F2_yield_linear_accumulation.
# With only an F2_diss source active (G overridden so e_sol and F_rad have no yield, hence no F2
# production by F-radical recombination), the dissolved F2 grows linearly:
#   C(t) = S t,  S = G * dose_rate / (100 eV * N_A) = 0.005 * 1e6 / (100*1.602176634e-19) / 6.02214076e23.
# At t = 10 s, C = 5.18214...e-3 mol/m^3.

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
  temperature = 673.15
  dose_rate = 1.0e6
  g_value_species = 'e_sol F_rad F2_diss'
  g_value_overrides = '0.0 0.0 0.005'
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
