# =====================================================================================================
# Training case 01: radiolytic F2 in a fluoride salt -- the G value and its units.
#
# The simplest possible chemistry input: one well-mixed cell of fluoride salt under a constant gamma
# dose rate. Only the direct F2 yield is switched on (the e_sol and F_rad G values are overridden to
# zero), so dissolved F2 has a constant source and no sink:
#
#   dC/dt = S,   S = G * D / (100 eV * N_A)       [mol/m^3/s]
#
# with G in molecules per 100 eV, D the volumetric dose rate [J/m^3/s = W/m^3] and N_A Avogadro's
# number. D = 1e6 W/m^3 is the order of the gamma heating in an MSR core (a few percent of a
# ~1e8 W/m^3 fission power density). The 'analytic' and 'rel_error' postprocessors check C = S t.
# =====================================================================================================

dose_rate = 1.0e6 # W/m^3
G_F2 = 0.005 # molecules / 100 eV

[Mesh]
  [cell]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  # One linear system per species (named <species>_sys by the action)
  linear_sys_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = fluoride
  temperature = 873.15
  dose_rate = ${dose_rate}
  g_value_species = 'e_sol F_rad F2_diss'
  g_value_overrides = '0.0 0.0 ${G_F2}'
[]

[Postprocessors]
  [c_F2_diss]
    type = ElementAverageValue
    variable = F2_diss
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [analytic]
    type = ParsedPostprocessor
    expression = '${G_F2} * ${dose_rate} / (100 * 1.602176634e-19) / 6.02214076e23 * t'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [rel_error]
    type = ParsedPostprocessor
    expression = 'if(analytic > 0, abs(c_F2_diss - analytic) / analytic, 0)'
    pp_names = 'c_F2_diss analytic'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys'
  scheme = implicit-euler
  end_time = 3600 # one hour
  dt = 60
  # Tight linear tolerance: the concentrations are small numbers in mol/m^3
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
