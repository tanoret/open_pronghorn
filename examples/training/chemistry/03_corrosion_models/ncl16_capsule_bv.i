# =====================================================================================================
# Training case 03, model 1 of 4: the NCL-16 conditions in a CLOSED, ISOTHERMAL Butler-Volmer capsule.
#
# Hastelloy N in purified fuel salt at the NCL-16 hot-leg temperature (704 C), with the wall
# Butler-Volmer law ([CorrosionPlating], topology = salt_only). i0 is seeded from the
# calibrated correlation for the NCL-16 features (hastelloy_n / fluoride_fuel / purified_baseline,
# nominal, dT = 166 C): 1.382 um/y on fresh salt. The salt depth V/A = 3.9e-4 m is the calibrated
# area-to-salt-mass ratio of the loop models (7.77 cm^2/g at 3.3 g/cm^3; ncl16_dridn.i).
#
# The applied overpotential is chosen so that the capsule's Nernst limit c_eq = c_ref exp(z f eta)
# is the SAME Cr equilibrium as the hot wall of the loop model (ncl16_loop_flow.i): 900 ppm at 704 C,
#   eta = ln(c_eq / c_ref) / (z f) = 0.170 V   (c_ref = 1 mol/m^3, f = F / (R T)).
# Without a cold leg the dissolved Cr has nowhere to go: the capsule is half saturated after 31 d and
# 99 % after 192 d, and the wall stops corroding at d = c_eq (V/A) M_Cr / rho_alloy = 0.13 um. It
# cannot produce the measured years-long hot-leg loss, and it has no cold-leg gain at all.
#
# Gold (3.37 y): salt Cr 900 ppm (= c_eq, c_eq_ratio = 1), loss 0.116 mg/cm^2, mean rate 0.039 um/y,
# against 2.9 mg/cm^2 and 1.016 um/y measured. This is a MODEL-FORM check (the wrong model for a
# loop), not a calibration. 'faraday_ratio' = 1 is a VERIFICATION (recession =
# dissolved inventory).
# =====================================================================================================

T_hot = 977.15 # K (704 C)
rho_salt = 3.3 # g/cm^3, LiF-BeF2-ThF4-UF4 fuel salt at ~700 C (estimate)
M_Cr = 51.9961 # g/mol
rho_alloy = 8.89 # g/cm^3, hastelloy_n in data/corrosion_database.json
VA = 3.9e-4 # m, salt volume / wetted area (calibrated area_to_salt_mass of ncl16_dridn.i)
Cr_eq_ppm = 900 # ppm, the hot-wall Cr equilibrium of the loop model
c_eq = ${fparse Cr_eq_ppm * rho_salt / M_Cr} # mol/m^3
eta = ${fparse log(c_eq) * 8.314462618 * T_hot / (2 * 96485.33212)} # V, Nernst limit = c_eq
years = 3.3652749258498744 # 29500 h
year = 31557600 # s

[Mesh]
  [salt]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 10
    xmax = ${VA}
  []
[]

[CorrosionPlating]
  topology = salt_only
  reaction_boundary = left
  elements = 'Cr'
  temperature = ${T_hot}
  reference_temperature = ${T_hot}
  material_class = hastelloy_n
  salt_class = fluoride_fuel
  redox_class = purified_baseline
  position_class = nominal
  flow_factor = 1.0
  delta_T_C = 166.0
  applied_overpotential = ${eta}
  default_salt_concentration = 0.0 # fresh (purified) salt
  reference_concentration = 1.0
  transient = true
  verbose = true
[]

[Postprocessors]
  [time_y]
    type = ParsedPostprocessor
    expression = 't / ${year}'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [corrosion_rate_um_y]
    type = SideAverageValue
    variable = corrosion_rate_um_y
    boundary = left
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [recession_um]
    type = SideAverageValue
    variable = recession_um
    boundary = left
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Cr_mean]
    type = ElementAverageValue
    variable = c_Cr
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [salt_Cr_ppm]
    type = ParsedPostprocessor
    expression = 'c_Cr_mean * ${M_Cr} / ${rho_salt}'
    pp_names = 'c_Cr_mean'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_eq_ratio]
    type = ParsedPostprocessor
    expression = 'c_Cr_mean / ${c_eq}'
    pp_names = 'c_Cr_mean'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Hot-wall mass loss [mg/cm^2] = recession x rho_alloy, and the exposure-averaged rate [um/y]
  [mass_loss_mg_cm2]
    type = ParsedPostprocessor
    expression = 'recession_um * 1e-4 * ${rho_alloy} * 1e3'
    pp_names = 'recession_um'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [mean_rate_um_y]
    type = ParsedPostprocessor
    expression = 'if(time_y > 0, recession_um / time_y, 0)'
    pp_names = 'recession_um time_y'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [faraday_ratio]
    type = ParsedPostprocessor
    expression = 'if(c_Cr_mean > 1e-12, recession_um / (c_Cr_mean * ${VA} * ${M_Cr} / (${rho_alloy} * 1e6) * 1e6), 1)'
    pp_names = 'recession_um c_Cr_mean'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  end_time = ${fparse years * year}
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 600
    growth_factor = 1.5
  []
  dtmax = ${fparse 0.05 * year}
  nl_abs_tol = 1.0e-13
  nl_rel_tol = 1.0e-10
[]

[Outputs]
  csv = true
[]
