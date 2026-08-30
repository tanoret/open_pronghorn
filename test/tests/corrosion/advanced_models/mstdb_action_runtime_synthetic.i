# Runtime regression for the finite-element MSTDB-TC Action path. The thermodynamic records are
# compact synthetic parser fixtures, not MSTDB-TC data and not scientific validation evidence.
#
# The salt starts at c_ref for Cr, Fe, and Ni. One short implicit step executes all three generated
# Butler-Volmer boundaries. Conservation converts each domain-average concentration increment back
# to an equivalent penetration rate, so this test detects a total-rate seed copied to every element
# instead of the calibrated source-resolved shares. The gold values are the one-EDGE2, one-step
# solution derived from the frozen synthetic C++ endpoint and production Butler-Volmer equation.
# The endpoint plans 0.44786428456393895, 0.11329080418926715, and 0.44141692471334704
# um/y (sum 1.002572013466553). The one-second implicit result differs by less than 1 ppm because
# the generated boundaries begin loading the salt immediately after the exactly c=c_ref initial state.
# The Action-generated corrosion_rate_um_y must report that full initial sum, and recession_um must
# integrate the same summed runtime rate rather than only the default chromium share.

L = 1.0e-3
rho = 8.89
seconds_per_year = 31557600

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
    xmax = ${L}
  []
[]

[CorrosionPlating]
  topology = salt_only
  reaction_boundary = left
  kinetics_model = mstdb_tc_standard_state
  database = corrosion_database.json
  advanced_database = advanced_models_synthetic.json
  fluoride_database = synthetic_fluoride_chemsage_fixture.dat
  chloride_database = synthetic_chloride_chemsage_fixture.dat
  elements = 'Cr Fe Ni'
  temperature = 923.15
  reference_temperature = 923.15
  reference_cold_temperature = 823.15
  reference_exposure_time = 31557600
  area_to_salt_mass = 0.25
  inventory_coupling_factor = 0.5
  material_class = hastelloy_n
  salt_class = fluoride_fuel
  redox_class = purified_baseline
  flow_factor = 1.0
  applied_overpotential = 0.1
  default_salt_concentration = 1.0
  reference_concentration = 1.0
  transient = true
[]

[Postprocessors]
  [c_cr]
    type = ElementAverageValue
    variable = c_Cr
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_fe]
    type = ElementAverageValue
    variable = c_Fe
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_ni]
    type = ElementAverageValue
    variable = c_Ni
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [dc_cr_dt]
    type = ChangeOverTimePostprocessor
    postprocessor = c_cr
    divide_by_dt = true
    execute_on = TIMESTEP_END
  []
  [dc_fe_dt]
    type = ChangeOverTimePostprocessor
    postprocessor = c_fe
    divide_by_dt = true
    execute_on = TIMESTEP_END
  []
  [dc_ni_dt]
    type = ChangeOverTimePostprocessor
    postprocessor = c_ni
    divide_by_dt = true
    execute_on = TIMESTEP_END
  []

  # Equivalent generated boundary currents, i = z F L dc_avg/dt [A/m^2].
  [current_cr]
    type = ParsedPostprocessor
    pp_names = dc_cr_dt
    expression = 'dc_cr_dt * ${L} * 2 * 96485.33212'
    execute_on = TIMESTEP_END
  []
  [current_fe]
    type = ParsedPostprocessor
    pp_names = dc_fe_dt
    expression = 'dc_fe_dt * ${L} * 2 * 96485.33212'
    execute_on = TIMESTEP_END
  []
  [current_ni]
    type = ParsedPostprocessor
    pp_names = dc_ni_dt
    expression = 'dc_ni_dt * ${L} * 2 * 96485.33212'
    execute_on = TIMESTEP_END
  []

  # J = L dc_avg/dt for this one-dimensional unit-area control volume. Faradaic mass
  # conversion gives rate = J M seconds_per_year / rho.
  [rate_cr]
    type = ParsedPostprocessor
    pp_names = dc_cr_dt
    expression = 'dc_cr_dt * ${L} * 51.9961 / ${rho} * ${seconds_per_year}'
    execute_on = TIMESTEP_END
  []
  [rate_fe]
    type = ParsedPostprocessor
    pp_names = dc_fe_dt
    expression = 'dc_fe_dt * ${L} * 55.845 / ${rho} * ${seconds_per_year}'
    execute_on = TIMESTEP_END
  []
  [rate_ni]
    type = ParsedPostprocessor
    pp_names = dc_ni_dt
    expression = 'dc_ni_dt * ${L} * 58.6934 / ${rho} * ${seconds_per_year}'
    execute_on = TIMESTEP_END
  []
  [rate_sum]
    type = LinearCombinationPostprocessor
    pp_names = 'rate_cr rate_fe rate_ni'
    pp_coefs = '1 1 1'
    execute_on = TIMESTEP_END
  []
  [action_rate_um_y]
    type = SideAverageValue
    variable = corrosion_rate_um_y
    boundary = left
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [action_recession_um]
    type = SideAverageValue
    variable = recession_um
    boundary = left
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 1.0
  num_steps = 1
  automatic_scaling = true
  nl_abs_tol = 1.0e-14
  nl_rel_tol = 1.0e-12
[]

[Outputs]
  csv = true
  execute_on = 'INITIAL TIMESTEP_END'
[]
