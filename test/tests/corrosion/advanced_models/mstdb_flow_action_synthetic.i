# Linear-FV Action smoke test for the static MSTDB-TC kinetics seed. The thermodynamic records are
# compact synthetic parser fixtures, not MSTDB-TC data and not scientific validation evidence.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  solve = false
  linear_sys_names = 'c_Cr_sys c_Fe_sys c_Ni_sys'
[]

[CorrosionPlatingFlow]
  kinetics_model = mstdb_tc_standard_state
  database = corrosion_database.json
  advanced_database = advanced_models_synthetic.json
  fluoride_database = MSTDB-TC_V3.1_Fluorides_No_Func.dat
  chloride_database = MSTDB-TC_V3.1_Chlorides_No_Func.dat
  elements = 'Cr Fe Ni'
  reaction_boundary = left
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
  reference_concentration = 1.0
  initial_concentration = 1.0
  temperature_dependent_kinetics = false
  diffusivity = 1.0e-9
  verbose = true
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
[]
