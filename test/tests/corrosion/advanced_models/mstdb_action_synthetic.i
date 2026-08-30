# Action-level smoke test for the static MSTDB-TC kinetics seed. The thermodynamic records are
# compact synthetic parser fixtures, not MSTDB-TC data and not scientific validation evidence.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 2
  []
[]

[Problem]
  solve = false
[]

[CorrosionPlating]
  topology = salt_only
  reaction_boundary = left
  kinetics_model = mstdb_tc_standard_state
  database = ../../../../data/corrosion_database.json
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
  transient = true
  verbose = true
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
[]
