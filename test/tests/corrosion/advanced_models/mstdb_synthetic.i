# This fixture exercises the MOOSE wrappers and strict file provenance with compact, synthetic
# ChemSage records. The records are deliberately not MSTDB-TC thermodynamic data and must not be
# used for scientific validation.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  solve = false
[]

[UserObjects]
  [mstdb_endpoint]
    type = MSTDBTCCorrosionUserObject
    corrosion_database = corrosion_database.json
    advanced_database = advanced_models_synthetic.json
    fluoride_database = MSTDB-TC_V3.1_Fluorides_No_Func.dat
    chloride_database = MSTDB-TC_V3.1_Chlorides_No_Func.dat
    hot_temperature = 923.15
    cold_temperature = 823.15
    exposure_time = 31557600
    flow_factor = 1.0
    area_to_salt_mass = 0.25
    inventory_coupling_factor = 0.5
    material_class = hastelloy_n
    salt_class = fluoride_fuel
    redox_class = purified_baseline
    execute_on = INITIAL
  []
[]

[Postprocessors]
  [front_rate_um_y]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = front_rate_um_y
    execute_on = INITIAL
  []
  [corrosion_rate_um_y]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = corrosion_rate_um_y
    execute_on = INITIAL
  []
  [front_depth_um]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = front_depth_um
    execute_on = INITIAL
  []
  [mass_loss_mg_cm2]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = mass_loss_mg_cm2
    execute_on = INITIAL
  []
  [mass_gain_mg_cm2]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = mass_gain_mg_cm2
    execute_on = INITIAL
  []
  [igc_depth_um]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = igc_depth_um
    execute_on = INITIAL
  []
  [source_fraction_cr]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = source_fraction_cr
    execute_on = INITIAL
  []
  [source_fraction_fe]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = source_fraction_fe
    execute_on = INITIAL
  []
  [source_fraction_ni]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = source_fraction_ni
    execute_on = INITIAL
  []
  [dissolved_inventory_cr_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = dissolved_inventory_cr_ppm
    execute_on = INITIAL
  []
  [dissolved_inventory_fe_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = dissolved_inventory_fe_ppm
    execute_on = INITIAL
  []
  [dissolved_inventory_ni_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = mstdb_endpoint
    quantity = dissolved_inventory_ni_ppm
    execute_on = INITIAL
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
[]

[Outputs]
  csv = true
  execute_on = INITIAL
[]
