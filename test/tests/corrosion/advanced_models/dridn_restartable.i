# Short deterministic DRIDN wrapper test. The duration is intentionally small because this checks
# MOOSE lifecycle, restartable state, units, and bookkeeping rather than long-exposure validation.

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
  [dridn]
    type = DRIDNCorrosionUserObject
    advanced_database = advanced_models_synthetic.json
    mass_fractions = '0.17 0.66 0.12'
    log_exchange_offsets = '0 -0.5 -1'
    affinity_baseline = '10 5 0'
    cold_capture_fraction = '0.2 0.4 0.6'
    initial_dissolved_ppm = '1 1 1'
    cr_fraction_ratio = 2
    density = 8
    flow_factor = 1.21
    selectivity_scale = 0.02
    product_activity_floor = 1
    initial_redox_shift = 0.5
    log_charge_base_no_redox = -3
    mass_transfer_rate = 100
    inventory_capacity = 1000
    area_to_salt_mass = 0.25
    inventory_scale_factor = 1
    inventory_coupling_factor = 0.5
    deposit_area_factor = 1.5
    mass_loss_fraction = 0.4
    cr_diffusion = 1e-14
    front_damage_multiplier = 1.2
    grain_boundary_length_multiplier = 0.8
    inventory_scale = explicit
    deposition_closure = flinak
    transient_redox = true
    stress_interfacial_activation = false
    fluoride_impurity_interfacial_activation = false
    chloride_salt = false
    charge_transfer_mode = legacy_irreversible
  []
[]

[Postprocessors]
  [./accepted_internal_steps]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = accepted_internal_steps
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./average_corrosion_rate_um_y]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = average_corrosion_rate_um_y
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./dissolved_inventory_cr_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = dissolved_inventory_cr_ppm
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./dissolved_inventory_fe_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = dissolved_inventory_fe_ppm
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./dissolved_inventory_ni_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = dissolved_inventory_ni_ppm
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./elapsed_time_y]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = elapsed_time_y
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./front_depth_um]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = front_depth_um
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./igc_depth_um]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = igc_depth_um
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./instantaneous_front_rate_um_y]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = instantaneous_front_rate_um_y
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mass_balance_relative_error]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = mass_balance_relative_error
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mass_gain_mg_cm2]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = mass_gain_mg_cm2
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mass_loss_mg_cm2]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = mass_loss_mg_cm2
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mass_recession_um]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = mass_recession_um
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./redox_shift]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = redox_shift
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./rejected_internal_steps]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = rejected_internal_steps
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./surface_availability_cr]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = surface_availability_cr
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./surface_availability_fe]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = surface_availability_fe
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./surface_availability_ni]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = surface_availability_ni
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 2
  dt = 157788 # 0.005 Julian years
[]

[Outputs]
  csv = true
  execute_on = 'INITIAL TIMESTEP_END'
[]
