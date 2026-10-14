# =====================================================================================================
# Training case 13: long-term corrosion of a Hastelloy N thermal-convection loop (DRIDN model).
#
# The Dynamic Redox-Inventory-Depletion Network (DRIDN) integrates the coupled Cr/Fe/Ni dissolution,
# the salt inventory, cold-leg capture (plating) and the depletion of the alloy surface over years
# of exposure. The physical context below is ORNL thermal-convection loop NCL-16 (ORNL-TM-4188,
# validation target M-014): purified fluoride fuel salt, 704 C hot leg, 166 C loop temperature
# span, 3.37 years. The measured salt chromium increase is ~+500 ppm.
#
# Units: MOOSE time is seconds; the model integrates in years internally.
# =====================================================================================================

years = 3.3652749258498744
n_steps = 40

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
    # Hastelloy N: Cr / Fe / Ni mass fractions and density [g/cm^3]
    mass_fractions = '0.07 0.05 0.71'
    density = 8.89
    cr_fraction_ratio = 1.0
    # Calibrated element kinetics and thermodynamic driving forces at the loop conditions
    log_exchange_offsets = '0.0 -1.0237419112804005 -2.3170915670911314'
    affinity_baseline = '18.112896721230936 8.835324825760875 0.3361703045970774'
    log_charge_base_no_redox = 0.398484197301191
    selectivity_scale = 0.004817896232653154
    initial_redox_shift = 0.0
    transient_redox = false
    # Loop: circulation, mass transfer and cold-leg capture of the dissolved metals
    flow_factor = 1.0
    mass_transfer_rate = 1004.9168083798145 # um/y
    cold_capture_fraction = '0.3143650162916868 0.665729274913501 0.7649127999143572'
    deposition_closure = fuel_like
    deposit_area_factor = 0.7064534573603768
    # Salt inventory: wetted area per salt mass of the loop and the dissolved-metal capacity
    inventory_scale = loop
    inventory_scale_factor = 1.0
    area_to_salt_mass = 7.767901069021746 # cm^2/g
    inventory_capacity = 839.1043065935307 # ppm
    inventory_coupling_factor = 1.0
    initial_dissolved_ppm = '0.9988232516525638 0.9988232516525638 0.9988232516525638'
    product_activity_floor = 0.9988232516525638
    # Alloy-side depletion and attack morphology
    cr_diffusion = 2.2516959233504542e-13 # cm^2/s
    mass_loss_fraction = 0.6752149812832213
    front_damage_multiplier = 1.0
    grain_boundary_length_multiplier = 0.5864505148098346
    stress_interfacial_activation = false
    fluoride_impurity_interfacial_activation = false
    chloride_salt = false
  []
[]

[Postprocessors]
  [time_y]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = elapsed_time_y
  []
  [salt_Cr_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = dissolved_inventory_cr_ppm
  []
  [salt_Fe_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = dissolved_inventory_fe_ppm
  []
  [salt_Ni_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = dissolved_inventory_ni_ppm
  []
  [corrosion_rate_um_y]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = instantaneous_front_rate_um_y
  []
  [front_depth_um]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = front_depth_um
  []
  [igc_depth_um]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = igc_depth_um
  []
  [mass_loss_mg_cm2]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = mass_loss_mg_cm2
  []
  [mass_gain_mg_cm2]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = mass_gain_mg_cm2
  []
  [surface_Cr_availability]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = surface_availability_cr
  []
  [mass_balance_error]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = mass_balance_relative_error
  []
[]

[Executioner]
  type = Transient
  end_time = '${fparse years * 365.25 * 86400}'
  dt = '${fparse years * 365.25 * 86400 / n_steps}'
[]

[Outputs]
  csv = true
[]
