# =====================================================================================================
# Training case 03, model 3 of 4: long-term corrosion of loop NCL-16 with the DRIDN model.
#
# The Dynamic Redox-Inventory-Depletion Network (DRIDN) integrates the coupled Cr/Fe/Ni dissolution,
# the salt inventory, cold-leg capture (plating), a bulk-capture sink and the depletion of the alloy
# surface over years of exposure. The physical case is ORNL thermal-convection loop NCL-16
# (ORNL-FL-04, ORNL-TM-4188): standard and Ti-modified Hastelloy N in purified LiF-BeF2-UF4, 704 C
# hot leg, 538 C cold leg, 3.37 years (29500 h), before the FeF2 additions.
#
# This is a CALIBRATION (in-sample), not a validation: the NCL-16 measurements M-011..M-015 are direct
# targets of the DRIDN fit, and the parameters below are the fitted values. The model reproduces them
# (gold CSV, final time) as
#
#   quantity                  measured     DRIDN
#   hot-leg loss  (M-011)     2.9 mg/cm2   3.13  'mass_loss_mg_cm2'
#   cold-leg gain (M-012)     1.7 mg/cm2   1.61  'mass_gain_mg_cm2'
#   rate          (M-013)     1.016 um/y   1.045 'mean_rate_um_y' (time-averaged recession rate)
#   salt Cr       (M-014)     +500 ppm     +533  'salt_Cr_ppm'
#   salt Fe       (M-015)     -100 ppm     +106  'salt_Fe_ppm'   <-- WRONG SIGN, see below
#
# Known shortcomings (read before using the numbers):
#   * Fe: the measured salt Fe DECREASED by 100 ppm (Fe2+ reduced by Cr from the alloy). The model
#     starts the salt at ~1 ppm Fe, so it can only add Fe (+106 ppm). The calibration matched M-015
#     only through a separate post-hoc stoichiometric Fe-decrease formula, not through this state.
#   * Ni: 378 ppm dissolved Ni is a calibration artefact (no Ni target exists). Ni is the noblest
#     component, and purified fluoride loops keep it near the analytic floor (tens of ppm at most).
#   * Mass ledger: of the metal released by the hot leg (released_*_ppm), most goes to an
#     empirical bulk-capture sink (bulk_captured_*_ppm: 72 % of the Cr, 64 % of all metal), the cold
#     leg receives 20 % of the Cr and only 8 % stays dissolved. Loss minus gain is therefore NOT the
#     dissolved inventory.
#   * The Cr plateau (~533 ppm after ~1 y) is set by inventory_capacity and the bulk-capture closure;
#     it persists with the cold-leg capture switched off (529 ppm).
#   * area_to_salt_mass = 7.77 cm^2/g is a shared FITTED value used for every calibration loop; a
#     ~1.7 cm ID tube full of ~2.1 g/cm^3 salt has 4 / (D rho) ~ 1.1 cm^2/g.
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
    # Design data -- Hastelloy N: Cr / Fe / Ni mass fractions and density [g/cm^3]
    mass_fractions = '0.07 0.05 0.71'
    density = 8.89
    cr_fraction_ratio = 1.0
    # CALIBRATED: element kinetics and thermodynamic driving forces at the loop conditions
    log_exchange_offsets = '0.0 -1.0237419112804005 -2.3170915670911314'
    affinity_baseline = '18.112896721230936 8.835324825760875 0.3361703045970774'
    log_charge_base_no_redox = 0.398484197301191
    selectivity_scale = 0.004817896232653154
    initial_redox_shift = 0.0
    transient_redox = false
    # CALIBRATED (except flow_factor): mass transfer and cold-leg capture of the dissolved metals
    flow_factor = 1.0
    mass_transfer_rate = 1004.9168083798145 # um/y
    cold_capture_fraction = '0.3143650162916868 0.665729274913501 0.7649127999143572'
    deposition_closure = fuel_like
    deposit_area_factor = 0.7064534573603768
    # CALIBRATED: salt inventory -- the wetted area per salt mass (a shared fitted value, not the
    # loop geometry), the dissolved-metal capacity and the activity floor
    inventory_scale = loop
    inventory_scale_factor = 1.0
    area_to_salt_mass = 7.767901069021746 # cm^2/g
    inventory_capacity = 839.1043065935307 # ppm
    inventory_coupling_factor = 1.0
    initial_dissolved_ppm = '0.9988232516525638 0.9988232516525638 0.9988232516525638'
    product_activity_floor = 0.9988232516525638
    # CALIBRATED: alloy-side depletion and attack morphology
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
  # Calibrated quantity: exposure-averaged recession rate (mass recession / time), cf. M-013
  [mean_rate_um_y]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = average_corrosion_rate_um_y
  []
  # Instantaneous attack-front rate (front advance, before the mass-loss fraction)
  [front_rate_um_y]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = instantaneous_front_rate_um_y
  []
  # Cumulative Cr released by the hot leg and the part removed by the bulk-capture closure [ppm]
  [released_Cr_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = cumulative_source_cr_ppm
  []
  [bulk_captured_Cr_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = bulk_captured_cr_ppm
  []
  # Cumulative Fe released by the hot leg and the part removed by the bulk-capture closure [ppm]
  [released_Fe_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = cumulative_source_fe_ppm
  []
  [bulk_captured_Fe_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = bulk_captured_fe_ppm
  []
  # Cumulative Ni released by the hot leg and the part removed by the bulk-capture closure [ppm]
  [released_Ni_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = cumulative_source_ni_ppm
  []
  [bulk_captured_Ni_ppm]
    type = AdvancedCorrosionModelPostprocessor
    model = dridn
    quantity = bulk_captured_ni_ppm
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
