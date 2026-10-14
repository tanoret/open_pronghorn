# =====================================================================================================
# Training case 11: the chromium-depletion layer -- selective dissolution of Cr from 316 stainless.
#
# Case 10 treated the wall as an infinite Cr source. In reality Cr must diffuse through the alloy to
# the corroding surface: the near-surface alloy becomes depleted in Cr (in molten fluorides this is
# the voided, Cr-depleted zone seen in metallography of loop specimens). Here the alloy is resolved
# as a second block, and the Butler-Volmer interface couples the two phases:
#
#   salt  (0.5 m reservoir):  dc/dt  = div(D_salt grad c),     flux in  = + i / (z F)
#   alloy (300 um coupon):    dcs/dt = div(D_alloy grad cs),   flux out = - i / (z F)
#
#   i = i0 [ (cs/cs0) exp(alpha z f eta) - (c/c_ref) exp(-alpha z f eta) ]     ('solid_activity')
#
# The anodic branch carries the alloy-surface Cr activity cs/cs0, so the rate falls as the surface
# depletes and becomes limited by solid-state diffusion.
#
# Case: 316 stainless in fuel salt at 650 C (ORNL-FL-12 of validation/corrosion, 21.6 um/y with the
# alloy at full activity). The alloy Cr diffusivity is an effective 1e-16 m^2/s lumping lattice and
# grain-boundary transport (lattice D_Cr ~ 1e-19 m^2/s at 650 C; grain boundaries dominate).
#
# Results (1 year): the rate drops from 23 to ~3 um/y within months, the surface Cr falls to ~48 %
# of the bulk and the depletion zone (cs < 95 % cs0, 'alloy_profile') reaches ~125 um ~ sqrt(4 D t).
#
# Checks:
#   * total Cr (salt + alloy) is conserved: 'cr_balance_error' ~ 1e-13
#   * the depletion layer grows as sqrt(D t): the integral depth 'depletion_depth_um' (Cr deficit /
#     surface deficit, 60 um after 1 year) over sqrt(D t) is 'depth_over_sqrt_Dt' ~ 1.06
#     (2/sqrt(pi) = 1.13 for a fixed surface concentration)
#   * the Faradaic recession equals the Cr lost by the alloy: 'faraday_ratio' = 1
#
# Exercise: L_salt=0.05 (a small capsule): the salt reaches the Nernst limit of case 10 within
# weeks, dissolution stops and Cr diffusing from the bulk *refills* the surface.
# =====================================================================================================

T = 923.15 # K (650 C)
D_alloy = 1.0e-16 # m^2/s effective Cr diffusivity in the alloy
L_salt = 0.5 # m, salt depth (volume / area)
L_alloy = 3.0e-4 # m, resolved alloy thickness
cs0 = ${fparse 0.17 * 8.00e6 / 51.9961} # mol/m^3, 17 wt% Cr in 316 (rho = 8.00 g/cm^3)
c_salt0 = 0.02 # mol/m^3, nearly clean salt
year = ${fparse 365.25 * 86400}

[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 1
    # salt reservoir (well mixed: few cells), then the alloy with 1 um cells
    dx = '${L_salt} ${L_alloy}'
    ix = '25 300'
    subdomain_id = '0 1'
  []
  [names]
    type = RenameBlockGenerator
    input = cmg
    old_block = '0 1'
    new_block = 'salt alloy'
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = names
    primary_block = salt
    paired_block = alloy
    new_boundary = interface
  []
[]

[CorrosionPlating]
  topology = two_block
  salt_block = 'salt'
  solid_block = 'alloy'
  interface_boundary = 'interface'
  elements = 'Cr'
  recession_element = Cr
  temperature = ${T}
  reference_temperature = ${T}
  material_class = stainless_316
  salt_class = fluoride_fuel
  redox_class = msre_or_fuel_baseline
  position_class = nominal
  flow_factor = 1.0
  delta_T_C = 110.0
  applied_overpotential = 0.1
  default_salt_concentration = ${c_salt0}
  reference_concentration = 1.0
  initial_condition_variables = 'cs_Cr'
  initial_condition_values = '${cs0}'
  solid_diffusivity = ${D_alloy}
  # Dissolution slows as the alloy surface depletes: anodic branch x cs / cs0
  solid_activity = true
  transient = true
  verbose = true
[]

[Postprocessors]
  [time_d]
    type = ParsedPostprocessor
    expression = 't / 86400'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [corrosion_rate_um_y]
    type = SideAverageValue
    variable = corrosion_rate_um_y
    boundary = interface
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [recession_um]
    type = SideAverageValue
    variable = recession_um
    boundary = interface
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Mean dissolved Cr [wppm] for a 3.3 g/cm^3 fuel salt
  [c_Cr_salt]
    type = ElementAverageValue
    variable = c_Cr
    block = salt
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [salt_Cr_wppm]
    type = ParsedPostprocessor
    expression = 'c_Cr_salt * 51.9961 / 3.3'
    pp_names = 'c_Cr_salt'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Surface Cr of the alloy relative to the bulk alloy
  [surface_Cr]
    type = ElementExtremeValue
    variable = cs_Cr
    block = alloy
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [surface_Cr_rel]
    type = ParsedPostprocessor
    expression = 'surface_Cr / ${cs0}'
    pp_names = 'surface_Cr'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Inventories per unit wetted area [mol/m^2]
  [cr_salt]
    type = ElementIntegralVariablePostprocessor
    variable = c_Cr
    block = salt
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_alloy]
    type = ElementIntegralVariablePostprocessor
    variable = cs_Cr
    block = alloy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_balance_error]
    type = ParsedPostprocessor
    expression = 'abs(cr_salt + cr_alloy - ${fparse c_salt0 * L_salt + cs0 * L_alloy}) '
                 '/ ${fparse cs0 * L_alloy}'
    pp_names = 'cr_salt cr_alloy'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Cr lost by the alloy expressed as a Faradaic recession [um], relative to recession_um
  [faraday_ratio]
    type = ParsedPostprocessor
    expression = 'dcr := ${fparse cs0 * L_alloy} - cr_alloy; '
                 'if(dcr > 0, recession_um / (dcr * 51.9961 / 8.00e6 * 1e6), 1)'
    pp_names = 'cr_alloy recession_um'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Integral depletion depth: the Cr deficit divided by the surface deficit
  [depletion_depth_um]
    type = ParsedPostprocessor
    expression = '(${fparse cs0 * L_alloy} - cr_alloy) / (${cs0} - surface_Cr + 1e-30) * 1e6'
    pp_names = 'cr_alloy surface_Cr'
    execute_on = 'TIMESTEP_END'
  []
  [depth_over_sqrt_Dt]
    type = ParsedPostprocessor
    expression = 'depletion_depth_um * 1e-6 / sqrt(${D_alloy} * t)'
    pp_names = 'depletion_depth_um'
    use_t = true
    execute_on = 'TIMESTEP_END'
  []
[]

[VectorPostprocessors]
  # Cr profile across the alloy (distance from the interface = x - L_salt)
  [alloy_profile]
    type = LineValueSampler
    variable = cs_Cr
    start_point = '${fparse L_salt + 5e-7} 0 0'
    end_point = '${fparse L_salt + 2.0e-4} 0 0'
    num_points = 100
    sort_by = x
    execute_on = 'FINAL'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  end_time = ${year}
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 100
    growth_factor = 1.3
  []
  dtmax = ${fparse year / 50}
  nl_abs_tol = 1.0e-10
  nl_rel_tol = 1.0e-10
[]

[Outputs]
  [out]
    type = CSV
    execute_vector_postprocessors_on = FINAL
  []
[]
