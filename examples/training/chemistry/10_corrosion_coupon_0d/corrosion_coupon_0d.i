# =====================================================================================================
# Training case 10: a stainless-steel coupon in a static fluoride-salt capsule.
#
# The wall (left boundary) is a Butler-Volmer electrode: chromium, the least noble alloy component,
# dissolves as Cr2+ into the salt,
#
#   i = i0 [ exp(alpha z f eta) - (c_Cr / c_ref) exp(-alpha z f eta) ],    f = F / (R T),
#
# and the Cr flux is i / (z F). The action seeds i0 from the calibrated corrosion correlation so that
# the initial rate reproduces the reference prediction for this alloy/salt/redox class. This is case
# ORNL-FL-01 of validation/corrosion (304L stainless in fuel salt, 688 C): model 16.82 um/y,
# measured 21.8 um/y in the loop.
#
# In a closed capsule the dissolved chromium has nowhere to go. As c_Cr builds up the cathodic
# (re-plating) branch grows and the net current vanishes at the Nernst equilibrium
#
#   c_eq = c_ref exp(z f eta)            (alpha_a = alpha_c = 0.5)
#
# = 11.2 mol/m^3 (~ 180 wppm in a 3.3 g/cm^3 salt) for eta = 0.1 V, approached with a time constant
# of ~ 13 days for a salt depth (volume/area) of 1 cm. The recession then stops at the Faraday value
#
#   d_final = (c_eq - c_0) L M_Cr / rho_alloy = 0.66 um,
#
# against 16.8 um/y x 120 d = 5.5 um if the initial rate were sustained.
#
# Lesson: corrosion in a static capsule is self-limiting; sustained attack needs a sink for the
# dissolved Cr -- the cold leg of a loop (case 13) or a more oxidizing salt (larger eta, case 12).
#
# Checks: 'rate_vs_reference' (initial rate / 16.82 um/y), 'c_eq_ratio' (wall concentration over the
# Nernst value at the end) and 'faraday_ratio' (recession over the dissolved Cr inventory).
# =====================================================================================================

T = 961.15 # K (688 C)
eta = 0.1 # V, metal-minus-salt potential difference
L = 1.0e-2 # m, salt depth = salt volume / wetted area
c0 = 1.0 # mol/m^3 initial dissolved Cr (~ 16 wppm), also the BV reference c_ref
M_Cr = 51.9961 # g/mol
rho = 8.00 # g/cm^3, stainless_304l in data/corrosion_database.json
day = 86400

[Mesh]
  [salt]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 20
    xmax = ${L}
    # Refine towards the corroding wall (left)
    bias_x = 1.2
  []
[]

[CorrosionPlating]
  topology = salt_only
  reaction_boundary = left
  elements = 'Cr'
  temperature = ${T}
  reference_temperature = ${T}
  # Calibrated correlation classes (data/corrosion_database.json)
  material_class = stainless_304l
  salt_class = fluoride_fuel
  redox_class = msre_or_fuel_baseline
  position_class = nominal
  flow_factor = 1.0
  delta_T_C = 100.0
  applied_overpotential = ${eta}
  default_salt_concentration = ${c0}
  transient = true
  verbose = true
[]

[Postprocessors]
  [time_d]
    type = ParsedPostprocessor
    expression = 't / ${day}'
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
  [c_Cr_wall]
    type = SideAverageValue
    variable = c_Cr
    boundary = left
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_Cr_mean]
    type = ElementAverageValue
    variable = c_Cr
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Initial rate relative to the validation-database prediction for ORNL-FL-01
  [rate_vs_reference]
    type = ParsedPostprocessor
    expression = 'corrosion_rate_um_y / 16.82090218683773'
    pp_names = 'corrosion_rate_um_y'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Wall concentration relative to the Nernst equilibrium c_ref exp(z f eta)
  [c_eq_ratio]
    type = ParsedPostprocessor
    expression = 'c_Cr_wall / (${c0} * exp(2 * 96485.33212 / (8.314462618 * ${T}) * ${eta}))'
    pp_names = 'c_Cr_wall'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Recession relative to the Faraday value from the dissolved Cr: (c_mean - c0) L M / rho
  [faraday_ratio]
    type = ParsedPostprocessor
    expression = 'if(c_Cr_mean > ${c0} + 1e-12, recession_um / '
                 '((c_Cr_mean - ${c0}) * ${L} * ${M_Cr} / (${rho} * 1e6) * 1e6), 1)'
    pp_names = 'recession_um c_Cr_mean'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  end_time = ${fparse 120 * day}
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 600
    growth_factor = 1.2
  []
  dtmax = ${fparse day / 2}
  nl_abs_tol = 1.0e-12
  nl_rel_tol = 1.0e-10
[]

[Outputs]
  csv = true
[]
