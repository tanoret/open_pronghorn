# =====================================================================================================
# Training case 12: which alloy, which salt, which redox state? -- comparing with loop data.
#
# One 0D coupon (as in case 10) whose initial Butler-Volmer rate is the calibrated corrosion
# correlation for the case features (alloy, salt, redox state, temperature, loop position, flow and
# thermal gradient). The driver run_scan.py reruns this input with command-line overrides for seven
# measured loop/capsule corrosion rates of validation/corrosion (data/validation_predictions.csv):
#
#   case         alloy             salt           redox                T [C]   measured [um/y]
#   ORNL-FL-04   hastelloy_n       fluoride_fuel  purified_baseline    704     1.0
#   FB-03        hastelloy_n       fluoroborate   purified_baseline    687     7.4
#   ORNL-FL-01   stainless_304l    fluoride_fuel  msre_or_fuel         688     21.8
#   ORNL-FL-12   stainless_316     fluoride_fuel  msre_or_fuel         650     25.4
#   MOD-F-07     gh3535            flinak         purified_baseline    700     23.3 (5.9-91.6)
#   ORNL-FL-15   stainless_316     flibe          impure_moisture      650     213
#   MOD-F-07s    hastelloy_n       flinak         stressed             700     217
#
# Run one case directly (defaults below are ORNL-FL-04):
#   open_pronghorn-opt -i coupon_rate.i
# or the whole scan with comparison table and figure:
#   python run_scan.py
# =====================================================================================================

T = 977.15 # K
material = hastelloy_n
salt = fluoride_fuel
redox = purified_baseline
position = nominal
flow = 1.0
dT = 166.0

[Mesh]
  [salt]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 2
    xmax = 1.0e-2
  []
[]

[CorrosionPlating]
  topology = salt_only
  reaction_boundary = left
  elements = 'Cr'
  temperature = ${T}
  reference_temperature = ${T}
  material_class = ${material}
  salt_class = ${salt}
  redox_class = ${redox}
  position_class = ${position}
  flow_factor = ${flow}
  delta_T_C = ${dT}
  applied_overpotential = 0.1
  default_salt_concentration = 1.0
  transient = true
[]

[Postprocessors]
  # Initial (fresh-salt) rate = the calibrated correlation for these features
  [corrosion_rate_um_y]
    type = SideAverageValue
    variable = corrosion_rate_um_y
    boundary = left
    execute_on = 'INITIAL'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 1
  num_steps = 1
  nl_abs_tol = 1.0e-12
[]

[Outputs]
  csv = true
[]
