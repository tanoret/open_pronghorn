# Parametrized 0D corrosion cell used to reproduce a single reference case with the mechanistic
# CorrosionPlating action. The case feature selectors (material_class, salt_class, redox_class,
# position_class, flow_factor, delta_T_C and the temperatures) are overridden on the command line by
# run_corrosion_validation.py. The solve is skipped: only the initial state is evaluated, where the
# Butler-Volmer boundary current Faradaically reproduces the calibrated dissolution rate, read out by
# the corrosion_rate_um_y auxiliary.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
    xmax = 1.0e-3
  []
[]

[Problem]
  solve = false
[]

[CorrosionPlating]
  topology = salt_only
  reaction_boundary = left
  elements = 'Cr'
  temperature = 923.15
  reference_temperature = 923.15
  material_class = hastelloy_n
  salt_class = fluoride_fuel
  redox_class = purified_baseline
  position_class = nominal
  flow_factor = 1.0
  delta_T_C = 0.0
  applied_overpotential = 0.1
  default_salt_concentration = 1.0
  transient = true
[]

[Postprocessors]
  [corrosion_rate_um_y]
    type = SideAverageValue
    variable = corrosion_rate_um_y
    boundary = left
    execute_on = 'INITIAL'
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1.0
[]

[Outputs]
  csv = true
  execute_on = 'INITIAL'
[]
