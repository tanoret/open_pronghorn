# Single-domain (salt-only) corrosion set up by the CorrosionPlating action. The metal wall is an
# external Butler-Volmer boundary. The action seeds the exchange current so that, at the applied
# overpotential and the initial (reference) concentration, the interfacial current reproduces the
# calibrated dissolution rate of the reference correlation. This case is ORNL-FL-01 (304L stainless,
# hot fluoride loop), whose reference corrosion rate is 16.761 um/y; the rate auxiliary reads that
# value at the initial state.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 10
    xmax = 1.0e-3
  []
[]

[CorrosionPlating]
  topology = salt_only
  reaction_boundary = left
  elements = 'Cr'
  temperature = 961.15
  reference_temperature = 961.15
  material_class = stainless_304l
  salt_class = fluoride_fuel
  redox_class = msre_or_fuel_baseline
  position_class = nominal
  flow_factor = 1.0
  delta_T_C = 100.0
  applied_overpotential = 0.1
  default_salt_concentration = 1.0
  transient = true
[]

[Postprocessors]
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
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 1.0e5
  num_steps = 2
  nl_abs_tol = 1.0e-12
[]

[Outputs]
  csv = true
[]
