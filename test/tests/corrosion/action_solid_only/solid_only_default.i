# Solid-only Action regression. The absent salt is held at c_ref and the Action's default
# counter-phase potential must realize phi_metal - phi_salt = applied_overpotential, not its negative.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 10
    xmax = 1.0e-3
  []
[]

[CorrosionPlating]
  topology = solid_only
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
  counter_concentration = 1
  default_salt_concentration = 1.0
  default_solid_concentration = 1.0e4
  transient = true
[]

[Postprocessors]
  [solid_concentration]
    type = ElementAverageValue
    variable = cs_Cr
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
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 1.0e5
  num_steps = 1
  automatic_scaling = true
  nl_abs_tol = 1.0e-12
[]

[Outputs]
  csv = true
[]
