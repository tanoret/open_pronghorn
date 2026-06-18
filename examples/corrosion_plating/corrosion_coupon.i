# =====================================================================================================
# Two-block molten salt corrosion coupon (Cr, Fe, Ni).
#
# A Hastelloy N coupon (the solid block) is immersed in initially clean fluoride salt (the salt
# block). The solid starts at the alloy composition (constant Cr, Fe and Ni molar concentrations) and
# the salt starts essentially free of dissolved metal. As corrosion proceeds the Butler-Volmer
# interface dissolves each metal into the salt: the salt concentrations rise from near zero and
# diffuse away from the wall into the bulk, while a chromium-depletion layer grows in the solid near
# the interface (chromium diffuses through the alloy toward the corroding surface).
#
# Geometry (1D): a 0.5 m salt reservoir (coarse, effectively well mixed) backs a 1 mm alloy coupon
# (finely resolved). The reservoir/area ratio is chosen so the dissolved chromium reaches the tens to
# hundreds of ppm expected over a couple of years at the calibrated rate -- a small bath would
# concentrate the corrosion product unrealistically fast.
#
# Run:   open_pronghorn-opt -i corrosion_coupon.i
# View:  the Exodus output shows c_Cr/c_Fe/c_Ni rising and diffusing into the salt, cs_Cr/cs_Fe/cs_Ni
#        depleting near the wall, and recession_um (the Faradaic-equivalent penetration).
# =====================================================================================================

[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 1
    # salt reservoir (0.5 m, 50 coarse cells) then the alloy coupon (1 mm, 100 cells of 10 micron)
    dx = '0.5 0.001'
    ix = '50 100'
    subdomain_id = '0 1'
  []
  [names]
    type = RenameBlockGenerator
    input = cmg
    old_block = '0 1'
    new_block = 'salt solid'
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = names
    primary_block = salt
    paired_block = solid
    new_boundary = interface
  []
[]

[CorrosionPlating]
  topology = two_block
  salt_block = 'salt'
  solid_block = 'solid'
  interface_boundary = 'interface'
  elements = 'Cr Fe Ni'
  recession_element = Cr
  temperature = 923.15
  reference_temperature = 923.15
  material_class = hastelloy_n
  salt_class = fluoride_fuel
  redox_class = oxidizing_fef2
  applied_overpotential = 0.1
  # The salt starts nearly clean (0.02 mol/m^3); the calibrated rate is referenced to 1 mol/m^3.
  default_salt_concentration = 0.02
  reference_concentration = 1.0
  # Constant alloy composition of the coupon (Hastelloy N: ~7 wt% Cr, ~4 wt% Fe, ~71 wt% Ni).
  initial_condition_variables = 'cs_Cr cs_Fe cs_Ni'
  initial_condition_values = '11968 6368 107540'
  # Representative solid-state diffusivity so the depletion layer is resolvable on the example mesh
  # (the calibrated chromium value is ~4e-18 m^2/s, which would confine depletion to ~1 cell).
  solid_diffusivity = 1.0e-15
  transient = true
[]

[Postprocessors]
  [salt_Cr_avg]
    type = ElementAverageValue
    variable = c_Cr
    block = salt
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [salt_Cr_ppm]
    type = ParsedPostprocessor
    pp_names = 'salt_Cr_avg'
    # ppm Cr by mass: c [mol/m^3] * 52 g/mol / (2000 kg/m^3 salt) -> mg/kg
    expression = 'salt_Cr_avg * 51.9961 / 2.0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [solid_Cr_min]
    type = ElementExtremeValue
    variable = cs_Cr
    block = solid
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_salt]
    type = ElementIntegralVariablePostprocessor
    variable = c_Cr
    block = salt
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_solid]
    type = ElementIntegralVariablePostprocessor
    variable = cs_Cr
    block = solid
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_total]
    type = LinearCombinationPostprocessor
    pp_names = 'cr_salt cr_solid'
    pp_coefs = '1.0 1.0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [recession_um]
    type = SideAverageValue
    variable = recession_um
    boundary = interface
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  # ~2 years of exposure in roughly month-long steps.
  end_time = 6.3e7
  dt = 1.5e6
  nl_abs_tol = 1.0e-10
[]

[Outputs]
  exodus = true
  csv = true
[]
