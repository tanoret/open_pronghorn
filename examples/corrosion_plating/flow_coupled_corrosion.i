# =====================================================================================================
# Flow-coupled molten salt wall corrosion (salt-only topology, one-way velocity coupling).
#
# Clean fluoride salt enters a heated channel and flows past a corroding alloy wall (the bottom
# boundary), which sheds chromium through a Butler-Volmer electrode reaction at the calibrated rate.
# The shed chromium is carried downstream by the flow and spread by dispersion, so a corrosion-product
# plume develops along the wall and grows toward the outlet. The plume reaches a steady shape within a
# few flow-through times, while the wall itself recedes slowly at the calibrated corrosion rate (a few
# micrometres per year), recorded by the recession field.
#
# The velocity is a prescribed developed channel profile (a functor) read by the corrosion advection
# kernel: this is the recommended one-way coupling (the flow drives the dissolved species, which do
# not feed back on the flow). For a resolved CFD velocity (e.g. the segregated Navier-Stokes solve in
# examples/molten_salt_corrosion/corrosion_channel.i) transfer vel_x/vel_y in with a MultiApp and point
# velocity_x/velocity_y at the transferred fields.
#
# Run:   open_pronghorn-opt -i flow_coupled_corrosion.i
# View:  the Exodus output shows the downstream chromium plume c_Cr and the wall field recession_um.
# =====================================================================================================

H = 0.01           # channel height [m]
L = 0.2            # channel length [m]
v_mean = 2.0e-4    # mean axial velocity [m/s] (slow loop flow)
D_eff = 1.0e-5     # effective dispersion coefficient [m^2/s] (turbulent dispersion >> molecular)
c_inlet = 0.02     # inlet (fresh) salt Cr concentration [mol/m^3]
c_ref = 1.0        # reference concentration at which the calibrated rate applies [mol/m^3]

[Mesh]
  [channel]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0.0
    xmax = ${L}
    ymin = 0.0
    ymax = ${H}
    nx = 80
    ny = 12
  []
[]

[Functions]
  [vel_profile]
    type = ParsedFunction
    expression = '6.0 * v_mean * (y / H) * (1.0 - y / H)'
    symbol_names = 'v_mean H'
    symbol_values = '${v_mean} ${H}'
  []
[]

[CorrosionPlating]
  topology = salt_only
  reaction_boundary = bottom
  elements = 'Cr'
  temperature = 923.15
  reference_temperature = 923.15
  material_class = stainless_316
  salt_class = fluoride_fuel
  redox_class = oxidizing_fef2
  applied_overpotential = 0.1
  salt_diffusivity = ${D_eff}
  default_salt_concentration = ${c_inlet}
  reference_concentration = ${c_ref}
  velocity_x = vel_profile
  transient = true
[]

[BCs]
  [inlet]
    type = ADDirichletBC
    variable = c_Cr
    boundary = left
    value = ${c_inlet}
  []
[]

[Postprocessors]
  [corrosion_rate_um_y]
    type = SideAverageValue
    variable = corrosion_rate_um_y
    boundary = bottom
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [recession_um]
    type = SideAverageValue
    variable = recession_um
    boundary = bottom
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_outlet_ppm]
    type = ParsedPostprocessor
    pp_names = 'cr_outlet'
    expression = 'cr_outlet * 51.9961 / 2.0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_outlet]
    type = SideAverageValue
    variable = c_Cr
    boundary = right
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  # The plume reaches steady state within a few flow-through times (~1000 s); the run then continues
  # for ~1 year so the wall recession accumulates at the steady corrosion rate.
  end_time = 3.15e7
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 200.0
    growth_factor = 2.0
    cutback_factor = 0.5
  []
  nl_abs_tol = 1.0e-11
[]

[Outputs]
  exodus = true
  csv = true
[]
