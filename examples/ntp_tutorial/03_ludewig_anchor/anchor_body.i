# Body of the stage-3 anchor inputs: physics, power, balances and executioner of one element,
# included after the mesh by anchor.i (axisymmetric) and by ../05_element_3d/anchor_3d.i.
# The including input sets the parameters listed in anchor.i and
#   power_shape_expression   the power shape in terms of 'axial_shape' (and x, y, z)

!include ../common/hydrogen_nasa.i
!include ../common/element_physics.i

[Functions]
  [axial_shape]
    type = PiecewiseConstant
    axis = y
    # each slice [z_j, z_j+1) takes its own factor: the breakpoints are the slice starts
    direction = left
    x = '0 0.06078333 0.12156667 0.18235 0.24313333 0.30391667 0.3647 0.42548333 0.48626667
         0.54705 0.60783333 0.66861667'
    y = '1.28 1.22 1.15 1.09 1.03 0.98 0.94 0.91 0.88 0.86 0.84 0.82'
  []
  # the power shape: the axial shape, times the transverse tilt of 3-D inputs
  [shape]
    type = ParsedFunction
    symbol_names = 'axial_shape volume_ratio'
    symbol_values = 'axial_shape volume_ratio'
    expression = 'volume_ratio * (${power_shape_expression})'
  []
  [q_fission]
    type = ParsedFunction
    symbol_names = 'shape'
    symbol_values = 'shape'
    # the shape factors average 1.0 over the 12 intervals; the power ramps up over ramp_time
    expression = '${q_avg} * shape / 1.0 * min(1, t / ${ramp_time})'
  []
[]

[Postprocessors]
  # exact bed volume over the bed volume of the mesh: the power density is scaled by it, so that
  # the element power is exact on a faceted (3-D) mesh too (1 on the axisymmetric mesh)
  [volume_ratio]
    type = ParsedPostprocessor
    pp_names = bed_volume
    expression = '${fparse element_fraction * pi * L / 3 * (r3_in^2 + r3_in * r3_out + r3_out^2
                                                              - r2_in^2 - r2_in * r2_out - r2_out^2)} / bed_volume'
    execute_on = initial
  []
  # the inlet flow is imposed by the boundary condition, so it is reported directly
  [mdot_in]
    type = ConstantPostprocessor
    value = ${mdot}
  []
  [mdot_out]
    type = NTPDarcyFluxIntegral
    variable = pressure
    boundary = outlet
    rho_mobility = rho_mobility
  []
  [H_in]
    type = ConstantPostprocessor
    value = ${fparse mdot * h_in}
  []
  [H_out]
    type = NTPDarcyFluxIntegral
    variable = pressure
    boundary = outlet
    rho_mobility = rho_mobility
    carried_quantity = h
  []
  [power]
    type = ElementIntegralFunctorPostprocessor
    functor = q_fission
    block = fuel_bed
  []
  # first-law residual (H_out - H_in - Q) / Q
  [energy_error]
    type = ParsedPostprocessor
    pp_names = 'H_in H_out power'
    expression = 'if(power > 0, (H_out - H_in - power) / power, 0)'
  []
  [mass_error]
    type = ParsedPostprocessor
    pp_names = 'mdot_in mdot_out'
    expression = '(mdot_out - mdot_in) / mdot_in'
  []
  [h_out_mixed]
    type = ParsedPostprocessor
    pp_names = 'H_out mdot_out'
    expression = 'H_out / mdot_out'
  []
  [T_fuel_max]
    type = ElementExtremeValue
    variable = T_solid
    block = fuel_bed
  []
  [T_hot_frit_max]
    type = ElementExtremeValue
    variable = T_solid
    block = hot_frit
  []
  [p_in]
    type = SideAverageValue
    variable = pressure
    boundary = inlet
  []
  [bed_volume]
    type = VolumePostprocessor
    block = fuel_bed
    execute_on = initial
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu NONZERO mumps'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  nl_max_its = 12
  line_search = none
  automatic_scaling = true
  compute_scaling_once = false
  steady_state_detection = true
  steady_state_tolerance = 1e-9
  steady_state_start_time = ${ramp_time}
  end_time = 1e3
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-5
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 6
  []
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]
