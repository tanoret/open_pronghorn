# Balances and figures of merit of the segregated element model (exact face fluxes of the
# Rhie-Chow object). The 3-D mesh is a fraction 'segregated_inlet_fraction' of the element; the
# flows are scaled to the whole element.

[Postprocessors]
  [mdot_in]
    type = NTPRhieChowFluxIntegral
    boundary = inlet
    rhie_chow_user_object = ins_rhie_chow_interpolator
    execute_on = 'timestep_end'
    outputs = none
  []
  [mdot_out]
    type = NTPRhieChowFluxIntegral
    boundary = outlet
    rhie_chow_user_object = ins_rhie_chow_interpolator
    execute_on = 'timestep_end'
    outputs = none
  []
  [H_in]
    type = NTPRhieChowFluxIntegral
    boundary = inlet
    rhie_chow_user_object = ins_rhie_chow_interpolator
    carried_quantity = h
    execute_on = 'timestep_end'
    outputs = none
  []
  [H_out]
    type = NTPRhieChowFluxIntegral
    boundary = outlet
    rhie_chow_user_object = ins_rhie_chow_interpolator
    carried_quantity = h
    execute_on = 'timestep_end'
    outputs = none
  []
  # the fission power the gas receives: with the cell volumes of the finite-volume system, as the
  # source kernel deposits it (the element quadrature differs on the non-planar cells of the 3-D
  # mesh by a few 1e-4)
  [power_mesh]
    type = NTPFVCellIntegral
    functor = ${power_density}
    block = fuel_bed
    execute_on = 'timestep_end'
    outputs = none
  []
  # element totals
  [mdot]
    type = ParsedPostprocessor
    pp_names = mdot_out
    expression = 'mdot_out / ${segregated_inlet_fraction}'
    execute_on = 'timestep_end'
  []
  [power]
    type = ParsedPostprocessor
    pp_names = power_mesh
    expression = 'power_mesh / ${segregated_inlet_fraction}'
    execute_on = 'timestep_end'
  []
  # (H_out - H_in - Q) / Q with outward fluxes
  [energy_error]
    type = ParsedPostprocessor
    pp_names = 'H_in H_out power_mesh'
    expression = 'if(power_mesh > 0, (H_out + H_in - power_mesh) / power_mesh, 0)'
    execute_on = 'timestep_end'
  []
  [mass_error]
    type = ParsedPostprocessor
    pp_names = 'mdot_in mdot_out'
    expression = '(mdot_out + mdot_in) / mdot_out'
    execute_on = 'timestep_end'
  []
  [T_mixed]
    type = NTPMixedTemperature
    fp = h2
    enthalpy_flow = H_out
    mass_flow = mdot_out
    pressure = ${p_out}
    execute_on = 'timestep_end'
  []
  [T_fuel_max]
    type = ElementExtremeValue
    variable = T_solid
    block = fuel_bed
    execute_on = 'timestep_end'
  []
  [p_in]
    type = SideAverageValue
    variable = pressure
    boundary = inlet
    execute_on = 'timestep_end'
  []
[]
