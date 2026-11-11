# Balances and figures of merit of the segregated element model (exact face fluxes of the
# Rhie-Chow object). The 3-D mesh is a fraction 'segregated_inlet_fraction' of the element; the
# flows are scaled to the whole element.

# the static outlet temperature is computed; it is written with segregated_static_outputs = 'csv'
segregated_static_outputs = none

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
  # The energy equation carries the enthalpy without the kinetic energy of the gas: the first law
  # conserves h + |u|^2 / 2, so T_mixed is the mixed stagnation temperature. The static one removes
  # the kinetic energy flow through the outlet, from the exact outlet face fluxes (about 15 K at
  # 700 m/s).
  [K_out]
    type = NTPRhieChowFluxIntegral
    boundary = outlet
    rhie_chow_user_object = ins_rhie_chow_interpolator
    normal_kinetic_energy = true
    rho = rho
    execute_on = 'timestep_end'
    outputs = none
  []
  [H_out_static]
    type = ParsedPostprocessor
    pp_names = 'H_out K_out'
    expression = 'H_out - K_out'
    execute_on = 'timestep_end'
    outputs = none
  []
  [T_mixed_static]
    type = NTPMixedTemperature
    fp = h2
    enthalpy_flow = H_out_static
    mass_flow = mdot_out
    pressure = ${p_out}
    execute_on = 'timestep_end'
    outputs = ${segregated_static_outputs}
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
