# Orifice designed during the PIMPLE (or SIMPLE) solve of the segregated element: add to
# anchor_pimple.i instead of the profile of the Newton anchor,
#
#   open_pronghorn-opt -i anchor_pimple.i orifice_lagged.i
#
# After every time step NTPOrificeDesignAux moves the cold-frit multiplier of each axial slice
# towards a uniform slice gas outlet temperature, at a fixed geometric mean exp(orifice_level). The
# segregated linear systems cannot hold the orifice as an unknown (as the Newton model does with
# orifice_implicit.i): the design lags the flow by one step, so the power ramp must be slow enough
# for the orifice to follow it.

orifice_level = 4
orifice_exponent = 0.3
orifice_max_change = 1.1
# the orifice is designed during the ramp and frozen once the power is reached: the lagged
# design-flow loop is not stable at a fixed power (it overshoots), and a fixed orifice is what an
# element would have
orifice_freeze_time = ${t_ramp}
n_slices = ${nz}

cold_frit_multiplier = chi_design

[MeshDivisions]
  [slices]
    type = CartesianGridDivision
    bottom_left = '${fparse -r5_out} 0 ${fparse -r5_out}'
    top_right = '${r5_out} ${L} ${r5_out}'
    nx = 1
    ny = ${n_slices}
    nz = 1
  []
[]

[AuxVariables]
  [chi_design]
    type = MooseVariableFVReal
    block = cold_frit
    initial_condition = ${fparse exp(orifice_level)}
  []
[]

[AuxKernels]
  [orifice_design]
    type = NTPOrificeDesignAux
    variable = chi_design
    block = cold_frit
    mesh_division = slices
    slice_temperatures = slice_outlet
    slice_temperature_vector = T_fluid
    inlet_temperature = ${T_in}
    level = ${orifice_level}
    exponent = ${orifice_exponent}
    max_change_factor = ${orifice_max_change}
  []
[]

[VectorPostprocessors]
  # gas temperature leaving each slice, averaged over the hot frit of the slice
  [slice_outlet]
    type = MeshDivisionFunctorReductionVectorPostprocessor
    mesh_division = slices
    functors = T_fluid
    reduction = average
    block = hot_frit
    execute_on = 'initial timestep_end'
  []
[]

[Postprocessors]
  [chi_max]
    type = ElementExtremeValue
    variable = chi_design
    block = cold_frit
  []
  [chi_min]
    type = ElementExtremeValue
    variable = chi_design
    value_type = min
    block = cold_frit
  []
[]

[Controls]
  [freeze_orifice]
    type = TimePeriod
    disable_objects = 'AuxKernels::orifice_design'
    start_time = ${orifice_freeze_time}
  []
[]
