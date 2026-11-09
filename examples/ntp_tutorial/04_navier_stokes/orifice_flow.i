# Orifice designed on the slice mass flows while the segregated element is solved: add to
# anchor_simple.i (SIMPLE, fixed-point continuation) or anchor_pimple.i (PIMPLE pseudo-transient)
# instead of the profile of the Newton anchor,
#
#   open_pronghorn-opt -i anchor_simple.i orifice_flow.i
#   open_pronghorn-opt -i anchor_pimple.i orifice_flow.i
#   open_pronghorn-opt -i anchor_simple.i orifice_flow.i simple_fixed_point=300 orifice_max_updates=200
#                                                  the orifice frozen after 200 solves (it drifts)
#
# After every solve (fixed-point iteration with SIMPLE, time step with PIMPLE)
# NTPOrificeFlowDesignAux moves the cold-frit multiplier chi_j of each axial slice so that the
# slice carries its target share w_j of the element flow, at a fixed geometric mean
# exp(orifice_level):
#
#   ln chi_j <- ln chi_j + alpha ln(m_j / (w_j mdot)),   then shifted to a mean orifice_level.
#
# The slice flows m_j are the exact Rhie-Chow fluxes through the plenum / cold-frit interface
# (NTPRhieChowSliceFlow). The flow through a slice falls monotonically with its drag, at the cold
# density of the frit, so this loop does not pass through the density of the hot gas: the
# temperature-driven design of orifice_lagged.i, which does, oscillates slice by slice.
#
# The target shares are the slice power shares: the same enthalpy rise per unit mass in every
# slice, so a uniform slice outlet temperature, up to the axial conduction in the bed and the frits.
# With orifice_share_exponent > 0 the shares are corrected further from the slice outlet
# temperatures, w_j <- w_j (rise_j / mean rise)^gamma, whenever the flows match their targets to
# orifice_share_tolerance. On this element the power shares already give a slice outlet spread of
# 0.1 % of the rise; the temperature correction re-introduces the density feedback and, with
# gamma >= 0.15, a slow oscillation of the inlet slice (README.md).

orifice_level = 3
orifice_flow_exponent = 0.7
orifice_max_change = 1.5
orifice_share_exponent = 0
orifice_max_share_change = 1.01
orifice_share_tolerance = 0.01
# the design is frozen after this many updates (to test the designed orifice as a fixed one)
orifice_max_updates = 1000000
n_slices = ${nz}
# azimuthal sectors of the design (3-D only; 1 = rings)
orifice_sectors = 1
# print the flow mismatch, the profile and the shares after every update
orifice_verbose = true

cold_frit_multiplier = chi_design

[MeshDivisions]
  # the columns of the design: axial slices (orifice_sectors = 1, the rings of the orifice), or
  # slices times azimuthal sectors of the half element (3-D)
  [slices]
    type = NTPElementColumnsDivision
    axial_min = 0
    axial_max = ${L}
    n_axial = ${n_slices}
    n_azimuthal = ${orifice_sectors}
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
    type = NTPOrificeFlowDesignAux
    variable = chi_design
    block = cold_frit
    mesh_division = slices
    slice_flows = slice_flow
    slice_powers = slice_power
    slice_power_vector = ${power_density}
    slice_temperatures = slice_outlet
    slice_temperature_vector = T_fluid
    inlet_temperature = ${T_in}
    level = ${orifice_level}
    flow_exponent = ${orifice_flow_exponent}
    max_change_factor = ${orifice_max_change}
    share_exponent = ${orifice_share_exponent}
    max_share_change = ${orifice_max_share_change}
    share_tolerance = ${orifice_share_tolerance}
    max_updates = ${orifice_max_updates}
    verbose = ${orifice_verbose}
  []
[]

[VectorPostprocessors]
  # mass flow entering each slice's cold frit from the plenum
  [slice_flow]
    type = NTPRhieChowSliceFlow
    rhie_chow_user_object = ins_rhie_chow_interpolator
    mesh_division = slices
    from_block = plenum
    to_block = cold_frit
    execute_on = 'timestep_end'
  []
  # fission power of each slice
  [slice_power]
    type = MeshDivisionFunctorReductionVectorPostprocessor
    mesh_division = slices
    functors = ${power_density}
    reduction = integral
    block = fuel_bed
    execute_on = 'initial timestep_end'
  []
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
  [T_slice_max]
    type = VectorPostprocessorReductionValue
    vectorpostprocessor = slice_outlet
    vector_name = T_fluid
    value_type = max
  []
  [T_slice_min]
    type = VectorPostprocessorReductionValue
    vectorpostprocessor = slice_outlet
    vector_name = T_fluid
    value_type = min
  []
  # spread of the slice outlet temperatures relative to the mean rise
  [T_slice_spread]
    type = ParsedPostprocessor
    pp_names = 'T_slice_max T_slice_min T_mixed'
    expression = '(T_slice_max - T_slice_min) / max(T_mixed - ${T_in}, 1)'
  []
[]

