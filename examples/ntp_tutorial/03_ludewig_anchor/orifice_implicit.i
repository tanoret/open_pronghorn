# Orifice designed implicitly, together with the flow: add to anchor.i with continuation.i,
#
#   open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i
#
# The cold-frit drag of axial slice j is multiplied by chi_j = exp(lambda_j), lambda being a scalar
# unknown with one component per slice. Its equations (NTPSliceUniformityConstraint) make the mean
# gas temperature of every slice's hot frit equal to the mean over all slices, which is a flow
# proportional to the slice fission power, and fix the mean of lambda to 'orifice_level' (the
# logarithm of the geometric mean of chi). Every Newton step therefore solves the orifice that
# matches the current power, so the power ramp never meets the flow-excursion fold of a fixed
# orifice (with uniform frits it sits near 8.5 % of the design power).
#
# The level gives every slice room to be opened as well as throttled. Referencing one slice
# instead (chi = 1 there) fails when that slice is not the most starved one: the starved slices
# then need their drag removed altogether, where their temperature no longer depends on chi.
# A level too high before the design is active throttles the element uniformly, which moves the
# fold of the still uniform frits to a lower power (3 % at orifice_level = 3.5).
orifice_level = 2

# The design equations are switched on gradually once the gas is heated (before that, every
# orifice gives the same uniform 30 K and the orifice cannot be determined).
design_on = 0.06
design_full = 0.2

# the division block has its own 'nz' parameter, so the slice count is passed under another name
n_slices = ${nz}

[Functions]
  [design_activation]
    type = ParsedFunction
    expression = 'min(1, max(0, (t - ${design_on}) / (${design_full} - ${design_on})))'
  []
[]

[MeshDivisions]
  [slices]
    type = CartesianGridDivision
    # covers the RZ strip and the 3-D element alike
    bottom_left = '${fparse -r5_out} 0 ${fparse -r5_out}'
    top_right = '${r5_out} ${L} ${r5_out}'
    nx = 1
    ny = ${n_slices}
    nz = 1
  []
[]

# one component per slice; MOOSE names scalar orders (TWELFTH for the 12-slice anchor mesh)
log_chi_order = TWELFTH

[Variables]
  [log_chi]
    family = SCALAR
    order = ${log_chi_order}
    initial_condition = 0
  []
[]

[FVKernels]
  [mass_flux]
    orifice_log_multiplier = log_chi
    mesh_division = slices
    orifice_blocks = cold_frit
  []
  [gas_advection]
    orifice_log_multiplier = log_chi
    mesh_division = slices
    orifice_blocks = cold_frit
  []
  [orifice_design]
    type = NTPSliceUniformityConstraint
    variable = T_fluid
    lambda = log_chi
    mesh_division = slices
    field = T_fluid
    block = hot_frit
    reference = mean
    mean_log_multiplier = ${orifice_level}
    # the mismatch is measured relative to the current mean temperature rise of the slices
    field_scale = 1
    field_scale_pp = mean_rise
    activation = design_activation
  []
[]

[FunctorMaterials]
  [mass_flux]
    orifice_log_multiplier = log_chi
    mesh_division = slices
    orifice_blocks = cold_frit
  []
[]

[VectorPostprocessors]
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
  [T_hot_frit_mean]
    type = ElementAverageValue
    variable = T_fluid
    block = hot_frit
    execute_on = 'initial linear timestep_end'
  []
  [mean_rise]
    type = ParsedPostprocessor
    pp_names = T_hot_frit_mean
    expression = 'T_hot_frit_mean - ${T_in}'
    execute_on = 'initial linear timestep_end'
  []
  # orifice multiplier of the exit slice, the most throttled one
  [chi_exit]
    type = ParsedPostprocessor
    pp_names = log_chi_exit
    expression = 'exp(log_chi_exit)'
  []
  [log_chi_exit]
    type = ScalarVariable
    variable = log_chi
    component = ${fparse nz - 1}
  []
[]

[Problem]
  # log_chi has no ScalarKernel: its equations come from the FV constraint
  kernel_coverage_check = false
[]
