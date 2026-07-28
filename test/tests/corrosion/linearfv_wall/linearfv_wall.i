# Linear finite-volume corrosion validation: a single well-mixed salt cell with a Butler-Volmer wall
# reaction (no flow, no diffusion), so the cell concentration rises only by the wall source. The
# CorrosionPlatingFlow action seeds the exchange current from the calibrated correlation, so the
# wall flux at the initial reference concentration reproduces the reference dissolution rate. This is
# the ORNL-FL-01 case (304L stainless, hot fluoride loop, reference rate 16.821 um/y) -- the same
# case validated through the finite-element path in test/tests/corrosion/action_salt_only.
#
# The corrosion current and the equivalent penetration rate are reported as postprocessors. The
# initial penetration rate equals the reference value, confirming the linear-FV Butler-Volmer wall
# boundary is consistent with the finite-element implementation.

L = 1.0e-3            # cell size [m]
M = 51.9961           # Cr molar mass [g/mol]
rho = 8.0             # 304L density [g/cm^3]
c_ref = 1.0           # reference concentration [mol/m^3]

[Mesh]
  [cell]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
    xmax = ${L}
  []
[]

[Problem]
  linear_sys_names = 'c_Cr_sys'
[]

[CorrosionPlatingFlow]
  elements = 'Cr'
  reaction_boundary = 'left'
  temperature = 961.15
  reference_temperature = 961.15
  material_class = stainless_304l
  salt_class = fluoride_fuel
  redox_class = msre_or_fuel_baseline
  flow_factor = 1.0
  delta_T_C = 100.0
  applied_overpotential = 0.1
  reference_concentration = ${c_ref}
  initial_concentration = ${c_ref}
  temperature_dependent_kinetics = false
  time_derivative = true
  # A diffusion kernel must be present for the wall flux boundary condition to be assembled; in a
  # single cell it has no internal faces, so this leaves the cell well mixed.
  diffusivity = 1.0e-9
[]

[LinearFVBCs]
  # The far face is sealed (no diffusive flux), so the cell concentration rises only by the wall.
  [sealed]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = c_Cr
    boundary = 'right'
    functor = 0
  []
[]

[Postprocessors]
  [c_Cr]
    type = ElementAverageValue
    variable = c_Cr
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Molar flux into the cell from the wall, inferred from the well-mixed mass balance
  # J = L * dc/dt [mol/m^2/s] (the wall face area is unity in 1D).
  [dc_dt]
    type = ChangeOverTimePostprocessor
    postprocessor = c_Cr
    change_with_respect_to_initial = false
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Equivalent penetration rate [um/y] = corrosion_current_to_um_y(J * z F), with the Faradaic
  # conversion J*M/(z F rho) expressed in um/y. Evaluated at the first step (c ~ c_ref).
  [corrosion_rate_um_y]
    type = ParsedPostprocessor
    pp_names = 'dc_dt'
    # Well-mixed mass balance: J = L * (dc/dt) [mol/m^2/s] (unit wall area in 1D). The Faradaic
    # penetration rate then reduces to rate_um_y = J * M / rho * SEC_PER_YEAR (the z F factors of the
    # current<->molar-flux and current->rate conversions cancel). dt = 10 s.
    expression = '(${L} * dc_dt / 10.0) * ${M} / ${rho} * 3.15576e7'
    execute_on = 'TIMESTEP_END'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 20
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'c_Cr_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  dt = 10.0
  num_steps = 1
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
