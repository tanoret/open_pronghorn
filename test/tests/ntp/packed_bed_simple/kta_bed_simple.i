# The KTA packed bed of ../packed_bed solved with the segregated SIMPLE algorithm: superficial
# velocity momentum with the drag of NTPPackedBedSegregatedFunctorMaterial (LinearFVMomentumFriction,
# force mu D u with D = (A + B|u|)/mu). The bed is 1-D, so the converged pressure drop must equal the
# analytic KTA 3102.3 drop and the velocity the imposed one. (The velocity is uniform, so SIMPLE's
# normalized momentum residual, which divides by the deviation of A u from its mean, stays at a
# round-off ratio: the run continues to the iteration limit and the pressure drop is checked.)
rho = 1.1726
mu = 1.83538e-5
eps = 0.385
d = 6.35e-3
L = 1.524
u_s = 0.5
model = kta

[Mesh]
  [g]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 1
    xmax = ${L}
    ymax = 0.1
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system'
[]

[Physics/NavierStokes/FlowSegregated/flow]
  compressibility = incompressible
  density = rho
  dynamic_viscosity = mu
  initial_velocity = '${u_s} 0 0'
  initial_pressure = 0
  inlet_boundaries = left
  momentum_inlet_types = fixed-velocity
  momentum_inlet_functors = '${u_s} 0'
  wall_boundaries = 'top bottom'
  momentum_wall_types = 'symmetry symmetry'
  outlet_boundaries = right
  momentum_outlet_types = fixed-pressure
  pressure_functors = 0
  friction_types = darcy
  friction_coeffs = darcy_drag
  momentum_advection_interpolation = upwind
[]

[FunctorMaterials]
  [gas]
    type = GenericFunctorMaterial
    prop_names = 'rho mu k cp'
    prop_values = '${rho} ${mu} 0.026 1006'
  []
  [bed]
    type = NTPPackedBedSegregatedFunctorMaterial
    porosity = ${eps}
    particle_diameter = ${d}
    pressure_drop_model = ${model}
    solid_conductivity = 1
  []
[]

[Postprocessors]
  [dp]
    type = PressureDrop
    pressure = pressure
    upstream_boundary = left
    downstream_boundary = right
    boundary = 'left right'
  []
  [dp_kta]
    type = ParsedPostprocessor
    expression = 'Re := ${rho} * ${u_s} * ${d} / (${mu} * (1 - ${eps})); (320 / Re + 6 / Re^0.1) * ${L} * (1 - ${eps}) / (${eps}^3 * ${d}) * 0.5 * ${rho} * ${u_s}^2'
  []
  [dp_ergun]
    type = ParsedPostprocessor
    expression = '(150 * ${mu} * (1 - ${eps})^2 / (${eps}^3 * ${d}^2) * ${u_s} + 1.75 * ${rho} * (1 - ${eps}) / (${eps}^3 * ${d}) * ${u_s}^2) * ${L}'
  []
  [error_kta]
    type = ParsedPostprocessor
    pp_names = 'dp dp_kta'
    expression = 'dp / dp_kta - 1'
  []
  [error_ergun]
    type = ParsedPostprocessor
    pp_names = 'dp dp_ergun'
    expression = 'dp / dp_ergun - 1'
  []
[]

[Executioner]
  type = SIMPLE
  rhie_chow_user_object = ins_rhie_chow_interpolator
  momentum_systems = 'u_system v_system'
  pressure_system = pressure_system
  momentum_equation_relaxation = 0.8
  pressure_variable_relaxation = 0.3
  num_iterations = 200
  pressure_absolute_tolerance = 1e-11
  momentum_absolute_tolerance = 1e-11
  momentum_petsc_options_iname = '-pc_type'
  momentum_petsc_options_value = 'lu'
  pressure_petsc_options_iname = '-pc_type'
  pressure_petsc_options_value = 'lu'
  momentum_l_abs_tol = 1e-14
  pressure_l_abs_tol = 1e-14
  momentum_l_tol = 0
  pressure_l_tol = 0
  continue_on_max_its = true
[]

[Outputs]
  csv = true
[]
