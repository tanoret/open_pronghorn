# A 1-D isothermal packed bed (TAMU air conditions) solved in the Darcy-Forchheimer form,
# -div(rho M grad p) = 0 with the NTPPackedBedFunctorMaterial mobility. The inlet pressure
# must equal the analytic KTA 3102.3 pressure drop, and the superficial speed the imposed one.
# Changing pressure_drop_model to ergun is checked against the analytic Ergun drop.
rho = 1.1726
mu = 1.83538e-5
eps = 0.385
d = 6.35e-3
L = 1.524
u_s = 0.5
[Mesh]
  [g]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 20
    xmax = ${L}
  []
[]
[Variables]
  [pressure]
    type = MooseVariableFVReal
    initial_condition = 0
  []
[]
[FVKernels]
  [mass]
    type = FVAnisotropicDiffusion
    variable = pressure
    coeff = rho_mobility
  []
[]
[FVBCs]
  [in]
    type = FVFunctorNeumannBC
    variable = pressure
    boundary = left
    functor = ${fparse rho * u_s}
  []
  [out]
    type = FVDirichletBC
    variable = pressure
    boundary = right
    value = 0
  []
[]
[FunctorMaterials]
  [props]
    type = ADGenericFunctorMaterial
    prop_names = 'rho mu k cp'
    prop_values = '${rho} ${mu} 0.026 1006'
  []
  [bed]
    type = NTPPackedBedFunctorMaterial
    porosity = ${eps}
    particle_diameter = ${d}
    solid_conductivity = 1
  []
[]
[Postprocessors]
  [p_in]
    type = SideAverageValue
    variable = pressure
    boundary = left
  []
  [u_in]
    type = SideAverageFunctorPostprocessor
    functor = superficial_speed
    boundary = left
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
    pp_names = 'p_in dp_kta'
    expression = 'p_in / dp_kta - 1'
  []
  [error_ergun]
    type = ParsedPostprocessor
    pp_names = 'p_in dp_ergun'
    expression = 'p_in / dp_ergun - 1'
  []
[]
[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
[]
[Outputs]
  csv = true
[]
