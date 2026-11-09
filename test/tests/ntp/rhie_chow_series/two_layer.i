# Two porous layers in series, solved with SIMPLE: a 4-cell layer of Darcy coefficient D_A
# (length L_A) and a 1-cell layer with D_B = 1000 D_A (length L_B), at a uniform superficial
# velocity u. The exact pressure drop is mu u (D_A L_A + D_B L_B). With the series face coupling and
# the series-corrected cell pressure gradient of NTPRhieChowMassFlux it is reproduced to 1e-9 and
# the cell velocities are uniform. The linear interpolation of the base RhieChowMassFlux
# (face_coupling = base, series_pressure_gradient = false) gives the interface face the mobility of
# the permeable layer: the drop is 44 % low.
mu = 1e-5
rho = 1
u = 0.1
D_A = 1e8
D_B = 1e11
L_A = 0.4
L_B = 0.1

[Mesh]
  [g]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${L_A} ${L_B}'
    ix = '4 1'
    dy = 0.05
    iy = 1
    subdomain_id = '0 1'
  []
  [names]
    type = RenameBlockGenerator
    input = g
    old_block = '0 1'
    new_block = 'layer_a layer_b'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system'
[]

[UserObjects]
  [ins_rhie_chow_interpolator]
    type = NTPRhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = flow_p_diffusion
  []
[]

[Physics/NavierStokes/FlowSegregated/flow]
  compressibility = incompressible
  density = ${rho}
  dynamic_viscosity = ${mu}
  initial_velocity = '${u} 0 0'
  initial_pressure = 0
  inlet_boundaries = left
  momentum_inlet_types = fixed-velocity
  momentum_inlet_functors = '${u} 0'
  wall_boundaries = 'top bottom'
  momentum_wall_types = 'symmetry symmetry'
  outlet_boundaries = right
  momentum_outlet_types = fixed-pressure
  pressure_functors = 0
  friction_blocks = 'layer_a; layer_b'
  friction_types = 'darcy; darcy'
  friction_coeffs = 'D_a; D_b'
  momentum_advection_interpolation = upwind
[]

[FunctorMaterials]
  [D_a]
    type = GenericVectorFunctorMaterial
    prop_names = D_a
    prop_values = '${D_A} ${D_A} ${D_A}'
    block = layer_a
  []
  [D_b]
    type = GenericVectorFunctorMaterial
    prop_names = D_b
    prop_values = '${D_B} ${D_B} ${D_B}'
    block = layer_b
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
  [dp_exact]
    type = ParsedPostprocessor
    expression = '${mu} * ${u} * (${D_A} * ${L_A} + ${D_B} * ${L_B})'
  []
  [error]
    type = ParsedPostprocessor
    pp_names = 'dp dp_exact'
    expression = 'dp / dp_exact - 1'
  []
  [mdot_out]
    type = NTPRhieChowFluxIntegral
    boundary = right
    rhie_chow_user_object = ins_rhie_chow_interpolator
  []
[]

[Executioner]
  type = SIMPLE
  rhie_chow_user_object = ins_rhie_chow_interpolator
  momentum_systems = 'u_system v_system'
  pressure_system = pressure_system
  momentum_equation_relaxation = 0.8
  pressure_variable_relaxation = 0.3
  num_iterations = 300
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
