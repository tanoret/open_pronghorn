# Stage 1: one packed bed, isothermal gas, solved with the segregated SIMPLE algorithm.
#
# The same bed as packed_bed.i, with the momentum equation for the superficial velocity u and the
# packed-bed drag as a Darcy friction mu D u, D = (A + B|u|) / mu, of
# NTPPackedBedSegregatedFunctorMaterial: the formulation of the segregated element model. The
# velocity is imposed at the inlet and the pressure is zero at the outlet; the pressure drop is the
# one this solver path predicts. validate.py runs it once per experimental point.
#
#   open_pronghorn-opt -i packed_bed_simple.i eps=0.385 d=6.35e-3 L=1.524 u=0.5
#
# (The velocity is uniform, so SIMPLE's normalized momentum residual, which divides by the
# deviation of A u from its mean, stays at a round-off ratio: the run continues to the iteration
# limit, where the pressure drop is converged to round-off.)

eps = 0.385          # bed porosity [-]
d = 6.35e-3          # particle diameter [m]
L = 1.524            # bed length [m]
u = 0.5              # superficial velocity [m/s]
rho = 1.1726         # gas density [kg/m^3]
mu = 1.83538e-5      # gas viscosity [Pa s]
model = kta          # kta or ergun
n_cells = 20

[Mesh]
  [g]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${n_cells}
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
