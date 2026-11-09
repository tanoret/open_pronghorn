# Stage 1: one packed bed, isothermal gas, solved with the element's Darcy-Forchheimer equation.
#
# The bed is a 1-D column of length L with superficial velocity u at the inlet and p = 0 at the
# outlet. The equation is the one the fuel element uses, -div(rho M grad p) = 0, with the mobility M
# of NTPPackedBedFunctorMaterial, so the inlet pressure is the pressure drop the correlation predicts
# along this solver path. validate.py runs it once per experimental point.
#
#   open_pronghorn-opt -i packed_bed.i eps=0.385 d=6.35e-3 L=1.524 u=0.5 rho=1.1726 mu=1.835e-5

eps = 0.385          # bed porosity [-]
d = 6.35e-3          # particle diameter [m]
L = 1.524            # bed length [m]
u = 0.5              # superficial velocity [m/s]
rho = 1.1726         # gas density [kg/m^3]
mu = 1.83538e-5      # gas viscosity [Pa s]
model = kta          # kta or ergun
n_cells = 20

[Mesh]
  [bed]
    type = GeneratedMeshGenerator
    dim = 1
    nx = ${n_cells}
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
    type = NTPDarcyMassFlux
    variable = pressure
    pressure = pressure
    rho_mobility = rho_mobility
  []
[]

[FVBCs]
  # an inflow mass flux rho u at the inlet, the outlet at zero gauge pressure
  [inlet]
    type = FVFunctorNeumannBC
    variable = pressure
    boundary = left
    functor = ${fparse rho * u}
  []
  [outlet]
    type = FVDirichletBC
    variable = pressure
    boundary = right
    value = 0
  []
[]

[FunctorMaterials]
  [gas]
    type = ADGenericFunctorMaterial
    prop_names = 'rho mu k cp'
    prop_values = '${rho} ${mu} 0.026 1006'
  []
  [bed]
    type = NTPPackedBedFunctorMaterial
    porosity = ${eps}
    particle_diameter = ${d}
    pressure_drop_model = ${model}
    solid_conductivity = 1
  []
[]

[Postprocessors]
  [dp]
    type = SideAverageValue
    variable = pressure
    boundary = left
  []
  [u_bed]
    type = ElementAverageValue
    variable = u_cell
  []
[]

[AuxVariables]
  [u_cell]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [u_cell]
    type = FunctorAux
    variable = u_cell
    functor = superficial_speed
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-14
[]

[Outputs]
  csv = true
[]
