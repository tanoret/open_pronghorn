# Stage 1: the hydrogen properties of NTPHydrogenFluidProperties on a temperature sweep.
#
# A 1-D mesh carries a linear temperature from T_lo to T_hi at pressure p; the gas density,
# enthalpy, isobaric specific heat, viscosity and conductivity are evaluated by the fluid
# properties at every cell and written to hydrogen_properties_out_props_*.csv.
# check_hydrogen.py compares them with the same model written in Python (common/hydrogen.py),
# which the post-processing scripts use.
#
#   open_pronghorn-opt -i hydrogen_properties.i

T_lo = 20
T_hi = 6000
p = 7e6
n = 299

!include ../common/hydrogen_nasa.i

[Mesh]
  [line]
    type = GeneratedMeshGenerator
    dim = 1
    nx = ${n}
    xmin = ${fparse T_lo - 0.5 * (T_hi - T_lo) / (n - 1)}
    xmax = ${fparse T_hi + 0.5 * (T_hi - T_lo) / (n - 1)}
  []
[]

[Problem]
  solve = false
[]

[AuxVariables]
  [T]
    type = MooseVariableFVReal
  []
  [rho]
    type = MooseVariableFVReal
  []
  [h]
    type = MooseVariableFVReal
  []
  [cp]
    type = MooseVariableFVReal
  []
  [mu]
    type = MooseVariableFVReal
  []
  [k]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  # the cell centers are the sample temperatures T_lo ... T_hi
  [T]
    type = FunctionAux
    variable = T
    function = x
    execute_on = initial
  []
  [rho]
    type = FunctorAux
    variable = rho
    functor = rho_h2
    execute_on = timestep_end
  []
  [h]
    type = FunctorAux
    variable = h
    functor = h_from_p_T
    execute_on = timestep_end
  []
  [cp]
    type = FunctorAux
    variable = cp
    functor = cp_h2
    execute_on = timestep_end
  []
  [mu]
    type = FunctorAux
    variable = mu
    functor = mu_h2
    execute_on = timestep_end
  []
  [k]
    type = FunctorAux
    variable = k
    functor = k_h2
    execute_on = timestep_end
  []
[]

[FunctorMaterials]
  [props]
    type = NTPGasPropertiesFunctorMaterial
    fp = h2
    pressure = ${p}
    temperature = T
    density_name = rho_h2
    viscosity_name = mu_h2
    specific_heat_name = cp_h2
    conductivity_name = k_h2
  []
  # h_from_p_T is the fluid-property enthalpy at (pressure, T_fluid); the material also needs an
  # enthalpy functor for its inverse T_from_p_h, which is not used here
  [enthalpy]
    type = LinearFVEnthalpyFunctorMaterial
    fp = h2
    pressure = ${p}
    T_fluid = T
    h = 0
  []
[]

[VectorPostprocessors]
  [props]
    type = ElementValueSampler
    variable = 'T rho h cp mu k'
    sort_by = x
    execute_on = timestep_end
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  [out]
    type = CSV
    execute_on = final
  []
[]
