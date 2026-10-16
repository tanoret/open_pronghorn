# Closed salt gap between a hot wall (left, 1000 K) and a cold wall (right, 900 K), set up by the
# CorrosionPlatingFlow action with the opt-in thermodynamic wall law. Each element relaxes toward
# its own equilibrium c_eq,k(T) = c_eq,k exp[-dH_k / R (1/T - 1/T_ref)]; with dH > 0 the hot wall
# dissolves and the cold wall plates, and the salt settles between the two wall equilibria.
# Cr is given a large equilibrium and Ni (the noble component) a small one, so Ni barely dissolves.
# The calibrated rate is split 90/10 between Cr and Ni (seed_fractions), so the total metal loss of
# fresh salt matches the calibrated rate instead of being counted once per element.

T_hot = 1000
T_cold = 900

[Mesh]
  [gap]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 50
    xmax = 0.1
  []
[]

[Functions]
  [T]
    type = ParsedFunction
    expression = '${T_hot} + (${T_cold} - ${T_hot}) * x / 0.1'
  []
[]

[Problem]
  linear_sys_names = 'c_Cr_sys c_Ni_sys'
[]

[CorrosionPlatingFlow]
  elements = 'Cr Ni'
  reaction_boundary = 'left right'
  temperature = T
  reference_temperature = 950
  material_class = hastelloy_n
  salt_class = fluoride_fuel
  redox_class = purified_baseline
  diffusivity = 1e-5
  equilibrium_concentrations = '10 0.01'
  dissolution_enthalpies = '50'
  seed_fractions = '0.9 0.1'
  initial_concentration = 0
  time_derivative = false
  verbose = true
[]

[Postprocessors]
  [Cr_hot]
    type = PointValue
    variable = c_Cr
    point = '0.001 0 0'
  []
  [Cr_cold]
    type = PointValue
    variable = c_Cr
    point = '0.099 0 0'
  []
  [Cr_mean]
    type = ElementAverageValue
    variable = c_Cr
  []
  [Ni_mean]
    type = ElementAverageValue
    variable = c_Ni
  []
  # Wall equilibria: the salt settles between them, below c_eq at the hot wall (dissolving) and
  # above c_eq at the cold wall (plating)
  [Cr_eq_hot]
    type = ParsedPostprocessor
    expression = '10 * exp(-50e3 / 8.314462618 * (1 / ${T_hot} - 1 / 950))'
  []
  [Cr_eq_cold]
    type = ParsedPostprocessor
    expression = '10 * exp(-50e3 / 8.314462618 * (1 / ${T_cold} - 1 / 950))'
  []
  # Net wall fluxes into the salt [mol/m^2/s]: dissolution (+) at the hot wall, plating (-) at the
  # cold wall, equal and opposite at steady state
  [Cr_flux_hot]
    type = ParsedPostprocessor
    expression = '1e-5 * (Cr_hot - Cr_cold) / 0.098'
    pp_names = 'Cr_hot Cr_cold'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 40
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Steady
  system_names = 'c_Cr_sys c_Ni_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  l_tol = 1e-12
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
