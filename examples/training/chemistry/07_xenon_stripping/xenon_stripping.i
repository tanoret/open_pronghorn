# =====================================================================================================
# Training case 07: removing the Xe-135 poison by off-gas sparging -- a molten salt reactor advantage.
#
# Xenon is insoluble in the salt. In an MSR it can be stripped continuously by sparging helium
# bubbles through the salt and venting the cover gas (MSRE: a spray ring in the pump bowl). The
# stripping adds a first-order removal r to the Xe-135 balance of case 06:
#
#   N_Xe = (Y_Xe F / N_A + lambda_I N_I) / (lambda_Xe + sigma_Xe phi + r)
#
# The stripped xenon goes to an off-gas inventory (Xe135_offgas, per unit salt volume) where it keeps
# decaying to Cs-135, so the off-gas system holds a decay heat and a Cs-135 source term.
#
# Each of the 5 cells of the mesh is an independent well-mixed salt sample with a different
# stripping rate r = 0, 1e-5, 1e-4, 1e-3, 1e-2 1/s (stripping time constants from infinity to
# 100 s), so a single run is a parameter scan. With phi = 3e17 n/m^2/s the fraction of the fission
# neutrons absorbed by Xe-135 is 2.1 % without stripping (about 2000 pcm of reactivity), 1.1 % for
# r = 1e-4 1/s and 0.19 % for r = 1e-3 1/s. Stripping must compete with decay + burnout,
# lambda_Xe + sigma_Xe phi = 1e-4 1/s, before it pays off. (MSRE achieved a Xe-135 poisoning of a
# few tenths of a percent with its pump-bowl spray, see ORNL-4865.)
# =====================================================================================================

P0 = 2.0e7 # W/m^3
phi0 = 3.0e17 # n/m^2/s
nu = 2.43
n_cells = 5

[Mesh]
  [scan]
    type = GeneratedMeshGenerator
    dim = 1
    nx = ${n_cells}
    xmax = ${n_cells}
  []
[]

[Functions]
  # Cell i (centroid x = i + 0.5): r = 0 for i = 0, else 10^(i - 6) 1/s
  [stripping_rate]
    type = ParsedFunction
    expression = 'if(x < 1, 0, 10^(floor(x) - 6))'
  []
[]

[Problem]
  linear_sys_names = 'I135_sys Xe135_sys Cs135_sys Xe135_offgas_sys Cs135_offgas_sys'
[]

[FissionProductInventory]
  nuclides = 'I135 Xe135 Cs135'
  power_density = ${P0}
  neutron_flux = ${phi0}
  # Noble gases are removed at the (spatially varying) stripping rate into the off-gas
  removal_groups = 'noble_gas'
  removal_rates = 'stripping_rate'
  removal_destinations = 'offgas'
[]

[AuxVariables]
  # Fraction of the fission neutrons absorbed by Xe-135, sigma_Xe N_Xe N_A phi / (nu F)
  [xe_poison_fraction]
    type = MooseLinearVariableFVReal
  []
  [r_strip]
    type = MooseLinearVariableFVReal
  []
[]

[AuxKernels]
  [xe_poison_fraction]
    type = ParsedAux
    variable = xe_poison_fraction
    expression = '2647600e-28 * Xe135 * 6.02214076e23 * ${phi0} '
                 '/ (${nu} * ${P0} / (200e6 * 1.602176634e-19))'
    coupled_variables = 'Xe135'
    execute_on = 'TIMESTEP_END'
  []
  [r_strip]
    type = FunctionAux
    variable = r_strip
    function = stripping_rate
    execute_on = 'INITIAL'
  []
[]

[VectorPostprocessors]
  # One row per stripping rate
  [scan]
    type = ElementValueSampler
    variable = 'r_strip I135 Xe135 Xe135_offgas Cs135_offgas xe_poison_fraction'
    sort_by = x
    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]
  [time_h]
    type = ParsedPostprocessor
    expression = 't / 3600'
    use_t = true
  []
  # Time histories without and with strong stripping
  [N_Xe135_no_strip]
    type = PointValue
    variable = Xe135
    point = '0.5 0 0'
  []
  [N_Xe135_strip_1e-3]
    type = PointValue
    variable = Xe135
    point = '3.5 0 0'
  []
  [poison_no_strip]
    type = PointValue
    variable = xe_poison_fraction
    point = '0.5 0 0'
  []
  [poison_strip_1e-3]
    type = PointValue
    variable = xe_poison_fraction
    point = '3.5 0 0'
  []
  # Analytic equilibrium for r = 1e-3 1/s
  [N_Xe135_strip_1e-3_analytic]
    type = ParsedPostprocessor
    expression = 'lI := log(2) / 23652; lX := log(2) / 32904; b := 0.99950473; '
                 'F := ${P0} / (200e6 * 1.602176634e-19); '
                 'NI0 := 0.0628187 * F / 6.02214076e23 / (lI + 80e-28 * ${phi0}); '
                 '((0.065385 - b * 0.0628187) * F / 6.02214076e23 + b * lI * NI0) '
                 '/ (lX + 2647600e-28 * ${phi0} + 1e-3)'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  dt = 600
  end_time = ${fparse 5 * 86400}
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  [out]
    type = CSV
    execute_vector_postprocessors_on = FINAL
  []
[]
