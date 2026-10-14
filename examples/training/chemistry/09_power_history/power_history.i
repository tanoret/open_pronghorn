# =====================================================================================================
# Training case 09: a month of reactor operation -- start-up, load follow, shutdown and decay.
#
# The fuel salt sees a realistic power history (relative power p(t), flux following the power):
#   day 0-1     start-up ramp 0 -> 100 %
#   day 1-10    full power
#   day 10-13   load follow at 50 %
#   day 13-20   full power
#   day 20      shutdown, followed by 10 days of decay
#
# Combining the physics of cases 05-07 in one input, this shows the different "memories" of the
# inventory:
#   * Xe-135 follows the power with a lag of hours and swings the reactivity on every power change.
#     With sparging at r = 1e-4 1/s (case 07) the equilibrium worth is ~ -1070 pcm and the
#     post-shutdown peak (2.5 h after the shutdown) only -1330 pcm; rerun with
#     'FissionProductInventory/removal_rates=1e-12' to see -4030 pcm at 8.5 h without stripping;
#   * I-131 (8 d) follows the power history of the last weeks;
#   * Cs-137 (30 y) integrates the energy produced: N_Cs = Y E / (E_f N_A) to within its decay
#     (0.1 % over a month) -- the classic burnup monitor ('Cs137_burnup_check').
# =====================================================================================================

P0 = 2.0e7 # W/m^3
phi0 = 3.0e17 # n/m^2/s
nu = 2.43
day = 86400

[Mesh]
  [cell]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Functions]
  [relative_power]
    type = PiecewiseLinear
    x = '0 ${day} ${fparse 10 * day} ${fparse 10 * day + 3600} ${fparse 13 * day} '
        '${fparse 13 * day + 3600} ${fparse 20 * day} ${fparse 20 * day + 60} ${fparse 30 * day}'
    y = '0 1 1 0.5 0.5 1 1 0 0'
  []
  [power_density]
    type = ParsedFunction
    expression = '${P0} * p'
    symbol_names = 'p'
    symbol_values = 'relative_power'
  []
  [flux]
    type = ParsedFunction
    expression = '${phi0} * p'
    symbol_names = 'p'
    symbol_values = 'relative_power'
  []
[]

[Problem]
  linear_sys_names = 'I135_sys Xe135_sys I131_sys Cs137_sys Xe135_offgas_sys'
[]

[FissionProductInventory]
  nuclides = 'I135 Xe135 I131 Cs137'
  power_density = power_density
  neutron_flux = flux
  removal_groups = 'noble_gas'
  removal_rates = '1e-4'
  removal_destinations = 'offgas'
[]

[Postprocessors]
  [time_d]
    type = ParsedPostprocessor
    expression = 't / ${day}'
    use_t = true
  []
  [power]
    type = FunctionValuePostprocessor
    function = relative_power
  []
  [N_I135]
    type = ElementAverageValue
    variable = I135
  []
  [N_Xe135]
    type = ElementAverageValue
    variable = Xe135
  []
  [N_I131]
    type = ElementAverageValue
    variable = I131
  []
  [N_Cs137]
    type = ElementAverageValue
    variable = Cs137
  []
  [N_Xe135_offgas]
    type = ElementAverageValue
    variable = Xe135_offgas
  []
  # Xe-135 reactivity worth at the full-power flux [pcm] (see case 06)
  [rho_Xe_pcm]
    type = ParsedPostprocessor
    expression = '-1e5 * 2647600e-28 * N_Xe135 * 6.02214076e23 * ${phi0} '
                 '/ (${nu} * ${P0} / (200e6 * 1.602176634e-19))'
    pp_names = 'N_Xe135'
  []
  # Energy produced per unit volume [J/m^3]
  [power_density_pp]
    type = FunctionValuePostprocessor
    function = power_density
    outputs = none
  []
  [energy]
    type = TimeIntegratedPostprocessor
    value = power_density_pp
  []
  [burnup_MWd_per_m3]
    type = ParsedPostprocessor
    expression = 'energy / (1e6 * ${day})'
    pp_names = 'energy'
  []
  # Cs-137 inventory relative to Y E / (E_f N_A)
  [Cs137_burnup_check]
    type = ParsedPostprocessor
    expression = 'if(energy > 0, N_Cs137 / (0.0618832 * energy / (200e6 * 1.602176634e-19) '
                 '/ 6.02214076e23), 1)'
    pp_names = 'N_Cs137 energy'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  end_time = ${fparse 30 * day}
  # 30-minute steps, 1 minute around the fast power changes (load changes and the shutdown)
  [TimeStepper]
    type = FunctionDT
    function = 'if((t > ${fparse 10 * day - 60} & t < ${fparse 10 * day + 7200}) | '
               '(t > ${fparse 13 * day - 60} & t < ${fparse 13 * day + 7200}) | '
               '(t > ${fparse 20 * day - 60} & t < ${fparse 20 * day + 7200}), 60, 1800)'
  []
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
