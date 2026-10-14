# =====================================================================================================
# Training case 06: iodine-135 / xenon-135 -- equilibrium poisoning and the post-shutdown xenon peak.
#
# Xe-135 has the largest thermal capture cross section of any nuclide (2.65e6 b). It is produced
# directly in fission and, mostly, by the decay of I-135 (6.6 h); it is removed by decay (9.1 h) and
# by neutron capture ("burnout"):
#
#   dN_I/dt  = Y_I F / N_A - (lambda_I + sigma_I phi) N_I
#   dN_Xe/dt = Y_Xe F / N_A + lambda_I N_I - (lambda_Xe + sigma_Xe phi) N_Xe
#
# The reactor runs at full power for 4 days (Xe reaches equilibrium) and then scrams: F and phi drop
# to zero. Burnout stops but the iodine keeps decaying into xenon, so Xe-135 *rises* for several hours
# before it decays away -- the "xenon pit" that can prevent a thermal reactor from restarting.
#
# Checks against the exact solutions (columns *_analytic):
#   equilibrium   N_I = Y_I F / (N_A (lambda_I + sigma_I phi)),
#                 N_Xe = (Y_Xe F / N_A + lambda_I N_I) / (lambda_Xe + sigma_Xe phi)
#   after scram   N_Xe(t) = N_Xe0 e^{-lambda_Xe t} + lambda_I N_I0 (e^{-lambda_I t} - e^{-lambda_Xe t})
#                                                    / (lambda_Xe - lambda_I)
# With phi = 3e17 n/m^2/s the peak comes 8.4 h after the scram at 1.89x the equilibrium inventory.
#
# The xenon worth is estimated as rho_Xe = -(sigma_Xe N_Xe N_A phi_0) / (nu F_0): the fraction of the
# fission neutrons that Xe-135 would absorb at full power (nu = 2.43). This is the reactivity the
# operator must overcome to restart.
# =====================================================================================================

P0 = 2.0e7 # W/m^3, full-power density
phi0 = 3.0e17 # n/m^2/s (3e13 n/cm^2/s), full-power one-group thermal flux
t_scram = ${fparse 4 * 86400}
t_end = ${fparse t_scram + 2 * 86400}
nu = 2.43

[Mesh]
  [cell]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Functions]
  # Relative power: 1 until the scram, then 0 (a step). The implicit time integration evaluates
  # the source at the end of each step, so the step that ends at the scram must still see full
  # power: compare with a 1 s tolerance on the scram time.
  [relative_power]
    type = ParsedFunction
    expression = 'if(t < ${t_scram} + 1, 1, 0)'
  []
  [power_density]
    type = ParsedFunction
    expression = '${P0} * relative_power'
    symbol_names = 'relative_power'
    symbol_values = 'relative_power'
  []
  # The flux follows the power (constant flux-to-power ratio)
  [flux]
    type = ParsedFunction
    expression = '${phi0} * relative_power'
    symbol_names = 'relative_power'
    symbol_values = 'relative_power'
  []
[]

[Problem]
  linear_sys_names = 'I135_sys Xe135_sys'
[]

[FissionProductInventory]
  nuclides = 'I135 Xe135'
  # Functors: here time-dependent functions
  power_density = power_density
  neutron_flux = flux
  # Start at zero inventory (fresh fuel)
[]

[Postprocessors]
  [time_h]
    type = ParsedPostprocessor
    expression = 't / 3600'
    use_t = true
  []
  [N_I135]
    type = ElementAverageValue
    variable = I135
  []
  [N_Xe135]
    type = ElementAverageValue
    variable = Xe135
  []
  # Xe-135 worth at full-power flux [pcm]
  [rho_Xe_pcm]
    type = ParsedPostprocessor
    expression = '-1e5 * 2647600e-28 * N_Xe135 * 6.02214076e23 * ${phi0} '
                 '/ (${nu} * ${P0} / (200e6 * 1.602176634e-19))'
    pp_names = 'N_Xe135'
  []
  # Exact equilibrium and post-scram solutions
  [N_I135_analytic]
    type = ParsedPostprocessor
    expression = 'lI := log(2) / 23652; F := ${P0} / (200e6 * 1.602176634e-19); '
                 'NI0 := 0.0628187 * F / 6.02214076e23 / (lI + 80e-28 * ${phi0}); '
                 'if(t < ${t_scram}, NI0, NI0 * exp(-lI * (t - ${t_scram})))'
    use_t = true
  []
  [N_Xe135_analytic]
    type = ParsedPostprocessor
    expression = 'lI := log(2) / 23652; lX := log(2) / 32904; b := 0.99950473; '
                 'F := ${P0} / (200e6 * 1.602176634e-19); '
                 'NI0 := 0.0628187 * F / 6.02214076e23 / (lI + 80e-28 * ${phi0}); '
                 'NX0 := ((0.065385 - b * 0.0628187) * F / 6.02214076e23 + b * lI * NI0) '
                 '/ (lX + 2647600e-28 * ${phi0}); s := t - ${t_scram}; '
                 'if(t < ${t_scram}, NX0, NX0 * exp(-lX * s) '
                 '+ b * lI * NI0 / (lX - lI) * (exp(-lI * s) - exp(-lX * s)))'
    use_t = true
  []
  # Largest Xe inventory so far, relative to the full-power equilibrium
  [N_Xe135_max]
    type = TimeExtremeValue
    postprocessor = N_Xe135
    value_type = max
    outputs = none
  []
  [xe_peak_ratio]
    type = ParsedPostprocessor
    expression = 'lI := log(2) / 23652; lX := log(2) / 32904; b := 0.99950473; '
                 'F := ${P0} / (200e6 * 1.602176634e-19); '
                 'NI0 := 0.0628187 * F / 6.02214076e23 / (lI + 80e-28 * ${phi0}); '
                 'NX0 := ((0.065385 - b * 0.0628187) * F / 6.02214076e23 + b * lI * NI0) '
                 '/ (lX + 2647600e-28 * ${phi0}); N_Xe135_max / NX0'
    pp_names = 'N_Xe135_max'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  end_time = ${t_end}
  # 10-minute steps, refined to 30 s for the first hour after the scram: the power step breaks the
  # smoothness that BDF2 relies on, so the steps right after it must be small. The scram time is
  # hit exactly at a step boundary.
  [TimeStepper]
    type = FunctionDT
    function = 'if(t >= ${t_scram} - 1 & t < ${t_scram} + 3600, 30, 600)'
  []
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
