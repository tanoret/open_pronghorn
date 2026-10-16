# =====================================================================================================
# Training case 03, model 2 of 4: NCL-16 as a closed hot/cold salt gap with the THERMODYNAMIC flow
# wall law ([CorrosionPlatingFlow] with equilibrium_concentrations).
#
# The loop is reduced to its essence: one salt volume wetting a hot wall (left, 704 C = 977.15 K) and
# a cold wall (right, 538 C = 811.15 K) of equal area. The gap width 7.8e-4 m makes V / A_total =
# 3.9e-4 m, the calibrated area-to-salt-mass ratio of the loop models (7.77 cm^2/g at 3.3 g/cm^3;
# ncl16_dridn.i). The salt is mixed across the gap (D = 1e-8 m^2/s, mixing time ~ 1 min, much shorter than
# the corrosion time scale): it plays the role of the circulating loop salt.
#
# Each wall drives each element k toward its own equilibrium (alpha_a = alpha_c = 0.5)
#   J_k = J_k,ref exp[Ea/R (1/T_ref - 1/T)] sqrt(c_eq,k(T)/c_eq,k) (1 - c_k / c_eq,k(T)),
#   c_eq,k(T) = c_eq,k exp[-dH_k / R (1/T - 1/T_ref)],   T_ref = 977.15 K (the hot leg).
# J_k,ref is the seed_fractions share (Cr 0.90, Fe 0.07, Ni 0.03: selective Cr leaching) of the
# calibrated 1.382 um/y of the NCL-16 features. With dH > 0 the hot wall dissolves and the cold wall
# plates, and the salt settles between the two wall equilibria.
#
# Equilibria (assumptions, chosen by hand -- a two-parameter CALIBRATION on M-014 and M-012):
#   Cr: c_eq = 900 ppm at 704 C, dH = 100 kJ/mol (72 ppm at 538 C). c_eq,hot is chosen so that the
#       loop-average salt Cr ends in the measured range (+500 ppm, M-014) and dH so that the cold-leg
#       gain is close to the measured 1.7 mg/cm^2 (M-012). The MSTDB-TC CrF2 saturation activities
#       (0.704 at 704 C, 0.483 at 538 C; validation/corrosion/calibration/results/advanced/
#       mstdb_thermochemistry_checks.json) correspond to only dH ~ 15 kJ/mol: see exercise 2.
#   Fe: c_eq = 40 ppm at 704 C, dH = 100 kJ/mol, starting from 150 ppm (an ASSUMED purified-salt Fe
#       level): Cr from the alloy reduces the Fe(II) of the salt and Fe plates out (M-015: -100 ppm).
#   Ni: c_eq = 1 ppm at 704 C (the noble component), dH = 100 kJ/mol.
#
# Gold (3.37 y): salt Cr 511 ppm (measured +500), salt Fe -127 ppm (measured -100), Ni 0.57 ppm;
# hot-wall loss 1.84 mg/cm^2 (measured 2.9), cold-wall gain 1.74 mg/cm^2 (measured 1.7), mean
# hot-wall recession 0.61 um/y (measured 1.016) and a steady hot-wall Cr rate of 0.54 um/y. The
# loss and the rate are PREDICTIONS of this calibrated model (they were not tuned). The hot wall
# dissolves and the cold wall plates Cr; Fe plates on both walls.
# 'inventory_balance_error' (7e-4) is a VERIFICATION: the dissolved metal = the time-integrated wall
# fluxes, i.e. the wall-flux postprocessors reproduce the boundary condition.
# =====================================================================================================

T_hot = 977.15 # K (704 C)
T_cold = 811.15 # K (538 C)
gap = 7.8e-4 # m, so that V / A_total = gap / 2 = 3.9e-4 m
D_mix = 1.0e-8 # m^2/s, salt mixing across the gap
rho_salt = 3.3 # g/cm^3, LiF-BeF2-ThF4-UF4 fuel salt at ~700 C (estimate)
rho_alloy = 8.89e6 # g/m^3, hastelloy_n
M_Cr = 51.9961
M_Fe = 55.845
M_Ni = 58.6934
ppm = ${fparse rho_salt}  # (mol/m^3) = ppm * rho_salt / M
Cr_eq_ppm = 900
Fe_eq_ppm = 40
Ni_eq_ppm = 1
Fe0_ppm = 150 # initial salt Fe
dH_Cr = 100 # kJ/mol
dH_Fe = 100
dH_Ni = 100
seed_Cr = 0.90
seed_Fe = 0.07
seed_Ni = 0.03
rate_ref = 1.38216 # um/y, calibrated rate of the NCL-16 features at T_ref (verbose output)
Ea = 54835.535 # J/mol, corrosion database Ea_corr_kJ_mol
years = 3.3652749258498744
year = 31557600 # s

c_eq_Cr = ${fparse Cr_eq_ppm * ppm / M_Cr}
c_eq_Fe = ${fparse Fe_eq_ppm * ppm / M_Fe}
c_eq_Ni = ${fparse Ni_eq_ppm * ppm / M_Ni}
# fresh-salt wall fluxes at T_ref [mol/m^2/s] (the action normalizes seed_fractions to sum to 1)
J_Cr = ${fparse seed_Cr / (seed_Cr + seed_Fe + seed_Ni) * rate_ref * 1e-6 / year * rho_alloy / M_Cr}
J_Fe = ${fparse seed_Fe / (seed_Cr + seed_Fe + seed_Ni) * rate_ref * 1e-6 / year * rho_alloy / M_Fe}
J_Ni = ${fparse seed_Ni / (seed_Cr + seed_Fe + seed_Ni) * rate_ref * 1e-6 / year * rho_alloy / M_Ni}

[Mesh]
  [gap]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 4
    xmax = ${gap}
  []
[]

[Functions]
  [T]
    type = ParsedFunction
    expression = '${T_hot} + (${T_cold} - ${T_hot}) * x / ${gap}'
  []
[]

[Problem]
  linear_sys_names = 'c_Cr_sys c_Fe_sys c_Ni_sys'
[]

# The dissolved metals are declared here (instead of by the action) so that Fe can start from the
# purified-salt level; the action releases into them (release_variables) and adds the wall BCs.
[Variables]
  [c_Cr]
    type = MooseLinearVariableFVReal
    solver_sys = c_Cr_sys
    initial_condition = 0
  []
  [c_Fe]
    type = MooseLinearVariableFVReal
    solver_sys = c_Fe_sys
    initial_condition = ${fparse Fe0_ppm * ppm / M_Fe}
  []
  [c_Ni]
    type = MooseLinearVariableFVReal
    solver_sys = c_Ni_sys
    initial_condition = 0
  []
[]

[LinearFVKernels]
  [Cr_time]
    type = LinearFVTimeDerivative
    variable = c_Cr
  []
  [Cr_diffusion]
    type = LinearFVDiffusion
    variable = c_Cr
    diffusion_coeff = ${D_mix}
  []
  [Fe_time]
    type = LinearFVTimeDerivative
    variable = c_Fe
  []
  [Fe_diffusion]
    type = LinearFVDiffusion
    variable = c_Fe
    diffusion_coeff = ${D_mix}
  []
  [Ni_time]
    type = LinearFVTimeDerivative
    variable = c_Ni
  []
  [Ni_diffusion]
    type = LinearFVDiffusion
    variable = c_Ni
    diffusion_coeff = ${D_mix}
  []
[]

[CorrosionPlatingFlow]
  elements = 'Cr Fe Ni'
  release_variables = 'c_Cr c_Fe c_Ni'
  reaction_boundary = 'left right'
  temperature = T
  reference_temperature = ${T_hot}
  material_class = hastelloy_n
  salt_class = fluoride_fuel
  redox_class = purified_baseline
  position_class = nominal
  delta_T_C = 166
  temperature_dependent_kinetics = true
  equilibrium_concentrations = '${c_eq_Cr} ${c_eq_Fe} ${c_eq_Ni}'
  dissolution_enthalpies = '${dH_Cr} ${dH_Fe} ${dH_Ni}'
  seed_fractions = '${seed_Cr} ${seed_Fe} ${seed_Ni}'
  verbose = true
[]

[Postprocessors]
  [time_y]
    type = ParsedPostprocessor
    expression = 't / ${year}'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Wall-cell concentrations (the value the BC sees) and the salt averages [ppm]
  [Cr_hot_cell]
    type = PointValue
    variable = c_Cr
    point = '${fparse gap / 8} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [Cr_cold_cell]
    type = PointValue
    variable = c_Cr
    point = '${fparse gap - gap / 8} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [Fe_hot_cell]
    type = PointValue
    variable = c_Fe
    point = '${fparse gap / 8} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [Fe_cold_cell]
    type = PointValue
    variable = c_Fe
    point = '${fparse gap - gap / 8} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [Ni_hot_cell]
    type = PointValue
    variable = c_Ni
    point = '${fparse gap / 8} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [Ni_cold_cell]
    type = PointValue
    variable = c_Ni
    point = '${fparse gap - gap / 8} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [c_Cr_avg]
    type = ElementAverageValue
    variable = c_Cr
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [c_Fe_avg]
    type = ElementAverageValue
    variable = c_Fe
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [c_Ni_avg]
    type = ElementAverageValue
    variable = c_Ni
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [salt_Cr_ppm]
    type = ParsedPostprocessor
    expression = 'c_Cr_avg * ${M_Cr} / ${ppm}'
    pp_names = 'c_Cr_avg'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [salt_Fe_ppm]
    type = ParsedPostprocessor
    expression = 'c_Fe_avg * ${M_Fe} / ${ppm}'
    pp_names = 'c_Fe_avg'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [salt_Ni_ppm]
    type = ParsedPostprocessor
    expression = 'c_Ni_avg * ${M_Ni} / ${ppm}'
    pp_names = 'c_Ni_avg'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [salt_Fe_change_ppm]
    type = ParsedPostprocessor
    expression = 'salt_Fe_ppm - ${Fe0_ppm}'
    pp_names = 'salt_Fe_ppm'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Wall equilibria [ppm]: the salt Cr settles between them
  [Cr_eq_hot_ppm]
    type = ParsedPostprocessor
    expression = '${Cr_eq_ppm}'
    execute_on = 'INITIAL'
  []
  [Cr_eq_cold_ppm]
    type = ParsedPostprocessor
    expression = '${Cr_eq_ppm} * exp(-${dH_Cr}e3 / 8.314462618 * (1 / ${T_cold} - 1 / ${T_hot}))'
    execute_on = 'INITIAL'
  []
  # ---- wall fluxes into the salt [mol/m^2/s], the same law as the boundary condition ----
  [J_Cr_hot]
    type = ParsedPostprocessor
    expression = '${J_Cr} * (1 - c / ${c_eq_Cr})'
    pp_names = 'Cr_hot_cell'
    pp_symbols = 'c'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [J_Fe_hot]
    type = ParsedPostprocessor
    expression = '${J_Fe} * (1 - c / ${c_eq_Fe})'
    pp_names = 'Fe_hot_cell'
    pp_symbols = 'c'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [J_Ni_hot]
    type = ParsedPostprocessor
    expression = '${J_Ni} * (1 - c / ${c_eq_Ni})'
    pp_names = 'Ni_hot_cell'
    pp_symbols = 'c'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [J_Cr_cold]
    type = ParsedPostprocessor
    expression = 'arr := exp(${Ea} / 8.314462618 * (1 / ${T_hot} - 1 / ${T_cold}));
                  eqf := exp(-${dH_Cr}e3 / 8.314462618 * (1 / ${T_cold} - 1 / ${T_hot}));
                  ${J_Cr} * arr * sqrt(eqf) * (1 - c / (${c_eq_Cr} * eqf))'
    pp_names = 'Cr_cold_cell'
    pp_symbols = 'c'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [J_Fe_cold]
    type = ParsedPostprocessor
    expression = 'arr := exp(${Ea} / 8.314462618 * (1 / ${T_hot} - 1 / ${T_cold}));
                  eqf := exp(-${dH_Fe}e3 / 8.314462618 * (1 / ${T_cold} - 1 / ${T_hot}));
                  ${J_Fe} * arr * sqrt(eqf) * (1 - c / (${c_eq_Fe} * eqf))'
    pp_names = 'Fe_cold_cell'
    pp_symbols = 'c'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [J_Ni_cold]
    type = ParsedPostprocessor
    expression = 'arr := exp(${Ea} / 8.314462618 * (1 / ${T_hot} - 1 / ${T_cold}));
                  eqf := exp(-${dH_Ni}e3 / 8.314462618 * (1 / ${T_cold} - 1 / ${T_hot}));
                  ${J_Ni} * arr * sqrt(eqf) * (1 - c / (${c_eq_Ni} * eqf))'
    pp_names = 'Ni_cold_cell'
    pp_symbols = 'c'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  # Areal mass flux of the hot / cold wall into the salt [mg/cm^2/y] (1 g/m^2 = 0.1 mg/cm^2)
  [hot_mass_rate]
    type = ParsedPostprocessor
    expression = '(Cr * ${M_Cr} + Fe * ${M_Fe} + Ni * ${M_Ni}) * 0.1 * ${year}'
    pp_names = 'J_Cr_hot J_Fe_hot J_Ni_hot'
    pp_symbols = 'Cr Fe Ni'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cold_mass_rate]
    type = ParsedPostprocessor
    expression = '(Cr * ${M_Cr} + Fe * ${M_Fe} + Ni * ${M_Ni}) * 0.1 * ${year}'
    pp_names = 'J_Cr_cold J_Fe_cold J_Ni_cold'
    pp_symbols = 'Cr Fe Ni'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [hot_Cr_mass_rate]
    type = ParsedPostprocessor
    expression = 'J_Cr_hot * ${M_Cr} * 0.1 * ${year}'
    pp_names = 'J_Cr_hot'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cold_Cr_mass_rate]
    type = ParsedPostprocessor
    expression = 'J_Cr_cold * ${M_Cr} * 0.1 * ${year}'
    pp_names = 'J_Cr_cold'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Time-integrated areal mass change [mg/cm^2]: hot-wall loss (+) and cold-wall gain (+)
  [hot_loss_rate_s]
    type = ParsedPostprocessor
    expression = 'hot_mass_rate / ${year}'
    pp_names = 'hot_mass_rate'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [cold_gain_rate_s]
    type = ParsedPostprocessor
    expression = '-cold_mass_rate / ${year}'
    pp_names = 'cold_mass_rate'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [mass_loss_mg_cm2]
    type = TimeIntegratedPostprocessor
    value = hot_loss_rate_s
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [mass_gain_mg_cm2]
    type = TimeIntegratedPostprocessor
    value = cold_gain_rate_s
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Exposure-averaged hot-wall recession rate [um/y] = loss / (rho_alloy t)
  [mean_rate_um_y]
    type = ParsedPostprocessor
    expression = 'if(time_y > 0, mass_loss_mg_cm2 / 1e3 / ${fparse rho_alloy / 1e6} * 1e4 / time_y, 0)'
    pp_names = 'mass_loss_mg_cm2 time_y'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Steady-state (final) hot-wall Cr recession rate [um/y]
  [hot_Cr_rate_um_y]
    type = ParsedPostprocessor
    expression = 'J_Cr_hot * ${M_Cr} / ${rho_alloy} * 1e6 * ${year}'
    pp_names = 'J_Cr_hot'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # VERIFICATION: dissolved metal (mol per m^2 of one wall, V/A_wall = gap) = integrated wall fluxes
  [net_flux]
    type = ParsedPostprocessor
    expression = 'a + b + c + d + e + f'
    pp_names = 'J_Cr_hot J_Fe_hot J_Ni_hot J_Cr_cold J_Fe_cold J_Ni_cold'
    pp_symbols = 'a b c d e f'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [net_released]
    type = TimeIntegratedPostprocessor
    value = net_flux
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [inventory_balance_error]
    type = ParsedPostprocessor
    expression = 'abs((c_Cr_avg + c_Fe_avg + c_Ni_avg - ${fparse Fe0_ppm * ppm / M_Fe}) * ${gap} - net_released)'
                 ' / (${fparse Fe0_ppm * ppm / M_Fe} * ${gap})'
    pp_names = 'c_Cr_avg c_Fe_avg c_Ni_avg net_released'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 10
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'c_Cr_sys c_Fe_sys c_Ni_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = bdf2
  end_time = ${fparse years * year}
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 3600
    growth_factor = 1.3
  []
  dtmax = ${fparse 0.005 * year}
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
