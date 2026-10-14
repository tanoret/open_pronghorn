# =====================================================================================================
# Training case 14: radiolysis-driven corrosion -- an irradiated chloride pool between 316 walls.
#
# Cases 03 and 10 come together. A molten chloride pool (600 C) sits between two 316 stainless
# walls. The walls corrode by Butler-Volmer kinetics releasing Cr(II) (case 10), and the gamma field
# drives the chlorine radical chain that oxidizes the dissolved chromium (case 03):
#
#   wall:   Cr(alloy) <-> Cr(II) + 2 e-         net rate ~ 1 - c_Cr(II),wall / c_eq
#   salt:   Cl2- + Cr(II)  -> Cr(III) + 2 Cl-   (k = 7.2e6)
#           Cl2- + Cr(III) -> Cr(II) + Cl2      (k = 1.4e6)
#
# The salt starts pre-equilibrated with the walls: Cr(II) = c_eq = c_ref exp(z f eta) = 14.3 mol/m^3
# (~ 470 wppm), so without irradiation the net wall current is zero -- the salt is passive (case 10).
# Under irradiation the radiolytic oxidant converts Cr(II) to Cr(III) (Cr(III)/Cr(II) -> 7.2/1.4,
# i.e. 84 % Cr(III) within minutes), the wall sees Cr(II) far below c_eq again and corrodes at ~ 90 %
# of the fresh-salt rate. Net effect: the alloy is the ultimate sink of the radiolytic oxidant
# (Cr(III) acts as a catalyst for Cl2- -> Cl2 + Cl-). After shutdown the Cr(III) remains (nothing
# reduces it) and the walls keep corroding until Cr(II) is back at c_eq: the salt remembers the
# irradiation.
#
# Power history (days): 0.5 d start-up ramp, 5 d at full power, shutdown, 5.5 d of decay.
# Mixing: natural convection stirs the pool; an eddy diffusivity D_mix = 1e-5 m^2/s (mixing time
# L^2 / D ~ 1 h) replaces the molecular value. Case 15 resolves the flow instead.
#
# Checks: 'c_eq_ratio' = 1 at t = 0 (passive salt); 'corrosion_rate_um_y', the wall loss inferred
# from the dissolved-Cr growth, (V/A) dc_Cr/dt M / rho, is ~ 0 before start-up and ~ 0.9 x the
# 39.6 um/y reference rate under irradiation (35 um/y), and stays ~ 33 um/y for days after shutdown.
# The Cl2 by-product accumulates (~ 64 mol/m^3 after 5 days) because this pool has no off-gas --
# case 04 shows the Henry-law sink. Exercise: dose_full=0 -> no corrosion at all.
# =====================================================================================================

dose_full = 1.0e4 # W/m^3 gamma dose rate at full power
T = 873.15 # K
Cl_bulk = 2.0e4 # mol/m^3
c_ref = 1.0 # mol/m^3 Butler-Volmer reference concentration of the calibrated rate
D_mix = 1.0e-5 # m^2/s eddy diffusivity
L = 0.2 # m pool width
eta = 0.1 # V
# Nernst equilibrium Cr(II) of the walls: initial (pre-equilibrated) salt
Cr_init = ${fparse c_ref * exp(2 * 96485.33212 / (8.314462618 * T) * eta)}
day = 86400

[Mesh]
  [pool]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 20
    xmax = ${L}
  []
[]

[Functions]
  [power_history]
    type = PiecewiseLinear
    x = '0 ${fparse 0.5 * day} ${fparse 5.5 * day} ${fparse 5.5 * day + 600} ${fparse 11 * day}'
    y = '0 1 1 0 0'
  []
  [dose_rate]
    type = ParsedFunction
    expression = '${dose_full} * p'
    symbol_names = 'p'
    symbol_values = 'power_history'
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys '
                     'Cr_II_sys Cr_III_sys Cr_I_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  metals = 'Cr'
  temperature = ${T}
  dose_rate = dose_rate
  diffusivity = ${D_mix}
  # Net oxidant yield: the fast e_sol/Cl recombination is lumped into a direct Cl2- yield (as in
  # examples/operating_msr_corrosion), which keeps a days-long transient affordable
  g_value_species = 'e_sol Cl_rad Cl2m_rad'
  g_value_overrides = '0.0 0.0 0.30'
  initial_condition_species = 'Cl_ion Cr_II'
  initial_condition_values = '${Cl_bulk} ${Cr_init}'
  time_derivative = true
[]

[CorrosionPlatingFlow]
  # Wall corrosion releases chromium as the radiolysis-tracked Cr(II)
  elements = 'Cr'
  release_variables = 'Cr_II'
  reaction_boundary = 'left right'
  temperature = ${T}
  reference_temperature = ${T}
  material_class = stainless_316
  salt_class = chloride
  redox_class = chloride_unspecified
  applied_overpotential = ${eta}
  reference_concentration = ${c_ref}
  temperature_dependent_kinetics = false
  verbose = true
[]

[Postprocessors]
  [time_d]
    type = ParsedPostprocessor
    expression = 't / ${day}'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [power]
    type = FunctionValuePostprocessor
    function = power_history
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cr_II_avg]
    type = ElementAverageValue
    variable = Cr_II
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cr_III_avg]
    type = ElementAverageValue
    variable = Cr_III
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cr_I_avg]
    type = ElementAverageValue
    variable = Cr_I
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cr_total]
    type = ParsedPostprocessor
    expression = 'Cr_II_avg + Cr_III_avg + Cr_I_avg'
    pp_names = 'Cr_II_avg Cr_III_avg Cr_I_avg'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cr_III_fraction]
    type = ParsedPostprocessor
    expression = 'Cr_III_avg / Cr_total'
    pp_names = 'Cr_III_avg Cr_total'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cr_II_wall]
    type = PointValue
    variable = Cr_II
    point = '${fparse L / 40} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_eq_ratio]
    type = ParsedPostprocessor
    expression = 'Cr_II_wall / ${Cr_init}'
    pp_names = 'Cr_II_wall'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Wall loss rate from the dissolved-Cr growth [um/y]: (V/A) dc/dt M / rho, V/A = L/2
  [dCr_total]
    type = ChangeOverTimePostprocessor
    postprocessor = Cr_total
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [dt_pp]
    type = TimestepSize
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [corrosion_rate_um_y]
    type = ParsedPostprocessor
    expression = 'if(dt_pp > 0, dCr_total / dt_pp * ${fparse L / 2} * 51.9961 / 8.00e6 * 1e6 '
                 '* ${fparse 365.25 * day}, 0)'
    pp_names = 'dCr_total dt_pp'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cl2m_rad_avg]
    type = ElementAverageValue
    variable = Cl2m_rad
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cl2_avg]
    type = ElementAverageValue
    variable = Cl2_diss
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 15
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys '
                 'Cr_II_sys Cr_III_sys Cr_I_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  end_time = ${fparse 11 * day}
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
    growth_factor = 1.3
  []
  dtmax = 1800
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
