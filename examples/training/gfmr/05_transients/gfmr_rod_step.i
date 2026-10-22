################################################################################
# GFMR training -- Stage 5: reactivity-insertion transient, point kinetics with
# circulating-fuel precursors and two-temperature feedback
#
# A control-rod step of rho_step = 0.3 $ (204 pcm, about 4 % of the total rod worth) at t = 1 s,
# from the nominal steady state, without scram. Solved as a MOOSE system of scalar ODEs (implicit
# BDF2):
#
#   n'      = (rho - beta_eff + rho_0) / Lambda * n + sum_i lambda_i c_i                 (power)
#   c_i'    = beta_i / Lambda * n - lambda_i c_i - (c_i - c_i,in) / t_core               (core)
#   c_i,in  = c_i(t - t_ex) exp(-lambda_i t_ex), lagged here by a first-order loop       (loop)
#             d_i' = (c_i - d_i) / t_ex - lambda_i d_i,  c_i,in = d_i
#   M_s c_s T_s' = (1 - f_g) P - mdot c_s (T_out - T_in) + UA (T_g - T_s)               (salt)
#   M_g c_g T_g' = f_g P - UA (T_g - T_s)                                                (graphite)
#   rho     = rho_step H(t - 1) + alpha_f (T_s - T_s0) + alpha_g (T_g - T_g0)
#
# The precursors that leave the core decay partly in the loop (t_ex = 6.0 s); the reactivity
# this costs at steady state is rho_0 (the drift reactivity), which the fission source must
# supply at t = 0 so that the initial state is critical. The heat exchanger holds T_in fixed (no
# loop thermal feedback: the most conservative boundary for a fast insertion).
#
# Data (all from the earlier stages):
#   beta_i, lambda_i, Lambda   OpenMC IFP and the delayed-neutron data (../02_neutronics/xs)
#   alpha_f = -2.7 pcm/K, alpha_g = -4.1 pcm/K   OpenMC +100 K perturbations
#   f_g = 7.87 %               OpenMC coupled n-gamma heating
#   UA, M_g c_g                stage-3 CFD (7.87 % of 150 MW crosses a 71 K graphite-salt
#                              difference) and the lattice graphite mass
#   t_core = 3.02 s            CFD core salt volume / flow;  t_ex = 6.04 s ASSUMED (stage 4)
#
# Run: open_pronghorn-opt -i gfmr_rod_step.i      (a few seconds, one rank)
################################################################################

P0 = 150e6
rho_step = ${fparse 0.3 * 0.00681035}

# kinetics (OpenMC IFP, ENDF/B-VIII.0 delayed data)
Lambda = 6.0511e-4
b1 = 2.44813e-4
b2 = 1.238017e-3
b3 = 1.147454e-3
b4 = 2.647081e-3
b5 = 1.081340e-3
b6 = 4.51647e-4
beta = ${fparse b1 + b2 + b3 + b4 + b5 + b6}
l1 = 0.0133410
l2 = 0.0327022
l3 = 0.120861
l4 = 0.303643
l5 = 0.853140
l6 = 2.86528

# circulation
t_core = 3.02036
t_ex = 6.04071

# feedback (OpenMC), 1/K
a_f = -2.7e-5
a_g = -4.1e-5

# thermal model (stage 3)
mdot = 2294.13
cp_s = 900.0
M_s = 6944.0
f_g = 0.0787
UA = 1.663e5
MC_g = 7.28e7
T_in = 613.8
# initial steady state of this model: the salt carries all 150 MW out (T_s0 is the core average,
# T_out = 2 T_s0 - T_in) and the graphite runs above it by its own heat over UA
T_s0 = ${fparse T_in + P0 / (2 * mdot * cp_s)}
T_g0 = ${fparse T_s0 + f_g * P0 / UA}

# steady precursor balance with drift: c_i = beta_i/Lambda n / (lambda_i + (1 - k_i)/t_core),
# k_i = d_i/c_i = (1/t_ex) / (1/t_ex + lambda_i) the fraction that survives the loop
k1 = ${fparse (1 / t_ex) / (1 / t_ex + l1)}
k2 = ${fparse (1 / t_ex) / (1 / t_ex + l2)}
k3 = ${fparse (1 / t_ex) / (1 / t_ex + l3)}
k4 = ${fparse (1 / t_ex) / (1 / t_ex + l4)}
k5 = ${fparse (1 / t_ex) / (1 / t_ex + l5)}
k6 = ${fparse (1 / t_ex) / (1 / t_ex + l6)}
c10 = ${fparse b1 / Lambda / (l1 + (1 - k1) / t_core)}
c20 = ${fparse b2 / Lambda / (l2 + (1 - k2) / t_core)}
c30 = ${fparse b3 / Lambda / (l3 + (1 - k3) / t_core)}
c40 = ${fparse b4 / Lambda / (l4 + (1 - k4) / t_core)}
c50 = ${fparse b5 / Lambda / (l5 + (1 - k5) / t_core)}
c60 = ${fparse b6 / Lambda / (l6 + (1 - k6) / t_core)}
# drift reactivity: the reactivity the steady state needs to stay critical with drift
rho_0 = ${fparse beta - Lambda * (l1 * c10 + l2 * c20 + l3 * c30 + l4 * c40 + l5 * c50 + l6 * c60)}

[Mesh]
  [g]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Variables]
  # relative power n = P / P0
  [n]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  []
  [c1]
    family = SCALAR
    order = FIRST
    initial_condition = ${c10}
  []
  [c2]
    family = SCALAR
    order = FIRST
    initial_condition = ${c20}
  []
  [c3]
    family = SCALAR
    order = FIRST
    initial_condition = ${c30}
  []
  [c4]
    family = SCALAR
    order = FIRST
    initial_condition = ${c40}
  []
  [c5]
    family = SCALAR
    order = FIRST
    initial_condition = ${c50}
  []
  [c6]
    family = SCALAR
    order = FIRST
    initial_condition = ${c60}
  []
  [d1]
    family = SCALAR
    order = FIRST
    initial_condition = ${fparse k1 * c10}
  []
  [d2]
    family = SCALAR
    order = FIRST
    initial_condition = ${fparse k2 * c20}
  []
  [d3]
    family = SCALAR
    order = FIRST
    initial_condition = ${fparse k3 * c30}
  []
  [d4]
    family = SCALAR
    order = FIRST
    initial_condition = ${fparse k4 * c40}
  []
  [d5]
    family = SCALAR
    order = FIRST
    initial_condition = ${fparse k5 * c50}
  []
  [d6]
    family = SCALAR
    order = FIRST
    initial_condition = ${fparse k6 * c60}
  []
  # core-average salt and graphite temperatures [degC]
  [Ts]
    family = SCALAR
    order = FIRST
    initial_condition = ${T_s0}
  []
  [Tg]
    family = SCALAR
    order = FIRST
    initial_condition = ${T_g0}
  []
[]

[Functions]
  [rho_ext]
    type = ParsedFunction
    expression = 'if(t >= 1, ${rho_step}, 0)'
  []
[]

[Postprocessors]
  [rho_ext]
    type = FunctionValuePostprocessor
    function = rho_ext
    execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
  []
  [power_MW]
    type = ScalarVariable
    variable = n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_salt]
    type = ScalarVariable
    variable = Ts
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_graphite]
    type = ScalarVariable
    variable = Tg
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [rho_feedback_pcm]
    type = ParsedPostprocessor
    expression = '1e5 * (${a_f} * (T_salt - ${T_s0}) + ${a_g} * (T_graphite - ${T_g0}))'
    pp_names = 'T_salt T_graphite'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [rho_net_dollar]
    type = ParsedPostprocessor
    expression = '(rho_ext + rho_feedback_pcm / 1e5) / ${beta}'
    pp_names = 'rho_ext rho_feedback_pcm'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_outlet]
    type = ParsedPostprocessor
    expression = '2 * T_salt - ${T_in}'
    pp_names = 'T_salt'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [drift_reactivity_pcm]
    type = Receiver
    default = ${fparse rho_0 * 1e5}
  []
[]

[ScalarKernels]
  # ---- power: n' = ((rho - beta + rho_0) n + Lambda sum lambda_i c_i) / Lambda
  [n_dt]
    type = ODETimeDerivative
    variable = n
  []
  [n_rhs]
    type = ParsedODEKernel
    variable = n
    coupled_variables = 'c1 c2 c3 c4 c5 c6 Ts Tg'
    postprocessors = 'rho_ext'
    expression = '-((rho_ext + ${a_f} * (Ts - ${T_s0}) + ${a_g} * (Tg - ${T_g0}) - ${beta} + ${rho_0}) / ${Lambda} * n
                   + ${l1} * c1 + ${l2} * c2 + ${l3} * c3 + ${l4} * c4 + ${l5} * c5 + ${l6} * c6)'
  []
  # ---- core precursors: production, decay, and exchange with the loop
  [c1_dt]
    type = ODETimeDerivative
    variable = c1
  []
  [c1_rhs]
    type = ParsedODEKernel
    variable = c1
    coupled_variables = 'n d1'
    expression = '-(${b1} / ${Lambda} * n - ${l1} * c1 - (c1 - d1) / ${t_core})'
  []
  [c2_dt]
    type = ODETimeDerivative
    variable = c2
  []
  [c2_rhs]
    type = ParsedODEKernel
    variable = c2
    coupled_variables = 'n d2'
    expression = '-(${b2} / ${Lambda} * n - ${l2} * c2 - (c2 - d2) / ${t_core})'
  []
  [c3_dt]
    type = ODETimeDerivative
    variable = c3
  []
  [c3_rhs]
    type = ParsedODEKernel
    variable = c3
    coupled_variables = 'n d3'
    expression = '-(${b3} / ${Lambda} * n - ${l3} * c3 - (c3 - d3) / ${t_core})'
  []
  [c4_dt]
    type = ODETimeDerivative
    variable = c4
  []
  [c4_rhs]
    type = ParsedODEKernel
    variable = c4
    coupled_variables = 'n d4'
    expression = '-(${b4} / ${Lambda} * n - ${l4} * c4 - (c4 - d4) / ${t_core})'
  []
  [c5_dt]
    type = ODETimeDerivative
    variable = c5
  []
  [c5_rhs]
    type = ParsedODEKernel
    variable = c5
    coupled_variables = 'n d5'
    expression = '-(${b5} / ${Lambda} * n - ${l5} * c5 - (c5 - d5) / ${t_core})'
  []
  [c6_dt]
    type = ODETimeDerivative
    variable = c6
  []
  [c6_rhs]
    type = ParsedODEKernel
    variable = c6
    coupled_variables = 'n d6'
    expression = '-(${b6} / ${Lambda} * n - ${l6} * c6 - (c6 - d6) / ${t_core})'
  []
  # ---- loop precursors (first-order lag of transit t_ex), decaying on the way
  [d1_dt]
    type = ODETimeDerivative
    variable = d1
  []
  [d1_rhs]
    type = ParsedODEKernel
    variable = d1
    coupled_variables = 'c1'
    expression = '-((c1 - d1) / ${t_ex} - ${l1} * d1)'
  []
  [d2_dt]
    type = ODETimeDerivative
    variable = d2
  []
  [d2_rhs]
    type = ParsedODEKernel
    variable = d2
    coupled_variables = 'c2'
    expression = '-((c2 - d2) / ${t_ex} - ${l2} * d2)'
  []
  [d3_dt]
    type = ODETimeDerivative
    variable = d3
  []
  [d3_rhs]
    type = ParsedODEKernel
    variable = d3
    coupled_variables = 'c3'
    expression = '-((c3 - d3) / ${t_ex} - ${l3} * d3)'
  []
  [d4_dt]
    type = ODETimeDerivative
    variable = d4
  []
  [d4_rhs]
    type = ParsedODEKernel
    variable = d4
    coupled_variables = 'c4'
    expression = '-((c4 - d4) / ${t_ex} - ${l4} * d4)'
  []
  [d5_dt]
    type = ODETimeDerivative
    variable = d5
  []
  [d5_rhs]
    type = ParsedODEKernel
    variable = d5
    coupled_variables = 'c5'
    expression = '-((c5 - d5) / ${t_ex} - ${l5} * d5)'
  []
  [d6_dt]
    type = ODETimeDerivative
    variable = d6
  []
  [d6_rhs]
    type = ParsedODEKernel
    variable = d6
    coupled_variables = 'c6'
    expression = '-((c6 - d6) / ${t_ex} - ${l6} * d6)'
  []
  # ---- salt: fission heat, advection out (T_out = 2 T_s - T_in), heat from the graphite
  [Ts_dt]
    type = ODETimeDerivative
    variable = Ts
  []
  [Ts_rhs]
    type = ParsedODEKernel
    variable = Ts
    coupled_variables = 'n Tg'
    expression = '-(((1 - ${f_g}) * ${P0} * n - ${mdot} * ${cp_s} * 2 * (Ts - ${T_in})
                    + ${UA} * (Tg - Ts)) / (${M_s} * ${cp_s}))'
  []
  # ---- graphite: gamma / neutron heating, conduction to the salt
  [Tg_dt]
    type = ODETimeDerivative
    variable = Tg
  []
  [Tg_rhs]
    type = ParsedODEKernel
    variable = Tg
    coupled_variables = 'n Ts'
    expression = '-((${f_g} * ${P0} * n - ${UA} * (Tg - Ts)) / ${MC_g})'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  start_time = 0
  end_time = 600
  dtmax = 5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    growth_factor = 1.2
    optimal_iterations = 5
  []
  nl_abs_tol = 1e-9
  nl_rel_tol = 1e-9
  line_search = none
[]

[Outputs]
  csv = true
[]
