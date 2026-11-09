# Steady power continuation for anchor.i:
#
#   open_pronghorn-opt -i anchor.i continuation.i [orifice_implicit.i]
#
# The time variable t is a continuation parameter, not physical time; every step is a steady
# Newton solve that starts from the previous one, and the time derivatives are switched off.
#   t in [-0.01, 0]       cold flow, no power
#   t in [0, t_design]    power held at s_design of nominal, so that the orifice design (if
#                         orifice_implicit.i is added) is activated where the element is stable
#   t in [t_design, 1]    power ramped linearly from s_design to 1
s_design = 0.05
t_design = 0.4

FVKernels/mass_time/enable = false
FVKernels/gas_time/enable = false
FVKernels/solid_time/enable = false

[Functions]
  [power_fraction]
    type = ParsedFunction
    expression = 'if(t < 0, 0,
                  if(t < ${t_design}, ${s_design} * min(1, t / 0.05),
                  ${s_design} + (1 - ${s_design}) * (t - ${t_design}) / (1 - ${t_design})))'
  []
  [q_fission]
    symbol_names = 'shape s'
    symbol_values = 'shape power_fraction'
    expression = '${q_avg} * shape * s'
  []
[]

[Postprocessors]
  # named differently from the function: a parsed expression would read a postprocessor of the
  # same name instead of the function
  [s_power]
    type = FunctionValuePostprocessor
    function = power_fraction
  []
[]

[Executioner]
  type = Transient
  start_time = -0.01
  end_time = 1
  nl_max_its = 40
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  line_search = bt
  steady_state_detection = false
  # at most 0.01 per step: the design phase then gives the orifice update about 40 iterations
  dtmax = 0.01
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.25
    cutback_factor = 0.5
    optimal_iterations = 10
  []
[]
