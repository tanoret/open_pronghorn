# Stage 3: is the designed orifice, once frozen, a stable state of the element in physical time?
#
#   open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i stability.i
#
# The run is the anchor (steady continuation, orifice designed with the flow) up to t = 1, the
# design point, held there by steady steps until t_freeze. From t_freeze the orifice is frozen (d lambda/dt = 0 instead of the design
# equations), the time derivatives of the gas mass, the gas energy and the particles are switched
# on, and the element is marched in physical time at the design power. For the first kick_time
# seconds the power of one axial slice (kick_slice) is raised by 'kick'. A stable element returns
# to the anchor; an unstable one (a flow excursion) runs away.
#
# The particle heat capacity of element_physics.i is a pseudo value that only shapes the
# continuation; rho_s and cp_s set a realistic one for the transient here.
kick = 0.10
kick_time = 0.1
kick_slice = 6
t_freeze = 1.05
t_transient = 1.5
rho_s = 1e4
cp_s = 250

[Functions]
  [q_fission]
    type = ParsedFunction
    symbol_names = 'shape s'
    symbol_values = 'shape power_fraction'
    expression = 'j := floor(y / ${fparse L / nz});
                  ${q_avg} * shape * s
                  * (1 + if(t > ${t_freeze} & t < ${fparse t_freeze + kick_time} & j = ${kick_slice}, ${kick}, 0))'
  []
  [power_fraction]
    type = ParsedFunction
    expression = 'if(t < 0, 0,
                  if(t < ${t_design}, ${s_design} * min(1, t / 0.05),
                  min(1, ${s_design} + (1 - ${s_design}) * (t - ${t_design}) / (1 - ${t_design}))))'
  []
[]

[ScalarKernels]
  # the orifice frozen once the transient starts
  [frozen_orifice]
    type = ODETimeDerivative
    variable = log_chi
    enable = false
  []
[]

[FunctorMaterials]
  [solid_inertia]
    prop_values = '${rho_s} ${cp_s}'
  []
[]

[Controls]
  [transient]
    type = TimePeriod
    enable_objects = 'ScalarKernels::frozen_orifice FVKernels::mass_time FVKernels::gas_time
                      FVKernels::solid_time'
    disable_objects = 'FVKernels::orifice_design'
    start_time = ${t_freeze}
    execute_on = 'initial timestep_begin'
  []
[]

[Executioner]
  end_time = ${fparse t_freeze + t_transient}
  # steps of at most 0.01 (as in continuation.i); in the transient they resolve the kick
  dtmax = 0.01
[]

[Outputs]
  file_base = stability_out
  exodus = false
[]
