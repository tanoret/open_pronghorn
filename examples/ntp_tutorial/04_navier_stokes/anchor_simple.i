# Stage 4: the anchor of stage 3 with the segregated Navier-Stokes model, solved to the design power
# with SIMPLE and a fixed-point continuation in the power.
#
#   open_pronghorn-opt -i anchor_simple.i orifice_flow.i     orifice designed on the slice flows
#
# Every fixed-point iteration is one SIMPLE solve (element_simple.i) of the flow, the gas energy
# and the particles at a fixed lagged density state, power and orifice. Between the solves the
# lagged state relaxes towards the solved one, the orifice design (orifice_flow.i) takes a step,
# and the power fraction steps from s_start to s_final over the first n_ramp iterations. It is
# then held for the remaining simple_fixed_point - n_ramp iterations, in which the state and the
# orifice converge. The CSV output has one line per fixed-point iteration (the column 'time' is
# always 1), the slice vectors at the end.
#
# Without orifice_flow.i the cold frit carries the orifice of the Newton anchor (orifice_profile.i),
# which starves the middle slices of this element (README.md).

!include ../common/ludewig_design1.i

nz = 12
n_ch = 1
n_hf = 1
n_bed = 1
n_cf = 1
n_pl = 1

mdot = ${fparse 21.054294673192363 / n_elements}
T_in = 30.0
p_out = 7.0e6
T_solid_init = 30.0
power_density = q_fission

# power continuation over the fixed-point iterations
s_start = 0.02
s_final = 1
n_ramp = 50
simple_fixed_point = 200
simple_iterations = 300

segregated_momentum_systems = 'u_system v_system'
segregated_velocities = 'vel_x vel_y'
segregated_inlet_velocity = '0 u_inlet'
segregated_inlet_fraction = 1
segregated_walls = 'axis outer_wall closed_ends'
segregated_wall_types = 'symmetry noslip noslip'
segregated_wall_energy_types = 'heatflux heatflux heatflux'
segregated_wall_energy_functors = '0 0 0'
segregated_initial_velocity = '0 0 0'

cold_frit_multiplier = chi_cf
!include orifice_profile.i

!include ../common/element_mesh.i
!include ../common/hydrogen_nasa.i
!include ../common/element_physics_segregated.i
!include ../common/element_segregated_rz.i
!include ../common/element_simple.i
!include ../common/element_segregated_postprocessors.i

[Functions]
  [shape]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '0 0.06078333 0.12156667 0.18235 0.24313333 0.30391667 0.3647 0.42548333 0.48626667
         0.54705 0.60783333 0.66861667'
    y = '1.28 1.22 1.15 1.09 1.03 0.98 0.94 0.91 0.88 0.86 0.84 0.82'
  []
  # power fraction of the solve after 'n_solves' fixed-point iterations
  [power_fraction]
    type = ParsedFunction
    symbol_names = 'k'
    symbol_values = 'n_solves'
    expression = '${s_start} + (${s_final} - ${s_start}) * min(1, k / ${n_ramp})'
  []
  [q_fission]
    type = ParsedFunction
    symbol_names = 'shape s'
    symbol_values = 'shape power_fraction'
    expression = '4.0e10 * shape * s'
  []
[]

[Postprocessors]
  [n_solves]
    type = NumFixedPointIterations
    execute_on = 'initial timestep_end'
  []
  [s_power]
    type = FunctionValuePostprocessor
    function = power_fraction
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  file_base = anchor_simple_out
  exodus = true
  # one line per fixed-point iteration (the continuation history); the slice vectors at the end
  [csv]
    type = CSV
    execute_on = 'multiapp_fixed_point_iteration_end'
    execute_vector_postprocessors_on = final
    new_row_detection_columns = all
  []
[]
