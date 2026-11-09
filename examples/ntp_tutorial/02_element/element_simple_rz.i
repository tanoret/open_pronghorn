# Stage 2: one fuel element, axisymmetric, with the segregated SIMPLE solver.
#
#   open_pronghorn-opt -i element_simple_rz.i                       cold flow, uniform frits
#   open_pronghorn-opt -i element_simple_rz.i s_power=0.05          5 % of the design power
#   open_pronghorn-opt -i element_simple_rz.i nz=24 n_pl=2 ...      a finer mesh
#
# s_power scales the design power density of the anchor (stage 3). Without power the energy
# equations are not solved (they would be at the inlet temperature).

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
s_power = 0
cold_frit_multiplier = 1
power_density = q_fission
simple_iterations = 300
simple_fixed_point = 30

segregated_momentum_systems = 'u_system v_system'
segregated_velocities = 'vel_x vel_y'
segregated_inlet_velocity = '0 u_inlet'
segregated_inlet_fraction = 1
segregated_walls = 'axis outer_wall closed_ends'
segregated_wall_types = 'symmetry noslip noslip'
segregated_wall_energy_types = 'heatflux heatflux heatflux'
segregated_wall_energy_functors = '0 0 0'
segregated_initial_velocity = '0 0 0'

!include ../common/element_mesh.i
!include ../common/hydrogen_nasa.i
!include ../common/element_physics_segregated.i
!include ../common/element_segregated_rz.i
!include ../common/element_simple.i
!include ../common/element_segregated_postprocessors.i

[Functions]
  # the inlet-peaked axial power shape of the anchor, 12 equal slices, mean 1
  [shape]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '0 0.06078333 0.12156667 0.18235 0.24313333 0.30391667 0.3647 0.42548333 0.48626667
         0.54705 0.60783333 0.66861667'
    y = '1.28 1.22 1.15 1.09 1.03 0.98 0.94 0.91 0.88 0.86 0.84 0.82'
  []
  [q_fission]
    type = ParsedFunction
    symbol_names = 'shape'
    symbol_values = 'shape'
    expression = '4.0e10 * ${s_power} * shape'
  []
[]

[Executioner]
  should_solve_energy = ${fparse s_power > 0}
  should_solve_solid_energy = ${fparse s_power > 0}
[]

[Outputs]
  csv = true
  exodus = true
[]
