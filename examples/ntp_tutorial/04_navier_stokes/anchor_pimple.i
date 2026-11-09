# Stage 4: the anchor of stage 3 with the segregated Navier-Stokes model (inertia and wall friction
# of the plenum and the chamber), marched to the design power with PIMPLE.
#
#   open_pronghorn-opt -i anchor_pimple.i orifice_flow.i     orifice designed on the slice flows
#   open_pronghorn-opt -i anchor_pimple.i                    orifice of the Newton anchor (stage 3)
#
# Without orifice_flow.i the cold-frit multiplier is the profile designed by the Newton anchor
# (stage 3, orifice_implicit.i), one value per axial slice (orifice_profile.i, written by
# extract_orifice.py); it does not fit this element (README.md). The power is ramped in pseudo-time
# from s_start to 1; each step converges the flow, the energy and the lagged state a little further
# (element_pimple.i). The run ends at the design power once the solution is steady.

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

# power ramp in pseudo-time, from s_start to s_final over t_ramp, then held
s_start = 0.02
s_final = 1
t_ramp = 20
t_end = 40
pimple_dt = 2e-3
pimple_dtmax = 0.05
pimple_end = ${t_end}
pimple_steady_state_detection = true
rho_cp_solid = 5e4

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
!include ../common/element_pimple.i
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
  [power_fraction]
    type = ParsedFunction
    expression = '${s_start} + (${s_final} - ${s_start}) * min(1, t / ${t_ramp})'
  []
  [q_fission]
    type = ParsedFunction
    symbol_names = 'shape s'
    symbol_values = 'shape power_fraction'
    expression = '4.0e10 * shape * s'
  []
[]

[Postprocessors]
  [s_power]
    type = FunctionValuePostprocessor
    function = power_fraction
  []
[]

[Executioner]
  steady_state_start_time = ${t_ramp}
[]

[Outputs]
  file_base = anchor_pimple_out
  csv = true
  exodus = true
[]
