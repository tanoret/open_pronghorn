# Stage 5: the stage-4 anchor of the segregated Navier-Stokes model in three dimensions, on half of
# the element (0 <= theta <= pi, symmetry planes at theta = 0 and pi), marched to the design power
# with PIMPLE and the orifice designed on the slice flows:
#
#   mpirun -n 8 open_pronghorn-opt -i anchor_pimple_3d.i ../04_navier_stokes/orifice_flow.i
#
# The 3-D counterpart of ../04_navier_stokes/anchor_pimple.i. The pseudo-transient also starts the
# flow: from rest, the first SIMPLE iterations of the 3-D element see only the gas viscosity in
# the free channels and diverge, while the time derivative bounds the momentum mobility by dt/rho.
#
# The mesh is faceted: the inlet flow is scaled by the exact over the meshed inlet area and the
# power density by the exact over the meshed bed volume, so the flow and the power of the half
# element are exact. The postprocessors report element totals.

!include ../common/ludewig_design1.i

nz = 12
n_ch = 1
n_hf = 1
n_bed = 1
n_cf = 1
n_pl = 1
n_theta = 8
r_axis = 4e-4

element_fraction = 0.5
# exact half-annulus inlet area over the faceted one
facet_area_ratio = ${fparse pi / (n_theta * sin(pi / n_theta))}
mdot = ${fparse 21.054294673192363 / n_elements * element_fraction * facet_area_ratio}
T_in = 30.0
p_out = 7.0e6
T_solid_init = 30.0
power_density = q_fission
cold_frit_multiplier = 1

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

segregated_momentum_systems = 'u_system v_system w_system'
segregated_velocities = 'vel_x vel_y vel_z'
segregated_inlet_velocity = '0 u_inlet 0'
segregated_inlet_fraction = ${element_fraction}
segregated_walls = 'axis outer_wall closed_ends symmetry_0 symmetry_pi'
segregated_wall_types = 'symmetry noslip noslip symmetry symmetry'
segregated_wall_energy_types = 'heatflux heatflux heatflux heatflux heatflux'
segregated_wall_energy_functors = '0 0 0 0 0'
segregated_initial_velocity = '0 0 0'

!include element_mesh_3d.i
!include ../common/hydrogen_nasa.i
!include ../common/element_physics_segregated.i
!include ../common/element_segregated_3d.i
!include ../common/element_pimple.i
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
  [power_fraction]
    type = ParsedFunction
    expression = '${s_start} + (${s_final} - ${s_start}) * min(1, t / ${t_ramp})'
  []
  [q_fission]
    type = ParsedFunction
    symbol_names = 'shape s volume_ratio'
    symbol_values = 'shape power_fraction volume_ratio'
    expression = '4.0e10 * shape * s * volume_ratio'
  []
[]

[Postprocessors]
  [s_power]
    type = FunctionValuePostprocessor
    function = power_fraction
  []
  # exact bed volume of the half element over the meshed one
  [volume_ratio]
    type = ParsedPostprocessor
    pp_names = bed_volume
    expression = '${fparse element_fraction * pi * L / 3 * (r3_in^2 + r3_in * r3_out + r3_out^2
                                                              - r2_in^2 - r2_in * r2_out - r2_out^2)} / bed_volume'
    execute_on = initial
  []
  # bed volume of the finite-volume cells, which the source kernel deposits the power in
  [bed_volume]
    type = NTPFVCellIntegral
    block = fuel_bed
    execute_on = initial
  []
[]

[Executioner]
  steady_state_start_time = ${t_ramp}
[]

[Outputs]
  file_base = anchor_pimple_3d_out
  csv = true
  exodus = true
[]
