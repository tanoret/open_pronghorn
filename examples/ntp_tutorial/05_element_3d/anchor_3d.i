# Stage 5: the 1 GW anchor of stage 3 in three dimensions, on half of the element
# (0 <= theta <= pi, symmetry planes at theta = 0 and pi), with the Newton Darcy-Forchheimer model.
#
#   open_pronghorn-opt -i anchor_3d.i ../03_ludewig_anchor/continuation.i \
#                         ../03_ludewig_anchor/orifice_implicit.i tolerance_3d.i   uniform power
#   ... tilt=0.05                                                              transverse power tilt
#
# The power density is the axial anchor shape times (1 + tilt x / r_bed(z)), r_bed the mean
# radius of the bed: an element next to the edge of the core sees a power gradient across it, which
# an axisymmetric model cannot represent. The tilt is symmetric about the plane z = 0, so the half
# element with symmetry planes holds the whole solution, and it does not change the element power.
# The mesh covers half of the element: the flows and the power reported are those of the half.

!include ../common/ludewig_design1.i

nz = 12
n_ch = 1
n_hf = 1
n_bed = 1
n_cf = 1
n_pl = 1
n_theta = 8
r_axis = 4e-4
tilt = 0

element_fraction = 0.5
mdot_bundle = 21.054294673192363
mdot = ${fparse mdot_bundle / n_elements * element_fraction}
T_in = 30.0
p_out = 7.0e6
T_solid_init = 30.0
h_in = ${fparse 8.314462618 / 2.01588e-3 * (2.34433112*T_in + 7.98052075e-03/2*T_in^2 - 1.94781510e-05/3*T_in^3 + 2.01572094e-08/4*T_in^4 - 7.37611761e-12/5*T_in^5 - 9.17935173e02)}
q_avg = 4.0e10
ramp_time = 2e-3
power_density = q_fission
cold_frit_multiplier = 1
# r_bed(z) = 0.5 (r2 + r3), linear in z between the inlet and the outlet ends
power_shape_expression = 'axial_shape * (1 + ${tilt} * x / (${fparse 0.5 * (r2_in + r3_in)}
                          + ${fparse 0.5 * (r2_out + r3_out - r2_in - r3_in) / L} * y))'

!include element_mesh_3d.i
!include ../03_ludewig_anchor/anchor_body.i
