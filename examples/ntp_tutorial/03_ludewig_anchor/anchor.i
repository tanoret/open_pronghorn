# Stage 3: one Ludewig Design 1 fuel element at the 1 GW, 7 MPa, 3000 K anchor state.
#
# The element flow is the first-law flow of the bundle divided among 37 elements:
#   mdot_bundle = Q / (h(3000 K) - h(30 K)) = 21.054 kg/s  (ideal-gas NASA hydrogen)
# The axial split of that flow between slices is NOT imposed: it follows from the porous-flow
# solution and the cold-frit drag multiplier chi_cf(z) (1 everywhere by default).
#
#   open_pronghorn-opt -i anchor.i                     uniform frits, power ramped to 1 GW
#                                                      in time (stops at the flow excursion)
#   open_pronghorn-opt -i anchor.i continuation.i      the same as a steady power continuation
#   open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i
#                                                      the orifice designed with the flow: the anchor
#
# Mesh: nz axial cells and the radial cells per region can be changed on the command line,
# e.g. nz=24 n_bed=4. The orifice of orifice_implicit.i has one component per axial cell: with
# nz != 12 also set its order, e.g. nz=24 log_chi_order=TWENTYFOURTH.

!include ../common/ludewig_design1.i

nz = 12
n_ch = 1
n_hf = 1
n_bed = 1
n_cf = 1
n_pl = 1

mdot_bundle = 21.054294673192363
mdot = ${fparse mdot_bundle / n_elements}
T_in = 30.0
p_out = 7.0e6
T_solid_init = 30.0
# NASA enthalpy of hydrogen at T_in (low-temperature branch of common/hydrogen.py)
h_in = ${fparse 8.314462618 / 2.01588e-3 * (2.34433112*T_in + 7.98052075e-03/2*T_in^2 - 1.94781510e-05/3*T_in^3 + 2.01572094e-08/4*T_in^4 - 7.37611761e-12/5*T_in^5 - 9.17935173e02)}

# Average fuel-bed power density and the inlet-peaked shape of the reference anchor (12 equal
# axial intervals, normalized to an arithmetic mean of one).
q_avg = 4.0e10
ramp_time = 2e-3
power_density = q_fission
cold_frit_multiplier = 1
element_fraction = 1
power_shape_expression = 'axial_shape'

!include ../common/element_mesh.i
!include anchor_body.i
