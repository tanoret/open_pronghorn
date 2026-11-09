# The 3-D half-element mesh, written to mesh_3d_in.e:
#   open_pronghorn-opt -i mesh_3d.i --mesh-only --allow-unused
!include ../common/ludewig_design1.i
nz = 12
n_ch = 2
n_hf = 1
n_bed = 2
n_cf = 1
n_pl = 1
n_theta = 12
r_axis = 4e-4
!include element_mesh_3d.i
