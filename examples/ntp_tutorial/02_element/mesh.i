# Stage 2: the mesh of one Ludewig Design 1 fuel element, written to mesh_in.e.
#
#   open_pronghorn-opt -i mesh.i --mesh-only --allow-unused
#   open_pronghorn-opt -i mesh.i --mesh-only --allow-unused nz=48 n_bed=8   (a finer mesh)
#
# (--allow-unused: the geometry file also defines the material parameters of the physics.)

!include ../common/ludewig_design1.i

nz = 12
n_ch = 1
n_hf = 1
n_bed = 1
n_cf = 1
n_pl = 1

!include ../common/element_mesh.i
