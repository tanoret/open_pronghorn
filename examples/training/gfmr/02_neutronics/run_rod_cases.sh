#!/bin/bash
# Sector diffusion at the three rod positions (nominal 36 cm, all out, all in).
# Usage (inside moose-dev-exec): ./run_rod_cases.sh [ranks]
n=${1:-24}
exe=${EXE:-../../../../open_pronghorn-opt}
mpirun -n $n $exe -i gfmr_diffusion_sector.i -w Outputs/file_base=sector_nominal > sector_nominal.log 2>&1
mpirun -n $n $exe -i gfmr_diffusion_sector.i -w Outputs/file_base=sector_rods_out \
  mesh_file=../01_mesh/gfmr_sector_3d_rods_out.e \
  Materials/inactive='absorber absorber_R' > sector_rods_out.log 2>&1
mpirun -n $n $exe -i gfmr_diffusion_sector.i -w Outputs/file_base=sector_rods_in \
  mesh_file=../01_mesh/gfmr_sector_3d_rods_in.e helium=helium > sector_rods_in.log 2>&1
for f in sector_nominal sector_rods_out sector_rods_in; do echo "$f $(tail -1 $f.csv)"; done
