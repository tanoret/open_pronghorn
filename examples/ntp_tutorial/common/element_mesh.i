# Axisymmetric (r, z) mesh of one conical particle-bed fuel element.
#
# x is the radius, y is the axial coordinate z. z = 0 is the narrow inlet end,
# where hydrogen enters the outer plenum; the chamber discharges at z = L.
#
#   r0 = 0 | chamber | r1 | hot frit | r2 | fuel bed | r3 | cold frit | r4 | plenum | r5
#
# Each interface radius varies linearly from r?_in (z = 0) to r?_out (z = L).
# The mesh is first built on a reference strip, region i occupying x in [i, i+1],
# and the nodes are then mapped onto the cone. The parameters below must be set
# by the including input:
#
#   L                         element length [m]
#   r1_in ... r5_in           interface radii at z = 0 [m]
#   r1_out ... r5_out         interface radii at z = L [m]
#   nz                        axial cells
#   n_ch n_hf n_bed n_cf n_pl radial cells per region
#
# Boundaries: inlet (plenum, z = 0), outlet (chamber, z = L), closed_ends (the rest of both ends),
# axis (r = 0) and outer_wall (r = r5). Every boundary face carries exactly one of these names.

[Mesh]
  [strip]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1 1 1 1 1'
    ix = '${n_ch} ${n_hf} ${n_bed} ${n_cf} ${n_pl}'
    dy = '${L}'
    iy = '${nz}'
    subdomain_id = '0 1 2 3 4'
  []
  # Boundaries are cut on the reference strip, where the region limits are fixed:
  # inlet = plenum at z = 0, outlet = chamber at z = L, closed_ends = the rest of both ends.
  [inlet]
    type = ParsedGenerateSideset
    input = strip
    combinatorial_geometry = 'y < 1e-9 & x > 4'
    new_sideset_name = inlet
  []
  [outlet]
    type = ParsedGenerateSideset
    input = inlet
    combinatorial_geometry = 'y > ${L} - 1e-9 & x < 1'
    new_sideset_name = outlet
  []
  [closed_ends]
    type = ParsedGenerateSideset
    input = outlet
    combinatorial_geometry = '(y < 1e-9 & x < 4) | (y > ${L} - 1e-9 & x > 1)'
    new_sideset_name = closed_ends
  []
  [cone]
    type = ParsedNodeTransformGenerator
    input = closed_ends
    constant_names = 'L r1a r1b r2a r2b r3a r3b r4a r4b r5a r5b'
    constant_expressions = '${L} ${r1_in} ${r1_out} ${r2_in} ${r2_out} ${r3_in} ${r3_out}
                            ${r4_in} ${r4_out} ${r5_in} ${r5_out}'
    # f = y/L is the axial fraction; within region i the node keeps its fraction x - i.
    x_function = 'f := y / L;
                  r1 := r1a + (r1b - r1a) * f; r2 := r2a + (r2b - r2a) * f;
                  r3 := r3a + (r3b - r3a) * f; r4 := r4a + (r4b - r4a) * f;
                  r5 := r5a + (r5b - r5a) * f;
                  if(x < 1, x * r1,
                  if(x < 2, r1 + (x - 1) * (r2 - r1),
                  if(x < 3, r2 + (x - 2) * (r3 - r2),
                  if(x < 4, r3 + (x - 3) * (r4 - r3),
                            r4 + (x - 4) * (r5 - r4)))))'
  []
  [names]
    type = RenameBlockGenerator
    input = cone
    old_block = '0 1 2 3 4'
    new_block = 'chamber hot_frit fuel_bed cold_frit plenum'
  []
  [axis]
    type = RenameBoundaryGenerator
    input = names
    old_boundary = 'left right'
    new_boundary = 'axis outer_wall'
  []
  # the ends are covered by inlet, outlet and closed_ends
  [single_names]
    type = BoundaryDeletionGenerator
    input = axis
    boundary_names = 'bottom top'
  []
  coord_type = RZ
  rz_coord_axis = Y
[]
