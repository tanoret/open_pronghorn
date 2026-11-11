# Three-dimensional mesh of one particle-bed fuel element: half of the element, 0 <= theta <= pi,
# with symmetry planes at theta = 0 and theta = pi.
#
# The reference strip of common/element_mesh.i (x = region index, y = axial coordinate z) is
# extruded in a third coordinate t in [0, 1] (n_theta layers), and the nodes are mapped onto the
# cones in cylindrical coordinates about the element axis y:
#
#   theta = pi t,   x = r(x_strip, y) cos(theta),   z = r(x_strip, y) sin(theta).
#
# The axis of the element, where the cells would collapse, is replaced by a thin cylinder of
# radius r_axis (boundary 'axis', a symmetry wall); its 0.4 mm radius is 2 % of the chamber
# radius at the inlet end. Parameters the including input must set: those of element_mesh.i, plus
#   n_theta   azimuthal cells over the half element
#   r_axis    radius of the axis cylinder [m]
#
# Boundaries: inlet, outlet, closed_ends, axis, outer_wall (as in 2-D), plus symmetry_0 and
# symmetry_pi (the two cut planes).

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
  [ends_only]
    type = BoundaryDeletionGenerator
    input = closed_ends
    boundary_names = 'bottom top'
  []
  [extrude]
    type = AdvancedExtruderGenerator
    input = ends_only
    direction = '0 0 1'
    heights = 1
    num_layers = ${n_theta}
    bottom_boundary = 100
    top_boundary = 101
  []
  [cut_planes]
    type = RenameBoundaryGenerator
    input = extrude
    old_boundary = '100 101 left right'
    new_boundary = 'symmetry_0 symmetry_pi axis outer_wall'
  []
  [cone]
    type = ParsedNodeTransformGenerator
    input = cut_planes
    constant_names = 'L r0 r1a r1b r2a r2b r3a r3b r4a r4b r5a r5b pi'
    constant_expressions = '${L} ${r_axis} ${r1_in} ${r1_out} ${r2_in} ${r2_out} ${r3_in} ${r3_out}
                            ${r4_in} ${r4_out} ${r5_in} ${r5_out} 3.141592653589793'
    x_function = 'f := y / L;
                  r1 := r1a + (r1b - r1a) * f; r2 := r2a + (r2b - r2a) * f;
                  r3 := r3a + (r3b - r3a) * f; r4 := r4a + (r4b - r4a) * f;
                  r5 := r5a + (r5b - r5a) * f;
                  r := if(x < 1, r0 + x * (r1 - r0),
                       if(x < 2, r1 + (x - 1) * (r2 - r1),
                       if(x < 3, r2 + (x - 2) * (r3 - r2),
                       if(x < 4, r3 + (x - 3) * (r4 - r3),
                                 r4 + (x - 4) * (r5 - r4)))));
                  r * cos(pi * z)'
    y_function = 'y'
    z_function = 'f := y / L;
                  r1 := r1a + (r1b - r1a) * f; r2 := r2a + (r2b - r2a) * f;
                  r3 := r3a + (r3b - r3a) * f; r4 := r4a + (r4b - r4a) * f;
                  r5 := r5a + (r5b - r5a) * f;
                  r := if(x < 1, r0 + x * (r1 - r0),
                       if(x < 2, r1 + (x - 1) * (r2 - r1),
                       if(x < 3, r2 + (x - 2) * (r3 - r2),
                       if(x < 4, r3 + (x - 3) * (r4 - r3),
                                 r4 + (x - 4) * (r5 - r4)))));
                  r * sin(pi * z)'
  []
  [names]
    type = RenameBlockGenerator
    input = cone
    old_block = '0 1 2 3 4'
    new_block = 'chamber hot_frit fuel_bed cold_frit plenum'
  []
  # the porous interfaces of the free channels (boundaries of the channel turbulence)
  [plenum_frit]
    type = SideSetsBetweenSubdomainsGenerator
    input = names
    primary_block = plenum
    paired_block = cold_frit
    new_boundary = plenum_frit
  []
  [chamber_frit]
    type = SideSetsBetweenSubdomainsGenerator
    input = plenum_frit
    primary_block = chamber
    paired_block = hot_frit
    new_boundary = chamber_frit
  []
[]
