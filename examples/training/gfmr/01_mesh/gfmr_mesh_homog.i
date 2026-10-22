################################################################################
# GFMR training -- Stage 1c: homogenized 3-D full-core mesh for neutronics and transients
#
# One hexagonal column per 14 cm unit cell (6 triangular sectors, salt + graphite smeared),
# with the 9 control-rod cells kept as their own blocks so that the rod can be inserted
# cell by cell, the 18 solid-graphite plugs, and the cylindrical reflector to R = 1.75 m.
#
# Blocks:  lattice_inner (rings 0-5), lattice_outer (rings 6-11), lattice_plug (no channel),
#          cr_cell_withdrawn / cr_cell_rodded (the 9 rod cells below / above the rod tip),
#          reflector
# Axial:   0-3.24 m (36 layers of 9 cm), 3.24-3.60 m rodded (4 layers of 9 cm)
#
# Run:  open_pronghorn-opt -i gfmr_mesh_homog.i --mesh-only gfmr_homog_3d.e
################################################################################

pitch = 0.14
r_reflector = 1.75
height = 3.6
cr_insertion = 0.36

[Mesh]
  [inner]
    type = SimpleHexagonGenerator
    hexagon_size = ${fparse pitch / 2}
    hexagon_size_style = apothem
    block_id = 1
    block_name = lattice_inner
  []
  [outer]
    type = SimpleHexagonGenerator
    hexagon_size = ${fparse pitch / 2}
    hexagon_size_style = apothem
    block_id = 2
    block_name = lattice_outer
  []
  [cr]
    type = SimpleHexagonGenerator
    hexagon_size = ${fparse pitch / 2}
    hexagon_size_style = apothem
    block_id = 3
    block_name = cr_cell_withdrawn
  []
  [plug]
    type = SimpleHexagonGenerator
    hexagon_size = ${fparse pitch / 2}
    hexagon_size_style = apothem
    block_id = 4
    block_name = lattice_plug
  []
  # 0 inner, 1 outer, 2 control rod cell, 3 graphite plug (same map as gfmr_mesh_2d.i)
  [core]
    type = PatternedHexMeshGenerator
    inputs = 'inner outer cr plug'
    pattern = '3 3 1 1 1 1 1 1 1 1 3 3;
              3 1 1 1 1 1 1 1 1 1 1 1 3;
             1 1 1 1 1 1 1 1 1 1 1 1 1 1;
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
           1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
         1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1;
        1 1 1 1 1 1 0 2 0 0 0 2 0 1 1 1 1 1 1;
       1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1;
      1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1;
     3 1 1 1 1 1 0 0 0 0 0 2 0 0 0 0 1 1 1 1 1 3;
    3 1 1 1 1 1 0 2 0 0 2 0 0 0 0 2 0 1 1 1 1 1 3;
     3 1 1 1 1 1 0 0 0 0 0 2 0 0 0 0 1 1 1 1 1 3;
      1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1;
       1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1;
        1 1 1 1 1 1 0 2 0 0 0 2 0 1 1 1 1 1 1;
         1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1;
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
           1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
             1 1 1 1 1 1 1 1 1 1 1 1 1 1;
              3 1 1 1 1 1 1 1 1 1 1 1 3;
               3 3 1 1 1 1 1 1 1 1 3 3'
    pattern_boundary = none
    external_boundary_id = 10000
    external_boundary_name = 'lattice_edge'
  []
  [reflector]
    type = PeripheralRingMeshGenerator
    input = core
    input_mesh_external_boundary = 'lattice_edge'
    peripheral_ring_radius = ${r_reflector}
    peripheral_layer_num = 4
    peripheral_ring_block_id = 40
    peripheral_ring_block_name = 'reflector'
    external_boundary_name = 'side'
  []
  [extrude]
    type = AdvancedExtruderGenerator
    input = reflector
    direction = '0 0 1'
    heights = '${fparse height - cr_insertion} ${cr_insertion}'
    num_layers = '36 4'
    # upper interval: the rod cells are rodded
    subdomain_swaps = '; 3 5'
    bottom_boundary = 100
    top_boundary = 101
  []
  [names]
    type = RenameBlockGenerator
    input = extrude
    old_block = '5'
    new_block = 'cr_cell_rodded'
  []
  [ends]
    type = RenameBoundaryGenerator
    input = names
    old_boundary = '100 101'
    new_boundary = 'bottom top'
  []
[]
