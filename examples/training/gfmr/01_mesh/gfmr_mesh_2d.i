################################################################################
# GFMR training -- Stage 1a: the 2-D core cross section with the Reactor-module mesh generators
#
# Geometry (Saltfoss GFMR reference specification, https://saltfossenergy.github.io/GFMR/):
#   * 14 cm hexagonal graphite unit cells, each with one 4 cm fuel channel;
#   * 379 channels: 82 inner (1.70 m/s), 288 outer (1.03 m/s), 9 control-rod channels
#     (Hf absorber r = 0.65 cm, He gap 0.1 cm, Hastelloy N guide tube 0.2 cm);
#   * 18 ring-11 corner positions without a channel (solid graphite);
#   * graphite reflector out to R = 1.75 m.
# The channel map is generated from ../common/gfmr.py by make_pattern.py (pattern below).
#
# Run:  open_pronghorn-opt -i gfmr_mesh_2d.i --mesh-only gfmr_2d.e
# Units: metres.
################################################################################

pitch = 0.14            # m, flat-to-flat of the unit cell
r_channel = 0.02        # m
r_abs = 0.0065          # m, Hf absorber
r_gap = 0.0075          # m, He gap outer radius
r_guide = 0.0095        # m, guide-tube outer radius
r_reflector = 1.75      # m

# Azimuthal resolution: 4 sectors per hexagon side = 24 cells around each channel
nsec = 4
# Near-wall salt layer (wall functions): a 3 mm ring next to the graphite puts the first cell
# centre at y = 1.5 mm, y+ = 26 (outer, 1.03 m/s) to 40 (inner, 1.70 m/s) at nominal flow
wall_layer = 0.003

[Mesh]
  # ---------------------------- unit cells ------------------------------
  # Fuel channel (inner / outer zone differ only in the block name, so that the CFD can prescribe
  # the zone mass flux and the chemistry can split the walls by zone).
  [cell_inner]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${nsec} ${nsec} ${nsec} ${nsec} ${nsec} ${nsec}'
    polygon_size = ${fparse pitch / 2}
    polygon_size_style = apothem
    ring_radii = '${fparse r_channel - wall_layer} ${r_channel}'
    ring_intervals = '3 1'
    ring_block_ids = '10 11 11'
    ring_block_names = 'salt_inner_tri salt_inner salt_inner'
    background_intervals = 3
    background_block_ids = '20'
    background_block_names = 'graphite'
    preserve_volumes = true
    flat_side_up = false
  []
  [cell_outer]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${nsec} ${nsec} ${nsec} ${nsec} ${nsec} ${nsec}'
    polygon_size = ${fparse pitch / 2}
    polygon_size_style = apothem
    ring_radii = '${fparse r_channel - wall_layer} ${r_channel}'
    ring_intervals = '3 1'
    ring_block_ids = '12 13 13'
    ring_block_names = 'salt_outer_tri salt_outer salt_outer'
    background_intervals = 3
    background_block_ids = '20'
    background_block_names = 'graphite'
    preserve_volumes = true
    flat_side_up = false
  []
  # Control-rod channel: absorber / helium gap / guide tube / salt annulus / graphite
  [cell_cr]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${nsec} ${nsec} ${nsec} ${nsec} ${nsec} ${nsec}'
    polygon_size = ${fparse pitch / 2}
    polygon_size_style = apothem
    ring_radii = '${r_abs} ${r_gap} ${r_guide} ${fparse r_guide + wall_layer} ${fparse r_channel - wall_layer} ${r_channel}'
    ring_intervals = '2 1 1 1 2 1'
    ring_block_ids = '30 31 32 33 14 14 14'
    ring_block_names = 'absorber_tri absorber helium guide_tube salt_cr salt_cr salt_cr'
    background_intervals = 3
    background_block_ids = '20'
    background_block_names = 'graphite'
    preserve_volumes = true
    flat_side_up = false
  []
  # Solid graphite position (no channel)
  [cell_graphite]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${nsec} ${nsec} ${nsec} ${nsec} ${nsec} ${nsec}'
    polygon_size = ${fparse pitch / 2}
    polygon_size_style = apothem
    ring_radii = '${r_channel}'
    ring_intervals = '3'
    ring_block_ids = '21 22'
    ring_block_names = 'graphite_plug_tri graphite_plug'
    background_intervals = 3
    background_block_ids = '20'
    background_block_names = 'graphite'
    preserve_volumes = true
    flat_side_up = false
  []

  # ---------------------------- core lattice -----------------------------
  # 0 inner channel, 1 outer channel, 2 control-rod channel, 3 solid graphite (see make_pattern.py)
  [core]
    type = PatternedHexMeshGenerator
    inputs = 'cell_inner cell_outer cell_cr cell_graphite'
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
    rotate_angle = 0
  []

  # ---------------------------- radial reflector --------------------------
  [reflector]
    type = PeripheralRingMeshGenerator
    input = core
    input_mesh_external_boundary = 'lattice_edge'
    peripheral_ring_radius = ${r_reflector}
    peripheral_layer_num = 6
    peripheral_radial_bias = 1.0
    peripheral_ring_block_id = 40
    peripheral_ring_block_name = 'reflector'
    external_boundary_name = 'reflector_outer'
  []
[]
