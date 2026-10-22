################################################################################
# GFMR training -- Stage 1b: 3-D CFD / conjugate-heat-transfer mesh of a 60-degree core sector
#
# The lattice of gfmr_mesh_2d.i has the D3 symmetry of its control rods: the three central rods
# sit at 60, 180 and 300 degrees, so the core has mirror lines every 60 degrees (0, 60, 120 ...)
# but is NOT six-fold rotationally symmetric. The 60-degree sector between the mirror lines at 0
# and 60 degrees is therefore an exact representation of the full core at nominal conditions
# (all rods at the same height); a 30-degree sector is not. HexagonMeshTrimmer cuts the sector
# along those two lines; the cut faces become the "symmetry" boundary.
#
# Axial layout (z = 0 at the bottom of the graphite):
#   0.00 - 3.24 m   graphite core, control rods withdrawn (guide tube filled with He below the tip)
#   3.24 - 3.60 m   rodded zone: the Hf absorbers are inserted 36 cm from the top (SPEC)
# The fuel channels are open at z = 0 (inlet) and z = 3.6 m (outlet); the lower and upper
# plena are not part of the GFMR specification and are not meshed.
#
# Run:  open_pronghorn-opt -i gfmr_mesh_3d.i --mesh-only gfmr_sector_3d_cfd.e
################################################################################

height = 3.6
cr_insertion = 0.36

# axial cells: 40 below the rod tip (8.1 cm), 6 in the rodded zone (6 cm)
nz_lower = 40
nz_rodded = 6

!include gfmr_mesh_2d.i

[Mesh]
  [sector]
    type = HexagonMeshTrimmer
    input = reflector
    # Retain one 60-degree sector (two 30-degree units) between symmetry lines 0 and 2
    center_trim_starting_index = 0
    center_trim_ending_index = 2
    center_trimming_section_boundary = 'symmetry'
  []

  # The guide tube interior below the rod tip is helium; above it the absorber and gap are present.
  [extrude]
    type = AdvancedExtruderGenerator
    input = sector
    direction = '0 0 1'
    heights = '${fparse height - cr_insertion} ${cr_insertion}'
    num_layers = '${nz_lower} ${nz_rodded}'
    # lower interval: absorber positions become helium (rod withdrawn)
    # block ids: 30 absorber_tri, 31 absorber, 32 helium (see gfmr_mesh_2d.i); 34 = helium_tri (new)
    subdomain_swaps = '30 34 31 32;'
    bottom_boundary = 100
    top_boundary = 101
  []
  [name_helium_tri]
    type = RenameBlockGenerator
    input = extrude
    old_block = '34'
    new_block = 'helium_tri'
  []
  [rename_ends]
    type = RenameBoundaryGenerator
    input = name_helium_tri
    old_boundary = '100 101'
    new_boundary = 'bottom top'
  []

  # ---------------------------- named surfaces ----------------------------
  # one inlet per velocity zone (the zones are fed at different velocities, SPEC)
  [inlet_inner]
    type = ParsedGenerateSideset
    input = rename_ends
    combinatorial_geometry = 'z < 1e-6'
    included_subdomains = 'salt_inner_tri salt_inner'
    new_sideset_name = 'inlet_inner'
    replace = true
  []
  [inlet_outer]
    type = ParsedGenerateSideset
    input = inlet_inner
    combinatorial_geometry = 'z < 1e-6'
    included_subdomains = 'salt_outer_tri salt_outer'
    new_sideset_name = 'inlet_outer'
    replace = true
  []
  [inlet_cr]
    type = ParsedGenerateSideset
    input = inlet_outer
    combinatorial_geometry = 'z < 1e-6'
    included_subdomains = 'salt_cr'
    new_sideset_name = 'inlet_cr'
    replace = true
  []
  [outlet]
    type = ParsedGenerateSideset
    input = inlet_cr
    combinatorial_geometry = 'z > ${fparse height - 1e-6}'
    included_subdomains = 'salt_inner_tri salt_inner salt_outer_tri salt_outer salt_cr'
    new_sideset_name = 'outlet'
    replace = true
  []
  [wall_inner]
    type = SideSetsBetweenSubdomainsGenerator
    input = outlet
    primary_block = 'salt_inner'
    paired_block = 'graphite'
    new_boundary = 'wall_inner'
  []
  [wall_outer]
    type = SideSetsBetweenSubdomainsGenerator
    input = wall_inner
    primary_block = 'salt_outer'
    paired_block = 'graphite'
    new_boundary = 'wall_outer'
  []
  [wall_cr]
    type = SideSetsBetweenSubdomainsGenerator
    input = wall_outer
    primary_block = 'salt_cr'
    paired_block = 'graphite'
    new_boundary = 'wall_cr'
  []
  [wall_guide]
    type = SideSetsBetweenSubdomainsGenerator
    input = wall_cr
    primary_block = 'salt_cr'
    paired_block = 'guide_tube'
    new_boundary = 'wall_guide'
  []
  # the 2-D cell generators leave unnamed interface sidesets (ids 1-11) on the channel walls;
  # linear FV allows one boundary per face, so they are removed
  [clean]
    type = BoundaryDeletionGenerator
    input = wall_guide
    boundary_names = '1 3 5 7 9 11'
  []
[]
