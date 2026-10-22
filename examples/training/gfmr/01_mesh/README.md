# Stage 1: GFMR meshes with the Reactor-module mesh generators

All meshes are built from the same channel map ([../common/gfmr.py](../common/gfmr.py),
reproduced by `make_pattern.py`) with MOOSE Reactor-module generators, in the style of the ABTR
training.

| Input | Product | Cells | Used by |
|---|---|---|---|
| `gfmr_mesh_2d.i` | 2-D full cross section | 70 k | inspection, the base of the 3-D sector |
| `gfmr_mesh_3d.i` | `gfmr_sector_3d_cfd.e`: 3-D 60° sector, explicit channels, zone inlets | 540 k | diffusion (stage 2), CFD (stage 3) |
| `gfmr_mesh_3d.i` (overrides, see below) | `gfmr_sector_3d_rods_in.e`, `_rods_out.e` | 540 k | rod worth (stage 2) |
| `gfmr_mesh_homog.i` | `gfmr_homog_3d.e`: homogenized hex columns, full core | 117 k | not used yet (for a future space-time transient) |

Run any of them with `open_pronghorn-opt -i <input> --mesh-only <output.e>`. The rod-position
meshes use command-line overrides of the extruder swaps:

```
open_pronghorn-opt -i gfmr_mesh_3d.i --mesh-only gfmr_sector_3d_rods_out.e \
  "Mesh/extrude/subdomain_swaps=30 34 31 32; 30 34 31 32"
open_pronghorn-opt -i gfmr_mesh_3d.i --mesh-only gfmr_sector_3d_rods_in.e \
  "Mesh/extrude/subdomain_swaps=;" Mesh/name_helium_tri/old_block=30 \
  Mesh/name_helium_tri/new_block=absorber_tri
```

## Construction

1. **Unit cells** (`PolygonConcentricCircleMeshGenerator`): 14 cm hexagons, 4 cell sectors per
   side, flats facing the six neighbours. The four cell types are:
   - inner and outer fuel channel: 4 cm salt circle with a 3 mm wall ring, in graphite;
   - control-rod channel: Hf absorber (r 6.5 mm), He gap (1 mm), Hastelloy N guide tube (2 mm),
     and a salt annulus with wall rings at both walls;
   - solid graphite plug.
2. **Lattice** (`PatternedHexMeshGenerator`): the 23-row pattern of 397 positions (82 inner, 288
   outer, 9 rod, 18 plugs).
3. **Reflector** (`PeripheralRingMeshGenerator`): graphite out to R = 1.75 m.
4. **Sector** (`HexagonMeshTrimmer`): the 60° wedge between two mirror lines (next section).
5. **Extrusion** (`AdvancedExtruderGenerator`): 3.6 m in 46 layers. The top 36 cm holds the
   inserted absorber, swapped to helium below the rod tip (`subdomain_swaps`).
6. **Surfaces**: zone inlets, outlet, the four salt walls (inner / outer / rod channel / guide
   tube), the symmetry planes and the reflector side. The unnamed interface sidesets left by the
   cell generators are removed (linear FV allows one boundary per face).

## Symmetry: why 60° and not 30°

The three central control rods sit at 60°, 180° and 300°. The core is therefore **not**
six-fold rotationally symmetric, but it has mirror lines every 60°. The 60° wedge between two
mirror lines is an exact representation of the full core when all rods are at the same height.
A 30° wedge is not: it would place a mirror plane through the asymmetric central rod triangle.

## Near-wall resolution

The salt walls are resolved with wall functions, which need the first cell centre in the log
layer. The 3 mm wall ring puts it at y = 1.5 mm, which gives the following first-cell y⁺ at
nominal flow:

| Zone | Velocity | y⁺ of first cell |
|---|---|---|
| inner | 1.70 m/s | 40 |
| outer | 1.03 m/s | 26 |
| rod annulus | 1.74 m/s | ≈ 40 (two walls) |

Without the ring (three uniform rings over the 2 cm radius) y⁺ is 60–90.

## Exercises

1. Plot `gfmr_sector_3d_cfd.e` and identify the nine-rod pattern in the 60° wedge (rods at the
   wedge apex and on the ring-4 corner).
2. Change `nsec` in `gfmr_mesh_2d.i` from 4 to 6 and compare the cell count.
3. Explain why `BoundaryDeletionGenerator` is needed before a linear-FV solve.
