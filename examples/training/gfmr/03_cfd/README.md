# Stage 3: 3-D conjugate CFD of the GFMR core

`gfmr_core_cfd.i` solves the steady turbulent flow and heat transfer of the 60° core sector. It
runs on 540 k cells with linear finite volumes and SIMPLE, and is driven by the diffusion power
map of stage 2.

## Model

| Part | Model |
|---|---|
| Salt flow | incompressible, standard k-ε, equilibrium (Newton) wall functions |
| Salt energy | advection + conduction k + μ_t c_p / Pr_t (Pr_t = 0.85) |
| Graphite | conduction, k = 25 W/m/K (SPEC) |
| Hastelloy N guide tube | conduction with the graphite k = 25 W/m/K (simplification: Hastelloy N is ~21 W/m/K near 700 °C; the 2 mm tube carries no net heat) |
| Salt-solid interfaces | `LinearFVTurbulentWallCHTBC` on both sides (Jayatilleke wall function) |
| Power | diffusion power map: 92.1 % in the salt as mapped, 7.9 % in the graphite shaped like the E > 5.53 keV flux (groups 1-2) |
| Inlets | SPEC zone velocities 1.70 / 1.03 / 1.74 m/s at 613.8 °C |
| Outlet | p = 0, outflow |
| Symmetry planes | slip / adiabatic |
| Outer graphite surfaces | adiabatic |
| Absorber and helium | removed (stagnant inside the guide tube) |

FUNaK properties are taken at 650 °C from the SPEC correlations: ρ = 4073.5 kg/m³,
μ = 14.9 mPa s, c_p = 900 J/kg/K. The conductivity k = 0.8 W/m/K is assumed, giving Pr = 16.8.

### Why a new wall condition

At Pr ≈ 17 most of the thermal resistance between the salt and the graphite sits in the
conductive sublayer, a fraction of a millimetre thick. A wall-function mesh does not resolve
that layer. Conducting through the near-wall cell with k + μ_t c_p / Pr_t overestimates the heat
transfer about 46-fold for the inner channel (1.8×10⁵ against 3.8×10³ W/m²/K at y⁺ = 40).

`LinearFVTurbulentWallCHTBC` uses the Jayatilleke P-function, which carries that resistance. It
is verified in `test/tests/kEpsilon/salt_pipe_heat_transfer`: a GFMR inner channel gives
Nu = 180 against 194 from Gnielinski, a ratio of 0.93, inside the ±10–20 % uncertainty of the
correlation. See `doc/content/source/turbulence/linearfvbcs/LinearFVTurbulentWallCHTBC.md`.

## Results

600 SIMPLE iterations, 24 ranks, 19 min. Momentum residuals reach 1e-6. The energy residual
levels off near 3e-5 while the interface heat flows are steady to 1e-5.

### Balances (verification)

| Check | Result |
|---|---|
| Energy balance (enthalpy rise − heat input) / heat input | 4.5×10⁻⁵ |
| Salt / graphite heat input of the sector | 23.03 / 1.97 MW (25.0 MW = 150 / 6) |
| Graphite heat leaving through the four walls | 1.97 MW (fluxes balanced to 1e-5) |
| Core mass flow (6 × sector) | 2294 kg/s (SPEC 2300) |
| Mixed outlet temperature | 686.4 °C (heat balance 686.4 °C) |

### Comparison with the specification

The channel pressure drops come from the CFD, with the specified velocities imposed at the inlets:

| Zone | CFD Δp | SPEC Δp | Ratio |
|---|---|---|---|
| inner (1.70 m/s) | 14.4 kPa | 14.3 kPa | 1.01 |
| outer (1.03 m/s) | 6.2 kPa | 6.0 kPa | 1.03 |
| rod annulus (1.74 m/s) | 35.4 kPa | 34 kPa | 1.04 |

The SPEC lists 800 °C as the graphite temperature (the temperature of the neutronics cross
sections). The CFD gives a moderator volume average of 721 °C, so the neutronics reference is
about 80 K hot for the graphite. With the −4.1 pcm/K graphite coefficient that is worth about
+330 pcm, which would bring the core close to critical (exercise 4).

### Temperatures

| Quantity | Value |
|---|---|
| Inlet | 613.8 °C |
| Mixed outlet: inner / outer / rod channels | 685.4 / 686.8 / 685.3 °C |
| Hottest salt cell (outer zone) | 719.5 °C |
| Graphite moderator: average / maximum | 721 / 830 °C |
| Radial reflector maximum | 859 °C (outer edge, r = 1.72 m) |
| Guide tube maximum | 690 °C |

**The zoned velocities do what the design intends.** The inner channels make 31 % of the power
with 22 % of the channels. At 1.70 m/s they leave at the same mixed temperature as the outer
channels at 1.03 m/s, within 1.5 K.

**The hottest graphite is in the radial reflector**, at 859 °C near its outer surface
(r = 1.72 m, mid-height). The moderator peaks at 830 °C at the lattice edge (r = 1.60 m), and the
central moderator at 818 °C. The reflector is heated but has no coolant, and its outer surface is
adiabatic here. Its heat must conduct inward to the outermost channels across tens of centimetres
of graphite.

**This edge result rests on two modelling choices and should be read as a sensitivity, not a
prediction.** First, the graphite heat (7.9 % of the power) is spread in proportion to the
fast flux, whereas gamma heating has a different and longer spatial distribution. Second, the
reflector's outer surface has no heat sink; a real vessel or downcomer would cool it.

### Limits of this calculation

- One steady state at nominal flow. There is no loop, so the inlet temperature is imposed.
- Standard k-ε with wall functions. The salt-pipe test puts the wall heat transfer at 0.93 ×
  Gnielinski; the graphite temperatures carry that uncertainty (about 3 K on a 30–40 K
  film drop).
- The power map is frozen: no temperature feedback on the neutronics (one-way coupling).
- The salt conductivity (0.8 W/m/K) is assumed.
- The graphite heat source follows the E > 5.53 keV flux (groups 1-2), not a gamma transport solution; the
  OpenMC heating tally (`heating_mesh`, 10 cm mesh) would give the real distribution.
- The outer reflector surface is adiabatic.

## How to run

```bash
# needs ../02_neutronics/sector_power.e (the diffusion power map on the CFD mesh):
#   cd ../02_neutronics && mpirun -n 24 open_pronghorn-opt -i gfmr_diffusion_sector.i -w \
#     Outputs/file_base=sector_power mesh_file=../01_mesh/gfmr_sector_3d_cfd.e
mpirun -n 24 open_pronghorn-opt -i gfmr_core_cfd.i     # ~19 min, 600 SIMPLE iterations
```

## Exercises

1. Why do the inner and outer channels need different velocities to reach similar outlet
   temperatures? Use the zone powers of stage 2.
2. Compare the computed channel pressure drops with the SPEC values.
3. Replace `LinearFVTurbulentWallCHTBC` with a plain Dirichlet/Robin CHT pair and explain the
   change in the graphite temperature.
4. The neutronics assumed graphite at 800 °C; the CFD gives 721 °C. Estimate the reactivity
   effect with the graphite coefficient of stage 2, and the change in k-eff.
5. Put a convective condition on `reflector_outer` and see how the reflector maximum moves.
