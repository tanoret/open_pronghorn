# GFMR training: theory manual

This is the reference for the GFMR training. It explains the physics and the numerical models
behind the five stages: the mesh, the neutronics, the conjugate CFD, the chemistry and the
transients. Each equation was checked against the inputs, the Python scripts and the C++ source
the stages use. Where a comment and the code disagree, **the code wins**. All numbers are quoted
from the stage READMEs and from `02_neutronics/xs/openmc_summary.json`. They are **code results**:
the GFMR specification is a reference problem without published results, so nothing here is
validated against measurements.

## 0. Purpose, notation and model map

### 0.1 The reactor in one paragraph

The GFMR (Saltfoss Energy ApS, DOI 10.5281/zenodo.23348200) is a 150 MWth graphite-moderated
reactor. Its fuel salt, 26.4 UF4 - 24.7 KF - 48.9 NaF mol% ("FUNaK", 2.12 % enriched), flows
upward through 379 channels of 4 cm diameter. The channels sit in a hexagonal lattice of 14 cm
graphite cells, surrounded by a graphite reflector out to R = 1.75 m. The core is 3.6 m high.
Nine hafnium control rods move inside Hastelloy N guide tubes placed in fuel channels. At nominal
conditions the rods are inserted 36 cm from the top. The salt enters at 613.8 °C and leaves at
686.2 °C (72.5 K rise at 2300 kg/s). This assumes the specified 650 °C fuel temperature is the
core-average salt temperature (an assumption, see the main README).

### 0.2 Notation

| symbol | meaning | unit |
|---|---|---|
| $p$ | lattice pitch (flat-to-flat of the unit cell), 0.14 | m |
| $(i, j)$ | axial hexagonal lattice coordinates | - |
| $g, g'$ | energy group index, $1 \dots G$ with $G = 8$, fast first | - |
| $\phi_g$ | scalar flux of group $g$ | 1/(m$^2$ s) |
| $D_g$ | diffusion coefficient | m |
| $\Sigma_{a,g}$, $\Sigma_{tr,g}$ | macroscopic absorption, transport cross section | 1/m |
| $\Sigma_{s,g'\to g}$ | scattering from $g'$ to $g$ (here the nu-scatter matrix) | 1/m |
| $\nu\Sigma_{f,g}$, $\kappa\Sigma_{f,g}$ | neutron production, recoverable energy release per unit flux | 1/m, J/m |
| $\chi_g$ | fission spectrum (fraction of fission neutrons born in $g$) | - |
| $k$, $\lambda$ | multiplication factor; MOOSE eigenvalue ($\lambda = 1/k$) | - |
| $\rho$ (neutronics) | reactivity, $(k-1)/k$ | pcm = 1e-5 |
| $\beta_i$, $\beta$ | delayed-neutron fraction of precursor group $i$; total | - |
| $\lambda_i$ | precursor decay constant (context distinguishes it from the eigenvalue) | 1/s |
| $C_i$ | precursor concentration | 1/m$^3$ |
| $\Lambda$ | prompt generation time | s |
| $v_g$ | group neutron speed | m/s |
| $\rho$ (flow), $\mu$, $c_p$, $k_f$ | salt density, viscosity, heat capacity, conductivity | kg/m$^3$, Pa s, J/(kg K), W/(m K) |
| $\mathbf u$, $p$ | velocity, pressure (context distinguishes $p$ from the pitch) | m/s, Pa |
| $k_t$, $\varepsilon$ | turbulent kinetic energy and its dissipation rate | m$^2$/s$^2$, m$^2$/s$^3$ |
| $\mu_t$ | turbulent (eddy) viscosity | Pa s |
| $y$, $y^+$ | distance of the near-wall cell centre from the wall, its wall unit | m, - |
| $u_\tau$ | friction velocity | m/s |
| $Pr$, $Pr_t$ | molecular and turbulent Prandtl numbers | - |
| $h_{wf}$ | wall-function heat-transfer coefficient | W/(m$^2$ K) |
| $k_s$, $d_S$ | solid conductivity, solid cell-centre-to-face distance | W/(m K), m |
| $q'''$ | volumetric heat source | W/m$^3$ |

The flow symbol $\rho$ (density) and the neutronics symbol $\rho$ (reactivity) never appear in
the same equation. The turbulent kinetic energy is written $k_t$ so that it is not confused with
$k$ (multiplication factor) or $k_f$ (conductivity); the inputs call it `TKE` and $\varepsilon$
`TKED`.

### 0.3 Model map

| stage | physics | MOOSE objects / tools | files |
|---|---|---|---|
| 1 | Reactor-module meshing | `PolygonConcentricCircleMeshGenerator`, `PatternedHexMeshGenerator`, `PeripheralRingMeshGenerator`, `HexagonMeshTrimmer`, `AdvancedExtruderGenerator`, `BoundaryDeletionGenerator` | `01_mesh/*.i`, `common/gfmr.py` |
| 2 | continuous-energy Monte Carlo; 8-group diffusion eigenvalue | OpenMC 0.16 (ENDF/B-VIII.0); `ArrayDiffusion`, `ArrayReaction`, `ArrayVacuumBC`, `Eigenvalue` | `02_neutronics/openmc/*.py`, `make_diffusion_input.py`, `gfmr_diffusion_*.i` |
| 3 | steady RANS k-ε, conjugate heat transfer | linear FV, `SIMPLE`, `RhieChowMassFlux`, `LinearFVTKESourceSink`, `LinearFVTurbulentWallCHTBC` | `03_cfd/gfmr_core_cfd.i`, `src/turbulence/linearfvbcs/LinearFVTurbulentWallCHTBC.C` |
| 4 | radiolysis, redox, corrosion, plating | `[MoltenSaltRadiolysis]`, `[CorrosionPlatingFlow]`, `[FissionProductInventory]` on the frozen CFD; 0-D loop balance | `04_chemistry/gfmr_chemistry.i`, `loop_balance.py` |
| 5 | point kinetics with circulating precursors, two-temperature feedback | scalar ODEs: `ODETimeDerivative`, `ParsedODEKernel`, `Transient` (BDF2) | `05_transients/gfmr_rod_step.i` |

### 0.4 Verification, calibration, code-to-code

The labels of `examples/training/chemistry/THEORY.md` section 6 are used here:

| label | meaning in this training |
|---|---|
| **VERIFICATION** | the code solves its own equations correctly (balances, invariants, a known solution) |
| **CALIBRATION** | a parameter was chosen to reproduce the quantity it is then compared with |
| **CODE-TO-CODE** | two independent methods for the same model (OpenMC against diffusion; CFD against a correlation) |
| **VALIDATION** | comparison with independent measurements; **none exists for the GFMR** |

## 1. Stage 1: mesh

### 1.1 The hexagonal lattice and the axial coordinates of `common/gfmr.py`

Every lattice position is named by two integers $(i, j)$, the **axial coordinates** of a
hexagonal grid. ("Axial" is the standard name of this hexagonal coordinate system. It has nothing
to do with the vertical $z$ axis.) The centre of the cell $(i, j)$ is

$$x = p\left(i + \frac{j}{2}\right), \qquad y = p\,\frac{\sqrt 3}{2}\,j ,$$

so the $i$ axis runs along $+x$ and the $j$ axis runs at 60° to it. The six neighbours of a cell
are at $(\pm1, 0)$, $(0, \pm1)$, $(1, -1)$ and $(-1, 1)$, at distance $p$ in the directions 0°,
60°, ..., 300°. The cells are hexagons of apothem $p/2$, and each flat faces one neighbour. The
mesh inputs build them with `flat_side_up = false`, so a vertex points along $+y$. The **ring** of a
position is its hexagonal distance from the centre,

$$r(i,j) = \max\bigl(|i|,\ |j|,\ |i + j|\bigr),$$

and ring $r \ge 1$ holds $6r$ positions. Rings 0 to 11 hold $1 + 3\cdot 11\cdot 12 = 397$
positions. `channel_type(i, j)` assigns them:

| rings | positions | content |
|---|---|---|
| 0-5 | 91 | inner zone: 82 inner channels + 9 control-rod channels |
| 6-10 | 240 | outer zone channels |
| 11 | 66 | 48 outer channels; the corner and its two ring neighbours at each of the 6 corners (18) are solid graphite plugs |

The rods are at the three ring-1 positions $(1,-1)$, $(-1,0)$, $(0,1)$ and at the six corners of
ring 4. This gives 82 + 288 + 9 = 379 channels, as specified. The map is not tabulated in the
specification. It was decoded from the published top view, and the decoded map reproduces the
three stated counts exactly. `make_pattern.py` writes it as the 23-row pattern of
`PatternedHexMeshGenerator` (rows $j = +11 \dots -11$, entries in $+x$).

The vertical coordinate has $z = 0$ at the bottom of the graphite (the salt inlet) and
$z = 3.6$ m at the top (the outlet). The rods enter from the top: the rodded zone is
$3.24 \le z \le 3.60$ m. The specification has no plena or axial reflectors, so none are meshed;
the top and bottom of the graphite are vacuum boundaries in the neutronics (section 2.4).

### 1.2 Symmetry: why a 60° sector is exact and a 30° sector is not

The lattice itself, the zone map and the ring-11 corner cuts have the full symmetry of the
hexagon (the dihedral group $D_6$: six rotations by multiples of 60° and six mirror lines, every
30°). The six ring-4 rods sit at the corners of a hexagon and are also $D_6$-symmetric. **The
three central rods break this.** With the centre formula above:

| rod | $(x, y)/p$ | angle |
|---|---|---|
| $(0, 1)$ | $(0.5,\ 0.866)$ | 60° |
| $(-1, 0)$ | $(-1,\ 0)$ | 180° |
| $(1, -1)$ | $(0.5,\ -0.866)$ | 300° |

An equilateral triangle has the symmetry group $D_3$: rotations by 120° and three mirror lines,
each through one vertex and the midpoint of the opposite side. Here the lines run through 60°/240°,
180°/0° and 300°/120°. Seen as half-lines from the centre, there is a mirror every 60°: at 0°,
60°, 120° and so on. A wedge between two adjacent mirror half-lines, e.g. $0° \le \theta \le 60°$,
is a fundamental domain of $D_3$: the reflections across its two edges generate the whole core.
Zero normal gradient (zero net current, zero heat flux, zero normal velocity) on the two cut planes
is therefore **exact**, as long as the solution keeps the symmetry. The rods must all be at the
same height, and no instability may break the symmetry.

A 30° wedge would need a mirror line at 30°. A reflection at 30° sends the rod at 60° to 0°,
where there is no rod. A 30° model would therefore represent a core with a rod at 0° as well:
six central rods instead of three. In `gfmr_mesh_3d.i`, `HexagonMeshTrimmer` keeps the sector
between trimming lines 0 and 2 (two 30° units) and names the cut faces `symmetry`. The model
represents 1/6 of the core: the neutronics normalizes the sector to 25 MW (150 MW / 6), and the
CFD core mass flow is $6\times$ the sector flow.

An asymmetric rod pattern (one rod stuck, one rod dropped) breaks the symmetry and needs the
full core. That is why the transients plan uses the full-core homogenized mesh
(`gfmr_homog_3d.e`).

### 1.3 Finite-volume mesh requirements

**One boundary per face.** In MOOSE's linear finite-volume (FV) system each boundary face
carries exactly one boundary condition per variable. The `LinearFV*BC` object is looked up from
the face's boundary id, so a face that belongs to two sidesets is ambiguous. The 2-D cell
generators leave unnamed interface sidesets (ids 1, 3, 5, 7, 9, 11) on the channel walls. These
coincide with the named walls `wall_inner`, `wall_outer`, `wall_cr`, `wall_guide`, so
`BoundaryDeletionGenerator` removes them before the CFD reads the mesh. The finite-element
diffusion solve does not care: its interior faces are not boundaries of any kernel, and its
reflective symmetry planes need no object at all (natural condition).

**Volume preservation.** `preserve_volumes = true` moves the polygon nodes so that the
discretized rings have the area of the true circles. Without it, the salt volume (and the mass
flux at a given inlet velocity) would be too small by the polygon factor.

**Near-wall resolution for wall functions.** The salt walls use wall functions (section 3.2).
These assume the first cell centre lies in the logarithmic layer, conventionally
$y^+ = \rho u_\tau y/\mu \gtrsim 30$. A low-Reynolds model would instead need $y^+ \lesssim 1$
and a much finer mesh. Each salt circle therefore ends with a 3 mm ring at the wall
(`wall_layer = 0.003`), so the first cell centre is at $y = 1.5$ mm. The first-cell $y^+$ at
nominal flow (stage-1 README) is:

| zone | velocity | $y^+$ |
|---|---|---|
| inner | 1.70 m/s | 45 |
| outer | 1.03 m/s | 29 |
| rod annulus | 1.74 m/s | ≈ 40 |

Without the ring (three uniform rings over the 2 cm radius), $y^+$ is 60-90. The outer zone, at 29,
sits just below 30: there the thermal wall condition blends toward the sublayer value
(section 3.3). The axial cells are much longer than the radial ones (46 layers over 3.6 m). This is
acceptable because the developed channel flow has small axial gradients; it is not acceptable for
the inlet entry region, which is not resolved.

**Mesh sizes.** The 2-D section has 70 k cells. The 3-D sector (46 layers: 40 below the rod tip,
6 in the rodded 36 cm) has 540 k cells. The homogenized full core (one hexagonal column per cell)
has 117 k cells. The diffusion solve uses the same 540 k-cell sector as the CFD
(3.6×10⁶ unknowns), so the power density lands on the CFD cells without interpolation between
meshes.

## 2. Stage 2: neutronics

### 2.1 The two-step method

Stage 2 is the classical two-step scheme. First, a continuous-energy Monte Carlo transport
calculation of the heterogeneous core (OpenMC 0.16, ENDF/B-VIII.0, $2\times10^5$ particles ×
200 active batches = $4\times10^7$ histories per case) gives the reference eigenvalues, rod
worths, temperature coefficients and kinetics parameters. It also tallies the few-group constants.
Second, a deterministic multigroup diffusion solve in MOOSE uses those constants on the CFD mesh.
The Monte Carlo is the reference; the diffusion is the fast, mesh-conforming model that feeds the
CFD and, later, the transients.

The 8 groups are the CASMO-8 structure (`openmc.mgxs.EnergyGroups('CASMO-8')`). The edges, fast
first, are 20 MeV, 0.821 MeV, 5.53 keV, 4 eV, 0.625 eV, 0.28 eV, 0.14 eV, 0.058 eV and 0. Five of
the eight groups are below 4 eV, which resolves the thermalization in the graphite.

### 2.2 The multigroup diffusion equation in the form the inputs solve

The steady multigroup diffusion equation (Duderstadt & Hamilton 1976, ch. 7) for group $g$ is

$$-\nabla\cdot\bigl(D_g\nabla\phi_g\bigr) + \Sigma_{t,g}\phi_g
 - \sum_{g'=1}^{G}\Sigma_{s,g'\to g}\phi_{g'}
 = \frac{1}{k}\,\chi_g\sum_{g'=1}^{G}\nu\Sigma_{f,g'}\phi_{g'} ,$$

where $\Sigma_{t,g}$ is the total cross section. The within-group scattering $\Sigma_{s,g\to g}$
appears on both sides of the balance and cancels:
$\Sigma_{t,g} - \Sigma_{s,g\to g} = \Sigma_{a,g} + \sum_{g'\ne g}\Sigma_{s,g\to g'}$ (absorption plus
out-scatter). `make_diffusion_input.py` writes the equation in matrix form for the array
variable $\boldsymbol\phi = (\phi_1,\dots,\phi_8)^T$:

$$-\nabla\cdot\bigl(\mathbf D\nabla\boldsymbol\phi\bigr) + \mathbf R\,\boldsymbol\phi
 = \frac{1}{k}\,\mathbf F\,\boldsymbol\phi ,$$

$$\mathbf R = \operatorname{diag}\Bigl(\Sigma_{a,g} + \sum_{g'\ne g}\Sigma_{s,g\to g'}\Bigr)
 - \mathbf S_\text{off}^T, \qquad
 (\mathbf S_\text{off})_{g g'} = \Sigma_{s,g\to g'}\ (g\ne g'),\ 0\ (g = g'),$$

$$\mathbf F_{g g'} = \chi_g\,\nu\Sigma_{f,g'} .$$

Here $\mathbf D = \operatorname{diag}(D_g)$. The **removal matrix** $\mathbf R$ has the removal
cross section on its diagonal and the negative in-scatter off the diagonal. Row $g$ of
$\mathbf R\boldsymbol\phi$ is the net loss from group $g$. In the committed `REM_core`, the first row
has only a diagonal entry (nothing scatters up into the fastest group). The thermal groups 4-8
have entries above the diagonal: this is upscatter in the graphite. The **fission matrix**
$\mathbf F$ is an outer product: rank one, nonzero only in the rows $g = 1$-3 where $\chi_g > 0$.
The eigenvalue spectrum $\chi_g$ is the total spectrum, which mixes the prompt and delayed spectra
by the IFP fractions,

$$\chi_g = (1-\beta)\,\chi_{p,g} + \beta\,\chi_{d,g}, \qquad
\chi_{d,g} = \frac{\sum_i \beta_i\,\chi_{d,i,g}}{\beta}.$$

The OpenMC constants are in cm units; the script converts them to SI ($D$ in m, $\Sigma$ in 1/m,
$\kappa\Sigma_f$ in J/m). The diffusion coefficient is $D_g = 1/(3\Sigma_{tr,g})$, from the
OpenMC `transport` cross section with the P0 correction. It is capped at 0.5 m in the helium gaps,
where $\Sigma_{tr}\to 0$: a 1 mm void region with infinite $D$ would make the operator singular.

**The MOOSE eigenvalue convention.** The `Eigenvalue` executioner solves the generalized problem

$$\mathbf A\,\boldsymbol\phi = \lambda\,\mathbf B\,\boldsymbol\phi ,$$

where $\mathbf A$ is assembled from the kernels without a tag (here `ArrayDiffusion` with $\mathbf D$
and `ArrayReaction` with $\mathbf R$), and $\mathbf B$ from the kernels tagged
`extra_vector_tags = 'eigen'`. An `ArrayReaction` adds $+\mathbf M\boldsymbol\phi$ to the residual.
A kernel with $\mathbf M = \mathbf F$ would therefore give $\mathbf B = -\mathbf F$ and
$\lambda = -1/k$. The inputs pass $\mathbf M = -\mathbf F$ instead (`F_total_neg` in the R-Z model,
`FN_<region>` in the sector model), so that

$$\lambda = \frac{1}{k}, \qquad k_\text{eff} = 1/\lambda$$

(the `k_eff` postprocessor). The fundamental mode is the eigenpair with the smallest $\lambda$, i.e.
the largest $k$. The executioner runs `free_power_iterations` (6 in the sector model) of the power
method $\mathbf A\boldsymbol\phi^{(n+1)} = \lambda^{(n)}\mathbf B\boldsymbol\phi^{(n)}$ to approach the
fundamental mode. It then solves the nonlinear eigenproblem by Newton (PJFNK, with BoomerAMG as the
preconditioner). An eigenvector is defined only up to a constant. `normalization = total_power`
scales it so that the kappa-fission power of the sector,

$$P = \int_{V_\text{salt}}\sum_g \kappa\Sigma_{f,g}\,\phi_g\,dV ,$$

equals 25 MW (150 MW / 6).

The spatial discretization is continuous Galerkin finite elements (CFEM, first-order Lagrange)
on the hexahedral and prismatic cells. The 8 groups are coupled in one array variable. On the
540 k-cell sector this is 3.6×10⁶ unknowns, solved in 6 min on 24 ranks.

### 2.3 Homogenization and cell-domain cross sections

A multigroup constant is a reaction-rate-preserving average of the continuous-energy cross section
$\sigma(\mathbf r, E)$ over a region $V$ and a group $g$:

$$\Sigma_{x,g} = \frac{\int_V\int_g \Sigma_x(\mathbf r,E)\,\phi(\mathbf r,E)\,dE\,dV}
                      {\int_V\int_g \phi(\mathbf r,E)\,dE\,dV} .$$

It is exact for the reaction rate of the reference flux, by construction. It is approximate as
soon as the flux shape changes (another rod position, another temperature, a perturbation): the
weighting spectrum is that of the nominal case. `run_reference.py` tallies this with
`openmc.mgxs.Library` in **cell domains** (`domain_type = 'cell'`). Each OpenMC cell (inner salt,
inner graphite, outer salt, outer graphite, rod-channel salt, guide tube, Hf rod, helium gap,
plug, reflector, ...) is one domain. A cell domain integrates over **all lattice instances** of
that cell, so a constant is a core average for its region type, not a local value. The 12 regions
of `gfmr_8g.json` are `reflector_lattice`, `plug_graphite`, `outer_salt`, `outer_graphite`,
`inner_salt`, `inner_graphite`, `cr_salt`, `cr_graphite`, `cr_guide`, `cr_rod`, `cr_he_low` and
`cr_gap`.

Two uses follow:

- **Heterogeneous sector model** (`gfmr_diffusion_sector.i`): each mesh block takes the constants
  of its own region. The graphite blocks of the mesh are not split by zone, so they all use
  `outer_graphite` (88 % of the moderator). The input comment states that the inner, outer and
  rod-channel graphite constants differ by < 0.5 %. In the committed `xs/gfmr_xs.i`, the
  diagonal removal entries of `inner_graphite` and `cr_graphite` differ from `outer_graphite` by up
  to about 0.9 % and 2.4 % respectively (fast groups); the $D_g$ values differ by < 0.2 %.
- **Homogenized R-Z model** (`gfmr_diffusion_rz.i`): one `core` region is the flux-volume
  weighting of all lattice regions,
  $\Sigma_{x,g}^\text{core} = \sum_m \Sigma_{x,m,g}\,\Phi_{m,g}\big/\sum_m\Phi_{m,g}$, with
  $\Phi_{m,g}$ the flux-volume integral of region $m$. The core $D$ is
  $1/(3\bar\Sigma_{tr})$ with the flux-weighted transport cross section (weighting $D$ directly
  would over-weight the low-density regions).

Flux-volume weighting preserves reaction rates but not the **leakage** between regions. The
diffusion coefficient of a region that is thin compared with a mean free path (a 3 mm wall ring,
a 1 mm gap, the 2 mm guide tube) is not physical. The comparison of section 2.6 shows that the
residual bias does not depend on the rods.

### 2.4 Boundary conditions: Marshak vacuum and reflection

At a free surface no neutrons come back in. In the $P_1$ (diffusion) approximation, the incoming
partial current is

$$J^-(\mathbf r) = \frac{\phi}{4} + \frac{D}{2}\,\nabla\phi\cdot\mathbf n ,$$

with $\mathbf n$ the outward normal. The **Marshak condition** $J^- = 0$ gives

$$-D\,\nabla\phi\cdot\mathbf n = \frac{1}{2}\,\phi .$$

`ArrayVacuumBC` adds $\alpha\,\phi_g$ to the weak form on the boundary, so `VAC_ALPHA` = 1/2 in
every group. It is applied on `top`, `bottom` and `reflector_outer`. (The equivalent extrapolated-
boundary view: the flux extrapolates to zero at a distance $2D$ beyond the surface, close to the
transport value $0.7104\,\lambda_{tr} = 2.13\,D$; Duderstadt & Hamilton 1976, ch. 4.) The two
symmetry planes need no object: the natural condition of the weak form is
$-D\nabla\phi\cdot\mathbf n = 0$, a reflective (zero net current) plane, which is the exact
condition for a mirror symmetry (section 1.2).

The bare top and bottom are a property of the specification as modelled, not of a real plant. With
no axial reflector or plena, the flux at the top of the core is low, and the leakage is large
(OpenMC leakage fraction 0.11721 ± 0.00005 in the graphite +100 K run log, `xs/openmc_run.log`).

### 2.5 Why diffusion fails near a black absorber, and what the rod factor means

Diffusion theory rests on Fick's law, $\mathbf J = -D\nabla\phi$, which follows from the $P_1$
expansion. It assumes an angular flux that is nearly isotropic, i.e. a region that is many mean
free paths thick, with weak absorption ($\Sigma_a \ll \Sigma_{tr}$) and far from sources, interfaces
and boundaries. A hafnium rod violates all three for thermal neutrons. The diagonal entries of the committed
`REM_cr_rod` (already scaled by 0.602, see below) are 121-244 1/m in groups 4-8, about 200-405 1/m
before scaling: a mean free path of 2.5-5 mm against a 6.5 mm rod radius. Neutrons entering the rod are absorbed within a mean free path of
the surface. The angular flux at the rod surface is strongly peaked inward (no neutrons come back
out), which is exactly where the $P_1$ closure is worst. Diffusion lets the flux inside the rod
couple to the outside through $D\nabla\phi$ as if the medium were scattering. With flux-volume
weighted constants this **overestimates** the absorption rate, and hence the rod worth
(Duderstadt & Hamilton 1976 discuss the classical remedies: blackness theory, an internal
boundary condition, or equivalence cross sections).

The stage-2 run with the raw constants gives a total rod worth of 5868 pcm against 4844 pcm from
OpenMC (21 % high). `make_diffusion_input.py --rod-factor 0.602` multiplies the **whole** Hf
removal matrix $\mathbf R_\text{cr\_rod}$ by 0.602. This reproduces the OpenMC total worth: 4845 pcm
against 4844 ± 21 pcm.

What the factor **is**: a one-parameter **calibration** of the absorber constants on one
quantity, the worth of the full insertion of all nine rods in this geometry, at nominal
temperatures, with this mesh and these 8 groups. It is an equivalence-type correction of the kind
used in lattice codes, determined here by matching an integral.

What it **is not**:

- not a prediction: the 4845 pcm agreement is in-sample by construction;
- not a physical property of hafnium: it compensates a method error, so it changes with the mesh,
  the group structure, the rod geometry and the spectrum;
- not a guarantee for partial insertions or the differential worth: the 36 cm insertion is worth
  46 pcm in diffusion against 16 ± 21 pcm in OpenMC. The OpenMC value is not significant, so this
  is neither a confirmation nor a refutation, but it is not covered by the calibration;
- not a correction of the scattering inside the rod: the factor scales absorption and in-scatter
  together, because it multiplies the whole matrix.

### 2.6 Comparison with the Monte Carlo reference

All values from the stage-2 README, on the final mesh with the 3 mm wall ring:

| quantity | diffusion | OpenMC | class |
|---|---|---|---|
| $k_\text{eff}$ nominal | 0.99624 (+115 pcm) | 0.99510 ± 0.00014 | CODE-TO-CODE |
| $k_\text{eff}$ rods out / rods in | 0.99670 / 0.95078 (+146 / +144 pcm) | 0.99525 / 0.94948 | CODE-TO-CODE |
| total rod worth | 4845 pcm | 4844 ± 21 pcm | **CALIBRATION** |
| power split inner / outer / rod | 30.8 / 66.6 / 2.6 % | 31.2 / 66.1 / 2.7 % | CODE-TO-CODE |
| axial peaking | 1.557 | 1.557 | CODE-TO-CODE |
| worth of the 36 cm insertion | 46 pcm | 16 ± 21 pcm | CODE-TO-CODE (not significant) |
| peak salt power density | 212 MW/m³ | - | model result |

The eigenvalue bias of about +130 pcm is the same with the rods in and out. It therefore comes from
the homogenization and the diffusion approximation in the lattice, not from the rods. The coarser
mesh without the wall ring gives the same $k$ within 1 pcm, so the bias is not a discretization
error either. Because the rod factor is fitted to a reactivity **difference**, it is insensitive to
a common bias.

### 2.7 Reactivity coefficients from perturbation runs

The reactivity of a state with multiplication factor $k$ is

$$\rho = \frac{k-1}{k}, \qquad \rho\,[\text{pcm}] = \frac{k-1}{k}\times10^{5}$$

(`rho_pcm` in `postprocess.py`). Reactivities of two states subtract; $k$ values do not. A worth or
a coefficient is the difference of two **independent** eigenvalue calculations:

$$w_\text{rods} = \rho(k_\text{out}) - \rho(k_\text{in}), \qquad
\alpha_X = \frac{\rho(k_{X+\Delta T}) - \rho(k_\text{nominal})}{\Delta T},\ \Delta T = 100\ \text{K}.$$

The perturbed runs change one thing: for the fuel, the salt temperature (Doppler broadening of
U-238, with the density from the SPEC law $\rho = 4.808 - 0.00113\,T$ g/cm³); for the graphite, the
moderator temperature (the thermal scattering law `c_Graphite`, interpolated in temperature). The
statistical uncertainty of a difference of independent runs is
$\sigma_{\Delta\rho} = \sqrt{\sigma_1^2 + \sigma_2^2}$, about 20 pcm here. A coefficient from
$\Delta T = 100$ K is therefore known to about 0.2 pcm/K. A smaller $\Delta T$ would drown the
signal in noise, and a larger one would leave the linear range.

The OpenMC results (`xs/openmc_summary.json`, 4×10⁷ histories per case):

| case | $k_\text{eff}$ |
|---|---|
| nominal (rods at 36 cm, fuel 650 °C, graphite 800 °C) | 0.99510 ± 0.00014 |
| rods withdrawn | 0.99525 ± 0.00015 |
| rods fully inserted (3.6 m) | 0.94948 ± 0.00013 |
| fuel +100 K | 0.99241 ± 0.00014 |
| graphite +100 K | 0.99103 ± 0.00015 |

| derived quantity | value |
|---|---|
| total rod worth (9 rods) | 4844 ± 21 pcm (7.1 \$) |
| worth of the nominal 36 cm insertion | 16 ± 21 pcm (not significant) |
| fuel temperature coefficient $\alpha_f$ | −2.7 ± 0.2 pcm/K |
| graphite temperature coefficient $\alpha_g$ | −4.1 ± 0.2 pcm/K |

The nominal reactivity is $\rho = (0.99510 - 1)/0.99510 \approx -490$ pcm: the core as specified
is subcritical (section 3.7 discusses the graphite temperature). In a circulating-fuel reactor the
fuel coefficient acts within the salt transit time and the graphite coefficient over the graphite
heat-up time (minutes). Both are negative here.

### 2.8 Kinetics parameters by iterated fission probability

The physical delayed-neutron fraction $\beta_0$ is the fraction of fission neutrons that are
delayed. Delayed neutrons are born with a softer spectrum than prompt ones, so in a thermal core
they are more likely to cause fission. The kinetics equations need the **effective** fraction, the
fraction weighted by the adjoint (importance) flux $\phi^\dagger$:

$$\beta_\text{eff} = \frac{\langle \phi^\dagger, \mathbf F_d\,\phi\rangle}
                         {\langle \phi^\dagger, \mathbf F\,\phi\rangle}, \qquad
\Lambda = \frac{\langle \phi^\dagger, \tfrac{1}{v}\,\phi\rangle}
               {\langle \phi^\dagger, \mathbf F\,\phi\rangle} ,$$

with $\mathbf F_d$ the delayed part of the fission operator and $\langle\cdot,\cdot\rangle$ the
integral over space and energy. $\Lambda$ is the prompt generation time. A Monte Carlo code does not
compute $\phi^\dagger$ directly. The **iterated fission probability** (IFP) method uses the fact that
the importance of a neutron is proportional to the number of its descendants in a later generation
(the asymptotic progeny). Each neutron carries the properties of its ancestors ($1/v$ and delayed
or prompt status) over a fixed number of generations, and the score of an ancestor is weighted by
the progeny of its descendants (Hurwitz 1964 for the IFP interpretation of the adjoint; Nauchi &
Kameyama 2010 and Kiedrowski 2011 for Monte Carlo IFP estimators of $\beta_\text{eff}$ and
$\Lambda$). OpenMC implements this with the tally scores `ifp-beta-numerator`,
`ifp-time-numerator` and `ifp-denominator`; the GFMR uses `ifp_n_generation = 10`.

| quantity | value |
|---|---|
| $\beta_\text{eff}$ (IFP) | 681 ± 5 pcm |
| $\Lambda$ (IFP) | 0.605 ms |

The six group fractions $\beta_{\text{eff},i}$ and the decay constants $\lambda_i$
(`BETA`, `LAMBDA` in `xs/gfmr_xs.i`, the decay constants from the inner-salt `decay-rate` MGXS)
are:

| group $i$ | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| $\beta_{\text{eff},i}$ (pcm) | 24.5 | 123.8 | 114.7 | 264.7 | 108.1 | 45.2 |
| $\lambda_i$ (1/s) | 0.01334 | 0.03270 | 0.1209 | 0.3036 | 0.8531 | 2.865 |

The generation time is long (0.6 ms) because a thermal neutron in graphite lives long. These
values are for **stationary** fuel; the effect of the salt flow on the delayed neutrons is the
subject of section 5.2.

### 2.9 Power deposition

The `kappa-fission` score is the recoverable energy per fission: fragments, betas, gammas and
neutrons. In a fuel salt the fragments and betas stay in the salt, but the gammas and the neutron
slowing-down heat the graphite. A coupled neutron-photon run (`heating` case) splits the total:
92.1 % in the salt, 7.9 % in the graphite, and negligible amounts in the rods and guide tubes
(`heating_fraction_by_material`). The CFD uses this split (section 3.6).

## 3. Stage 3: conjugate CFD

### 3.1 RANS and the standard k-ε model

The salt is incompressible, with constant properties taken at 650 °C from the SPEC correlations:
$\rho$ = 4073.5 kg/m³, $\mu$ = 14.9 mPa s, $c_p$ = 900 J/(kg K), and the **assumed**
$k_f$ = 0.8 W/(m K). This gives $Pr = c_p\mu/k_f$ = 16.8. The channel Reynolds numbers
$Re = \rho u D_h/\mu$ are about 18 600 (inner), 11 300 (outer) and 10 000 (rod annulus, hydraulic
diameter $D_h = D - 2r_\text{guide}$ = 2.1 cm); `common/gfmr.py` prints them. The flow is fully
turbulent, but not by a large margin.

The Reynolds-averaged Navier-Stokes (RANS) equations govern the mean velocity $\mathbf u$ and
pressure $p$. Averaging the Navier-Stokes equations produces the Reynolds stress
$-\rho\,\overline{\mathbf u'\mathbf u'}$, which the Boussinesq hypothesis models with an eddy
viscosity $\mu_t$:

$$\nabla\cdot\mathbf u = 0,$$

$$\nabla\cdot(\rho\,\mathbf u\mathbf u) = -\nabla p
 + \nabla\cdot\Bigl[(\mu+\mu_t)\bigl(\nabla\mathbf u + \nabla\mathbf u^T\bigr)\Bigr].$$

(The steady model has no time derivative; gravity is absorbed in the pressure because the density
is constant.) The **standard k-ε model** (Launder & Spalding 1974) closes $\mu_t$ with two transport
equations, for the turbulent kinetic energy $k_t$ and its dissipation rate $\varepsilon$:

$$\nabla\cdot(\rho\,\mathbf u\,k_t) = \nabla\cdot\Bigl[\Bigl(\mu + \frac{\mu_t}{\sigma_k}\Bigr)\nabla k_t\Bigr]
 + G_k - \rho\varepsilon ,$$

$$\nabla\cdot(\rho\,\mathbf u\,\varepsilon) = \nabla\cdot\Bigl[\Bigl(\mu + \frac{\mu_t}{\sigma_\varepsilon}\Bigr)\nabla\varepsilon\Bigr]
 + C_{1\varepsilon}\frac{\varepsilon}{k_t}G_k - C_{2\varepsilon}\rho\frac{\varepsilon^2}{k_t} ,$$

$$\mu_t = \rho\,C_\mu\frac{k_t^2}{\varepsilon}, \qquad
G_k = \mu_t\,\bigl(\nabla\mathbf u + \nabla\mathbf u^T\bigr):\nabla\mathbf u ,$$

where $G_k$ is the production of turbulent energy by the mean shear. The constants are those of
Launder & Spalding, used unchanged in `gfmr_core_cfd.i`:

| $C_\mu$ | $\sigma_k$ | $\sigma_\varepsilon$ | $C_{1\varepsilon}$ | $C_{2\varepsilon}$ |
|---|---|---|---|---|
| 0.09 | 1.0 | 1.3 | 1.44 | 1.92 |

At the inlets $k_t = \tfrac32 (I u_\text{inner})^2$ with intensity $I$ = 5 %, and
$\varepsilon = C_\mu^{3/4}k_t^{3/2}/\ell$ with the mixing length $\ell = 0.07\,D$. The same values
are used in all three zones. The inlet profile is flat; the channels are 90 diameters long, so the
inlet state is forgotten well before the outlet. `C_pl = 1e10` effectively disables the production
limiter.

The salt energy equation is advection-diffusion with the turbulent heat flux modelled by gradient
diffusion with a turbulent Prandtl number $Pr_t$ = 0.85:

$$\nabla\cdot(\rho c_p\,\mathbf u\,T) = \nabla\cdot\bigl(k_\text{eff}\nabla T\bigr) + q'''_\text{salt},
\qquad k_\text{eff} = k_f + \frac{\mu_t c_p}{Pr_t}.$$

The graphite and the guide tube solve steady conduction,
$-\nabla\cdot(k_s\nabla T_s) = q'''_\text{graphite}$, with $k_s$ = 25 W/(m K) (SPEC). The guide tube
uses the graphite conductivity in the input (`k_graphite` on the whole `solid` block list), not a
Hastelloy N value.

### 3.2 Wall functions and the friction velocity

The near-wall cell (centre at $y$ = 1.5 mm) cannot resolve the viscous sublayer. A wall function
assumes the universal **logarithmic law of the wall** between the cell and the wall:

$$u^+ = \frac{u_\parallel}{u_\tau} = \frac{1}{\kappa}\ln\bigl(E\,y^+\bigr), \qquad
y^+ = \frac{\rho\,u_\tau\,y}{\mu}, \qquad u_\tau = \sqrt{\tau_w/\rho} ,$$

where $u_\parallel$ is the cell velocity parallel to the wall, $\tau_w$ the wall shear stress, and
$\kappa$ = 0.4187, $E$ = 9.793 (`NS::von_karman_constant`, `NS::E_turb_constant`). With the
**equilibrium Newton** treatment (`wall_treatment = eq_newton`), $u_\tau$ solves

$$f(u_\tau) = \frac{u_\tau}{\kappa}\ln\Bigl(\frac{E\,u_\tau\,y}{\nu}\Bigr) - u_\parallel = 0, \qquad
u_\tau^{(n+1)} = u_\tau^{(n)} - \frac{f}{f'},\quad
f' = \frac{1}{\kappa}\Bigl[1 + \ln\Bigl(\frac{E\,u_\tau y}{\nu}\Bigr)\Bigr],$$

with $\nu = \mu/\rho$ (`NS::findUStar`, relative tolerance 1e-6, usually 3-4 iterations). The wall
shear $\tau_w = \rho u_\tau^2$ enters the momentum equation through a wall viscosity
(`LinearFVTurbulentViscosityWallFunctionBC`). In the near-wall cells (log-layer branch), the $k_t$ production and
dissipation take the equilibrium forms of Launder & Spalding (1974),
$G_k = \tau_w\,C_\mu^{1/4}k_t^{1/2}/(\kappa y)$ and $\varepsilon = C_\mu^{3/4}k_t^{3/2}/(\kappa y)$
(`LinearFVTKESourceSink`, `LinearFVTKEDSourceSink`).
The equilibrium treatment assumes that production balances dissipation in the wall cell, which
holds in a developed channel flow.

### 3.3 The thermal law of the wall and the Jayatilleke P-function

The thermal analogue of the log law, for a wall at temperature $T_w$ and a cell at $T_P$, is

$$T^+ = \frac{(T_w - T_P)\,\rho c_p u_\tau}{q_w''} =
Pr_t\Bigl[\frac{1}{\kappa}\ln\bigl(E\,y^+\bigr) + P\Bigl(\frac{Pr}{Pr_t}\Bigr)\Bigr],$$

with $q_w''$ the wall heat flux. The **P-function** of Jayatilleke (1969) is

$$P = 9.24\Bigl[\Bigl(\frac{Pr}{Pr_t}\Bigr)^{3/4} - 1\Bigr]\Bigl[1 + 0.28\,e^{-0.007\,Pr/Pr_t}\Bigr].$$

The physical picture: very close to the wall, heat moves only by molecular conduction, across a
**conductive sublayer** whose thickness in wall units scales roughly as $Pr^{-1/3}$ of the viscous
sublayer. Outside it, the turbulent transport takes over. The log term carries the turbulent part
of the resistance; $P$ is the extra resistance of the sublayer relative to a fluid with
$Pr = Pr_t$ (for which $P = 0$). The wall-function heat-transfer coefficient is

$$h_{wf} = \frac{q_w''}{T_w - T_P} = \frac{\rho c_p u_\tau}{T^+}\quad (y^+\ge 30),\qquad
h_{wf} = \frac{k_f}{y}\quad (y^+\le 5),$$

blended linearly in $y^+$ between 5 and 30, with $T^+$ evaluated at $y^+ = 30$ in the blend
(`wallFunctionHTC()` in `LinearFVTurbulentWallCHTBC.C`). The flow quantities are lagged by one
iteration, as in the momentum wall function.

**Why a high-Pr salt needs it.** For the GFMR salt, $Pr/Pr_t = 16.8/0.85 \approx 19.8$, and the
formula gives $P \approx 96$. At $y^+$ = 30-45 the log term $\ln(Ey^+)/\kappa$ is only about 14.
**The sublayer carries about 7/8 of the total thermal resistance between the wall and the cell.**
A plain conduction treatment of the near-wall cell, with the face flux
$k_\text{eff}(T_P - T_w)/y$ and $k_\text{eff} = k_f + \mu_t c_p/Pr_t$ evaluated at the cell
centre, uses the eddy conductivity of the log layer across the whole distance $y$ and omits the
sublayer. For the GFMR inner channel (first cell at $y^+ = 40$, $\mu_t = \kappa y^+ \mu$) it gives
$k_\text{eff}/y \approx 1.8\times10^5$ W/m²/K against $h_{wf} \approx 3.8\times10^3$ W/m²/K, a 46-fold overestimate. For a
gas or water ($Pr$ of order 1), $P$ is small and the difference is minor. That is why
standard CFD often gets away without it, and why a molten salt does not.

**Verification of the wall law.** `test/tests/kEpsilon/salt_pipe_heat_transfer` runs a GFMR inner
channel (FUNaK, 4 cm, 1.7 m/s, $Re = 1.85\times10^4$, $Pr$ = 16.8) inside a uniformly heated
graphite sleeve, with the first cell at $y^+$ = 33. It gives $Nu$ = 180 against 194 from the
Gnielinski correlation (Gnielinski 1976),

$$Nu = \frac{(f/8)(Re - 1000)\,Pr}{1 + 12.7\sqrt{f/8}\,\bigl(Pr^{2/3} - 1\bigr)} ,$$

with $f$ the Darcy friction factor. The ratio 0.93 is inside the ±10-20 % uncertainty usually
quoted for the correlation. This is a **CODE-TO-CODE** check against a correlation, not a
validation with salt data.

### 3.4 Two-sided conjugate coupling: `LinearFVTurbulentWallCHTBC`

The salt and the solid are separate linear systems (`energy_system`, `solid_energy_system`),
solved one after the other in each SIMPLE iteration. The interface `cht_interfaces` of the
executioner exchanges boundary data between them. The same BC object is placed on both sides of
each wall, selected with `side`.

Between the near-wall fluid cell P and the near-wall solid cell S there are **two resistances in
series**: the wall function $1/h_{wf}$ (P to wall) and the half solid cell $d_S/k_s$ (wall to S).
The overall coefficient is

$$U = \Bigl(\frac{1}{h_{wf}} + \frac{d_S}{k_s}\Bigr)^{-1}.$$

- **Fluid side.** Outward flux $h_{wf}(T_P - T_w)$, with the solid surface temperature $T_w$
  handed over from the solid (`interface_temperature_solid_<boundary>`). The fluid hands $T_P$
  to the solid.
- **Solid side.** Outward flux $U(T_S - T_P)$, with $T_P$ from the fluid
  (`interface_temperature_fluid_<boundary>`). The solid hands back the surface temperature

$$T_w = T_S + U\,\frac{d_S}{k_s}\,(T_P - T_S),$$

which is the temperature at the junction of the two series resistances.

Each side is **implicit in its own cell temperature** (the coefficient $h_{wf}$ or $U$ goes to the
matrix diagonal) and lags only the other side's temperature. At convergence the two lagged values
equal the current ones, and both fluxes are $U(T_S - T_P) = h_{wf}(T_w - T_P)$, so the heat flow is
conserved exactly across the interface.

**Why it is stable without a relaxation coefficient.** The classical Dirichlet-Neumann partition
(fluid gets the solid temperature, solid gets the fluid flux) leaves the solid with only a flux
condition on its cooled surfaces. In the GFMR graphite, which is cooled only by the salt (the outer
surfaces are adiabatic), that makes the solid problem a pure Neumann problem whose temperature level
is undetermined. The fixed-point iteration then oscillates or diverges unless it is relaxed, and the
relaxation coefficient becomes a numerical parameter. Robin-Robin schemes with an arbitrary
coefficient avoid this but need a coefficient that is tuned for the convergence rate. Here each side
has a Robin condition whose coefficient is the **physical** one ($h_{wf}$ or $U$), and whose
reference temperature is a temperature the other side actually has. Neither sub-problem loses its
temperature reference, and no free coefficient enters the converged answer. The input therefore runs
with `max_cht_fpi = 1` (one exchange per SIMPLE iteration) and all CHT relaxation factors equal to 1.
The residual risk is the usual one of a lagged coupling: it converges at the rate of the outer
SIMPLE loop, not in one step.

### 3.5 The SIMPLE segregated solution

The linear FV system is solved with the SIMPLE algorithm (Patankar & Spalding 1972; Patankar 1980).
Each iteration:

1. **Momentum predictor.** Assemble and solve the three momentum components with the lagged
   pressure, mass fluxes and $\mu_t$, under-relaxed by `momentum_equation_relaxation` = 0.7.
2. **Pressure correction.** Form $\mathbf H/\mathbf A$ (the momentum operator without the pressure
   gradient, divided by its diagonal) and its face flux by **Rhie-Chow interpolation**
   (`RhieChowMassFlux`, Rhie & Chow 1983). Rhie-Chow removes the checkerboard pressure mode of a
   collocated grid. Solve the pressure Poisson equation
   $\nabla\cdot(\mathbf A^{-1}\nabla p) = \nabla\cdot(\mathbf H/\mathbf A)$, relax the pressure by
   `pressure_variable_relaxation` = 0.3, and correct the face mass fluxes so that they are
   divergence-free.
3. **Turbulence.** Solve $k_t$ then $\varepsilon$ (equation and field relaxation 0.5), update
   $\mu_t$ (`kEpsilonViscosityAux`).
4. **Energy.** Solve the salt energy (relaxation 0.95) and the solid energy, with the CHT exchange
   of section 3.4 between them.

The run uses 600 iterations (24 ranks, about 20 min). The momentum residuals reach 1e-6. The energy
residual levels off near 3e-5, while the interface heat flows are steady to 1e-5. The advected
quantities use first-order upwind interpolation, which is diffusive but robust; the nonorthogonal
correction is off.

### 3.6 Heat sources and the energy-balance verification

The diffusion power map (`sector_power.e`, the same 540 k-cell mesh) provides the salt power
density, rescaled so that the salt receives 92.1 % of the 25 MW sector power. The graphite receives
7.9 % (the OpenMC coupled n-γ fraction), distributed in proportion to the fast flux (groups 1-2,
$E >$ 5.53 keV). This shape is an assumption: gamma heating has a different and longer range.

The global balance for the sector is

$$\epsilon_E = \frac{\dot H_\text{out} - \dot H_\text{in} - (Q_\text{salt} + Q_\text{graphite})}
                   {Q_\text{salt} + Q_\text{graphite}}, \qquad
\dot H = \int_A \rho c_p T\,\mathbf u\cdot\mathbf n\,dA ,$$

(`energy_balance`; `enthalpy_in` is negative because the flow enters). It closes only if the
graphite heat leaves through the interfaces into the salt, so it tests the CHT exchange as well as
the energy equations. Results (stage-3 README):

| check | result | class |
|---|---|---|
| energy balance $\epsilon_E$ | 4.5×10⁻⁵ | VERIFICATION |
| salt / graphite heat input of the sector | 23.03 / 1.97 MW (25.0 MW = 150/6) | VERIFICATION |
| graphite heat leaving through the four walls | 1.97 MW (fluxes balanced to 1e-5) | VERIFICATION |
| core mass flow (6 × sector) | 2294 kg/s (SPEC 2300) | consistency with SPEC velocities |
| mixed outlet temperature | 686.4 °C (heat balance 686.4 °C) | VERIFICATION |
| channel Δp inner / outer / rod | 14.4 / 6.2 / 35.4 kPa vs SPEC 14.3 / 6.0 / 34 kPa | CODE-TO-CODE |

The core mass flow is 2294 kg/s, not 2300, because the specified zone velocities at the 650 °C
density give 2294 kg/s; the velocities are imposed, not the flow rate.

### 3.7 Results to interpret

The zoned velocities equalize the outlet temperatures: the inner channels make 31 % of the power with
22 % of the channels, and at 1.70 m/s they leave at the same mixed temperature as the outer channels
at 1.03 m/s, within 1.5 K (685.4 / 686.8 / 685.3 °C inner / outer / rod). The hottest salt cell is at
719.5 °C. The graphite moderator averages 721 °C with a maximum of 830 °C; the radial reflector
reaches 859 °C at its adiabatic outer surface (a sensitivity to the two modelling choices of
section 3.6 and the adiabatic boundary, not a prediction).

The neutronics assumed graphite at 800 °C. The CFD moderator average, 721 °C, is about 80 K lower.
With $\alpha_g$ = −4.1 pcm/K this is worth about +330 pcm, which would bring the core close to
critical. The coupling is one-way here (frozen power map, no feedback), so this is an estimate,
not a coupled result.

## 4. Stage 4: chemistry

Stage 4 is described in `04_chemistry/README.md`. This section states the physics it uses. The shared models (equations, databases, numerical treatment, limitations) are
documented in `examples/training/chemistry/THEORY.md` and are not repeated here.

### 4.1 Radiolysis in a fluoride fuel salt

In a fuel salt the dose is deposited mostly by the fission fragments in the salt itself, not by
gammas from outside. The dose rate is therefore close to the fission power density,
$\dot D \approx (0.85\text{-}0.9)\,q'''$ (chemistry THEORY section 2.1). The stage-4 plan uses
about $0.87\,q'''$ from the stage-2 power density. The radiolytic source of a species $i$ is
$S_i = G_i\,\dot D/(100\ \text{eV}\cdot N_A)$, with $G_i$ the yield in molecules per 100 eV and
$N_A$ Avogadro's number (chemistry THEORY section 2.2). The primary products in a fluoride are
solvated electrons and fluorine radicals. At reactor temperature they recombine within about
1e-10 s, so the net radiolytic effect on the salt is small, unlike in a frozen salt. The database
yields are gamma yields; fission-fragment yields in molten fluorides are not available (chemistry
THEORY section 2.9).

### 4.2 Fission oxidation and redox buffering

The net oxidant of a fuel salt is fission. A fissioned UF4 frees four fluorides; the fission
products bind fewer of them, so each fission leaves an excess of oxidizing equivalents (about 0.8
per fission in the ORNL MSRE experience used by the chemistry training, section 2.8). Without a
reductant the free fluorine potential rises, and the salt oxidizes the structural metals. The
classical buffer is a small U(III)/U(IV) ratio: U(III) is oxidized to U(IV) first, and the fluorine
potential stays low while U(III) lasts. A plant restores the ratio by adding a reductant (Be in the
MSRE). The redox potential set by $\text{U(III)}/\text{U(IV)}$ fixes the equilibrium concentration
of CrF2 in the salt, and hence the driving force for corrosion.

### 4.3 Element-resolved corrosion and plating

Hastelloy N corrodes mainly by selective dissolution of Cr, the least noble major element:
$\text{Cr} + 2\,\text{UF}_4 \rightleftharpoons \text{CrF}_2 + 2\,\text{UF}_3$. The dissolution
equilibrium depends on temperature. In a loop with a temperature difference, the hot walls dissolve
Cr until the salt reaches their equilibrium, and the cold walls, where the equilibrium is lower,
plate it out. The salt settles between the two equilibria and carries metal from hot to cold at a
steady rate set by the wall mass transfer and the kinetics. `[CorrosionPlatingFlow]` represents
this with a Butler-Volmer wall law per element whose equilibrium concentration follows a van 't Hoff
dependence on the wall temperature (chemistry THEORY section 4.4):

$$c_{\text{eq},k}(T) = c_{\text{eq},k}\exp\Bigl[-\frac{\Delta H_k}{R}\Bigl(\frac1T - \frac1{T_\text{ref}}\Bigr)\Bigr],$$

with $k$ the element (Cr, Fe, Ni), $\Delta H_k$ the dissolution enthalpy and $R$ the gas constant.
Noble fission products (Mo, Ru, Tc, ...) also deposit on metal and graphite surfaces (chemistry
THEORY section 4.7).

### 4.4 GFMR-specific application

**Once-through core, loop balance.** The salt crosses the GFMR core in $t_\text{core} = 3.02$ s,
while the chemistry evolves over minutes (I-135, Xe-135) to days (redox). The 3-D core model
needs inlet values, which come from a 0-D loop balance (`loop_balance.py`) with an ASSUMED loop:
ex-core volume $2\times$ the core salt ($t_\text{loop} = 9.06$ s), U(III)/U held at 1 %, and 10 %
Xe stripping per pass. For a species with a core source $S$ and loop-average removal rate $r$,

$$\bar c = \frac{S\,V_\text{core}/V_\text{loop}}{\lambda + r},\qquad
c_\text{in} = \bar c - \tfrac12 S\,t_\text{core}\,\frac{V_\text{core}}{V_\text{loop}} .$$

**Verification.** The 3-D core rise of I-135 and Xe-135 (flow-weighted, last minus first layer)
matches the 0-D core-pass increment to 4 % and 5 % (VERIFICATION).

**The numerical floor.** A source-free species (`F_ion`) leaves the 3-D core with a mixed relative
error of $1.1\times10^{-5}$, and up to 0.6 % in single inlet-corner cells. The per-pass redox change
($\Delta[\mathrm{U(III)}]/[\mathrm{U(III)}] = 1.1\times10^{-7}$) and the guide-tube Cr exchange
($<10^{-6}$) are below this floor. Their 3-D fields therefore carry no information, and the redox
and corrosion results are taken from the 0-D loop. The fission products rise by
$1.6\times10^{-4}$ (I-135) per pass, above the floor.

**Results.** Fission oxidation is 0.54 mol eq/day for the whole loop. Holding U(III)/U at 1 % needs
the same amount of reductant: 43 g/day of U metal ($\mathrm{U + 3\,UF_4 \to 4\,UF_3}$) or
2.4 g/day of Be. The loop-average Xe-135 concentration is $8.47\times10^{-6}$ mol/m³. Stripping
($1.2\times10^{-2}$ s⁻¹) dominates burnout ($8.1\times10^{-5}$ s⁻¹) and decay
($2.1\times10^{-5}$ s⁻¹), so the xenon poison is small. It is therefore set mostly by the assumed
stripping efficiency.

## 5. Stage 5: transients

Stage 5 is described in `05_transients/README.md`. One scenario is done: an unprotected 0.3 $ rod step.

### 5.1 Point kinetics and space-time diffusion

The time-dependent multigroup diffusion equation with $I$ = 6 precursor groups is

$$\frac{1}{v_g}\frac{\partial\phi_g}{\partial t} = \nabla\cdot(D_g\nabla\phi_g) - \sum_{g'}R_{gg'}\phi_{g'}
 + (1-\beta)\,\chi_{p,g}\sum_{g'}\nu\Sigma_{f,g'}\phi_{g'} + \sum_{i=1}^{I}\chi_{d,g}\lambda_i C_i ,$$

$$\frac{\partial C_i}{\partial t} + \nabla\cdot(\mathbf u\,C_i) = \beta_i\sum_{g'}\nu\Sigma_{f,g'}\phi_{g'} - \lambda_i C_i
 \qquad (+\ \text{turbulent diffusion of } C_i).$$

The advection term $\nabla\cdot(\mathbf u C_i)$ is what distinguishes a circulating-fuel reactor:
the precursors move with the salt. `xs/gfmr_xs.i` already contains the $1/v_g$ (`INV_V_*`), the
prompt and delayed spectra (`CHI_P`, `CHI_D`), the prompt fission matrix (`F_prompt`), `BETA` and
`LAMBDA` needed for these equations. Solving them in space and time is **space-time kinetics**.

**Point kinetics** factors the flux into an amplitude and a fixed shape,
$\phi(\mathbf r,E,t) = n(t)\,\psi(\mathbf r,E)$, and projects the equations onto the adjoint of the
shape (Keepin 1965; Ott & Neuhold 1985):

$$\frac{dn}{dt} = \frac{\rho(t) - \beta_\text{eff}}{\Lambda}\,n + \sum_i \lambda_i c_i, \qquad
\frac{dc_i}{dt} = \frac{\beta_{\text{eff},i}}{\Lambda}\,n - \lambda_i c_i\ (-\ \text{loop loss}),$$

with $n$ the amplitude (power) and $c_i$ the weighted precursor amplitudes. It is accurate when the
shape barely changes (small, uniform perturbations). It is poor when the perturbation is local, as
for a single rod, or when the shape changes strongly, as for the axial shift of the power when the
rods (at the top of the core) move. Space-time diffusion captures these at a much higher cost. With
$\beta_\text{eff}$ = 681 pcm and $\Lambda$ = 0.605 ms, the prompt-critical threshold of the
stationary-fuel core is a reactivity of 681 pcm (1 \$). The total rod worth, 4844 pcm, is 7.1 \$.

### 5.2 Precursor drift in a circulating-fuel reactor

In a stationary-fuel core all delayed neutrons are born where their precursors were made. In a
circulating-fuel reactor the salt carries the precursors (a) downstream inside the core, toward the
outlet where the importance is lower, and (b) out of the core into the external loop, where they
decay without contributing to the chain reaction. Both effects reduce the effective delayed
fraction. At steady state the core must be supercritical by the lost delayed fraction to remain
critical: the **reactivity loss from precursor drift**.

A lumped (0-D) estimate treats the core as well mixed, with a core residence time $\tau_c$ and a
loop transit time $\tau_\ell$. The precursors that leave the core and survive the loop
(a fraction $e^{-\lambda_i\tau_\ell}$) re-enter it. The steady precursor balance of the core,
$0 = \beta_i F - \lambda_i C_i - C_i/\tau_c + C_i e^{-\lambda_i\tau_\ell}/\tau_c$ with $F$ the fission
neutron production rate, gives the fraction of group $i$ that decays inside the core:

$$\frac{\beta_{i,\text{circ}}}{\beta_i} =
\frac{\lambda_i}{\lambda_i + \dfrac{1}{\tau_c}\bigl(1 - e^{-\lambda_i\tau_\ell}\bigr)} ,\qquad
\Delta\rho_\text{drift} \approx \sum_i \beta_{\text{eff},i}\Bigl(1 - \frac{\beta_{i,\text{circ}}}{\beta_i}\Bigr).$$

Long-lived groups ($\lambda_i\tau_\ell \ll 1$) circulate many times and decay where they happen to
be: the ratio tends to $\tau_c/(\tau_c + \tau_\ell)$, the fraction of time the salt spends in the
core. Short-lived groups ($\lambda_i\tau_c \gg 1$) decay in the core, and the ratio tends to 1. The lumped form ignores the importance gradient along
the channel. The space-time model with the advection term captures that gradient directly, with the
salt velocity of stage 3. The GFMR core residence time follows from the channel length (3.6 m) and
the zone velocities, but **the loop transit time is not in the specification** and must be assumed.
With $\tau_c = 3.02$ s and an ASSUMED ex-core transit $\tau_\ell = 6.04$ s, the plug-flow loop of the
formula above gives a drift loss of 325 pcm. The first-order-lag loop of the stage-5 ODEs, in which
a fraction $(1/\tau_\ell)/(1/\tau_\ell + \lambda_i)$ survives the loop, gives 295 pcm (0.43 $). The
30 pcm difference is the effect of the loop's residence-time distribution. The loss also makes the transients depend on
flow: when the flow falls, fewer precursors leave the core, and the recovered delayed neutrons add
reactivity.

### 5.3 Temperature feedback

The reactivity feedback is

$$\rho_\text{fb}(t) = \alpha_f\,\bigl(\bar T_f(t) - \bar T_{f,0}\bigr) + \alpha_g\,\bigl(\bar T_g(t) - \bar T_{g,0}\bigr),$$

with $\bar T_f$, $\bar T_g$ the (importance-weighted) mean salt and graphite temperatures, and the
OpenMC coefficients $\alpha_f$ = −2.7 ± 0.2 pcm/K and $\alpha_g$ = −4.1 ± 0.2 pcm/K. The two act on
different time scales. The fuel feedback is prompt: the fission heat is deposited in the salt
itself. The graphite feedback is delayed by the heat capacity of the moderator: 1700 kg/m³ (SPEC)
times the **assumed** 1900 J/(kg K), heated by 7.9 % of the power and by conduction from the salt.
In a space-time model the coefficients are replaced by temperature-dependent cross sections, or
applied locally; the planned model uses porous-medium salt and graphite energy equations on the
homogenized full-core mesh.

### 5.4 Scenarios

| scenario | initiator | dominant physics |
|---|---|---|
| rod withdrawal | rods move out from a deep insertion (the nominal 36 cm insertion is worth < 50 pcm, so the full-insertion worth of 4844 pcm is the relevant one) | reactivity insertion vs prompt fuel feedback; the axial shape shift |
| loss of flow | pump coast-down | salt heats up (negative fuel feedback); less precursor drift (positive reactivity) |
| overcooling | colder salt at the inlet | positive reactivity from the negative coefficients; power rise |
| loss of heat sink | no heat removal in the loop | inlet temperature rises with the loop transit time; power falls by feedback |

The last three need the lumped loop with a transit time, which is not part of the specification.

### 5.5 GFMR-specific application

`gfmr_rod_step.i` solves the point-kinetics form of 5.1 with six precursor groups in the core and
six in the loop, and two temperature nodes (salt, graphite). It uses the OpenMC kinetics and
feedback coefficients, and the thermal constants of the stage-3 CFD: graphite-salt conductance
$UA = 1.66\times10^5$ W/K, graphite time constant 438 s.

**Drift and prompt criticality.** Drift leaves $\beta - \rho_0 = 386$ pcm of delayed neutrons in
the core. A 0.3 $ step (204 pcm, in units of the static $\beta$) is therefore 0.53 $ of the
effective delayed fraction, and the core becomes prompt critical at 0.57 $ of the static $\beta$.

**Verification.** The null transient holds the steady state to 8 digits over 300 s. The final
power after the step (1.3789 at 600 s) matches the analytic feedback balance
$\rho_\text{step} + \alpha_f\Delta T_s + \alpha_g\Delta T_g = 0$, which gives 1.3797.

**The transient.** The power jumps to about 2.1× in the first second, as the prompt jump with the
effective delayed fraction predicts. The fuel coefficient turns it over at 2.78× nominal (417 MW)
at 15 s; the core-average outlet then peaks at 798 °C. The graphite, heating on its 438 s time
constant, brings the power down to 1.38× over about 10 minutes, with the salt 13.8 K and the
graphite 40.8 K above their initial values.

## 6. References

- J. J. Duderstadt, L. J. Hamilton, *Nuclear Reactor Analysis*, Wiley (1976). Multigroup diffusion,
  Marshak and extrapolated boundary conditions, homogenization, control-rod treatment.
- R. E. Marshak, "Note on the spherical harmonic method as applied to the Milne problem for a
  sphere", Physical Review 71, 443 (1947). The Marshak vacuum condition.
- G. R. Keepin, *Physics of Nuclear Kinetics*, Addison-Wesley (1965). Delayed neutrons, point kinetics.
- K. O. Ott, R. J. Neuhold, *Introductory Nuclear Reactor Dynamics*, American Nuclear Society
  (1985). Point kinetics, feedback, circulating-fuel precursor balance.
- H. Hurwitz, "Physical interpretation of the adjoint flux: iterated fission probability", in
  *Naval Reactors Physics Handbook*, vol. 1 (1964). The IFP interpretation of the adjoint.
- Y. Nauchi, T. Kameyama, "Development of calculation technique for iterated fission probability and
  reactor kinetic parameters using continuous-energy Monte Carlo method", Journal of Nuclear Science
  and Technology 47, 977 (2010).
- B. C. Kiedrowski, F. B. Brown, P. P. H. Wilson, "Adjoint-weighted tallies for k-eigenvalue
  calculations with continuous-energy Monte Carlo", Nuclear Science and Engineering 168, 226 (2011).
- P. K. Romano et al., "OpenMC: A state-of-the-art Monte Carlo code for research and development",
  Annals of Nuclear Energy 82, 90 (2015). The OpenMC IFP tallies (`ifp-*` scores) are described in
  the OpenMC 0.16 user documentation.
- D. A. Brown et al., "ENDF/B-VIII.0: The 8th major release of the nuclear reaction data library",
  Nuclear Data Sheets 148, 1 (2018).
- B. E. Launder, D. B. Spalding, "The numerical computation of turbulent flows", Computer Methods in
  Applied Mechanics and Engineering 3, 269 (1974). Standard k-ε model and wall functions.
- C. L. V. Jayatilleke, "The influence of Prandtl number and surface roughness on the resistance of
  the laminar sub-layer to momentum and heat transfer", Progress in Heat and Mass Transfer 1, 193
  (1969). The P-function.
- V. Gnielinski, "New equations for heat and mass transfer in turbulent pipe and channel flow",
  International Chemical Engineering 16, 359 (1976).
- S. V. Patankar, D. B. Spalding, "A calculation procedure for heat, mass and momentum transfer in
  three-dimensional parabolic flows", International Journal of Heat and Mass Transfer 15, 1787
  (1972); S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere (1980). SIMPLE.
- C. M. Rhie, W. L. Chow, "Numerical study of the turbulent flow past an airfoil with trailing edge
  separation", AIAA Journal 21, 1525 (1983). Rhie-Chow interpolation.
- Saltfoss Energy ApS, "Generic FUNaK-Fuelled Thermal Molten Salt Reactor (GFMR) Model",
  https://saltfossenergy.github.io/GFMR/, DOI 10.5281/zenodo.23348200.
- In the repository: `examples/training/chemistry/THEORY.md`;
  `doc/content/source/turbulence/linearfvbcs/LinearFVTurbulentWallCHTBC.md`;
  `test/tests/kEpsilon/salt_pipe_heat_transfer`; the stage READMEs of this training.
