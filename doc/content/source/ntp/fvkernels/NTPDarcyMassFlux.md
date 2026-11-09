# NTPDarcyMassFlux

`NTPDarcyMassFlux` adds the gas mass conservation $\nabla \cdot \dot{\mathbf{m}} = 0$ to the pressure
variable of a particle-bed element. The face mass flux, positive along the face normal, is

\begin{equation}
\dot{m}_f = -(\rho M)_f \left[ \nabla p + \nabla J \right]_f \cdot \mathbf{n}_f ,
\end{equation}

with $\rho M$ the diagonal mass mobility (`rho_mobility`, e.g. from
[NTPPackedBedFunctorMaterial](NTPPackedBedFunctorMaterial.md) or
[NTPDuctFlowFunctorMaterial](NTPDuctFlowFunctorMaterial.md)). Without a momentum flux $J$, the kernel is
an anisotropic diffusion of the pressure with coefficient $\rho M$, as in
[FVAnisotropicDiffusion](FVAnisotropicDiffusion.md). It differs in how the face mobility is evaluated
and interpolated, as described below. The same
face flux is shared with [NTPDarcyAdvection](NTPDarcyAdvection.md),
[NTPDarcyMassFluxFunctorMaterial](NTPDarcyMassFluxFunctorMaterial.md) and
[NTPDarcyFluxIntegral](NTPDarcyFluxIntegral.md), so mass and energy see identical fluxes.

- +Face mobility+ (`face_mobility = true`, default). The mobility of each side is evaluated at the
  face with the compact two-point face pressure gradient. Otherwise it is evaluated in the cell with the
  reconstructed cell gradient. The face evaluation is exact for the flow across a porous layer one cell
  thick, bounded by porosity jumps.
- +Face interpolation+ (`coeff_interp_method`). The default `series` combines the two sides as series
  resistances of the half segments between the cell centers,
  $(\rho M)_f = 1 / \left[ (1 - g_C)/(\rho M)_C + g_C/(\rho M)_N \right]$. The MOOSE harmonic average
  uses the linear interpolation weights instead, which attribute each half segment to the wrong
  cell. The two agree on uniform meshes only. The `harmonic` and `average` options are available for
  comparison.
- +Momentum flux+ (`momentum_flux`, `momentum_blocks`). With $J = \rho u_z^2$, the flux relation
  in a slender free channel becomes the steady one-dimensional momentum balance
  $d(\rho u_z^2)/dz + dp/dz = -(\text{friction})\, u_z$. This balance carries the acceleration of
  heated gas and the pressure recovery of a dead-ended manifold. $J$ acts only on internal faces
  whose two cells are both in `momentum_blocks`, through its axial ($y$) jump between the cell
  centers.
- +Orifice solved with the flow+ (`orifice_log_multiplier`, `mesh_division`, `orifice_blocks`).
  The scalar variable $\lambda$ has one component per slice of the mesh division. The mobility of a cell in
  `orifice_blocks` and slice $j$ is divided by $\exp(\lambda_j)$ before the face interpolation. The
  kernel reads the scalar so that its derivatives enter the Jacobian. It is used with
  [NTPSliceUniformityConstraint](NTPSliceUniformityConstraint.md).

On a boundary face, the mobility is that of the cell inside the domain.

!syntax parameters /FVKernels/NTPDarcyMassFlux

!syntax inputs /FVKernels/NTPDarcyMassFlux

!syntax children /FVKernels/NTPDarcyMassFlux
