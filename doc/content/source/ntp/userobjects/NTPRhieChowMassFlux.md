# NTPRhieChowMassFlux

`NTPRhieChowMassFlux` is a [RhieChowMassFlux](RhieChowMassFlux.md) for the segregated
(SIMPLE/PIMPLE) solvers in flows with steps of many orders of magnitude in the momentum diagonal,
as between the free channels and the packed beds of a particle-bed fuel element. The cell momentum
balance is $\mathbf{u}_C = -M_C (\nabla p_C + \mathbf{f}_C)$, with $M = V/A$ the cell mobility and
$\mathbf{f} = \mathbf{H}/V$ the force density of the off-diagonal and source terms. The base class
interpolates $\rho M$ and $\rho M \mathbf{f}$ linearly to the faces. Across a step in $M$, this gives
the face the mobility of the more permeable side, so a porous layer one cell thick loses its flow
resistance. With `face_coupling = series` (default), the internal-face coupling is replaced by the
series form of the two half segments between the cell centers:

\begin{equation}
\dot{m}_f = -(\rho M)_f \left[ (\nabla p)_f + \mathbf{f}_f \right] \cdot \mathbf{n}, \quad
(\rho M)_f = \frac{1}{(1 - g_C)/(\rho M)_C + g_C/(\rho M)_N}, \quad
\mathbf{f}_f = (1 - g_C)\, \mathbf{f}_C + g_C\, \mathbf{f}_N .
\end{equation}

This is the exact face flux of two uniform half segments in series, and it reduces to the base
interpolation for smooth fields. `face_coupling = base` keeps the base-class coupling.

- `series_pressure_gradient` (default true). The Green-Gauss cell pressure gradient of the momentum
  equations uses the linear average of the two cells as face pressure. Across a mobility jump,
  this is far from the face pressure of flux continuity, and the permeable cell gets a large
  spurious gradient and velocity. With this option, the gradient is corrected with
  $p_f = (w_C p_C + w_N p_N)/(w_C + w_N)$, $w = \rho M / d$ being the half-segment conductances.
  The correction is applied once per new gradient: before each momentum assembly, and at the end of a
  solve for the output. The cell velocity is then recomputed.
- `isotropic_zero_flux_boundaries` (default true). On boundaries with a
  [LinearFVExtrapolatedPressureBC](LinearFVExtrapolatedPressureBC.md) or a
  [LinearFVPressureSymmetryBC](LinearFVPressureSymmetryBC.md), the face coefficient is
  replaced by its normal component, made isotropic. Otherwise, the tangential part of the anisotropic
  pressure diffusion carries a mass flux through walls and symmetry planes that are not aligned with
  the axes (the conical outer wall, the cut planes of a 3-D sector). All other boundary faces keep
  the base treatment.
- `normal_face_coupling` (default false). Reduces the coefficient of internal faces to its normal
  part, so the pressure diffusion has no tangential term. An option for meshes not aligned with the
  axes; the element inputs do not need it.

By default the object runs on `NONLINEAR` (before every linear system assembly) and `TIMESTEP_END`. It updates
the coupling fields of the pressure assembly after each computation of $\mathbf{H}/A$. Declare it under the name
of the flow Physics' Rhie-Chow object (by default `ins_rhie_chow_interpolator`); the Physics then
does not add its own.

!syntax parameters /UserObjects/NTPRhieChowMassFlux

!syntax inputs /UserObjects/NTPRhieChowMassFlux

!syntax children /UserObjects/NTPRhieChowMassFlux
