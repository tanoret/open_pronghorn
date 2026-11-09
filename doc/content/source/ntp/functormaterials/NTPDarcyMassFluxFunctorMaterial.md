# NTPDarcyMassFluxFunctorMaterial

`NTPDarcyMassFluxFunctorMaterial` reconstructs the cell superficial mass flux
$\mathbf{G} = \rho \mathbf{u}$ [kg/(m$^2$ s)] of the Darcy-type flow of the Newton element model from
its face mass fluxes $\dot{m}_f$,

\begin{equation}
\mathbf{G}_C = \frac{1}{V_C} \sum_f \dot{m}_f A_f \left( \mathbf{x}_f - \mathbf{x}_C \right),
\end{equation}

and provides the three components as AD functors (`mass_flux_x_name`, `mass_flux_y_name`,
`mass_flux_z_name`; defaults `G_x`, `G_y`, `G_z`). The face fluxes are computed with the same
calculator as [NTPDarcyMassFlux](NTPDarcyMassFlux.md). The parameters `rho_mobility`,
`face_mobility`, `coeff_interp_method` and the orifice parameters (`orifice_log_multiplier`,
`mesh_division`, `orifice_blocks`) have the same meaning as in that kernel and should be set to the
same values.

The interphase heat transfer depends on $\mathbf{G}$ (`mass_flux_x`/`mass_flux_y` of
[NTPPackedBedFunctorMaterial](NTPPackedBedFunctorMaterial.md)). Through this AD functor, the
dependence of $\mathbf{G}$ on the pressure and the temperatures enters the Newton Jacobian. An auxiliary
variable would instead be a constant, lagged by one residual evaluation.

On the domain boundary, a face carries the kernel flux where the pressure is imposed (Dirichlet).
On `inflow_boundaries`, it carries the imposed inflow mass flux `inflow_mass_fluxes` (one functor
per boundary, positive into the domain). Other boundary faces carry no flux. The functors can only
be evaluated on cells or faces. On a face, they return the value of the cell on the requested
side, or of the face's element if no side is requested.

!syntax parameters /FunctorMaterials/NTPDarcyMassFluxFunctorMaterial

!syntax inputs /FunctorMaterials/NTPDarcyMassFluxFunctorMaterial

!syntax children /FunctorMaterials/NTPDarcyMassFluxFunctorMaterial
