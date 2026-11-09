# NTPRhieChowMassFluxAux

`NTPRhieChowMassFluxAux` computes one component (`component` = `x`, `y` or `z`) of the cell
superficial mass flux $\mathbf{G} = \rho \mathbf{u}$ [kg/(m$^2$ s)] of the segregated solvers. It
reconstructs $\mathbf{G}$ from the face mass fluxes $\dot{m}_f$ of the Rhie-Chow object
`rhie_chow_user_object`, the same fluxes the discrete balances use:

\begin{equation}
\mathbf{G}_C = \frac{1}{V_C} \sum_f \dot{m}_f A_f \left( \mathbf{x}_f - \mathbf{x}_C \right).
\end{equation}

This is the same reconstruction that [NTPDarcyMassFluxFunctorMaterial](NTPDarcyMassFluxFunctorMaterial.md)
performs for the Newton Darcy model. In porous layers one cell thick, the cell velocity of the
segregated solvers comes from a pressure gradient reconstructed across porosity jumps and can be
far from the flux through the layer. The reconstructed $\mathbf{G}$ is the flow to use for the
interphase heat transfer and the inertial drag, through `mass_flux_x`/`mass_flux_y`/`mass_flux_z` of
[NTPPackedBedSegregatedFunctorMaterial](NTPPackedBedSegregatedFunctorMaterial.md). The variable
must be elemental.

!syntax parameters /AuxKernels/NTPRhieChowMassFluxAux

!syntax inputs /AuxKernels/NTPRhieChowMassFluxAux

!syntax children /AuxKernels/NTPRhieChowMassFluxAux
