# NTPRhieChowFluxIntegral

`NTPRhieChowFluxIntegral` integrates the face mass fluxes stored by the Rhie-Chow object of the
segregated solvers (`rhie_chow_user_object`, a [RhieChowMassFlux](RhieChowMassFlux.md) such as
[NTPRhieChowMassFlux](NTPRhieChowMassFlux.md)) over a boundary, positive out of the domain,
optionally times a carried quantity $q$ (`carried_quantity`, default 1):

\begin{equation}
\sum_{f \in \partial \Omega} \dot{m}_f\, q_f\, A_f .
\end{equation}

$q_f$ is evaluated at the boundary face from the inside cell. For a variable, this is its boundary
value: the Dirichlet value of an inflow or the extrapolated value of an outflow. This is what the
advection kernels transport. Because these are the fluxes the discrete balances see, the mass flow
rate ($q = 1$) and, with the enthalpy as carried quantity, the enthalpy flow rate of a boundary
are exact.

!syntax parameters /Postprocessors/NTPRhieChowFluxIntegral

!syntax inputs /Postprocessors/NTPRhieChowFluxIntegral

!syntax children /Postprocessors/NTPRhieChowFluxIntegral
