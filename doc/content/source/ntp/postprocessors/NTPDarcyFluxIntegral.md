# NTPDarcyFluxIntegral

`NTPDarcyFluxIntegral` integrates the Darcy mass flux over a boundary, optionally times a carried
quantity $q$ (`carried_quantity`, default 1):

\begin{equation}
\int_{\partial \Omega} -(\rho M) \nabla p \cdot \mathbf{n}\; q \, dA .
\end{equation}

With $q = 1$ the result is the mass flow rate through the boundary [kg/s]. With $q = h$ it is the
enthalpy flow rate [W]. `variable` is the finite-volume pressure variable and `rho_mobility` the
diagonal mass mobility. The mass flux uses the same face pressure gradient (including the
boundary-condition treatment) and the same one-sided mobility as
[NTPDarcyAdvection](NTPDarcyAdvection.md) on a boundary face. $q$ is taken from the cell inside the domain,
which is the upwind value of an outflow. The result is therefore exactly the flow that the discrete
conservation equations see leave the domain.

Use it on boundaries where the pressure is prescribed, such as an outlet. On a boundary with an
imposed mass flux, the reconstructed boundary gradient is not that flux, and the imposed value
itself should be reported instead.

!syntax parameters /Postprocessors/NTPDarcyFluxIntegral

!syntax inputs /Postprocessors/NTPDarcyFluxIntegral

!syntax children /Postprocessors/NTPDarcyFluxIntegral
