# NTPDuctFrictionSegregatedFunctorMaterial

`NTPDuctFrictionSegregatedFunctorMaterial` gives the wall friction of a free channel (the plenum or
the chamber of a particle-bed element) for the segregated (SIMPLE/PIMPLE) solvers. The friction is
the Darcy coefficient vector `drag_name` of
[LinearFVMomentumFriction](LinearFVMomentumFriction.md), whose force per unit volume is
$\mu D \mathbf{u}$.

The channels are resolved by only a few cells across their gap, which cannot resolve the turbulent
wall layer. The wall shear is therefore modeled as the friction of a duct of hydraulic diameter
$D_h$, the same law as [NTPDuctFlowFunctorMaterial](NTPDuctFlowFunctorMaterial.md) uses in the
Newton path:

\begin{equation}
D_z = \frac{1}{\mu} \left( \frac{32 \mu}{D_h^2} + \frac{f \rho |u_z|}{2 D_h} \right), \qquad D_r = 0.
\end{equation}

The friction acts only on the flow along the channel axis ($y$). The radial velocity of a channel
cell is the cross flow that feeds or drains the adjacent porous layer, not duct flow. `laminar_term = false`
drops the Hagen-Poiseuille term. By default, $|u_z|$ is the absolute value of the axial velocity
functor `v`, lagged by one iteration. With the optional `speed` functor, it is that functor instead.
An under-relaxed speed, updated once per outer
iteration, keeps the Picard iteration on the turbulent term, whose coefficient grows with the
solution, from oscillating.

!syntax parameters /FunctorMaterials/NTPDuctFrictionSegregatedFunctorMaterial

!syntax inputs /FunctorMaterials/NTPDuctFrictionSegregatedFunctorMaterial

!syntax children /FunctorMaterials/NTPDuctFrictionSegregatedFunctorMaterial
