# NTPDuctFlowFunctorMaterial

`NTPDuctFlowFunctorMaterial` gives the mobility of a free channel (the plenum or the chamber of a
particle-bed element) for a Darcy-type description of the flow, $\mathbf{u} = -M \nabla p$, with $y$
the channel axis (RZ, or 3-D with $x$ and $z$ transverse). It is the duct counterpart of
[NTPPackedBedFunctorMaterial](NTPPackedBedFunctorMaterial.md) in the Newton (AD) path.

Along the axis, the friction of a duct of hydraulic diameter $D_h$ is the laminar Hagen-Poiseuille
term plus the Darcy-Weisbach term with friction factor $f$:

\begin{equation}
-\frac{\partial p}{\partial z} = \left( \frac{32 \mu}{D_h^2} + \frac{f \rho |u_z|}{2 D_h} \right) u_z,
\qquad
M_z = \frac{2}{a + \sqrt{a^2 + 4 b |\partial p / \partial z|}}, \quad a = \frac{32 \mu}{D_h^2}, \quad
b = \frac{f \rho}{2 D_h}.
\end{equation}

The model neglects inertia (acceleration). Across the channel, the mobility is
`radial_mobility_ratio` (default 100) times $M_z$: a free channel offers almost no resistance to the
radial redistribution that feeds or drains the porous layers. The object outputs the mass mobility
vector $\rho M$ (`rho_mobility_name`), the velocity $-M \nabla p$ (`velocity_name`) and the axial
speed $M_z |\partial p / \partial z|$ (`speed_name`). `hydraulic_diameter` is a functor and
`friction_factor` is a constant.

!syntax parameters /FunctorMaterials/NTPDuctFlowFunctorMaterial

!syntax inputs /FunctorMaterials/NTPDuctFlowFunctorMaterial

!syntax children /FunctorMaterials/NTPDuctFlowFunctorMaterial
