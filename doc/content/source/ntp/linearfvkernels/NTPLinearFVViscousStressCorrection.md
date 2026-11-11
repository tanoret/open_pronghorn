# NTPLinearFVViscousStressCorrection

`NTPLinearFVViscousStressCorrection` adds the part of the Newtonian viscous stress of a compressible
gas that the Laplacian momentum flux $\mu \nabla \mathbf{u}$ omits, on the internal faces between
cells of its blocks:

\begin{equation}
\boldsymbol{\tau} = \mu \left(\nabla \mathbf{u} + \nabla \mathbf{u}^T
  - \tfrac{2}{3} (\nabla \cdot \mathbf{u}) \mathbf{I}\right),
\qquad
\boldsymbol{\tau}_\text{extra} \cdot \mathbf{n}
  = \mu \left(\nabla \mathbf{u}^T - \tfrac{2}{3} (\nabla \cdot \mathbf{u}) \mathbf{I}\right) \cdot \mathbf{n},
\end{equation}

for the momentum component $i$: $\mu (\partial_i u_j - \tfrac{2}{3} \delta_{ij} \nabla \cdot \mathbf{u}) n_j$.
The term is explicit: the cell gradients of the last iteration, interpolated linearly to the face.
In RZ the divergence includes $u_r / r$. The face viscosity is the linear average of the two cells
(`mu` is $\mu + \mu_t$ in a turbulent flow).

Faces with a cell outside the blocks and boundary faces get no contribution. The kernel is meant for
the free channels of a particle-bed element (the plenum and the chamber): across a porosity jump
the superficial velocity gradient is not a strain rate, and in the packed beds the drag exceeds the
viscous stress by many orders of magnitude. With `gradient_functors` the cell gradients are those
of other functors, e.g. the `<component>_channel` velocities of
[NTPChannelVelocity](NTPChannelVelocity.md), whose gradients ignore the porous neighbours.

`examples/ntp_tutorial/common/channel_turbulence.i` uses it.

!syntax parameters /LinearFVKernels/NTPLinearFVViscousStressCorrection

!syntax inputs /LinearFVKernels/NTPLinearFVViscousStressCorrection

!syntax children /LinearFVKernels/NTPLinearFVViscousStressCorrection
