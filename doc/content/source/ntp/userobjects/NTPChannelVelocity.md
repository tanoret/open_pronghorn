# NTPChannelVelocity

`NTPChannelVelocity` registers, for each component of `velocity`, a functor named
`<component><suffix>` (`vel_x_channel` by default): the velocity with a cell gradient that ignores
the porous neighbours of the free-channel blocks `channel_blocks`. The gradient is the Green-Gauss
reconstruction

\begin{equation}
\nabla \phi_C = \frac{1}{V_C} \sum_f (\phi_f - \phi_C) A_f \mathbf{n}_f ,
\end{equation}

exact for a constant field in Cartesian and RZ coordinates, with $\phi_f$ the linear interpolation
on a face between two channel cells, the functor's face value on a boundary, and $\phi_C$ on a face
towards a cell outside the channels: no normal gradient through a porous interface. Outside the
channel blocks the functor's own gradient is returned.

Between a free channel and a packed bed the superficial velocity drops by orders of magnitude over
one cell (20 m/s against 0.04 m/s at the chamber of the fuel element). The cell gradient of the
velocity reads that jump as a shear of thousands of 1/s, 50 to 160 times the strain of the channel
flow, and a $k$-$\epsilon$ model fed with it produces turbulence without bound. The channel closures
that need the strain rate, the $k$-$\epsilon$ production and turbulent viscosity and
[NTPLinearFVViscousStressCorrection](NTPLinearFVViscousStressCorrection.md), take these functors.

`examples/ntp_tutorial/common/channel_turbulence.i` uses it.

!syntax parameters /UserObjects/NTPChannelVelocity

!syntax inputs /UserObjects/NTPChannelVelocity

!syntax children /UserObjects/NTPChannelVelocity
