# CorrosionAdvection

Advection of a dissolved species by a prescribed velocity field, in non-conservative form
$(\psi,\; \vec{u}\cdot\nabla c)$. The velocity components are read as functors, so the field can be
supplied by a separate flow solve (for example the segregated Navier-Stokes velocity `vel_x`/`vel_y`)
as a one-way coupling: the flow drives the species, which do not feed back on the flow. For strictly
mass-conserving transport use the linear finite-volume advection objects of the flow framework; this
kernel is intended for the nonlinear Newton corrosion system where advection is a secondary transport
contribution.

!syntax parameters /Kernels/CorrosionAdvection

!syntax inputs /Kernels/CorrosionAdvection

!syntax children /Kernels/CorrosionAdvection
