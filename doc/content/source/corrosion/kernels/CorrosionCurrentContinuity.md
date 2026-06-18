# CorrosionCurrentContinuity

Current-continuity (charge-conservation) equation for the salt electric potential under the
dilute-solution / MacInnes closure. The ionic current density is

!equation
\vec{i} = -F\sum_k z_k D_k \nabla c_k - \kappa\nabla\phi, \qquad \kappa = \frac{F^2}{RT}\sum_k z_k^2 D_k c_k,

and $\nabla\cdot\vec{i} = 0$ gives the weak (grad-test) residual

!equation
(\nabla\psi,\; \kappa\nabla\phi) + \left(\nabla\psi,\; F\sum_k z_k D_k \nabla c_k\right).

The conductivity $\kappa$ is built from the coupled concentrations, so the potential responds to the
evolving ion field. Electroneutrality is handled separately; the otherwise-undetermined constant mode
of $\phi$ must be pinned with a reference (e.g. an `ADDirichletBC`).

!syntax parameters /Kernels/CorrosionCurrentContinuity

!syntax inputs /Kernels/CorrosionCurrentContinuity

!syntax children /Kernels/CorrosionCurrentContinuity
