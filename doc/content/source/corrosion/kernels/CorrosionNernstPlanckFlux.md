# CorrosionNernstPlanckFlux

Nernst-Planck diffusion and electromigration flux of one dissolved cation $M^{z+}$ in the molten
salt. The molar flux is $\vec{N} = -D\nabla c - \tfrac{zF}{RT}D\,c\,\nabla\phi$, giving the weak
(grad-test) residual

!equation
(\nabla\psi,\; D\nabla c) + \left(\nabla\psi,\; \frac{zF}{RT}D\,c\,\nabla\phi\right).

The migration term is added only when the `potential` variable is coupled; otherwise the kernel is
pure diffusion. Advection is added separately by [CorrosionAdvection](CorrosionAdvection.md), and the
Butler-Volmer flux enters as the natural boundary term handled by the interface kernel or boundary
condition.

!syntax parameters /Kernels/CorrosionNernstPlanckFlux

!syntax inputs /Kernels/CorrosionNernstPlanckFlux

!syntax children /Kernels/CorrosionNernstPlanckFlux
