# NTPPackedBedSegregatedFunctorMaterial

`NTPPackedBedSegregatedFunctorMaterial` provides the packed-bed closures of one region of a
particle-bed fuel element for the segregated (SIMPLE/PIMPLE, linear finite-volume) solvers, which
solve for the superficial velocity. All outputs are `Real` functors. The correlations are the ones
of [NTPPackedBedFunctorMaterial](NTPPackedBedFunctorMaterial.md) (KTA 3102.3 or Ergun pressure
drop, Achenbach or KTA 3102.2 Nusselt number). The Newton and the segregated methods therefore
solve the same closure equations.

- `drag_name` (vector) is the Darcy coefficient of
  [LinearFVMomentumFriction](LinearFVMomentumFriction.md), whose force per unit volume is
  $\mu D \mathbf{u}$. From the packed-bed law $-\nabla p = (A + B|\mathbf{u}|)\,\mathbf{u}$,

  \begin{equation}
  D = \frac{\chi (A + B |\mathbf{u}|)}{\mu},
  \end{equation}

  isotropic, with $\chi$ the product of `pressure_drop_multiplier` and the functor
  `orifice_multiplier`. $B$ depends on $|\mathbf{u}|$, which the segregated iterations lag by one
  iteration.
- `alpha_name` is the interphase coefficient $h a_s = \chi_h Nu\, k_f / d \cdot \chi_a 6 (1-\varepsilon)/d$
  [W/(m$^3$ K)].
- `kappa_name` (vector) is the particle conductivity
  $[(1-\varepsilon) k_s + \varepsilon k_f]\,(\chi_r, \chi_z, \chi_r)$, with $y$ the element axis.
- `speed_name` is the superficial speed used in the correlations.

By default the speed is $\sqrt{u^2 + v^2}$ from the velocity functors `u` and `v`. There is no
parameter for a third velocity component, so in 3-D use the mass-flux inputs to include it. With
`mass_flux_x` and `mass_flux_y` (and optionally `mass_flux_z`),
for example from [NTPRhieChowMassFluxAux](NTPRhieChowMassFluxAux.md), the speed is instead
$|\mathbf{G}| / \rho$. This speed is exact in porous layers one cell thick.

`examples/ntp_tutorial/01_closures/packed_bed_simple.i` uses this object.

!syntax parameters /FunctorMaterials/NTPPackedBedSegregatedFunctorMaterial

!syntax inputs /FunctorMaterials/NTPPackedBedSegregatedFunctorMaterial

!syntax children /FunctorMaterials/NTPPackedBedSegregatedFunctorMaterial
