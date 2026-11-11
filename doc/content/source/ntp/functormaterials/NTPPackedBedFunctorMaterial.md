# NTPPackedBedFunctorMaterial

`NTPPackedBedFunctorMaterial` provides the closures of one packed-bed region (a frit or the fuel
bed) of a particle-bed fuel element, for a Darcy-Forchheimer description of the flow
$\mathbf{u} = -M \nabla p$, with $\mathbf{u}$ the superficial velocity. It is used with the Newton
(AD) path, i.e. [NTPDarcyMassFlux](NTPDarcyMassFlux.md) and [NTPDarcyAdvection](NTPDarcyAdvection.md).

The pressure drop is written as $-\nabla p = (A + B |\mathbf{u}|)\, \mathbf{u}$, with
$Re_m = \rho |\mathbf{u}| d / (\mu (1 - \varepsilon))$ and, for `pressure_drop_model = kta` (KTA 3102.3)
or `ergun`,

\begin{equation}
A_\mathrm{KTA} = \frac{160 \mu (1-\varepsilon)^2}{\varepsilon^3 d^2}, \quad
B_\mathrm{KTA} = \frac{3 \rho (1-\varepsilon)}{\varepsilon^3 d} Re_m^{-0.1}, \qquad
A_\mathrm{Ergun} = \frac{150 \mu (1-\varepsilon)^2}{\varepsilon^3 d^2}, \quad
B_\mathrm{Ergun} = \frac{1.75 \rho (1-\varepsilon)}{\varepsilon^3 d}.
\end{equation}

Both coefficients are multiplied by `pressure_drop_multiplier` and by the functor
`orifice_multiplier`, which can vary along the element to model orifice shaping. An orifice that
is solved together with the flow is applied by the flux kernels instead (`orifice_log_multiplier`
of [NTPDarcyMassFlux](NTPDarcyMassFlux.md)). Since $|\mathbf{u}| = M |\nabla p|$, the mobility
follows from $B |\nabla p| M^2 + A M - 1 = 0$:

\begin{equation}
M = \frac{2}{A + \sqrt{A^2 + 4 B |\nabla p|}}.
\end{equation}

This expression is explicit for Ergun. For KTA, the weak $Re_m^{-0.1}$ dependence of $B$ is resolved
by twelve fixed-point iterations started from the Ergun value. A small offset in $|\nabla p|$ keeps
the derivative finite where the gradient vanishes. The object outputs the mass mobility vector
$\rho M$ (`rho_mobility_name`, isotropic), the velocity $-M \nabla p$ and the superficial speed.

The interphase coefficient is $\alpha = h a_s$ [W/(m$^3$ K)], with $a_s = \chi_a\, 6 (1-\varepsilon)/d$ and
$h = \chi_h Nu\, k_f / d$. The Nusselt number comes from Achenbach (1995) or KTA 3102.2 with
$Re = \rho |\mathbf{u}| d / \mu$ (both Reynolds numbers are floored by `minimum_reynolds`). It is recommended to take $|\mathbf{u}|$ from the cell mass flux reconstructed from the face fluxes
(`mass_flux_x`/`mass_flux_y`/`mass_flux_z`, e.g. from
[NTPDarcyMassFluxFunctorMaterial](NTPDarcyMassFluxFunctorMaterial.md)). Without these inputs the
speed comes from the cell pressure gradient, which can be wrong by a large factor in a porous layer
one cell thick. The particle conductivity vector is
$[(1-\varepsilon) k_s + \varepsilon k_f]\,(\chi_r, \chi_z, \chi_r)$, with $x$ the radius and $y$ the
element axis.

With `particle_emissivity` $e > 0$ the radiative conductivity of the bed (Zehner, Bauer and
Schlünder, simplified as in VDI Heat Atlas D6.3) is added to the particle conductivity,

\begin{equation}
k_\mathrm{rad} = \frac{4 \sigma e}{2 - e}\, d\, T_s^3,
\end{equation}

with $T_s$ the particle temperature (`solid_temperature`). For 0.5 mm particles at 3000 K and
$e = 0.8$ it is about 2 W/(m K). With `dispersion_coefficient` $C > 0$ the object also outputs the
thermal dispersion conductivity of the gas, `dispersion_name`,

\begin{equation}
k_\mathrm{disp} = C\, Pe\, k_f = C\, |\mathbf{G}|\, c_p\, d,
\end{equation}

(Yagi and Kunii; Bauer), $C \approx 0.1$ for the transverse and $0.5$ for the longitudinal
dispersion; the gas conduction of the inputs adds it to $\varepsilon k_f$. Both are off by default.

`examples/ntp_tutorial/01_closures/packed_bed.i` uses this object.

!syntax parameters /FunctorMaterials/NTPPackedBedFunctorMaterial

!syntax inputs /FunctorMaterials/NTPPackedBedFunctorMaterial

!syntax children /FunctorMaterials/NTPPackedBedFunctorMaterial
