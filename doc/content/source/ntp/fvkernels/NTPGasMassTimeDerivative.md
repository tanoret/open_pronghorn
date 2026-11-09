# NTPGasMassTimeDerivative

`NTPGasMassTimeDerivative` adds the mass storage of a gas in a porous medium to the pressure
equation,

\begin{equation}
\frac{\partial (\varepsilon \rho)}{\partial t} = \varepsilon \left( \frac{\partial \rho}{\partial p}
\frac{\partial p}{\partial t} + \frac{\partial \rho}{\partial T} \frac{\partial T}{\partial t} \right),
\qquad \rho = \rho(p, T),
\end{equation}

where the variable is the pressure $p$, $T$ is the gas temperature functor `T_fluid`, $\varepsilon$ the
`porosity` functor, and $\rho(p, T)$ and its derivatives come from the fluid properties `fp`. The porosity
is treated as constant in time. The kernel is used with the Darcy-Forchheimer mass equation of
[NTPDarcyMassFlux](NTPDarcyMassFlux.md).

In steady calculations, the term only damps the pseudo-transient iterations. It does so with the
physical compressibility, so pressure and mass flow relax together.

!syntax parameters /FVKernels/NTPGasMassTimeDerivative

!syntax inputs /FVKernels/NTPGasMassTimeDerivative

!syntax children /FVKernels/NTPGasMassTimeDerivative
