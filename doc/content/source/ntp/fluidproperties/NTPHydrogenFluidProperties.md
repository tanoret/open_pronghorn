# NTPHydrogenFluidProperties

`NTPHydrogenFluidProperties` models hot hydrogen for nuclear thermal propulsion as a thermally
ideal, calorically imperfect gas. The specific heat is the NASA 7-coefficient polynomial, the
viscosity follows Sutherland's law and the conductivity follows the Eucken relation:

\begin{equation}
\rho = \frac{p}{R T}, \qquad R = \frac{R_u}{M}, \qquad
\frac{c_p}{R} = a_1 + a_2 T + a_3 T^2 + a_4 T^3 + a_5 T^4,
\end{equation}

\begin{equation}
\frac{h}{R} = a_1 T + \frac{a_2}{2} T^2 + \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4 + \frac{a_5}{5} T^5 + a_6,
\qquad e = h - R T,
\end{equation}

\begin{equation}
\frac{s}{R} = a_1 \ln T + a_2 T + \frac{a_3}{2} T^2 + \frac{a_4}{3} T^3 + \frac{a_5}{4} T^4 + a_7
- \ln \frac{p}{10^5\ \mathrm{Pa}},
\end{equation}

\begin{equation}
\mu = \mu_0 \left( \frac{T}{T_0} \right)^{3/2} \frac{T_0 + S}{T + S}, \qquad
k = \mu \left( c_p + f R \right).
\end{equation}

The defaults are $M = 2.01588 \times 10^{-3}$ kg/mol, $R_u = 8.314462618$ J/(mol K), the
GRI-Mech 3.0 coefficients for H$_2$ (`coefficients_low` and `coefficients_high`, seven each, with
the branch switch at 1000 K), $\mu_0 = 8.76 \times 10^{-6}$ Pa s, $T_0 = 300$ K, $S = 72$ K and the
diatomic Eucken factor $f = 5/4$. The ratio of specific heats is $\gamma = c_p / (c_p - R)$, the speed of
sound $c = \sqrt{\gamma R T}$, and the thermal expansion coefficient $1/T$.

The object provides the $(p, T)$, $(v, e)$ and $(p, h)$ interfaces with analytic derivatives.
Only the density, the specific volume and the entropy depend on pressure; every other property
depends on temperature only. $T(h)$ and $T(e)$
come from a Newton iteration on the monotone polynomial. If the iteration fails to converge, the
object throws a MOOSE exception. Outside `[T_min, T_max]` (20 to 3500 K by default, the range of the fit) the
specific heat, viscosity and conductivity stay at their end-point values, their temperature
derivatives are zero, and $h$ continues linearly with the end-point $c_p$. This keeps $h$ monotone
and smooth.

The model does not include dissociation, ortho/para conversion or real-gas effects. In the AD
overloads that return partial derivatives, the property value carries the full AD derivative. The
returned partial derivatives are treated as constants, except those of the density.

The NTP tutorial (`examples/ntp_tutorial`, stage `01_closures` and `common/hydrogen_nasa.i`) uses
this object and compares it against the same equations in Python.

!syntax parameters /FluidProperties/NTPHydrogenFluidProperties

!syntax inputs /FluidProperties/NTPHydrogenFluidProperties

!syntax children /FluidProperties/NTPHydrogenFluidProperties
