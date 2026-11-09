# NTPLinearFVEnthalpyExchange

`NTPLinearFVEnthalpyExchange` adds the volumetric heat exchange $\alpha (T_s - T)$ between a porous
solid and a gas to a gas energy equation solved for the specific enthalpy $h$ (linear finite
volume, segregated solvers). The inputs are `alpha` [W/(m$^3$ K)], the solid temperature `T_solid`,
the gas temperature $T(h)$ (`T_fluid`) and the gas specific heat `cp`. The exchange is linearized
about the current state $(h_0, T_0, c_{p,0})$ with $dT/dh = 1/c_p$:

\begin{equation}
\alpha (T_s - T) \approx \alpha (T_s - T_0) - \frac{\alpha}{c_{p,0}} (h - h_0).
\end{equation}

The kernel adds $\alpha / c_{p,0}$ to the matrix diagonal of $h$ and $\alpha (T_s - T_0 + h_0 / c_{p,0})$ to the
right-hand side, both times the cell volume. At convergence ($h = h_0$), the source is exactly
$\alpha (T_s - T)$. [LinearFVVolumetricHeatTransfer](LinearFVVolumetricHeatTransfer.md) on the
enthalpy variable would instead give $\alpha (T_s - h)$, with inconsistent units.

+Partial elimination+ (optional, enabled by `solid_conductivity`). The segregated solvers solve the
gas and the particle energy one after the other. Where the exchange $\alpha V$ dominates the
particle balance, this Gauss-Seidel coupling contracts by only about $1 - G/(\alpha V)$ per
iteration, with $G$ the conduction conductance of the particles of the cell (in a packed bed $G$ is
about 1 % of $\alpha V$). With the neighbors held
fixed, the particles follow the gas, $\delta T_s = \beta\, \delta T$ with
$\beta = \alpha V / (\alpha V + G)$. Accounting for that response in the gas equation,

\begin{equation}
\alpha (T_s - T) \approx \alpha (T_s - T_0) - \alpha (1 - \beta)(T - T_0),
\end{equation}

makes the coupled pair converge at the rate $G/(\alpha V)$ instead. The diagonal becomes
$(1 - \beta)\, \alpha / c_{p,0}$, and the right-hand side
$\alpha (T_s - T_0 + (1 - \beta)\, h_0 / c_{p,0})$. $G$ is the two-point conductance of the
particle conductivity vector `solid_conductivity` to the neighboring cells in the kernel's blocks.
$1 - \beta$ is floored by `elimination_floor` (default 0.05), so a cell with little flow keeps a gas
diagonal. Whatever $\beta$, the source is still exactly $\alpha (T_s - T)$ at convergence: $G$
only sets the convergence rate.

!syntax parameters /LinearFVKernels/NTPLinearFVEnthalpyExchange

!syntax inputs /LinearFVKernels/NTPLinearFVEnthalpyExchange

!syntax children /LinearFVKernels/NTPLinearFVEnthalpyExchange
