# LinearFVTurbulentWallCHTBC

!syntax description /LinearFVBCs/LinearFVTurbulentWallCHTBC

`LinearFVTurbulentWallCHTBC` couples a turbulent fluid resolved with wall functions to a
conducting solid in the segregated linear finite-volume conjugate heat transfer (CHT) solve
(`cht_interfaces` of the [SIMPLE](SIMPLE.md) / [PIMPLE](PIMPLE.md) executioners). The same object
is placed on both sides of the interface, selected with [!param](/LinearFVBCs/LinearFVTurbulentWallCHTBC/side).

## Wall heat-transfer coefficient

The wall-function coefficient is the thermal law of the wall of Jayatilleke (1969):

\begin{equation}
T^+ = Pr_t \left( \frac{1}{\kappa} \ln(E y^+) + P\!\left(\frac{Pr}{Pr_t}\right) \right),
\qquad
P = 9.24 \left( \left(\frac{Pr}{Pr_t}\right)^{3/4} - 1 \right)
          \left( 1 + 0.28\, e^{-0.007\, Pr/Pr_t} \right),
\end{equation}

\begin{equation}
h_{wf} = \frac{\rho c_p u_\tau}{T^+} \quad (y^+ \ge 30), \qquad
h_{wf} = \frac{k}{y} \quad (y^+ \le 5),
\end{equation}

blended linearly in between, with $u_\tau$ from the same wall treatment as the momentum wall
function ([!param](/LinearFVBCs/LinearFVTurbulentWallCHTBC/wall_treatment)). The $P$-function
carries the thermal sublayer resistance, which grows strongly with the Prandtl number: for a
molten fluoride salt ($Pr \approx 15$-$20$) it is the dominant resistance, and a model that only
conducts through the near-wall cell with $k + \mu_t c_p / Pr_t$ overpredicts the wall heat transfer
by more than an order of magnitude (46-fold for the GFMR inner channel at $y^+ = 40$).

## Coupling

Between the near-wall fluid cell P and the near-wall solid cell S there are two resistances in
series, the wall function $1/h_{wf}$ and the half solid cell $d_S/k_S$:

- `side = fluid`: outward flux $h_{wf}(T_P - T_w)$, with the solid surface temperature $T_w$ from
  the CHT handler (`interface_temperature_solid_<boundary>`); the condition hands $T_P$ to the solid.
- `side = solid`: outward flux $U(T_S - T_P)$, $U = (1/h_{wf} + d_S/k_S)^{-1}$, with $T_P$ from the
  fluid (`interface_temperature_fluid_<boundary>`); the condition hands back the surface
  temperature $T_w = T_S + U d_S/k_S\,(T_P - T_S)$.

Each side is implicit in its own temperature and lags the other, so neither sub-problem loses its
temperature reference (a solid that is cooled only by the fluid stays well posed), no stabilizing
coefficient enters the answer, and both fluxes agree at convergence.

The near-wall cell must lie in the logarithmic layer ($y^+ \gtrsim 30$), as for any wall function.

## Verification

`test/tests/kEpsilon/salt_pipe_heat_transfer` runs the GFMR inner channel (FUNaK, 4 cm, 1.7 m/s,
$Re = 1.85 \times 10^4$, $Pr = 16.8$) inside a uniformly heated graphite sleeve. With the first cell
at $y^+ = 33$ the energy balance closes to $3.5 \times 10^{-5}$ and $Nu = 180$, against $194$ from the
Gnielinski correlation (ratio 0.93, inside its $\pm 10$-$20\,\%$ uncertainty); `nusselt.py`
evaluates the comparison.

## Example input syntax

!listing test/tests/kEpsilon/salt_pipe_heat_transfer/salt_pipe.i block=LinearFVBCs/fluid_wall_T

!syntax parameters /LinearFVBCs/LinearFVTurbulentWallCHTBC

!syntax inputs /LinearFVBCs/LinearFVTurbulentWallCHTBC

!syntax children /LinearFVBCs/LinearFVTurbulentWallCHTBC
