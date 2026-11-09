# NTPMixedTemperature

`NTPMixedTemperature` computes the mixed-mean (bulk) temperature of a stream from its enthalpy
flow rate $H$ (`enthalpy_flow`, [W]) and its mass flow rate $\dot{m}$ (`mass_flow`, [kg/s]), both
postprocessors, by inverting the fluid properties `fp`:

\begin{equation}
T_\mathrm{mix} = T\!\left(p,\, h = \frac{H}{\dot{m}}\right).
\end{equation}

The pressure $p$ (`pressure`, a postprocessor or a number, default $10^6$ Pa) only matters for a
fluid whose enthalpy depends on it. The postprocessor returns 0 when $\dot{m} = 0$. $H$ and
$\dot{m}$ typically come from [NTPRhieChowFluxIntegral](NTPRhieChowFluxIntegral.md) or
[NTPDarcyFluxIntegral](NTPDarcyFluxIntegral.md) with the enthalpy and with 1 as carried quantity.

!syntax parameters /Postprocessors/NTPMixedTemperature

!syntax inputs /Postprocessors/NTPMixedTemperature

!syntax children /Postprocessors/NTPMixedTemperature
