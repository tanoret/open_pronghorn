# NTPOrificeFlowDesignAux

`NTPOrificeFlowDesignAux` designs the axial orifice profile of a particle-bed element while the
element is solved by the segregated solvers, on the slice mass flows rather than on the slice gas
temperatures ([NTPOrificeDesignAux](NTPOrificeDesignAux.md)). The cold-frit drag multiplier $\chi$
has one value per axial slice of `mesh_division`. It is updated after every solve (every time step
with PIMPLE, every fixed-point iteration with SIMPLE) in two nested loops.

**Flow matching** (every update). With $\dot{m}_j$ the mass flow entering the cold frit of slice $j$
(`slice_flows`, e.g. [NTPRhieChowSliceFlow](NTPRhieChowSliceFlow.md)), $w_j$ its target share and
$\dot{m} = \sum_k \dot{m}_k$,

\begin{equation}
\ln \chi_j \leftarrow \ln \chi_j + \mathrm{clamp}\!\left( \alpha \ln \frac{\dot{m}_j}{w_j \dot{m}},\,
-\ln F,\, \ln F \right),
\end{equation}

after which all $\ln \chi_j$ are shifted so that their mean is `level`, i.e. the geometric mean of
$\chi$ is $\exp(\text{level})$, as for [NTPOrificeDesignAux](NTPOrificeDesignAux.md). $\alpha$ is
`flow_exponent` and $F$ is `max_change_factor`. The flow through a slice falls monotonically with
its drag, and the frit is passed by the cold gas: this loop does not go through the density of the
heated gas, which couples a temperature-driven update to the flow it changes.

**Flow shares** (when $\max_j |\ln(\dot{m}_j / w_j \dot{m})| <$ `share_tolerance`). The shares start
as the slice power shares $w_j = Q_j / \sum_k Q_k$ (`slice_powers`): with the same enthalpy rise per
unit mass, the gas leaves every slice at the same temperature. Axial conduction in the bed, the
heat exchange of the frits and the mixing in the chamber shift the slice outlet temperatures
$T_j$ (`slice_temperatures`) from that, and the shares are corrected with the rises
$\Delta T_j = T_j - T_\mathrm{in}$ over their mean $\overline{\Delta T}$:

\begin{equation}
w_j \leftarrow w_j \left( \frac{\Delta T_j}{\overline{\Delta T}} \right)^{\gamma},
\quad \text{then normalized to } \sum_j w_j = 1,
\end{equation}

each factor limited to [1/`max_share_change`, `max_share_change`]. $\gamma$ is `share_exponent`;
$\gamma = 0$ keeps the power shares. The shares are not corrected while $\overline{\Delta T}$ is below
`minimum_mean_rise`.

Updates happen only on `TIMESTEP_END`, at most once per time step and fixed-point iteration, and
stop after `max_updates` updates, which freezes the designed profile. Initially
$\chi = \exp(\text{level})$ in every slice. The auxiliary variable, which must be elemental, receives
$\chi$ of the cell's slice and can be passed as `orifice_multiplier` to
[NTPPackedBedSegregatedFunctorMaterial](NTPPackedBedSegregatedFunctorMaterial.md).
`examples/ntp_tutorial/04_navier_stokes/orifice_flow.i` uses it.

!syntax parameters /AuxKernels/NTPOrificeFlowDesignAux

!syntax inputs /AuxKernels/NTPOrificeFlowDesignAux

!syntax children /AuxKernels/NTPOrificeFlowDesignAux
