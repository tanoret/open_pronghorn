# NTPOrificeDesignAux

`NTPOrificeDesignAux` designs the axial orifice profile of a particle-bed element while the element
is solved by the segregated solvers. In those solvers, the orifice cannot be an unknown of the linear
systems. The cold-frit drag multiplier $\chi$ has one value per axial slice of `mesh_division`. It is
updated towards a uniform slice gas outlet temperature after every solve: every time step with
PIMPLE, every fixed-point iteration with SIMPLE. The update uses the temperature rise
$\Delta T_j = T_j - T_\mathrm{in}$ of each slice over the inlet temperature (`inlet_temperature`), and
the mean rise $\overline{\Delta T}$:

\begin{equation}
\ln \chi_j \leftarrow \ln \chi_j + \mathrm{clamp}\!\left( -\beta \ln \frac{\Delta T_j}{\overline{\Delta T}},\,
-\ln F,\, \ln F \right),
\end{equation}

after which all $\ln \chi_j$ are shifted so that their mean is `level`, i.e. the geometric mean of
$\chi$ is $\exp(\text{level})$. A slice hotter than the mean carries too little flow and is opened
(smaller $\chi$). The level gives every slice room to be opened as well as throttled. This is the
same reference as [NTPSliceUniformityConstraint](NTPSliceUniformityConstraint.md) with
`reference = mean`. $\beta$ is `exponent` (default 0.5). $F$ is `max_change_factor` (default 1.2), which
limits each change. $\Delta T_j$ is floored at $10^{-3}\, \overline{\Delta T}$.

The slice temperatures come from the vector `slice_temperature_vector` of the vector postprocessor
`slice_temperatures`, with one entry per slice. The profile is left unchanged while
$\overline{\Delta T}$ is below `minimum_mean_rise`, or while the vector does not have one entry per slice.
Updates happen only on `TIMESTEP_END`, at most once per time step and fixed-point iteration. Initially
$\chi = \exp(\text{level})$ in every slice. The auxiliary variable, which must be elemental, receives
$\chi$ of the cell's slice and can be passed as `orifice_multiplier` to
[NTPPackedBedSegregatedFunctorMaterial](NTPPackedBedSegregatedFunctorMaterial.md).
`examples/ntp_tutorial/04_navier_stokes/orifice_lagged.i` uses it, and is added to
`anchor_pimple.i` on the command line.

!syntax parameters /AuxKernels/NTPOrificeDesignAux

!syntax inputs /AuxKernels/NTPOrificeDesignAux

!syntax children /AuxKernels/NTPOrificeDesignAux
