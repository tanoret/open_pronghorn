# NTPRhieChowSliceFlow

`NTPRhieChowSliceFlow` computes the mass flow rate through the interface between two groups of
blocks, one value per division of a mesh division (`mesh_division`, e.g. the axial slices of a
particle-bed element). It sums the face mass fluxes stored by the Rhie-Chow object of the
segregated solvers (`rhie_chow_user_object`, a [RhieChowMassFlux](RhieChowMassFlux.md) such as
[NTPRhieChowMassFlux](NTPRhieChowMassFlux.md)) over the faces between `from_block` and `to_block`:

\begin{equation}
\dot{m}_j = \sum_{f \in \Gamma_j} \dot{m}_f\, A_f ,
\end{equation}

with $\dot{m}_f$ oriented from `from_block` into `to_block` and $\Gamma_j$ the interface faces whose
`to_block` cell lies in division $j$. Because these are the fluxes the discrete mass and energy
balances see, the flows are exact and add up to the flow through the whole interface. Between the
plenum and the cold frit of the fuel element they are the slice flows that
[NTPOrificeFlowDesignAux](NTPOrificeFlowDesignAux.md) designs the orifice on.

The vector `flow` holds the flows, the vector `division` the division indices.
`examples/ntp_tutorial/04_navier_stokes/orifice_flow.i` uses it.

!syntax parameters /VectorPostprocessors/NTPRhieChowSliceFlow

!syntax inputs /VectorPostprocessors/NTPRhieChowSliceFlow

!syntax children /VectorPostprocessors/NTPRhieChowSliceFlow
