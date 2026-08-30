# CorrosionTotalRateAux

`CorrosionTotalRateAux` is the multi-species counterpart of
[CorrosionRateAux](CorrosionRateAux.md). It evaluates the Butler-Volmer current for every configured
species, converts each current to an equivalent penetration rate using that species' valence and
molar mass, and sums those rates. In `mode = recession`, it integrates the summed rate in time.

The finite-element [CorrosionPlating](../actions/CorrosionPlatingAction.md) action uses this object
for `corrosion_rate_um_y` and `recession_um` when
`kinetics_model = mstdb_tc_standard_state`. The resulting diagnostic is the total Cr+Fe+Ni front
rate represented by the MSTDB-TC source fractions, rather than only the recession-controlling
element's share. The reduced-empirical Action path continues to use the legacy single-element
`CorrosionRateAux` behavior.

Supply either one coupled `concentrations` variable per species or one shared external
`concentration_functor`. The per-species electrochemical vectors must all have the same length.
Each phase potential independently accepts a coupled variable, a functor, or a fixed value; sources
for the same phase are mutually exclusive.

!syntax parameters /AuxKernels/CorrosionTotalRateAux

!syntax inputs /AuxKernels/CorrosionTotalRateAux

!syntax children /AuxKernels/CorrosionTotalRateAux
