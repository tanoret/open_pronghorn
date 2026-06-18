# CorrosionRateAux

Diagnostic Butler-Volmer output for an electrode boundary: the interfacial current density
[A/m^2] (`mode = current`), the equivalent penetration rate [um/y] (`mode = penetration_rate`), or
the time-integrated metal recession / plating thickness [um] (`mode = recession`). The same kinetics
as the Butler-Volmer objects are evaluated (non-AD), so the value tracks the corrosion or plating the
interface is driving. The [CorrosionPlating](CorrosionPlatingAction.md) action creates the
`corrosion_rate_um_y` and `recession_um` fields with this auxiliary kernel.

!syntax parameters /AuxKernels/CorrosionRateAux

!syntax inputs /AuxKernels/CorrosionRateAux

!syntax children /AuxKernels/CorrosionRateAux
