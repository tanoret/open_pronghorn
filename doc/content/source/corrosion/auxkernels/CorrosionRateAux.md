# CorrosionRateAux

Diagnostic Butler-Volmer output for an electrode boundary: the interfacial current density
[A/m^2] (`mode = current`), the equivalent penetration rate [um/y] (`mode = penetration_rate`), or
the time-integrated metal recession / plating thickness [um] (`mode = recession`). The same kinetics
as the Butler-Volmer objects are evaluated (non-AD), so the value tracks the corrosion or plating the
interface is driving. The [CorrosionPlating](../actions/CorrosionPlatingAction.md) action creates the
`corrosion_rate_um_y` and `recession_um` fields with this auxiliary kernel.

The salt-ion concentration must be supplied through exactly one of the coupled `concentration`
variable and `concentration_functor`. Each phase potential may independently come from one coupled
variable, one functor (`salt_potential_functor` or `metal_potential_functor`), or its fixed `*_value`
fallback; supplying multiple sources for one phase is an input error. Existing inputs that couple
`concentration` and use coupled or fixed phase potentials retain their original behavior. The
functor forms let a solid-only diagnostic evaluate the same external-salt concentration and
potential used by its [ButlerVolmerBC](../bcs/ButlerVolmerBC.md), including spatially or temporally
varying inputs.

For a source-resolved model whose front rate is the sum of several species' equivalent penetration
rates, use [CorrosionTotalRateAux](CorrosionTotalRateAux.md). `CorrosionRateAux` intentionally keeps
its original single-species semantics.

!syntax parameters /AuxKernels/CorrosionRateAux

!syntax inputs /AuxKernels/CorrosionRateAux

!syntax children /AuxKernels/CorrosionRateAux
