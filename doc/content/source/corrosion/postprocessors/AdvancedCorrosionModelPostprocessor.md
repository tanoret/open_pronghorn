# AdvancedCorrosionModelPostprocessor

`AdvancedCorrosionModelPostprocessor` exposes one scalar selected by `quantity` from either an
[MSTDBTCCorrosionUserObject](../userobjects/MSTDBTCCorrosionUserObject.md) or a
[DRIDNCorrosionUserObject](../userobjects/DRIDNCorrosionUserObject.md). Separate instances can report
multiple quantities without duplicating transient state or advancing the model again.

The MSTDB-TC endpoint provides total front/corrosion rates and depths, mass loss/gain, IGC depth,
diffusivity, redox and Fe(II) diagnostics, and per-element affinity, source fraction, saturation,
cold capture, deposit fraction, and dissolved inventory. DRIDN additionally provides surface
availability, cumulative source, coupon deposit, bulk capture, elapsed time, adaptive-step counts,
and its salt-side mass-balance error.

!syntax parameters /Postprocessors/AdvancedCorrosionModelPostprocessor

!syntax inputs /Postprocessors/AdvancedCorrosionModelPostprocessor

!syntax children /Postprocessors/AdvancedCorrosionModelPostprocessor
