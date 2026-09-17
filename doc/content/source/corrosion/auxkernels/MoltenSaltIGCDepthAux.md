# MoltenSaltIGCDepthAux

!syntax description /AuxKernels/MoltenSaltIGCDepthAux

This elemental diagnostic applies the reduced-order model's calibrated IGC-depth relation to an
already accumulated uniform corrosion-depth variable. It does not multiply the current corrosion rate
by the total reactor time; the input field is expected to contain the integrated corrosion-rate
history. For example:

```
[AuxKernels]
  [igc_depth]
    type = MoltenSaltIGCDepthAux
    variable = igc_depth_um
    corrosion_depth = corrosion_depth_um
    salt_class = fluoride_fuel
    redox_class = purified_baseline
  []
[]
```

The output and coupled `corrosion_depth` variables must be elemental. Standard variable coupling
allows MOOSE to execute the AuxKernel that produces `corrosion_depth` before this diagnostic at each
execution point. The calibrated morphology equation and damage parameters remain in
`MoltenSaltCorrosionModel`.

!syntax parameters /AuxKernels/MoltenSaltIGCDepthAux

!syntax inputs /AuxKernels/MoltenSaltIGCDepthAux

!syntax children /AuxKernels/MoltenSaltIGCDepthAux
