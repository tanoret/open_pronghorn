# MoltenSaltCorrosionRateAux

!syntax description /AuxKernels/MoltenSaltCorrosionRateAux

This elemental diagnostic evaluates the effective corrosion correlation with the local value of
the `temperature` functor. All other corrosion features remain scalar input parameters. For example:

```
[AuxKernels]
  [corrosion_rate]
    type = MoltenSaltCorrosionRateAux
    variable = corrosion_rate_um_y
    temperature = T_fluid
    boundary = 'shield_wall reflector_wall'
    material_class = hastelloy_n
    salt_class = fluoride_fuel
    redox_class = purified_baseline
    flow_factor = 0.75
    delta_T_C = 0
  []
[]
```

The output variable must be elemental. The named temperature may be any scalar functor, including
an FV variable. On a boundary, an FV temperature is sampled from the adjacent element; this object
does not reconstruct or impose a boundary-face temperature.

!syntax parameters /AuxKernels/MoltenSaltCorrosionRateAux

!syntax inputs /AuxKernels/MoltenSaltCorrosionRateAux

!syntax children /AuxKernels/MoltenSaltCorrosionRateAux
