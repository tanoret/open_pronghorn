# Radiolytic fluorine off-gassing over an MSR power cycle

This example exercises the full gas-handling capability over a semi-realistic reactor power
transient: **power ramp-up, steady full-power operation, and ramp-down/shutdown**.

## Scenario

A section of fluoride fuel salt (FLiBe-relevant, ~600 C) sits below a swept cover gas with a free
surface at the top. The reactor power history is prescribed as a volumetric dose rate
(`power_history`): it ramps from zero to full power over 0-100 s, holds full power to 400 s, ramps
back to zero over 400-500 s, and then stays shut down. The radiolytic F2 production follows this
power.

During operation the dissolved F2:

- is generated volumetrically by radiolysis (the net F2 G-value, scaled by the power history),
- transfers to the gas phase (bubbles) by interfacial mass transfer ([MSRGasExchange](../../doc/content/source/msr/linearfvkernels/MSRGasExchange.md)),
- degasses directly at the free surface to the cover gas ([MSRDegassingBC](../../doc/content/source/msr/linearfvbcs/MSRDegassingBC.md)),

while the evolved gas phase rises buoyantly, disperses, and vents through the top free surface.

## What it demonstrates

- a **time-dependent radiolytic source** driven by a reactor power history (`dose_rate` functor),
- **gas-phase tracking with buoyancy and dispersion** (`gas_buoyancy_velocity`, `gas_dispersivity`),
- **two off-gas pathways**: volumetric bubble exchange + venting, and direct surface degassing
  (`MSRDegassingBC`),
- the dissolved and gas-phase F2 inventories **building during operation and clearing after
  shutdown**, with the off-gas tracking the power.

## Output

The CSV reports the reactor power fraction and the dissolved / gas-phase / total F2 inventories;
the Exodus output shows the `F2_diss` and `F2_gas` fields building up while the power is on and
clearing once it is off. With the default settings the dissolved F2 inventory peaks at the end of
the steady-operation phase and decays by roughly two orders of magnitude within ~200 s of shutdown
as it vents and degasses.

## Notes

- The radical yields are set to zero so the model follows the net F2 yield (G = 0.005, the Davis
  2022 FLiBe value); this isolates the off-gas dynamics and keeps the chemistry inexpensive. To
  include the explicit radical chemistry, remove the `g_value_*` overrides (and expect a stiffer,
  more iteration-hungry solve).
- A stagnant pool is used (transport by diffusion, gas buoyancy and venting). For coupling to a
  solved flow field, see `examples/molten_salt_corrosion`.
