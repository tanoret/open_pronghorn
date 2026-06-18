# Complete flowing-MSR example: flow + power + radiolysis + corrosion

`flowing_msr_corrosion.i` couples, in a single linear finite-volume segregated (SIMPLE) solve, all of
the physics of an operating molten salt reactor fuel-salt channel:

| Physics | Model |
| :------ | :---- |
| Momentum + mass | incompressible Navier-Stokes with Rhie-Chow (`LinearWCNSFVMomentumFlux`, pressure projection) |
| Turbulence | k-epsilon with wall functions (turbulent viscosity `mu_t`) |
| Nuclear power + heat | cosine neutron-flux-shaped volumetric heat source, advected enthalpy, molecular + turbulent conduction |
| Radiolysis | F2 produced with the same neutron-flux shape, advected and diffused (molecular + turbulent), with the evolved gas rising buoyantly (`[MoltenSaltRadiolysis]`) |
| Corrosion | chromium, iron and nickel shed from the alloy walls by a temperature-dependent Butler-Volmer reaction (`[CorrosionPlatingFlow]`), transported up the core by the same turbulent flow |

## The coupling

A vertical fuel-salt channel carries salt upward through the active core between two Hastelloy N
walls. The neutron flux (cosine shape, peaked at the mid-plane) drives both the nuclear heat and the
radiolytic source. The energy equation produces the temperature field, which feeds back into:

- the chemistry (the radiolysis gas-exchange equilibrium), and
- the **Arrhenius corrosion kinetics** — the Butler-Volmer exchange current is anchored to the
  calibrated reference rate (7.4 um/y for this oxidizing Hastelloy N case) and scaled by
  `exp[(Ea/R)(1/T_ref - 1/T)]`, so the wall corrodes fastest where the salt is hottest (peaking near
  7 um/y at the core mid-plane; see the `corrosion_rate_um_y` field and `corrosion_rate_max`).

The corrosion products (`c_Cr`, `c_Fe`, `c_Ni`) are released at the walls by
[CorrosionLinearFVButlerVolmerBC](../../doc/content/source/corrosion/linearfvbcs/CorrosionLinearFVButlerVolmerBC.md)
and transported as passive scalars by the solved turbulent flow and the same effective (molecular +
turbulent) diffusivity as the radiolysis species, giving the near-wall corrosion-product plumes that
rise through the core. Because the radiolytically oxidizing environment is what drives the corrosion,
the corrosion uses the oxidizing redox class of the calibrated correlation.

The dissolved concentrations are small (sub-ppm) because this is a single pass of fast-flowing salt;
in a recirculating loop the inventory accumulates over many passes to the tens-to-hundreds of ppm seen
in the closed-cell coupon and loop models (`examples/corrosion_plating`, `validation/corrosion`).

```
open_pronghorn-opt -i flowing_msr_corrosion.i
```

The Exodus output shows the velocity, temperature, turbulent viscosity, the F2 radiolysis fields, the
chromium/iron/nickel corrosion plumes, and the temperature-dependent wall corrosion-rate field.

## Architecture note

Corrosion is solved here on the **linear finite-volume** basis so that the corrosion products are
passive scalars in the same segregated flow solve as the radiolysis and energy — the natural way to
couple corrosion into a flowing reactor. The complementary **finite-element** corrosion framework
(`[CorrosionPlating]`) adds the full electromigration and electric-potential (electrophoresis) solve
and two-block (salt + solid) interface kinetics, for detailed corrosion/plating studies where the
flow is prescribed or one-way coupled. Both share the same calibrated kinetics and database.
