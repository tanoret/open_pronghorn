# Stage 1: closures

Before the fuel-element model is assembled, this stage checks the two closures every later stage
rests on: the hydrogen properties and the packed-bed pressure drop. Each check runs MOOSE through
the same objects the element model uses.

All commands run inside the MOOSE environment (`moose-dev-shell`, or `moose-dev-exec`), from this
directory. The scripts find the executable at `../../../open_pronghorn-opt`; pass `--exe` for a
different one.

## 1.1 Hydrogen properties

```bash
python3 check_hydrogen.py
```

`hydrogen_properties.i` sweeps the temperature from 20 to 6000 K at 7 MPa. At each point it
evaluates the density, enthalpy, isobaric specific heat, viscosity and conductivity of
`NTPHydrogenFluidProperties` (included from `../common/hydrogen_nasa.i`). The model is:

- a thermal ideal gas, `rho = p / (R T)` with `R = 8.314462618 / 2.01588e-3` J/(kg K);
- the GRI-Mech 3.0 NASA 7-coefficient fit for `cp(T)` and `h(T)`, switching at 1000 K and valid up
  to 3500 K; above it `cp`, `mu` and `k` are held at their 3500 K values and `h` continues linearly
  (the sweep to 6000 K also checks this continuation);
- Sutherland viscosity and Eucken conductivity, `k = mu (cp + 5/4 R)`.

There is no dissociation and no ortho/para conversion. The script compares MOOSE with the same
equations written in Python (`../common/hydrogen.py`, which later post-processing uses). It also
prints the first-law flow of the 1 GW anchor:

| quantity | value |
|---|---|
| largest relative difference, MOOSE vs Python (all five properties) | about 7e-14 |
| h(3000 K) - h(30 K) | 4.749625e7 J/kg |
| bundle flow 1 GW / Δh | 21.054295 kg/s |
| element flow (37 elements) | 0.569035 kg/s |

## 1.2 Packed-bed pressure drop

```bash
python3 validate.py              # KTA 3102.3 (default); --model ergun for Ergun
```

The fuel bed and both frits are packed beds. Their pressure gradient follows the KTA 3102.3
correlation

```
-dp/dx = (A + B |u|) u,
A = 160 mu (1 - eps)^2 / (eps^3 d^2),
B = 3 rho (1 - eps) / (eps^3 d) Re_m^-0.1,      Re_m = rho |u| d / (mu (1 - eps)),
```

with `u` the superficial velocity. Equivalently, the modified friction factor is
`f Re_m = 160 + 3 Re_m^0.9`.

The element model solves this law in two ways, and both are exercised here on a 1-D bed:

- `packed_bed.i`: the Newton path. The pressure equation is `-div(rho M grad p) = 0`, with the
  mobility `M` solved from the law in closed form (`NTPPackedBedFunctorMaterial`,
  `NTPDarcyMassFlux`).
- `packed_bed_simple.i`: the segregated SIMPLE path. The momentum equation is solved for `u`, with
  the law as a Darcy friction `mu D u`, `D = (A + B|u|)/mu`
  (`NTPPackedBedSegregatedFunctorMaterial`).

`validate.py` runs both inputs at each of the 46 points of the Texas A&M air data set. These are
three particle sizes in a horizontal bed, `D/d` from 6.3 to 19 and `Re_m` from 260 to 7900. The
data are Table A-5 of G. Zou and R. Hu, *SAM Code Validation on Frictional Pressure Drop through
Pebble Beds*, ANL/NSE-20/5 (2020). The script compares `f Re_m` with the measurements and with the
correlation evaluated in closed form.

| check | result |
|---|---|
| Newton path vs closed-form KTA | max relative difference 8e-12 |
| SIMPLE path vs closed-form KTA | max relative difference 2e-13 |
| KTA vs Texas A&M data (46 points) | mean abs. error 5.32 %, max 13.73 % |

Both solver paths implement the correlation exactly. The correlation agrees with the data within
its stated 15 % uncertainty. Its largest errors are on the 12.7 mm spheres at high `Re_m`, where
the bed-to-particle diameter ratio `D/d = 9.5` puts a large fraction of the bed near the wall.

![Texas A&M validation](validation.png)

The script writes `validation_points.csv` (all points, all predictions), `validation_summary.txt`
and `validation.png`.

`data/missouri_air_pressure_drop.csv` (Table B-2 of the same report) holds a second data set: an
air bed at three diameter ratios. Comparing with it needs the bed height, porosity and particle
diameter of those experiments, which are given in the report and not repeated here.

## What is not validated

The interphase heat transfer (Achenbach correlation in the bed, KTA in the frits) and the
effective particle-bed conductivity enter the element model as published correlations. No public
heated-hydrogen packed-bed data are used in this tutorial. Their uncertainty carries over to the
fuel temperatures of stages 3–5; the mixed outlet temperature does not depend on them.
