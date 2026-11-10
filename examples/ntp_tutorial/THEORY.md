# Theory: thermal-hydraulics of a particle-bed fuel element

This note derives the equations that the tutorial solves and states every closure it uses. It
also describes the two solution methods: a Newton Darcy-Forchheimer model, and a segregated
SIMPLE/PIMPLE Navier-Stokes model. The tutorial stages refer to the numbered equations.

## 1. The element

A particle-bed reactor (PBR) fuel element is a pair of coaxial porous cylinders (frits) holding a
packed bed of coated fuel particles between them (Ludewig et al., *Prog. Nucl. Energy* 30(1),
1996). Hydrogen enters the outer annular **plenum** at the narrow end. It flows radially inward
through the **cold frit**, the **fuel bed** and the **hot frit** into the central **chamber**, and
leaves axially towards the nozzle. The regions are cones, so their radii vary linearly along the
element (`common/ludewig_design1.i`). The model is axisymmetric, in `(r, z)`, with `z = 0` at the
inlet end.

| region | porosity | particle diameter | particle conductivity |
|---|---|---|---|
| plenum, chamber | 1 | — | — |
| cold frit (Be) | 0.32 | 80 µm | 190 W/(m K) |
| fuel bed | 0.40 | 500 µm | 22 W/(m K), ×0.65 axial, ×0.85 radial |
| hot frit (TaC) | 0.28 | 120 µm | 33 W/(m K) |

Design 1 has 37 elements at 1 GW and a 7 MPa chamber pressure. The gas leaves at 3000 K; the
element flow is the first-law flow `mdot = Q / (37 [h(3000 K) − h(30 K)]) = 0.569035 kg/s`.

## 2. Hydrogen

The gas is a thermal ideal gas with a temperature-dependent specific heat (`NTPHydrogenFluidProperties`):

```
rho = p / (R T),                     R = 8.314462618 / 2.01588e-3 J/(kg K)            (2.1)
cp / R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4,   h = ∫ cp dT   (NASA 7-coefficient,
         GRI-Mech 3.0, switch at 1000 K, valid to 3500 K)                                 (2.2)
mu = 8.76e-6 (T / 300)^1.5 (300 + 72) / (T + 72)   (Sutherland)                          (2.3)
k  = mu (cp + 5/4 R)                                (Eucken)                              (2.4)
```

There is no dissociation and no ortho/para conversion. Stage 1 checks the MOOSE implementation
against the same equations in Python (`common/hydrogen.py`).

## 3. Governing equations

The velocity `u` is the **superficial** velocity, the volume flow per unit total area. The
interstitial velocity is `u / eps`. The porous regions have local thermal non-equilibrium: the
gas temperature `T` and the particle temperature `T_s` are separate unknowns.

Mass (steady):

```
div(rho u) = 0                                                                            (3.1)
```

Momentum. In the packed beds the drag (4.1) is many orders of magnitude larger than the inertia
and viscous terms, so

```
-grad p = (A + B |u|) u                       in the frits and the bed                    (3.2)
```

In the free channels (plenum, chamber) the momentum balance keeps inertia, viscous stress and a
wall-friction term:

```
div(rho u u) = -grad p + div(tau) - F_wall    in the plenum and the chamber               (3.3)
```

The two solution methods (Section 6) differ only in how they treat (3.3).

Gas energy, with the specific enthalpy `h(T)` of (2.2):

```
div(rho u h) = div(eps k grad T) + alpha (T_s - T)                                        (3.4)
```

Particle energy (porous regions):

```
0 = div(kappa_s grad T_s) + q''' - alpha (T_s - T)                                        (3.5)
```

`q'''` is the fission power density, nonzero in the fuel bed only.

## 4. Closures

### 4.1 Packed-bed pressure drop (KTA 3102.3)

```
A = 160 mu (1 - eps)^2 / (eps^3 d^2),
B = 3 rho (1 - eps) / (eps^3 d) Re_m^(-0.1),       Re_m = rho |u| d / (mu (1 - eps))      (4.1)
```

The equivalent friction factor is `f Re_m = 160 + 3 Re_m^0.9`. The Ergun law (`A` with 150,
`B = 1.75 rho (1 - eps) / (eps^3 d)`) is available for comparison. Stage 1 validates (4.1)
against air data.

### 4.2 Interphase heat transfer

`alpha = h a_s`, with the specific surface `a_s = 6 (1 - eps) / d`. The Nusselt number
`Nu = h d / k` comes from the Achenbach correlation in the fuel bed and from KTA 3102.2 in the
frits, both at the particle Reynolds number `Re = rho |u| d / mu`.

### 4.3 Effective particle-bed conductivity

```
kappa_s = [(1 - eps) k_s + eps k] (chi_r, chi_z)                                          (4.2)
```

`chi_r` and `chi_z` are the bed anisotropy multipliers of the table in Section 1 (1 in the frits).

### 4.4 Free-channel wall friction

The plenum and chamber are a few millimetres wide and are resolved by one or a few cells across.
That is too few to resolve a turbulent wall layer, so the wall shear is the friction of a duct of
hydraulic diameter `D_h` (twice the gap of the plenum annulus, the diameter of the chamber):

```
F_wall = (32 mu / D_h^2 + f rho |u| / (2 D_h)) u_z e_z,     f = 0.02                      (4.3)
```

### 4.5 Orificing

Uniform frits do not distribute the flow in proportion to the power. The split follows from the
pressure differences that the plenum and chamber flows leave across the porous layers, and those
differ from the axial power shape. The cold frit therefore carries an axial drag multiplier
`chi(z)`, the orifice, constant on each of `n` axial slices:

```
-grad p = chi(z) (A + B |u|) u                in the cold frit                            (4.4)
```

The tutorial designs `chi` together with the flow, so that the gas leaves every slice at the
same temperature. With `<T>_j` the gas temperature leaving slice `j` (averaged over its hot frit)
and `lambda_j = ln chi_j`, the design conditions are

```
<T>_j = (1/n) sum_k <T>_k   (n - 1 independent conditions),    (1/n) sum_k lambda_k = lambda_0   (4.5)
```

The level `lambda_0` (`orifice_level`) gives every slice room to be opened as well as throttled.
With one slice held at `chi = 1` instead, a starved slice may need its drag removed altogether.
There its temperature no longer depends on `chi`, and the design is singular.

- **Newton model.** The `lambda_j` are unknowns of the Newton system, and (4.5) are their
  equations (`NTPSliceUniformityConstraint`). The design equations are switched on gradually as
  the gas heats.
- **Segregated model, on the slice flows** (`NTPOrificeFlowDesignAux`, stage 4). The segregated
  linear systems cannot hold the `lambda_j` as unknowns. They are updated after every solve
  (fixed-point iteration or time step) towards target slice flows instead of (4.5) directly:

  ```
  lambda_j <- lambda_j + alpha ln(m_j / (w_j mdot)),   then shifted to a mean lambda_0        (4.6)
  ```

  `m_j` is the mass flow entering the cold frit of slice `j`, the sum of the Rhie-Chow face fluxes
  through its plenum interface (`NTPRhieChowSliceFlow`). The shares `w_j = Q_j / Q` are the slice
  fission powers: with the same enthalpy rise per unit mass, every slice reaches the same outlet
  temperature, (4.5) up to the axial conduction in the bed. `m_j` falls monotonically with
  `lambda_j`, and the cold frit sees the gas at its inlet density, so (4.6) does not close a loop
  through the density of the heated gas. Optionally, the shares are corrected from the outlet
  temperatures, `w_j <- w_j (rise_j / mean rise)^gamma`, which does.
- **Segregated model, on the slice temperatures** (`NTPOrificeDesignAux`). The `chi_j` are updated
  after every time step towards (4.5), damped: `chi_j <- chi_j (mean rise / rise_j)^beta`. The
  update closes its loop through the gas density and has no stable fixed point (stage 4).

## 5. Boundary conditions

| boundary | flow | gas energy | particles |
|---|---|---|---|
| inlet (plenum annulus, z = 0) | mass flow `mdot` | `T = 30 K` | — |
| outlet (chamber, z = L) | `p = 7 MPa` | outflow | — |
| axis (r = 0) | symmetry | adiabatic | — |
| outer wall, closed ends | no slip / no flux | adiabatic | adiabatic |

## 6. Solution methods

### 6.1 Newton, Darcy-Forchheimer (stage 3)

In every region the momentum balance is reduced to a mobility law `u = -M grad p`. In the porous
regions, `|u| = M |grad p|` turns (3.2) into `B |grad p| M^2 + A M - 1 = 0`, whose positive root
is

```
M = 2 / (A + sqrt(A^2 + 4 B |grad p|))                                                    (6.1)
```

With Ergun, `B` does not depend on `|u|` and (6.1) is explicit. The weak `Re_m^-0.1` dependence of
the KTA `B` is a fixed point. It starts from Ergun and alternates (6.1) with `|u| = M |grad p|`
(12 iterations, converged to round-off).

In the free channels `M` follows from (4.3), with inertia neglected. Mass (3.1) becomes a pressure
equation, `div(rho M grad p) = 0`. The gas energy uses the same face mass fluxes in convective
form, `sum_f mdot_f (h_upwind - h_P)`, so mass and energy see identical fluxes. The face mobility
combines the two half cells **in series** (`NTPSeriesInterpolation`):

```
M_f = 1 / ((1 - g_C) / M_C + g_C / M_N),       g_C = |x_N - x_f| / |x_N - x_C|            (6.2)
```

This is the exact mobility of two uniform half segments. MOOSE's harmonic average uses the
linear-interpolation weights instead, `1 / (g_C / M_C + (1 - g_C) / M_N)`. The two agree on uniform
meshes only. Next to a 2 mm frit cell, the MOOSE weights charge the frit with the resistance of
half the neighbouring 8.75 mm bed cell.

The speed in the interphase heat-transfer coefficient is the cell mass flux reconstructed from the
face fluxes,

```
G_C = (1 / V_C) sum_f mdot_f A_f (x_f - x_C)                                              (6.3)
```

This is exact for a uniform flux. In a one-cell porous layer, the velocity from the cell pressure
gradient can be wrong by a factor of five. (6.3) is evaluated as an AD functor of the solution
(`NTPDarcyMassFluxFunctorMaterial`), so the Newton Jacobian contains its derivatives.

All equations are solved together by Newton's method. The steady state is followed in the power
by **continuation**: a sequence of steady solves, each starting from the previous one, as the
power rises from a few per cent to 100 %.

### 6.2 Segregated SIMPLE / PIMPLE Navier-Stokes (stages 2, 4 and 5)

The momentum equation is solved for the superficial velocity in all five regions. The packed-bed
law (3.2) and the duct friction (4.3) enter as a Darcy friction `mu D u`, `D = (A + B|u|)/mu`.
The plenum and chamber resolve inertia and viscous stress (3.3). Pressure and velocity are coupled
by the Rhie-Chow interpolation of the SIMPLE algorithm (steady) or of PIMPLE (pseudo-transient).
Seven points are specific to the particle-bed element:

1. **Series face coupling** (`NTPRhieChowMassFlux`). From the momentum matrix, each cell has a
   mobility `M = V/A` and a force density `f = H/V`: `u_C = -M_C (grad p + f)_C`. The standard
   Rhie-Chow interpolation averages `rho M` linearly to the faces. Next to a free channel, whose
   `M` is ten orders of magnitude larger, a one-cell frit then loses its resistance. The face flux
   is instead the series form (6.2) with the segment-averaged force density,
   `mdot_f = -(rho M)_f [(grad p)_f + f_f] . n`. The cell pressure gradient of the momentum
   equations (Green-Gauss) is corrected consistently, with the face pressure of flux continuity,
   `p_f = (w_C p_C + w_N p_N) / (w_C + w_N)` and `w = rho M / d` the half-segment conductances,
   instead of the linear average. The linear average gives a free-channel cell next to a frit a
   spurious pressure gradient of the frit's magnitude. On walls and other boundaries without a
   normal pressure gradient, the face coefficient is made isotropic. The tangential part of the
   anisotropic pressure diffusion would otherwise carry mass through walls not aligned with the
   axes, the conical outer wall for instance.
2. **Lagged density.** The segregated pressure equation has no compressibility term. With
   `rho(p, T)` updated inside the pressure–velocity iterations, the factor 1e10 between the free and
   packed momentum diagonals amplifies pressure corrections without bound. The flow properties are
   therefore evaluated at a lagged state `(p_rho, T_rho)`. That state moves towards the solved
   state by 30 % per outer (fixed-point) iteration with SIMPLE, or per time step with PIMPLE. The
   converged state is that of the fully coupled equations.
3. **Lagged friction speed.** The turbulent duct friction (4.3) grows with the axial speed `|u_z|`.
   The radial velocity of a channel cell is the cross flow into the frit, not duct flow. The speed
   is under-relaxed and updated once per outer iteration; updated every iteration, the friction
   diverges.
4. **Enthalpy exchange with partial elimination** (`NTPLinearFVEnthalpyExchange`). The energy
   equation is solved for `h`. The exchange `alpha (T_s - T(h))` is linearised with
   `dT/dh = 1/cp`. The gas and particle energies are solved one after the other. Where the exchange
   `alpha V` dominates the particle conduction conductance `G` of a cell, that sequence contracts
   by only `1 - G/(alpha V)`, about 0.998 per iteration. The gas equation therefore accounts for the
   response of the particles of the cell, `dT_s = beta dT` with `beta = alpha V / (alpha V + G)`:

   ```
   alpha (T_s - T) ~ alpha (T_s - T0) - alpha (1 - beta) (T - T0)
   ```

   This is exact at convergence for any `beta`. The pair then converges at the rate `G/(alpha V)`.
   The particle conduction takes the harmonic mean of the two cells' conductivities at a face, as
   the Newton model does. Next to the fuel bed the cold frit conducts ten times better (an
   effective 129 against 13 W/(m K)); taken on the face directly, that conductivity lets the bed
   heat the cold-frit gas to twice the temperature of the Newton model.
5. **Cell mass flux.** The cell mass flux for the closures (`|u|` in (4.1) and in `Nu`) is
   reconstructed from the face fluxes (6.3), here from those of the Rhie-Chow object
   (`NTPRhieChowMassFluxAux`).
6. **Start of the flow.** With all face fluxes zero, the first momentum matrix of the free
   channels holds only the viscous term, and the first pressure correction overshoots by orders
   of magnitude. SIMPLE starts from a plug-flow estimate of the velocity
   (`segregated_plug_velocity`); PIMPLE may start from rest, the time derivative bounding the
   mobility by `dt / rho`. In RZ both starts converge; in 3-D only the flowing one.
7. **Orifice columns in 3-D.** The flow design (4.6) acts on columns, an axial slice times an
   azimuthal sector (`NTPElementColumnsDivision`). With one multiplier per slice (rings), the
   columns of a slice are parallel heated channels behind one throttle, and the flow excursion
   of the frozen orifice (stage 4) grows between them around the element.

## 7. Figures of merit

- `T_mixed`: mixed-mean outlet gas temperature, from the outlet enthalpy flow,
  `T(h = H_out / mdot)`.
- `T_fuel,max`: largest particle temperature in the fuel bed.
- `p_in`: inlet pressure; `p_in - p_out` is the element pressure drop.
- `chi_max`: largest orifice multiplier.
- Energy balance error: `(H_out - H_in - Q) / Q`, from the exact face fluxes.
