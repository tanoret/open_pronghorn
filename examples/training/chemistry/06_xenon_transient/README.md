# Case 06: I-135 / Xe-135 -- equilibrium poisoning and the post-shutdown xenon peak

**Module 2 (fission products) | ~20 min | 1 cell, serial, 9 s**

## Objective

Model a decay chain with neutron capture and time-dependent power and flux, and reproduce the
classic xenon transient after a scram.

## Physics

Xe-135 has the largest thermal capture cross section of any nuclide (2.65e6 b). It comes from
fission directly and, mostly, from the decay of I-135 (6.6 h). It is removed by decay (9.1 h) and
by capture ("burnout"):

$$\frac{dN_I}{dt} = \frac{Y_I F}{N_A} - (\lambda_I + \sigma_I\phi)N_I, \qquad
  \frac{dN_{Xe}}{dt} = \frac{Y_{Xe} F}{N_A} + \lambda_I N_I - (\lambda_{Xe} + \sigma_{Xe}\phi)N_{Xe}$$

The reactor runs 4 days at full power, so Xe reaches equilibrium, and then scrams: $F = \phi = 0$.
Burnout stops but the iodine keeps feeding the xenon, so **Xe-135 rises for hours**. This is the
"xenon pit" that can stop a thermal reactor from restarting.

Reactivity worth: $\rho_{Xe} \approx -\sigma_{Xe} N_{Xe} N_A \phi_0 / (\nu F_0)$, the fraction of
fission neutrons Xe-135 would absorb at full power.

## The input, step by step

[xenon_transient.i](xenon_transient.i):

- `power_density = power_density` and `neutron_flux = flux` are **functors**: here
  `PiecewiseLinear` functions that drop to zero at the scram.
- The action links I135 -> Xe135 automatically from the database (decay branching 0.9995) and
  orders the systems parent-before-daughter. One Gauss-Seidel sweep per step is then exact.
- `N_Xe135_analytic` holds the equilibrium value before the scram and the exact two-exponential
  solution after it.

## Run

```bash
open_pronghorn-opt -i xenon_transient.i
python plot.py
```

## Expected results

![](../figures/case06_xenon_transient.png)

- Equilibrium Xe-135 matches the analytic value to $10^{-4}$ (the approach before 90 h is the
  build-up from fresh fuel; the analytic column is the equilibrium value only).
- After the scram the model follows the exact solution to $2\times10^{-4}$.
- Peak: **8.3 h after the scram, at 1.89x** the equilibrium inventory. The worth goes from
  -2130 pcm to about -4000 pcm.

## Exercises

1. Run at a higher flux (`phi0=1e18`). How do the equilibrium worth, peak ratio and peak time
   change?
2. Why is there no xenon peak in a fast reactor (case 16)?
3. Estimate the peak time analytically: set $dN_{Xe}/dt = 0$ after the scram.

<details><summary>Solutions</summary>

1. Burnout dominates more, so the equilibrium Xe is lower relative to its production, but the
   iodine inventory is the same. The peak is stronger and later: **4.5x at 10.2 h**, with the
   worth going from -2490 to -11150 pcm. This is why high-flux thermal reactors are the most
   sensitive.
2. Fast spectra have negligible Xe-135 capture (the 2.65e6 b is thermal), so burnout is
   negligible at power. There is no "burnout lost" at shutdown, and the inventory only decays.
3. The peak is where $b\lambda_I N_I = \lambda_{Xe} N_{Xe}$. With
   $A = b\lambda_I N_{I0}/(\lambda_{Xe} - \lambda_I)$ and the post-scram solution,

   $$t_{peak} = \frac{1}{\lambda_{Xe} - \lambda_I}
     \ln\frac{\lambda_{Xe}(N_{Xe0} - A)}{-A\,\lambda_I} = 8.37\ \mathrm{h}$$

   for the equilibrium inventories of the base case. The run gives 8.33 h: the 0.04 h difference
   is the first-order (implicit Euler) error of the post-scram time steps.

</details>
