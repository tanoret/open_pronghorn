# Fission-Product Inventory in Circulating Fuel

In a molten salt reactor the fuel is dissolved in the salt, so fission products are born in a
moving fluid and are subject to chemical separation processes that do not exist in solid fuel. The
[FissionProductInventory](FissionProductInventoryAction.md) action models the inventory $N_i$
[mol/m$^3$] of each tracked nuclide with the transport form of the Bateman equations

\begin{equation}
\frac{\partial N_i}{\partial t} + \nabla\cdot(\mathbf{u} N_i) - \nabla\cdot(D \nabla N_i) =
\frac{Y_i F}{N_A} + \sum_p \left(\lambda_p b_{pi} + \sigma_p \phi\, \delta_{pi}\right) N_p
- \left(\lambda_i + \sigma_i \phi + r_i\right) N_i .
\end{equation}

The fission rate density follows from the power density, $F = q'''/E_f$ with $E_f \approx 200$ MeV.

## Chemical groups and removal

Fission products in fluoride and chloride salts fall into four behavioural groups:

| Group | Examples | Behaviour in the salt | Removal |
| :---- | :------- | :-------------------- | :------ |
| Noble gases | Kr-85, Xe-133, Xe-135 | insoluble, migrate to bubbles and the cover gas | off-gas sparging, $r = \epsilon \dot{V}_\text{He}/V$ |
| Noble metals | Mo-99, Tc-99, Ru-103, Ru-106, Nb-95, Te-132 | reduced metals, not soluble; colloids | plating on metal walls ($k_m A/V$), entrainment into off-gas |
| Halogens | I-131, I-135 | stable fluoride/iodide anions | stay in the salt |
| Salt seekers | Sr-90, Zr-95, Cs-137, Sm-149 | stable fluorides | stay in the salt (chemical processing) |

In the MSRE most of the noble metals were found on the metal surfaces of the loop and in the
off-gas system rather than in the salt (ORNL-4865), which is the motivation for tracking the plated
and off-gas inventories explicitly.

## Xenon poisoning and stripping

For the 135 chain with the precursor tellurium lumped into iodine, the equilibrium inventories are

\begin{equation}
N_I = \frac{Y_I F}{N_A(\lambda_I + \sigma_I\phi)}, \qquad
N_{Xe} = \frac{Y_{Xe} F/N_A + b\, \lambda_I N_I}{\lambda_{Xe} + \sigma_{Xe}\phi + r_\text{strip}} .
\end{equation}

The capture cross section of Xe-135 (2.65$\times 10^6$ b) makes burnout dominant at high flux; a
stripping rate $r_\text{strip}$ comparable to $\lambda_{Xe} + \sigma_{Xe}\phi$ reduces the poison
substantially, and after shutdown ($\phi = 0$) the decay of the iodine inventory produces the
xenon peak.

## Numerics

The equations are discretized with the linear finite-volume method; each inventory is a separate
linear system. The loss terms are implicit and the feeds use the latest parent iterate, so solving
the systems parent before daughter makes each implicit-Euler step exact for a decay chain in a single
sweep. Verification against the analytic saturation, Bateman and xenon-equilibrium solutions is in
`test/tests/msr/fission_products`.
