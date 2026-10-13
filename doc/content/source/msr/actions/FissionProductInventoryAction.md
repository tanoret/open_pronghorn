# FissionProductInventory

The `[FissionProductInventory]` block assembles a fission-product inventory (depletion) problem in a
circulating fuel salt from the built-in nuclide database. For every nuclide listed in `nuclides` it
creates a `MooseLinearVariableFVReal` molar inventory $N_i$ [mol/m$^3$] in its own linear system
`<nuclide>_sys` and the kernels of the Bateman equations with molten salt reactor removal terms

\begin{equation}
\frac{\partial N_i}{\partial t} + \nabla\cdot(\mathbf{u} N_i) - \nabla\cdot(D \nabla N_i) =
\frac{Y_i F}{N_A} + \sum_p \left(\lambda_p b_{pi} + \sigma_p \phi\, \delta_{pi}\right) N_p
- \left(\lambda_i + \sigma_i \phi + \sum_k r_{k,i}\right) N_i ,
\end{equation}

where $F$ is the fission rate density, $Y_i$ the effective fission yield, $b_{pi}$ the decay branching
of a tracked parent $p$ to $i$, $\delta_{pi}=1$ when $i$ is the capture product of $p$, $\sigma$ the
one-group capture cross section, $\phi$ the one-group neutron flux and $r_{k,i}$ the first-order
removal rates of $i$'s chemical group. See [theory/fissionProductInventory.md] for the derivation and
the MSR removal pathways.

The objects created are:

- `LinearFVTimeDerivative` (unless `time_derivative = false`),
- `LinearFVScalarAdvection` (with `rhie_chow_user_object`) or `LinearFVAdvection` (with a constant
  `velocity`) and `LinearFVDiffusion` (with `diffusivity`) for the salt inventories,
- `LinearFVSource` for the fission source $Y_i F / N_A$,
- [FissionProductFeed](FissionProductFeed.md) for every decay or capture link between tracked
  nuclides and for the transfer from the salt to a removed inventory,
- [FissionProductLoss](FissionProductLoss.md) for decay, capture and removal,
- [MSRDegassingBC](MSRDegassingBC.md) (with zero back pressure) for noble-metal plating on
  `plating_boundaries`.

## Nuclide database

Half-lives, decay branchings, 2200 m/s capture cross sections, fission yields (U-235 thermal and
fast, independent and cumulative) and the chemical group of each nuclide are read from
`data/fission_product_database.json`, generated from ENDF/B-VIII.0 by
`scripts/fission_products/build_fission_product_database.py`. Short-lived isomers that are not
tracked are lumped into their ground state (for example Xe-135m into Xe-135). A different database
can be selected with `database`.

## Effective yields

With the default `yield_type = effective` the direct fission source of each nuclide is its
cumulative yield minus the fission production that reaches it through *tracked* parents by decay,

\begin{equation}
Y_i = \mathrm{CY}_i - \sum_{p\ \text{tracked}} b_{pi}\, \mathrm{CY}_p ,
\end{equation}

so that untracked short-lived precursors (e.g. Te-135) are lumped into their first tracked
descendant and no production is double counted. Tracking `I135 Xe135` thus gives Xe-135 only its
direct yield and the Te-135 production feeds I-135. The `verbose` output prints the yields used.

## Removal pathways

Each entry of `removal_groups` (`noble_gas`, `noble_metal`, `halogen`, `salt_seeking`) defines a
first-order removal pathway with rate `removal_rates` [1/s] and destination `removal_destinations`:

- `offgas`: sparging of noble gases by helium bubbles, entrainment of noble-metal particles in the
  off-gas,
- `plated`: deposition of noble metals on the metal surfaces of the loop; for a well-mixed volume
  the rate is $k_m A/V$,
- `processed`: chemical clean-up (not tracked).

With `track_removed_inventories = true` (default) the action also creates `<nuclide>_offgas` and
`<nuclide>_plated` inventories (per unit salt volume) for every removed nuclide and its tracked
daughters. These inventories are out of the neutron flux but keep decaying, so plated Mo-99 becomes
plated Tc-99 and the total chain inventory is conserved:
$\frac{d}{dt}\sum (N + N_\text{plated} + N_\text{offgas}) = \sum Y F/N_A - \text{decay out of the chain}$.

Resolved wall plating is imposed with `plating_boundaries` and `plating_mass_transfer_coefficient`
$k_m$ [m/s]: the noble metals leave the salt through the flux $J = k_m N$.

## Solver ordering

The Bateman system is lower triangular. The action orders the inventories parent before daughter
(printed with `verbose`); when the systems are listed in that order in `[Problem] linear_sys_names`
and `Executioner/system_names`, one Gauss-Seidel sweep per time step is exact for the coupled chain.
The inventories are small numbers, so set a tight `l_tol` (e.g. `1e-14`).

## Example

!listing test/tests/msr/fission_products/xe_equilibrium.i block=FissionProductInventory

!syntax parameters /FissionProductInventory

!syntax inputs /FissionProductInventory
