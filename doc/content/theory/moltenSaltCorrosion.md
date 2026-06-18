# Molten salt corrosion and plating theory

The corrosion and plating framework models the electrochemical loss of structural metal into a molten
salt (and its redeposition) as a mechanistic, spatial finite-element problem. It couples species
transport, the electric potential, and electrode kinetics, and is parameterized from a validated
effective Butler-Volmer correlation.

## Species transport (Nernst-Planck)

Each dissolved cation $M^{z_k+}$ obeys mass conservation with a Nernst-Planck flux,

!equation
\frac{\partial c_k}{\partial t} + \nabla\cdot\vec{N}_k = 0, \qquad
\vec{N}_k = \underbrace{-D_k\nabla c_k}_{\text{diffusion}} \;\underbrace{-\; \frac{z_k F}{RT}D_k c_k\nabla\phi}_{\text{electromigration}} \;+\; \underbrace{\vec{u}\,c_k}_{\text{advection}} ,

where $c_k$ is the molar concentration [mol/m^3], $D_k$ the diffusivity, $z_k$ the charge number,
$\phi$ the electric potential, and $\vec{u}$ the salt velocity. The migration term is the
electrophoresis driving force. In the solid alloy the controlling metal (chromium) diffuses by
solid-state diffusion with an Arrhenius diffusivity $D_s(T) = D_{s,\text{ref}}\exp[(E_a/R)(1/T_{ref}-1/T)]$.

## Electric potential (current continuity)

Under the dilute-solution / MacInnes closure the ionic current density is

!equation
\vec{i} = F\sum_k z_k \vec{N}_k = -F\sum_k z_k D_k \nabla c_k - \kappa\nabla\phi, \qquad
\kappa = \frac{F^2}{RT}\sum_k z_k^2 D_k c_k ,

with $\kappa$ the ionic conductivity. Charge conservation $\nabla\cdot\vec{i} = 0$ closes the system
for $\phi$. Electroneutrality fixes the relative ion populations; the otherwise-undetermined constant
mode of $\phi$ is pinned with an applied reference potential. For a supporting electrolyte $\kappa$ is
dominated by inert ions and is nearly constant, and the minor corroding cations migrate as dilute
tracers.

## Electrode kinetics (Butler-Volmer)

At the metal/salt interface the reaction $M(s)\rightleftharpoons M^{z+} + z\,e^-$ proceeds at the
Butler-Volmer rate, written in concentration-explicit form:

!equation
i_{BV} = i_0\left[\exp\!\left(\frac{\alpha_a z F\,\eta}{RT}\right) - \frac{c}{c_{ref}}\exp\!\left(-\frac{\alpha_c z F\,\eta}{RT}\right)\right], \qquad
\eta = (\phi_{metal} - \phi_{salt}) - E_0 .

The anodic (dissolution) branch carries the solid metal activity, taken constant ($=1$), so the
corrosion rate is set by the overpotential and stays bounded even for fresh, nearly clean salt; the
cathodic (plating) branch carries the dissolved-ion activity $c/c_{ref}$, so deposition accelerates as
the salt loads. The reaction reaches equilibrium ($i_{BV}=0$) when
$c = c_{ref}\exp[(\alpha_a+\alpha_c)\,zF\eta/RT]$. A positive (anodic) current dissolves the metal into
the salt (corrosion); a negative (cathodic) current plates it out (deposition). The molar flux follows
Faraday's law, $J = i_{BV}/(zF)$, and the metal recession rate is $v = M\,i_{BV}/(zF\rho)$. The same
kinetics are applied as an interface condition between a salt block and a solid block, or as a boundary
condition when only one phase is modeled.

## Calibration and validation

The exchange current $i_0$, the effective overpotential and the solid chromium diffusivity are seeded
from the calibrated effective Butler-Volmer correlation
[MoltenSaltCorrosionModel](source/corrosion/base/MoltenSaltCorrosionModel.md), a C++ port of the
validated reference model fit to 76 molten salt corrosion/plating cases and 43 measurement targets
(median direct-target factor error 1.11; 92.9 % within a factor of two). The port is reproduced term
for term by the unit tests, and the mechanistic framework reproduces the reference dissolution rate of
all 76 cases end-to-end; see `validation/corrosion`.
