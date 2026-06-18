# MSRReaction

`MSRReaction` adds the contribution of a single mass-action chemical reaction to one molten salt
radiolysis species in a linear finite-volume system. The forward rate of a reaction is

\begin{equation}
r = k(T) \prod_{j \in \text{reactants}} [S_j]^{o_j},
\end{equation}

where $[S_j]$ are the species concentrations (mol/m$^3$), $o_j$ the kinetic orders, and $k(T)$ the
temperature-dependent rate constant evaluated with the Arrhenius forms of the
[molten salt radiolysis database](MoltenSaltRadiolysisAction.md):

\begin{equation}
k(T) =
\begin{cases}
A\,\exp\!\left(-E_a/(R T)\right), & A,\,E_a \text{ supplied},\\[4pt]
k_\text{ref}\,\exp\!\left[\dfrac{E_a}{R}\left(\dfrac{1}{T_\text{ref}}-\dfrac{1}{T}\right)\right], & k_\text{ref},\,T_\text{ref},\,E_a \text{ supplied},\\[4pt]
k_\text{ref}, & \text{only } k_\text{ref} \text{ supplied}.
\end{cases}
\end{equation}

A given species participates in a reaction either as a reactant (it is consumed) or as a product
(it is produced). Following the linear finite-volume convention used by `LinearFVReaction` (an
implicit matrix term $c\,u$) and `LinearFVSource` (an explicit right-hand side term), the two cases
are assembled differently and selected with the `mode` parameter.

## Consumption (implicit sink)

When the species $i$ is a reactant with stoichiometric coefficient $\nu_i$ and kinetic order $o_i$,
its consumption rate is $\nu_i\,r$. One power of the variable is kept implicit so that the
contribution sits on the matrix diagonal:

\begin{equation}
\underbrace{\nu_i\, k(T)\, [i]^{\,o_i-1} \prod_{m \ne i} [S_m]^{o_m}}_{\text{matrix coefficient}} \; [i] \; = \; \nu_i\, r .
\end{equation}

The partner reactant concentrations $[S_m]$ are evaluated as functors. Keeping one power of $[i]$
implicit guarantees a non-negative diagonal, which preserves the M-matrix property and the
convergence of the segregated fixed-point iteration. Lagged concentrations are clamped at zero
inside the rate products.

## Production (explicit source)

When the species $i$ is a product with stoichiometric coefficient $\nu_i$, the production rate
$\nu_i\,r$ is added to the right-hand side using the lagged reactant concentrations:

\begin{equation}
\text{RHS} \mathrel{+}= \nu_i\, k(T) \prod_{j \in \text{reactants}} [S_j]^{o_j}.
\end{equation}

## Coupling and linearization

The reactant concentrations of other species are read as functors at the current element and state.
Within the segregated multi-system fixed-point solve, these return the latest available iterate, so
the mass-action coupling is linearized in a Gauss-Seidel sense. For weakly coupled networks a few
fixed-point iterations per time step are sufficient; stiffer networks require more iterations or
smaller time steps.

This object is normally created by the [MoltenSaltRadiolysis action](MoltenSaltRadiolysisAction.md),
which instantiates one `MSRReaction` per (reaction, species) pair from the built-in database, but it
can also be used directly in a `[LinearFVKernels]` block.

!syntax parameters /LinearFVKernels/MSRReaction

!syntax inputs /LinearFVKernels/MSRReaction

!syntax children /LinearFVKernels/MSRReaction
