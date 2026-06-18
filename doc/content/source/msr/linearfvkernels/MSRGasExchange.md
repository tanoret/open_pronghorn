# MSRGasExchange

`MSRGasExchange` adds the gas-liquid exchange of a dissolved diatomic (Cl$_2$, F$_2$) between the
molten salt and a well-mixed headspace to a linear finite-volume system. It implements the lumped
two-volume model of the validated 0D radiolysis model:

\begin{equation}
R = k_L a \left( C_\text{liq} - k_H\, p_\text{gas} \right), \qquad p_\text{gas} = C_\text{gas}\, R_g\, T,
\end{equation}

where $k_L a$ is the overall mass-transfer coefficient (1/s), $k_H$ the Henry coefficient
(mol/(m$^3$ Pa)), and $R_g$ the universal gas constant. The gas inventory is carried as a gas-phase
concentration $C_\text{gas} = n_\text{gas}/V_\text{gas}$, so the ideal-gas pressure reduces to the
purely local expression $p_\text{gas} = C_\text{gas} R_g T$ and no auxiliary postprocessor is
required. The dissolved and gas equations are

\begin{equation}
\frac{\mathrm{d}C_\text{liq}}{\mathrm{d}t} = -R, \qquad
\frac{\mathrm{d}C_\text{gas}}{\mathrm{d}t} = \frac{V_\text{liq}}{V_\text{gas}}\, R,
\end{equation}

which conserve the total amount of the diatomic. Each instance assembles one of the two equations,
selected with `mode`:

- `mode = liquid` (variable $C_\text{liq}$): implicit sink $k_L a$ on the diagonal and an explicit
  equilibrium-restoring source $k_L a\, k_H R_g T\, C_\text{gas}$ on the right-hand side.
- `mode = gas` (variable $C_\text{gas}$): implicit sink $(V_\text{liq}/V_\text{gas})\,k_L a\, k_H R_g T$
  on the diagonal and an explicit gain $(V_\text{liq}/V_\text{gas})\,k_L a\, C_\text{liq}$ on the
  right-hand side.

The variable's own term is implicit (a positive diagonal) and the partner concentration is lagged on
the right-hand side, matching the segregated fixed-point linearization. At steady state $R = 0$ and
the dissolved concentration satisfies the Henry-law partition $C_\text{liq} = k_H\, p_\text{gas}$.

This object is normally created by the [MoltenSaltRadiolysis action](MoltenSaltRadiolysisAction.md)
for each gas listed in `gas_species`, but it can also be used directly in a `[LinearFVKernels]` block.

## Gas-phase transport

The gas-phase concentration created for each diatomic is itself transported: the action adds buoyant
advection (the flow velocity plus an upward rise/slip velocity, set with `gas_buoyancy_velocity`) and
dispersion (`gas_dispersivity`) for the gas-phase species, so the evolved gas rises and spreads and
can be vented at a free surface with a `LinearFVAdvectionDiffusionOutflowBC`. Direct escape of the
*dissolved* species at a free surface to a swept cover gas is handled separately by
[MSRDegassingBC](MSRDegassingBC.md). Together these provide the volumetric (bubble) and surface
(degassing) off-gas pathways.

!syntax parameters /LinearFVKernels/MSRGasExchange

!syntax inputs /LinearFVKernels/MSRGasExchange

!syntax children /LinearFVKernels/MSRGasExchange
