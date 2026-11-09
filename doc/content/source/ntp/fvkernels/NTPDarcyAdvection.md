# NTPDarcyAdvection

`NTPDarcyAdvection` adds the upwind advection of a functor $q$ (`advected_quantity`) with the face
mass flux $\dot{m}_f$ of a Darcy-Forchheimer pressure field. This is the same flux that
[NTPDarcyMassFlux](NTPDarcyMassFlux.md) uses for the mass equation, and the shared parameters
(`pressure`, `rho_mobility`, `face_mobility`, `coeff_interp_method`, momentum flux and orifice)
have the same meaning. It is used with `advected_quantity` set to the specific enthalpy (or an
enthalpy rise) for the gas energy equation.

The term is assembled in convective form,
$\rho \mathbf{u} \cdot \nabla q = \nabla \cdot (\rho \mathbf{u} q) - q \nabla \cdot (\rho \mathbf{u})$.
Each cell $P$ receives

\begin{equation}
\sum_f \dot{m}_f A_f \left( q_\mathrm{up} - q_P \right),
\end{equation}

with $q_\mathrm{up}$ the upwind value. At a converged flow ($\sum_f \dot{m}_f A_f = 0$), this equals the
conservative divergence exactly. During the iterations, it keeps a positive diagonal, equal to the sum
of $|\dot{m}_f| A_f$ over the inflow faces, even where the current pressure iterate does not
conserve mass. The conservative form has a zero diagonal in a cell with inflow and no outflow, which
makes the Newton Jacobian singular.

On a boundary, the flux is computed from the boundary pressure gradient. An outflow face
contributes nothing (it carries the cell value). An inflow face on one of `inflow_boundaries` brings
in `boundary_value`, which is required when `inflow_boundaries` is set. Only differences of $q$
enter, so $q$ may be shifted by a constant if `boundary_value` is shifted the same way. An example
is the enthalpy rise $h - h_\mathrm{in}$ with `boundary_value = 0`.

!syntax parameters /FVKernels/NTPDarcyAdvection

!syntax inputs /FVKernels/NTPDarcyAdvection

!syntax children /FVKernels/NTPDarcyAdvection
