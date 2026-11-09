# NTPSliceUniformityConstraint

`NTPSliceUniformityConstraint` assembles equations for a scalar unknown $\lambda$ (`lambda`), with one
component per slice of a [mesh division](syntax/MeshDivisions/index.md) (`mesh_division`). The
equations make the slice averages of a `field` $q$ uniform. A typical use is the logarithm of
the orifice multiplier of each axial slice of a fuel element, solved together with the flow.
$\langle q \rangle_j$ is the volume average of $q$ over slice $j$ within the kernel's blocks.

With `reference = slice` (default) and reference slice $r$ (`reference_slice`), the residuals are

\begin{equation}
R_j = w \frac{\langle q \rangle_j - \langle q \rangle_r}{s} + (1 - w + \delta)\, \lambda_j \quad (j \ne r),
\qquad R_r = \lambda_r .
\end{equation}

With `reference = mean`, the slice averages are equalized with their volume average $q_\mathrm{mean}$,
and the mean of the components is fixed to $\lambda_\mathrm{mean}$ (`mean_log_multiplier`):

\begin{equation}
R_j = w \frac{\langle q \rangle_j - q_\mathrm{mean}}{s} + (1 - w + \delta)(\lambda_j - \lambda_\mathrm{mean})
\quad (j < n - 1), \qquad R_{n-1} = \frac{1}{n} \sum_k \lambda_k - \lambda_\mathrm{mean}.
\end{equation}

With the mean reference, no slice is pinned. With a fixed reference, some slice might need its drag
removed altogether ($\lambda_j \to -\infty$, where the slice average no longer depends on
$\lambda_j$). With the mean reference, that slice keeps a finite drag, and the level
$\exp(\lambda_\mathrm{mean})$ sets the margin. The order of `lambda` must equal the number of slices.

The scale $s$ (`field_scale`, or `field_scale_pp` floored by `field_scale`) makes the equations
dimensionless, so their conditioning does not depend on the power level. The regularization $\delta$
(`regularization`, default $10^{-3}$) keeps $\lambda$ defined, and pins it to 0, where the field does
not depend on it, e.g. a uniform field at zero power. The activation $w(t)$ (`activation`, clipped to
$[0, 1]$) switches the design on smoothly once the field differs between slices: before the gas is
heated, the orifice has no effect on the field. Every slice must contain cells of the kernel's
blocks.

Nothing is added to the equation of the kernel's `variable`. The Jacobian with respect to the field
comes from automatic differentiation. The scalar must enter the field equations for the system to be
well posed, here through `orifice_log_multiplier` of [NTPDarcyMassFlux](NTPDarcyMassFlux.md) and
[NTPDarcyAdvection](NTPDarcyAdvection.md). `examples/ntp_tutorial/03_ludewig_anchor/orifice_implicit.i`
uses it.

!syntax parameters /FVKernels/NTPSliceUniformityConstraint

!syntax inputs /FVKernels/NTPSliceUniformityConstraint

!syntax children /FVKernels/NTPSliceUniformityConstraint
