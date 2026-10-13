# FissionProductFeed

`FissionProductFeed` adds the production of a nuclide from a parent inventory $N_p$ to a linear
finite-volume system,

\begin{equation}
S = \left(\lambda_p b + \sigma_p \phi + r\right) N_p ,
\end{equation}

where $\lambda_p b$ is the partial decay constant of the parent to this nuclide, $\sigma_p \phi$ the
parent capture rate (when this nuclide is the capture product) and $r$ a first-order transfer rate,
e.g. the rate at which a salt inventory is moved to a plated or off-gas inventory. The parent is
evaluated at its latest iterate and goes on the right-hand side; ordering the parent's system before
this one makes the coupling implicit within one fixed-point sweep.

It is normally created by the [FissionProductInventory action](FissionProductInventoryAction.md).

!syntax parameters /LinearFVKernels/FissionProductFeed

!syntax inputs /LinearFVKernels/FissionProductFeed

!syntax children /LinearFVKernels/FissionProductFeed
