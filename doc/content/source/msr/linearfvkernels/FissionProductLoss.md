# FissionProductLoss

`FissionProductLoss` adds the implicit first-order loss of a nuclide inventory by radioactive decay,
neutron capture and removal to a linear finite-volume system,

\begin{equation}
L = \left(\lambda + \sigma \phi + \sum_k r_k\right) N ,
\end{equation}

where $\lambda$ is the decay constant, $\sigma$ the one-group capture cross section (barn),
$\phi$ the one-group neutron flux and $r_k$ the removal-rate functors (off-gas sparging,
plating, processing). The loss coefficient goes on the matrix diagonal.

It is normally created by the [FissionProductInventory action](FissionProductInventoryAction.md).

!syntax parameters /LinearFVKernels/FissionProductLoss

!syntax inputs /LinearFVKernels/FissionProductLoss

!syntax children /LinearFVKernels/FissionProductLoss
