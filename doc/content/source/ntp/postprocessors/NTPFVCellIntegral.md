# NTPFVCellIntegral

`NTPFVCellIntegral` integrates a functor over blocks with the cell volumes of the finite-volume
system:

\begin{equation}
I = \sum_{C} \phi(\mathbf{x}_C)\, V_C ,
\end{equation}

where $\phi(\mathbf{x}_C)$ is the functor at the cell centroid and $V_C$ is the cell volume of the
finite-volume element information, times the coordinate factor (the volume the linear
finite-volume kernels use). This is exactly the integral a `LinearFVSource` deposits, so a balance
against it closes to the solver tolerance. With the default `functor = 1` it is the finite-volume
volume of the blocks.

[ElementIntegralFunctorPostprocessor](ElementIntegralFunctorPostprocessor.md) and
[VolumePostprocessor](VolumePostprocessor.md) integrate with the element quadrature instead. On
cells with non-planar faces, such as the hexahedra of the 3-D fuel element mapped onto cones
(`examples/ntp_tutorial/05_element_3d`), the first-order quadrature volume differs from the
finite-volume one by a few 1e-4. That shows as an energy balance error of the same size in the
segregated model, although the solution conserves energy to round-off.
`examples/ntp_tutorial/common/element_segregated_postprocessors.i` uses it for the fission power,
and the 3-D segregated inputs use it for the bed volume that normalizes the power density.

!syntax parameters /Postprocessors/NTPFVCellIntegral

!syntax inputs /Postprocessors/NTPFVCellIntegral

!syntax children /Postprocessors/NTPFVCellIntegral
