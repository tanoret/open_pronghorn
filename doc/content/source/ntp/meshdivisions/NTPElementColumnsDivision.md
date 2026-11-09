# NTPElementColumnsDivision

`NTPElementColumnsDivision` divides a cylindrical fuel element into columns: `n_axial` equal
slices along the element axis between `axial_min` and `axial_max`, times `n_azimuthal` equal
sectors of the angle range [`theta_min`, `theta_max`] about the axis. Slice $j$ (from `axial_min`)
and sector $k$ (from `theta_min`) form the division

\begin{equation}
i = j\, n_\theta + k .
\end{equation}

The axial coordinate is measured along `axis_direction` from `center`, the angle from
`azimuthal_start` ($e_1$) towards $e_1 \times$ `axis_direction`. With the default axis $y$ and start
$x$, $\theta = \mathrm{atan2}(z, x)$, so the default range $[0, \pi]$ is the half element of
`examples/ntp_tutorial/05_element_3d`. Cell centroids slightly outside the axial range (a faceted
mesh) or the angle range (round-off on a symmetry plane) belong to the nearest slice or sector,
so every cell has a division.

With one sector the division is the axial slicing of the orifice design of
`examples/ntp_tutorial/04_navier_stokes/orifice_flow.i`, where it feeds
[NTPRhieChowSliceFlow](NTPRhieChowSliceFlow.md) and
[NTPOrificeFlowDesignAux](NTPOrificeFlowDesignAux.md). With several sectors (`orifice_sectors` in
that input) the design shapes the cold-frit flow around the element as well as along it, which the
3-D segregated element needs to reach the design power.

!syntax parameters /MeshDivisions/NTPElementColumnsDivision

!syntax inputs /MeshDivisions/NTPElementColumnsDivision

!syntax children /MeshDivisions/NTPElementColumnsDivision
