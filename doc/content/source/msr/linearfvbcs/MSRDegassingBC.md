# MSRDegassingBC

`MSRDegassingBC` is a surface degassing (outgassing) boundary condition for a dissolved molten salt
radiolysis species at a free surface or cover-gas interface. The dissolved species leaves the melt
across the surface at the convective mass-transfer rate

\begin{equation}
J = k_\text{surf}\,\left(C - k_H\, p_\text{cover}\right),
\end{equation}

where $k_\text{surf}$ is the surface mass-transfer coefficient [m/s], $k_H$ the Henry coefficient
[mol/(m$^3$ Pa)] and $p_\text{cover}$ the partial pressure of the species in the cover gas [Pa]. It
is the mass-transfer analogue of a convective heat-transfer boundary condition and represents direct
escape of a volatile radiolysis product (e.g. Cl$_2$, F$_2$) to a swept cover gas — the off-gas
pathway — as distinct from the volumetric bubble exchange handled by [MSRGasExchange](MSRGasExchange.md).

By default `cover_gas_pressure = 0`, i.e. a continuously swept clean cover gas, so the surface acts
as a first-order sink that removes all dissolved species above zero. A non-zero
`cover_gas_pressure` (constant or functor) imposes a finite equilibrium concentration
$C_\text{eq} = k_H\, p_\text{cover}$ below which the surface no longer degasses.

The flux already carries the full transfer coefficient, so the diffusion kernel does not multiply it
by the molecular diffusivity (`includesMaterialPropertyMultiplier` returns `true`). The species must
have a diffusion (dispersion) kernel for this boundary flux to be assembled; the
[MoltenSaltRadiolysis action](MoltenSaltRadiolysisAction.md) adds one whenever a `diffusivity` is
supplied.

For a well-mixed cell of surface area $A$ and volume $V$ with a clean cover gas, the dissolved
inventory then decays as $C(t) = C_0\exp(-k_\text{surf} A\,t / V)$; this is verified analytically in
`test/tests/msr/gas_transport`.

!syntax parameters /LinearFVBCs/MSRDegassingBC

!syntax inputs /LinearFVBCs/MSRDegassingBC

!syntax children /LinearFVBCs/MSRDegassingBC
