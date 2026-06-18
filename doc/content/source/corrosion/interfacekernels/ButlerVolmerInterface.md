# ButlerVolmerInterface

Butler-Volmer mass flux of a metal cation across a metal/salt interface modeled as two mesh blocks.
The this-side variable is the salt cation concentration $c$ and the neighbor variable is the
solid-side metal concentration $c_s$. The interfacial molar flux $J = i_{BV}/(zF)$ enters the salt as
a source and leaves the solid as a sink, so a single object conserves mass and handles both
dissolution (anodic, $i_{BV}>0$) and plating (cathodic, $i_{BV}<0$).

!equation
i_{BV} = i_0\left[\exp\!\left(\frac{\alpha_a z F\,\eta}{RT}\right) - \frac{c}{c_{ref}}\exp\!\left(-\frac{\alpha_c z F\,\eta}{RT}\right)\right], \qquad \eta = (\phi_{metal} - \phi_{salt}) - E_0

The phase potentials may be coupled variables (`phi_salt`, `phi_solid`; fully coupled, AD-exact) or
fixed constants (`salt_potential`, `metal_potential`) for a supporting-electrolyte problem. The
companion [ButlerVolmerPotentialInterface](ButlerVolmerPotentialInterface.md) couples the current into
the potential equations.

!syntax parameters /InterfaceKernels/ButlerVolmerInterface

!syntax inputs /InterfaceKernels/ButlerVolmerInterface

!syntax children /InterfaceKernels/ButlerVolmerInterface
