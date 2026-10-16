# CorrosionLinearFVButlerVolmerBC

Butler-Volmer corrosion / plating electrode reaction as a linear finite-volume wall boundary
condition, for a flowing molten salt solved on the segregated (SIMPLE/PIMPLE) basis. It is the
linear-FV counterpart of [ButlerVolmerBC](../bcs/ButlerVolmerBC.md), used by the
[CorrosionPlatingFlow](../actions/CorrosionPlatingFlowAction.md) action so the dissolved corrosion
products are passive scalars in the same flow solve as the radiolysis chemistry and energy.

With the concentration-explicit Butler-Volmer kinetics and a prescribed (metal minus salt)
overpotential, the molar flux of the cation into the salt is linear in the cell concentration,

!equation
J_{in} = A - B\,c, \qquad
A = \frac{i_0}{zF}\exp\!\left(\frac{\alpha_a z F \eta}{R T}\right), \qquad
B = \frac{i_0}{zF\,c_{ref}}\exp\!\left(-\frac{\alpha_c z F \eta}{R T}\right),

so it assembles directly into the linear system without lagging (the outward boundary flux is
$B c - A$, returned as the gradient matrix and right-hand-side contributions). The exchange current is
Arrhenius in the local temperature, $i_0(T) = i_0\,\exp[(E_a/R)(1/T_{ref} - 1/T)]$, so the wall
corrodes faster where the salt is hotter; the temperature, metal potential and salt potential are read
as functors, which lets the boundary couple one-way to a flow / energy / radiolysis solution. A
diffusion kernel must be present on the variable for the wall flux to be assembled.

## Temperature-dependent equilibrium

With a constant overpotential the zero-current concentration
$c_{eq} = c_{ref}\exp[(\alpha_a + \alpha_c) z F \eta / R T]$ *falls* with temperature, so a closed
loop would plate at its hot walls and dissolve at its cold walls, the opposite of what is observed.
When `equilibrium_concentration` is positive the overpotential is instead derived from a
temperature-dependent equilibrium,

!equation
c_{eq}(T) = c_{eq,ref}\,\exp\!\left[-\frac{\Delta H}{R}\left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right],
\qquad
\eta(T) = \frac{RT}{(\alpha_a + \alpha_c) z F}\ln\frac{c_{eq}(T)}{c_{ref}},

so the net flux vanishes at $c = c_{eq}(T)$. A positive dissolution enthalpy $\Delta H$
(`dissolution_enthalpy`) makes the equilibrium rise with temperature: hot walls dissolve and cold
walls plate. Different elements are given different equilibria, so a noble component (Ni) barely
dissolves while Cr does. The potential functors must not be given in this mode.

!syntax parameters /LinearFVBCs/CorrosionLinearFVButlerVolmerBC

!syntax inputs /LinearFVBCs/CorrosionLinearFVButlerVolmerBC

!syntax children /LinearFVBCs/CorrosionLinearFVButlerVolmerBC
