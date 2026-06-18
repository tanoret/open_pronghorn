# ButlerVolmerBC

Butler-Volmer electrode reaction applied as a boundary condition for a single-domain (salt-only or
solid-only) corrosion/plating problem, where the metal/salt interface is an external boundary. The
same kinetics as [ButlerVolmerInterface](ButlerVolmerInterface.md) are used; the absent counter-phase
is supplied through functors (`applied_potential`, `counter_concentration`).

The boundary can add the molar flux $J = i_{BV}/(zF)$ to a concentration variable
(`flux_type = species`) or the current $i_{BV}$ to a potential variable (`flux_type = charge`). When
the salt is modeled (`metal_domain = false`) the reacting concentration is the boundary concentration
variable; when the solid is modeled (`metal_domain = true`) the salt concentration is supplied as a
fixed functor.

!syntax parameters /BCs/ButlerVolmerBC

!syntax inputs /BCs/ButlerVolmerBC

!syntax children /BCs/ButlerVolmerBC
