# ButlerVolmerPotentialInterface

Butler-Volmer charge flux coupling the interfacial current into the current-continuity (potential)
equations on the two blocks. The this-side variable is the salt-phase potential $\phi_{salt}$ and the
neighbor variable is the metal-phase potential $\phi_{solid}$; the Butler-Volmer current $i_{BV}$
enters the salt charge balance as a source and leaves the metal charge balance as a sink, conserving
charge across the double layer. Add one instance per tracked element (the currents are additive). Use
together with [ButlerVolmerInterface](ButlerVolmerInterface.md) when the salt potential is solved and
the interfacial current must feed back on the field.

!syntax parameters /InterfaceKernels/ButlerVolmerPotentialInterface

!syntax inputs /InterfaceKernels/ButlerVolmerPotentialInterface

!syntax children /InterfaceKernels/ButlerVolmerPotentialInterface
