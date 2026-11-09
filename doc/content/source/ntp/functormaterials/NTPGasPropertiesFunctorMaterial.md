# NTPGasPropertiesFunctorMaterial

`NTPGasPropertiesFunctorMaterial` defines four functor properties from a
[SinglePhaseFluidProperties](SinglePhaseFluidProperties.md) object (`fp`), evaluated at a pressure
functor and a temperature functor:

\begin{equation}
\rho(p, T), \qquad \mu(p, T), \qquad c_p(p, T), \qquad k(p, T).
\end{equation}

`density_name`, `viscosity_name`, `specific_heat_name` and `conductivity_name` set the property
names (defaults `rho`, `mu`, `cp`, `k`). The object defines nothing else. Two instances can
therefore coexist on the same blocks under different names. In a segregated solve, for example,
one instance can give the flow properties at a lagged state and the other the energy properties
at the current gas temperature.

`NTPGasPropertiesFunctorMaterial` creates `Real` functors, which the linear finite-volume objects
use. Its AD counterpart, `ADNTPGasPropertiesFunctorMaterial`, has the same parameters and
creates AD functors.

!syntax parameters /FunctorMaterials/NTPGasPropertiesFunctorMaterial

!syntax inputs /FunctorMaterials/NTPGasPropertiesFunctorMaterial

!syntax children /FunctorMaterials/NTPGasPropertiesFunctorMaterial

## ADNTPGasPropertiesFunctorMaterial

`ADNTPGasPropertiesFunctorMaterial` is the AD version of the object above. Its pressure and
temperature functors, and the four properties it creates, are AD functors.

!syntax parameters /FunctorMaterials/ADNTPGasPropertiesFunctorMaterial

!syntax inputs /FunctorMaterials/ADNTPGasPropertiesFunctorMaterial

!syntax children /FunctorMaterials/ADNTPGasPropertiesFunctorMaterial
