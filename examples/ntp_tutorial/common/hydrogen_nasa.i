# Hydrogen gas properties used by every MOOSE input of the NTP tutorial.
#
# NTPHydrogenFluidProperties is the ideal-gas NASA-polynomial hydrogen of the reference channel model:
#   thermal ideal gas rho = p/(R T); NASA 7-coefficient cp(T), GRI-Mech 3.0 H2 fit,
#   switch at 1000 K; Sutherland mu(T); Eucken k(T) = mu (cp + 5/4 R).
# The same equations are written in Python in common/hydrogen.py, which the
# post-processing scripts use. Stage 1 checks the two against each other.
#
# Include this file with  !include ../common/hydrogen_nasa.i  from a stage directory.

[FluidProperties]
  [h2]
    type = NTPHydrogenFluidProperties
  []
[]
