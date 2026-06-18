# CorrosionSolidDiffusion

Solid-state diffusion of a metal species in the alloy (e.g. chromium depletion feeding the interface
dissolution). The residual is the standard $(\nabla\psi,\; D_s\nabla c_s)$. The solid-state
diffusivity is supplied as a constant; the [CorrosionPlating](CorrosionPlatingAction.md) action sets
it from the calibrated Arrhenius correlation $D_s(T)$ at the operating temperature.

!syntax parameters /Kernels/CorrosionSolidDiffusion

!syntax inputs /Kernels/CorrosionSolidDiffusion

!syntax children /Kernels/CorrosionSolidDiffusion
