# Merged after continuation.i for the 3-D inputs: the 3-D system has a higher round-off floor of
# the residual than the axisymmetric one, so the absolute tolerance of the continuation steps is
# raised (the balances stay at 1e-7).
[Executioner]
  nl_abs_tol = 1e-5
[]
