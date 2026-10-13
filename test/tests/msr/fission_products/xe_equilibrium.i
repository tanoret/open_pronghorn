# I-135 / Xe-135 equilibrium in a 0D salt volume with burnout and sparging.
# Steady state (exact for implicit Euler at any dt):
#   N_I  = Y_I F / (N_A lambda_I)
#   N_Xe = (Y_Xe F / N_A + b lambda_I N_I) / (lambda_Xe + sigma_Xe phi + r_strip)
# The effective yields are Y_I = CY_I135 (Te-135 untracked, lumped) and Y_Xe = CY_Xe135 - b CY_I135.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  linear_sys_names = 'I135_sys Xe135_sys Cs135_sys Xe135_offgas_sys Cs135_offgas_sys'
[]

[FissionProductInventory]
  nuclides = 'I135 Xe135 Cs135'
  fission_rate_density = 6.2e17
  neutron_flux = 1e18
  removal_groups = 'noble_gas'
  removal_rates = '1e-4'
  removal_destinations = 'offgas'
  verbose = true
[]

[Postprocessors]
  [N_I135]
    type = ElementAverageValue
    variable = I135
  []
  [N_Xe135]
    type = ElementAverageValue
    variable = Xe135
  []
  [N_Xe135_offgas]
    type = ElementAverageValue
    variable = Xe135_offgas
  []
[]

[Executioner]
  type = Transient
  system_names = 'I135_sys Xe135_sys Cs135_sys Xe135_offgas_sys Cs135_offgas_sys'
  dt = 3600
  num_steps = 200
  # The inventories are small numbers (mol/m^3); converge the linear solves tightly.
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
