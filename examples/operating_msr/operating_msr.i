# =====================================================================================================
# Operating MSR example: radiolytic fluorine off-gassing over a reactor power cycle.
#
# A section of fluoride fuel salt (FLiBe-relevant) sits below a swept cover gas. The reactor power is
# ramped up, held at full power, and ramped back down; the radiolytic F2 production follows the power
# history. The dissolved F2:
#   * is generated volumetrically by radiolysis (G-value source scaled by the power history),
#   * transfers to the gas phase (bubbles) by interfacial mass transfer (MSRGasExchange),
#   * degasses directly at the free surface to the cover gas (MSRDegassingBC),
# while the evolved gas phase rises buoyantly, disperses, and vents at the top free surface. The
# dissolved and gas-phase F2 inventories build up during operation and decay after shutdown, and the
# off-gas tracks the power.
#
# The radical yields are set to zero so that the model follows the net F2 yield (G = 0.005, the Davis
# 2022 FLiBe value); this isolates the off-gas dynamics and keeps the (otherwise stiff) chemistry
# inexpensive. A stagnant pool is used (transport by diffusion, gas buoyancy and venting); for flow
# coupling see examples/molten_salt_corrosion.
#
# Run:  open_pronghorn-opt -i operating_msr.i
# View: the Exodus output shows F2_diss and F2_gas building and clearing over the power cycle; the
#       CSV reports the power history and the dissolved / gas-phase / total F2 inventories.
# =====================================================================================================

dose_full = 1.0e5     # J/m^3/s volumetric dose at full power
kH_F2 = 1.0e-6        # Henry coefficient for F2 [mol/(m^3 Pa)]
T = 873.15            # K (FLiBe ~600 C)

[Mesh]
  [pool]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 20
    xmax = 0.4
    ymax = 1.0
  []
[]

[Functions]
  # Reactor power history (volumetric dose rate): ramp up 0->full over 0-100 s, hold to 400 s,
  # ramp down to 0 over 400-500 s, then shut down.
  [power_history]
    type = PiecewiseLinear
    x = '0    100        400        500   100000'
    y = '0    ${dose_full} ${dose_full} 0     0'
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys F2_gas_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = fluoride
  gas_species = 'F2'
  temperature = ${T}
  dose_rate = power_history          # radiolytic source follows the reactor power
  g_value_species = 'e_sol F_rad'    # follow the net F2 yield only (radicals lumped out)
  g_value_overrides = '0.0 0.0'
  diffusivity = 1.0e-9               # dissolved-species molecular diffusivity
  kLa = 0.02                         # interfacial transfer dissolved -> gas (bubbles)
  volume_ratio = 1.0
  gas_buoyancy_velocity = '0 0.05 0' # buoyant rise of the evolved gas
  gas_dispersivity = 1.0e-3
[]

[LinearFVBCs]
  # Evolved gas vents through the free surface.
  [vent_F2_gas]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = F2_gas
    boundary = top
    use_two_term_expansion = false
  []
  # Dissolved F2 degasses directly at the free surface to the swept (clean) cover gas.
  [degas_F2_diss]
    type = MSRDegassingBC
    variable = F2_diss
    boundary = top
    mass_transfer_coefficient = 5.0e-4
    henry_coefficient = ${kH_F2}
    cover_gas_pressure = 0.0
  []
[]

[Postprocessors]
  [reactor_power_fraction]
    type = FunctionValuePostprocessor
    function = power_history
    scale_factor = '${fparse 1.0 / dose_full}'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [F2_dissolved_total]
    type = ElementIntegralVariablePostprocessor
    variable = F2_diss
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [F2_gas_total]
    type = ElementIntegralVariablePostprocessor
    variable = F2_gas
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [F2_diss_max]
    type = ElementExtremeValue
    variable = F2_diss
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 10
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys F2_gas_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = 700.0
  dt = 10.0
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  exodus = true
  csv = true
[]
