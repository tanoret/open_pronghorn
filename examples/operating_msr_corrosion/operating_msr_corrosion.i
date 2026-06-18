# =====================================================================================================
# Coupled radiolysis + corrosion over an operating MSR power cycle.
#
# A molten chloride coolant pool sits between two structural-alloy (316 stainless) walls. The reactor
# power is ramped up, held at full power, and ramped down. Two physics are explicitly coupled through
# the dissolved chromium redox couple:
#
#   * CORROSION: the alloy walls dissolve through a Butler-Volmer reaction, releasing chromium into
#     the salt as Cr(II). The release rate is seeded from the validated corrosion correlation.
#
#   * RADIOLYSIS: irradiation produces the chlorine radical chain; the oxidant Cl2.- oxidizes the
#     dissolved Cr(II) -> Cr(III) (and solvated electrons reduce it), exactly the reactions in the
#     chemistry database. The dissolved Cr(II) shed by corrosion is therefore the SAME species the
#     radiolysis network acts on.
#
# The coupling is two-way:
#   - corrosion feeds Cr(II) into the salt (a source for the radiolytic oxidation, gettering the
#     radiolytic oxidant), and
#   - radiolysis sets the Cr(III)/Cr(II) redox state of the salt (its oxidizing power), which is what
#     drives corrosion in the first place (here represented by the oxidizing redox class).
#
# Over the cycle the salt chromium inventory grows steadily from corrosion, while its oxidation state
# Cr(III)/Cr(II) tracks the reactor power: it rises under irradiation (radiolytic oxidation of the
# corrosion product) and relaxes after shutdown as corrosion keeps supplying fresh Cr(II).
#
# Run:  open_pronghorn-opt -i operating_msr_corrosion.i
# View: the Exodus output shows Cl2m_rad (radiolytic oxidant), Cr_II / Cr_III (corrosion product and
#       its oxidized form); the CSV reports the power history, the chromium inventory and redox state.
# =====================================================================================================

dose_full = 1.0e5     # J/m^3/s volumetric dose at full power
T = 873.15            # K (~600 C chloride coolant)
Cl_bulk = 2.0e4       # mol/m^3 bulk chloride
Cr_init = 2.0         # mol/m^3 initial dissolved Cr(II) (legacy corrosion inventory)
D_species = 1.0e-9    # m^2/s molecular diffusivity

[Mesh]
  [pool]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 16
    xmax = 0.2
    ymax = 1.0
  []
[]

[Functions]
  # Reactor power history (volumetric dose rate): ramp up 0->full over 0-100 s, hold to 400 s, ramp
  # down to 0 over 400-500 s, then shut down.
  [power_history]
    type = PiecewiseLinear
    x = '0    100          400          500   100000'
    y = '0    ${dose_full} ${dose_full} 0     0'
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys '
                     'Cr_II_sys Cr_III_sys Cr_I_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  metals = 'Cr'                       # track Cr(II)/Cr(III)/Cr(I) and their radiolytic reactions
  temperature = ${T}
  dose_rate = power_history           # radiolytic source follows the reactor power
  diffusivity = ${D_species}
  # Lump the fast electron/atom kinetics into a net radiolytic-oxidant yield: produce the oxidant
  # Cl2.- directly proportional to dose (as the other examples lump radicals into the net product).
  # This keeps the explicit oxidant -> Cr(II) -> Cr(III) coupling while removing the stiff solvated-
  # electron chemistry, so the power-cycle transient is affordable.
  g_value_species = 'e_sol Cl_rad Cl2m_rad'
  g_value_overrides = '0.0 0.0 0.30'
  initial_condition_species = 'Cl_ion Cr_II'
  initial_condition_values = '${Cl_bulk} ${Cr_init}'
  time_derivative = true
[]

[CorrosionPlatingFlow]
  # Wall corrosion releases chromium into the salt as the radiolysis-tracked Cr(II) species, coupling
  # the two physics through the shared dissolved-chromium redox couple.
  elements = 'Cr'
  release_variables = 'Cr_II'
  reaction_boundary = 'left right'
  temperature = ${T}
  reference_temperature = ${T}
  material_class = stainless_316
  salt_class = chloride
  redox_class = chloride_unspecified
  applied_overpotential = 0.1
  reference_concentration = ${Cr_init}
  temperature_dependent_kinetics = false
[]

[Postprocessors]
  [reactor_power_fraction]
    type = FunctionValuePostprocessor
    function = power_history
    scale_factor = '${fparse 1.0 / dose_full}'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cr_II_total]
    type = ElementIntegralVariablePostprocessor
    variable = Cr_II
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cr_III_total]
    type = ElementIntegralVariablePostprocessor
    variable = Cr_III
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cr_I_total]
    type = ElementIntegralVariablePostprocessor
    variable = Cr_I
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cr_total]
    type = ParsedPostprocessor
    pp_names = 'Cr_II_total Cr_III_total Cr_I_total'
    expression = 'Cr_II_total + Cr_III_total + Cr_I_total'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Salt oxidation state: the coupling observable, set by the radiolytic oxidation of corrosion Cr.
  [redox_ratio_CrIII_CrII]
    type = ParsedPostprocessor
    pp_names = 'Cr_II_total Cr_III_total'
    expression = 'Cr_III_total / (Cr_II_total + 1e-30)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Cl2m_rad_max]
    type = ElementExtremeValue
    variable = Cl2m_rad
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 15
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys '
                 'Cr_II_sys Cr_III_sys Cr_I_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0.0
  end_time = 480.0
  dt = 5.0
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  exodus = true
  csv = true
[]
