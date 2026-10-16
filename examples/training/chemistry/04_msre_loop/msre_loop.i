# =====================================================================================================
# Training case 04: the MSRE fuel loop as a lumped multi-region model, checked against the MSRE record.
#
# Why 0D: the loop transit time (26 s) is far shorter than every chemistry time (Xe-135 9.1 h, Mo-99
# 66 h, corrosion years), so the salt is well mixed (case 05 shows that a 2D loop inventory equals the
# volume-averaged 0D one to ~1e-3). The regions enter through their rates: the core (flux, fission),
# the pump-bowl spray ring (stripping), the graphite (Xe and noble-metal sink, Xe burnout) and the
# 79.2 m^2 of Hastelloy N (noble-metal plating, Cr corrosion). Inventories are per m^3 of salt.
#
# Data: ORNL-4865 (Compere 1975): 7.4 MW, 2.01 m^3 salt (0.65 m^3 in the core, 88 % of the fissions),
# loop-averaged thermal flux 2.9e12 n/cm^2/s, 79.2 m^2 metal, 225 m^2 graphite. ORNL-TM-3464 (Engel &
# Steffy 1971): spray 50 of 1200 gpm, salt-to-graphite Xe mass transfer 0.01-0.18 ft/h, ORNL's own
# no-bubble prediction 1.08 % dk/k, MEASURED 0.26-0.28 % dk/k. ORNL-TM-3063 (McCoy 1970): salt Cr
# 40 -> ~100 ppm over a few weeks after the U-233 start, equivalent to 0.3 mil (7.6 um) Cr removal.
#
# 1. Xenon. Removal of the salt Xe: spray r_s = (50/1200) eps / 26 s (eps = 0.10), graphite uptake
#    r_g = k_g A_g / V (k_g = 0.06 ft/h) into Xe135_graphite, which decays and burns out in the core
#    flux phi_c = phi V / V_core x 0.88 = 7.9e16 n/m^2/s, and circulating bubbles r_b. Poisoning
#    rho = sigma (phi N_salt + phi_c N_graphite) / (nu F), without importance weighting.
#    A priori (r_b = 0): 1.033 % (CODE-TO-MODEL check against ORNL's 1.08 %); 97 % of it is the Xe in
#    the graphite. Measured 0.26-0.28 %: the a priori model over-predicts 4x -- the historical MSRE
#    finding that circulating helium bubbles strip xenon. CALIBRATION: r_b = 2.14e-3 1/s gives 0.270 %
#    (run with r_bubble=2.1414e-3). PREDICTION (not validated, TM-3464 gives no clean measured value):
#    the poisoning reaches 50 % / 90 % of equilibrium 14 h / 38 h after a step to full power.
#    Note: the action takes one removal rate per (group, destination); spray and bubbles are summed.
# 2. Noble metals Mo-99, Ru-103, Te-132: plating r_m = k_m A_metal / V (k_m = 1e-5 m/s), graphite
#    r_mg = k_mg A_g / V (k_mg = 2e-6 m/s, into <X>_graphite), off-gas 5e-7 1/s. VALIDATION against
#    ORNL-4865 Table 12.3 (salt 1-20 %, metal 5-30 % (Te 20-90 %), graphite 5-30 %): the model gives
#    salt 0.4-0.5 % (Ru-103 0.06 %), metal 63 %, graphite 36 %, off-gas 0.08 % -- salt too low and
#    metal too high: deposition at the mass-transfer limit with perfect sticking does not describe the
#    MSRE, as ORNL concluded.
# 3. Corrosion: CorrosionPlating salt_only on a salt column of depth V / A_metal = 2.54 cm with the
#    MSRE features of targets M-003/M-004 (hastelloy_n, fluoride_fuel, msre_or_fuel_baseline,
#    core_or_reactor, flow 2, 650 C): 3.00 um/y. Salt Cr rises by +26 ppm in 3 weeks and +50 ppm in
#    6 weeks (McCoy: +60 ppm in "a few weeks"); a CALIBRATION check, since M-003/M-004 are fit targets.
#
# VERIFICATION: xe_poison_ratio (= 1.0000), Mo_salt_ratio (= 1.0000), Ru_salt_fraction_equilibrium
# (analytic), rate_vs_correlation (= 1 at t = 0) and faraday_ratio (= 1).
# =====================================================================================================

# ---- MSRE data (ORNL-4865, ORNL-TM-3464) -------------------------------------------------------
P_th = 7.4e6 # W, isotopic full power (ORNL-4865; 8 MW nominal)
V = 2.01 # m^3 circulating fuel salt (71 ft^3)
V_core = 0.65 # m^3 fuel salt in the core (23 ft^3)
f_core = 0.88 # fraction of the fissions in the core channels
phi = 2.9e16 # n/m^2/s thermal flux averaged over the circulating salt (2.9e12 n/cm^2/s)
phi_core = ${fparse phi * V / V_core * f_core} # n/m^2/s mean flux in the core (graphite Xe burnout)
A_metal = 79.2 # m^2 Hastelloy N wetted by fuel salt (852 ft^2)
A_graphite = 225.0 # m^2 graphite surface
nu = 2.43 # neutrons per fission (U-235)
P = ${fparse P_th / V} # W/m^3 loop-averaged fission power density (3.68 MW/m^3)
Fd = ${fparse P / (200e6 * 1.602176634e-19)} # fissions/m^3/s

# ---- xenon removal pathways [1/s] ------------------------------------------------------------
eps_spray = 0.10 # a priori spray-ring stripping efficiency
r_spray = ${fparse 50 / 1200 * eps_spray / 26} # 50 of 1200 gpm through the spray, 26 s loop transit
k_g = ${fparse 0.06 * 0.3048 / 3600} # m/s salt-to-graphite Xe mass transfer (0.06 ft/h, TM-3464)
r_graphite = ${fparse k_g * A_graphite / V}
r_bubble = 0 # circulating-bubble stripping (0 = design model; 2.14e-3 = calibrated to 0.27 %)
r_xe_offgas = ${fparse r_spray + r_bubble} # one pathway per (group, destination): see the header

# ---- noble-metal deposition [1/s] -----------------------------------------------------------
k_metal = 1.0e-5 # m/s turbulent salt-to-metal mass transfer (pipes, heat exchanger)
k_mg = 2.0e-6 # m/s salt-to-graphite mass transfer (laminar core channels)
r_metal = ${fparse k_metal * A_metal / V}
r_mgraph = ${fparse k_mg * A_graphite / V}
r_nm_offgas = 5.0e-7 # 1/s noble-metal entrainment to the off-gas (order of magnitude)

# ---- corrosion (CorrosionPlating salt_only, MSRE features of targets M-003/M-004) ---------------
T = 923.15 # K (650 C)
L = ${fparse V / A_metal} # m salt volume per wetted metal area
c0 = ${fparse 40 * 2.27 / 51.9961} # mol/m^3 Cr, 40 ppm at the U-233 start (TM-3063)
rho_salt = 2.27 # g/cm^3 MSRE fuel salt
x_Cr = 0.07 # Cr weight fraction of Hastelloy N

day = 86400
N_A = 6.02214076e23
sx = 2647600e-28 # m^2 Xe-135 capture cross section

[Mesh]
  [loop]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
    xmax = ${L}
  []
[]

[Problem]
  nl_sys_names = 'nl0'
  linear_sys_names = 'I135_sys Xe135_sys Mo99_sys Ru103_sys Te132_sys Xe135_offgas_sys Mo99_offgas_sys '
                     'Ru103_offgas_sys Te132_offgas_sys Mo99_plated_sys Ru103_plated_sys Te132_plated_sys '
                     'Xe135_graphite_sys Mo99_graphite_sys Ru103_graphite_sys Te132_graphite_sys'
[]

[FissionProductInventory]
  nuclides = 'I135 Xe135 Mo99 Ru103 Te132'
  yield_set = U235_thermal
  power_density = ${P}
  neutron_flux = ${phi}
  removal_groups = 'noble_gas noble_gas noble_metal noble_metal noble_metal'
  removal_rates = '${r_xe_offgas} ${r_graphite} ${r_metal} ${r_mgraph} ${r_nm_offgas}'
  removal_destinations = 'offgas processed plated processed offgas'
  verbose = true
[]

# Graphite inventories [mol per m^3 of circulating salt], fed by the graphite removal pathways of the
# action (destination 'processed' removes them from the salt; these variables keep them).
[Variables]
  [Xe135_graphite]
    type = MooseLinearVariableFVReal
    solver_sys = Xe135_graphite_sys
  []
  [Mo99_graphite]
    type = MooseLinearVariableFVReal
    solver_sys = Mo99_graphite_sys
  []
  [Ru103_graphite]
    type = MooseLinearVariableFVReal
    solver_sys = Ru103_graphite_sys
  []
  [Te132_graphite]
    type = MooseLinearVariableFVReal
    solver_sys = Te132_graphite_sys
  []
[]

[LinearFVKernels]
  [Xe135_graphite_time]
    type = LinearFVTimeDerivative
    variable = Xe135_graphite
  []
  [Xe135_graphite_uptake]
    type = FissionProductFeed
    variable = Xe135_graphite
    parent = Xe135
    transfer_rate = ${r_graphite}
  []
  [Xe135_graphite_loss] # decay + burnout in the core graphite
    type = FissionProductLoss
    variable = Xe135_graphite
    decay_constant = ${fparse log(2) / 32904}
    capture_cross_section = 2647600
    neutron_flux = ${phi_core}
  []
  [Mo99_graphite_time]
    type = LinearFVTimeDerivative
    variable = Mo99_graphite
  []
  [Mo99_graphite_uptake]
    type = FissionProductFeed
    variable = Mo99_graphite
    parent = Mo99
    transfer_rate = ${r_mgraph}
  []
  [Mo99_graphite_loss]
    type = FissionProductLoss
    variable = Mo99_graphite
    decay_constant = ${fparse log(2) / 237513.6}
  []
  [Ru103_graphite_time]
    type = LinearFVTimeDerivative
    variable = Ru103_graphite
  []
  [Ru103_graphite_uptake]
    type = FissionProductFeed
    variable = Ru103_graphite
    parent = Ru103
    transfer_rate = ${r_mgraph}
  []
  [Ru103_graphite_loss]
    type = FissionProductLoss
    variable = Ru103_graphite
    decay_constant = ${fparse log(2) / 3390941.0}
  []
  [Te132_graphite_time]
    type = LinearFVTimeDerivative
    variable = Te132_graphite
  []
  [Te132_graphite_uptake]
    type = FissionProductFeed
    variable = Te132_graphite
    parent = Te132
    transfer_rate = ${r_mgraph}
  []
  [Te132_graphite_loss]
    type = FissionProductLoss
    variable = Te132_graphite
    decay_constant = ${fparse log(2) / 276825.6}
  []
[]

# Hastelloy N corrosion: the wall of the 1-cell salt column of depth L = V / A_metal.
[CorrosionPlating]
  topology = salt_only
  reaction_boundary = left
  elements = 'Cr'
  temperature = ${T}
  reference_temperature = ${T}
  material_class = hastelloy_n
  salt_class = fluoride_fuel
  redox_class = msre_or_fuel_baseline
  position_class = core_or_reactor
  flow_factor = 2.0
  delta_T_C = 0.0
  applied_overpotential = 0.1
  applied_potential = 0.1
  default_salt_concentration = ${c0}
  reference_concentration = ${c0}
  salt_diffusivity = ${fparse k_metal * L / 2} # wall-to-cell-centre resistance = 1 / k_metal
  transient = true
[]

[Postprocessors]
  [time_d]
    type = ParsedPostprocessor
    expression = 't / ${day}'
    use_t = true
  []
  # ---- xenon --------------------------------------------------------------------------------
  [Xe_salt]
    type = ElementAverageValue
    variable = Xe135
  []
  [Xe_graphite]
    type = ElementAverageValue
    variable = Xe135_graphite
  []
  [Xe_offgas]
    type = ElementAverageValue
    variable = Xe135_offgas
  []
  [xe_poison_pct] # sigma (phi N_salt + phi_core N_graphite) / (nu F), in % dk/k
    type = ParsedPostprocessor
    expression = '100 * ${sx} * ${N_A} * (${phi} * Xe_salt + ${phi_core} * Xe_graphite) / (${nu} * ${Fd})'
    pp_names = 'Xe_salt Xe_graphite'
  []
  [xe_poison_graphite_share]
    type = ParsedPostprocessor
    expression = '${phi_core} * Xe_graphite / (${phi} * Xe_salt + ${phi_core} * Xe_graphite)'
    pp_names = 'Xe_salt Xe_graphite'
  []
  [xe_poison_analytic_pct] # steady state with I-135 in equilibrium (cumulative Xe yield 0.065385)
    type = ParsedPostprocessor
    expression = 'lX := log(2) / 32904; Ns := 0.065385 * ${Fd} / ${N_A} '
                 '/ (lX + ${sx} * ${phi} + ${r_spray} + ${r_bubble} + ${r_graphite}); '
                 'Ng := ${r_graphite} * Ns / (lX + ${sx} * ${phi_core}); '
                 '100 * ${sx} * ${N_A} * (${phi} * Ns + ${phi_core} * Ng) / (${nu} * ${Fd})'
  []
  [xe_poison_ratio]
    type = ParsedPostprocessor
    expression = 'xe_poison_pct / xe_poison_analytic_pct'
    pp_names = 'xe_poison_pct xe_poison_analytic_pct'
  []
  # ---- noble metals: inventories and fractions -----------------------------------------------
  [Mo_s]
    type = ElementAverageValue
    variable = Mo99
  []
  [Mo_m]
    type = ElementAverageValue
    variable = Mo99_plated
  []
  [Mo_g]
    type = ElementAverageValue
    variable = Mo99_graphite
  []
  [Mo_o]
    type = ElementAverageValue
    variable = Mo99_offgas
  []
  [Ru_s]
    type = ElementAverageValue
    variable = Ru103
  []
  [Ru_m]
    type = ElementAverageValue
    variable = Ru103_plated
  []
  [Ru_g]
    type = ElementAverageValue
    variable = Ru103_graphite
  []
  [Ru_o]
    type = ElementAverageValue
    variable = Ru103_offgas
  []
  [Te_s]
    type = ElementAverageValue
    variable = Te132
  []
  [Te_m]
    type = ElementAverageValue
    variable = Te132_plated
  []
  [Te_g]
    type = ElementAverageValue
    variable = Te132_graphite
  []
  [Te_o]
    type = ElementAverageValue
    variable = Te132_offgas
  []
  [Mo_salt_fraction]
    type = ParsedPostprocessor
    expression = 'Mo_s / (Mo_s + Mo_m + Mo_g + Mo_o)'
    pp_names = 'Mo_s Mo_m Mo_g Mo_o'
  []
  [Mo_metal_fraction]
    type = ParsedPostprocessor
    expression = 'Mo_m / (Mo_s + Mo_m + Mo_g + Mo_o)'
    pp_names = 'Mo_s Mo_m Mo_g Mo_o'
  []
  [Mo_graphite_fraction]
    type = ParsedPostprocessor
    expression = 'Mo_g / (Mo_s + Mo_m + Mo_g + Mo_o)'
    pp_names = 'Mo_s Mo_m Mo_g Mo_o'
  []
  [Mo_offgas_fraction]
    type = ParsedPostprocessor
    expression = 'Mo_o / (Mo_s + Mo_m + Mo_g + Mo_o)'
    pp_names = 'Mo_s Mo_m Mo_g Mo_o'
  []
  [Ru_salt_fraction]
    type = ParsedPostprocessor
    expression = 'Ru_s / (Ru_s + Ru_m + Ru_g + Ru_o)'
    pp_names = 'Ru_s Ru_m Ru_g Ru_o'
  []
  [Ru_metal_fraction]
    type = ParsedPostprocessor
    expression = 'Ru_m / (Ru_s + Ru_m + Ru_g + Ru_o)'
    pp_names = 'Ru_s Ru_m Ru_g Ru_o'
  []
  [Ru_graphite_fraction]
    type = ParsedPostprocessor
    expression = 'Ru_g / (Ru_s + Ru_m + Ru_g + Ru_o)'
    pp_names = 'Ru_s Ru_m Ru_g Ru_o'
  []
  [Ru_offgas_fraction]
    type = ParsedPostprocessor
    expression = 'Ru_o / (Ru_s + Ru_m + Ru_g + Ru_o)'
    pp_names = 'Ru_s Ru_m Ru_g Ru_o'
  []
  [Te_salt_fraction]
    type = ParsedPostprocessor
    expression = 'Te_s / (Te_s + Te_m + Te_g + Te_o)'
    pp_names = 'Te_s Te_m Te_g Te_o'
  []
  [Te_metal_fraction]
    type = ParsedPostprocessor
    expression = 'Te_m / (Te_s + Te_m + Te_g + Te_o)'
    pp_names = 'Te_s Te_m Te_g Te_o'
  []
  [Te_graphite_fraction]
    type = ParsedPostprocessor
    expression = 'Te_g / (Te_s + Te_m + Te_g + Te_o)'
    pp_names = 'Te_s Te_m Te_g Te_o'
  []
  [Te_offgas_fraction]
    type = ParsedPostprocessor
    expression = 'Te_o / (Te_s + Te_m + Te_g + Te_o)'
    pp_names = 'Te_s Te_m Te_g Te_o'
  []
  [Mo_salt_analytic] # steady salt inventory, no tracked parent (cumulative yield 0.0610873)
    type = ParsedPostprocessor
    expression = '0.0610873 * ${Fd} / ${N_A} / (log(2) / 237513.6 + 8e-28 * ${phi} '
                 '+ ${r_metal} + ${r_mgraph} + ${r_nm_offgas})'
  []
  [Mo_salt_ratio]
    type = ParsedPostprocessor
    expression = 'Mo_s / Mo_salt_analytic'
    pp_names = 'Mo_s Mo_salt_analytic'
  []
  [Ru_salt_fraction_equilibrium] # analytic: lambda / (lambda + r_m + r_g + r_o), sigma phi ~ 0
    type = ParsedPostprocessor
    expression = 'l := log(2) / 3390941.0; l / (l + ${r_metal} + ${r_mgraph} + ${r_nm_offgas})'
  []
  # ---- corrosion ----------------------------------------------------------------------------
  [corrosion_rate_um_y]
    type = SideAverageValue
    variable = corrosion_rate_um_y
    boundary = left
  []
  [recession_um]
    type = SideAverageValue
    variable = recession_um
    boundary = left
  []
  [c_Cr]
    type = ElementAverageValue
    variable = c_Cr
  []
  [salt_Cr_ppm]
    type = ParsedPostprocessor
    expression = 'c_Cr * 51.9961 / ${rho_salt}'
    pp_names = 'c_Cr'
  []
  [salt_Cr_rise_ppm] # compare with M-004: 40 -> ~100 ppm over a few weeks
    type = ParsedPostprocessor
    expression = 'salt_Cr_ppm - 40'
    pp_names = 'salt_Cr_ppm'
  []
  [cr_removal_depth_um] # uniform removal of the Cr only (x_Cr = 7 wt%), as in TM-3063 (M-003)
    type = ParsedPostprocessor
    expression = 'recession_um / ${x_Cr}'
    pp_names = 'recession_um'
  []
  [rate_vs_correlation] # initial rate / correlation rate 3.00317 um/y printed by the action
    type = ParsedPostprocessor
    expression = 'corrosion_rate_um_y / 3.00317'
    pp_names = 'corrosion_rate_um_y'
  []
  [faraday_ratio] # salt Cr gain over the Cr removed from the metal (conservation)
    type = ParsedPostprocessor
    expression = 'if(recession_um > 0, (c_Cr - ${c0}) * ${L} * 51.9961 / 8.89e6 / (recession_um * 1e-6), 1)'
    pp_names = 'c_Cr recession_um'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  l_tol = 1e-14
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-10
  end_time = ${fparse 42 * day}
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 60
    growth_factor = 1.25
  []
  dtmax = 14400
[]

[Outputs]
  csv = true
[]
