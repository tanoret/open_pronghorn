# =====================================================================================================
# Training case 05 (capstone), step 2 of 2: fuel-salt chemistry of the MSFR primary circuit.
#
# Restarts the converged flow / k-epsilon / energy of msfr_flow.i (msfr_flow_restart.e: VTB EVOL
# geometry, LiF-ThF4-UF4, 3 GWth Griffin power shape) and adds, in the same SIMPLE solve:
#   * redox: the U(IV)/U(III) buffer (2.5 mol% UF4, U(III)/U target 1 %). Radiolysis is charge
#     consistent: the fission-fragment dose (0.87 x the fission power) makes e_sol / F pairs that
#     recombine in the hot salt, so the net radiolytic yields are zero (G_e = G_F = G_F2 = 0, no F2).
#     The salt is oxidized by FISSION instead (0.8 equivalents per fission, U(III) -> U(IV)) and
#     reduced by the redox control of the fuel-processing unit (pump region, tau 1 day);
#   * [CorrosionPlatingFlow] element-resolved wall law on the Hastelloy-N-class walls: per-element
#     equilibria c_eq,k(T) (Cr 4, Fe 0.05, Ni 1e-4 mol/m^3 at 973 K, dH = 50 kJ/mol) and selective Cr
#     leaching (seed fractions 0.90 / 0.08 / 0.02), so hot walls dissolve and cold walls plate;
#   * [FissionProductInventory] I-135 / Xe-135 and Mo-99 / Ru-103 / Te-132 (U235_fast yields as a
#     stand-in for U-233), helium bubbling in the pump region (EVOL 30 s gas-extraction time), noble
#     metals plating on the HX plates (a = 600 m^2/m^3, Chilton-Colburn k_m) and the cavity walls.
# Results (8 ranks, ~6 min):
#   VERIFICATION (conservation / 0D): redox control = fission oxidation to 5e-10, U total 1200.0,
#   Cr wall net / dissolved 5e-6, Cr_vs_0d 1.000, Xe135_vs_0d 1.023, Mo99_vs_0d 1.008, Mo-99 balance 3e-11.
#   Model results (no measured MSFR chemistry exists): U(III)/U 0.96 %, reductant 10.7 mol eq/day,
#   Cr 4.02 mol/m^3 = 49 wppm between the hot / cold leg equilibria 4.65 / 2.46, 40 % of the wall area
#   dissolves (hot upper cavity wall + hot leg), Cr loss there 1.1 mg/cm^2/y, 0.52 kg/y moved hot -> cold; Ni stays at
#   c_eq; 99.94 % of the Xe-135 stripped; Mo-99 99.5 % on the HX plates, 0.5 % on the cavity walls.
# =====================================================================================================

# --- fuel salt LiF-ThF4-UF4 (VTB / Brovchenko 2019, @ ~1000 K) ---
rho = 4284            # kg/m^3
cp = 1594             # J/kg/K
drho_dT = 0.882       # kg/m^3/K
mu = 0.0166           # Pa s
k_salt = 1.7          # W/m/K
alpha_b = ${fparse drho_dT / rho}
Pr_t = 0.9

# --- VTB operating point ---
T_HX = 873.15         # K heat-exchanger secondary (ambient) temperature
h_hx = ${fparse 600 * 20e3}   # W/m^3/K volumetric HX conductance (VTB)
F_hx = 4.0e3          # kg/m^4 Forchheimer coefficient of the HX (VTB, non-standard formulation)
pump_force = -24000.  # N/m^3 pump body force (downward); VTB uses -20000 with the mixing length, tuned here
                      # to the VTB mass flow rate (~19200 kg/s) with k-epsilon wall functions

# --- k-epsilon ---
C_mu = 0.09
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
walls = 'shield_wall reflector_wall'
wall_treatment = 'eq_newton'
fluid = 'fuel pump hx'

# --- chemistry ---
D_mol = 2.0e-9        # m^2/s molecular diffusivity of the dissolved species
Sc_t = 1.0            # turbulent Schmidt number (VTB value for the precursors)
# Dose: in a fuel salt the fission fragments deposit their kinetic energy (~84 % of the recoverable
# 200 MeV) in the salt itself, plus beta (~3-4 %) and a few % of the prompt/delayed gamma; the
# radiolytic dose rate is therefore ~0.85-0.9 x the local fission power density.
dose_fraction = 0.87
eV = 1.602176634e-19
N_A = 6.02214076e23
F_salt = 8.0e4        # mol/m^3 fluoride (LiF 77.5 / ThF4 20 / UF4 2.5 mol%, ~4.8e4 mol salt/m^3 at 4284 kg/m^3)
U_total = 1200.0      # mol/m^3 uranium (2.5 mol% UF4)
x_UIII = 0.01         # U(III)/U target held by the redox control (MSRE practice ~1 %)
U_III_0 = ${fparse x_UIII * U_total}
U_IV_0 = ${fparse U_total - U_III_0}
# Fission oxidizes the salt: the fission products of one UF4 bind on average ~3.2 F, so ~0.8 F per
# fission is freed and oxidizes U(III) to U(IV) (ORNL MSRE experience; 0.8 equivalents per fission).
eq_per_fission = 0.8
E_fission = 200.0     # MeV per fission (same as the FissionProductInventory action)
fission_eq_per_J = ${fparse eq_per_fission / (E_fission * 1e6 * eV * N_A)}   # mol eq / J
# Redox control: the MSFR design assumes a continuous redox control by the fuel-processing unit
# (a salt bypass contacted with a reductant, e.g. Th or U metal, and returned at the target
# U(III)/U). Modelled in the pump region as a relaxation of U(III) and U(IV) toward the target
# composition with a loop-averaged time constant tau_redox; it conserves uranium at steady state
# and its net reduction rate must equal the fission oxidation (redox_balance_error).
tau_redox = 8.64e4    # s, loop-averaged redox-control time constant (1 day)
# Off-gas: helium bubbling with a loop-average removal time of 30 s (EVOL reference design),
# placed in the pump region (bubble injection/separation close to the pump) -> local rate
V_salt = 23.66        # m^3 salt volume of the VTB RZ model (salt_volume postprocessor)
V_pump = 2.382        # m^3 pump-region volume
tau_bubbling = 30.0   # s
r_bubble = ${fparse V_salt / V_pump / tau_bubbling}
k_ctrl = ${fparse V_salt / V_pump / tau_redox}   # 1/s local redox-control rate in the pump region
# HX plates: Chilton-Colburn k_m = h/(rho cp) (Pr/Sc)^(2/3) with the VTB HX h and area density
a_hx = 600.0          # m^2/m^3
h_wall_hx = 20e3      # W/m^2/K
Pr_mol = ${fparse mu * cp / k_salt}
Sc_mol = ${fparse mu / (rho * D_mol)}
k_m_hx = ${fparse h_wall_hx / (rho * cp) * (Pr_mol / Sc_mol)^(2/3)}
r_hx_plate = ${fparse k_m_hx * a_hx}
# Cavity and duct walls: Sh = 0.023 Re^0.8 Sc^(1/3), u ~ 1 m/s, L ~ 1 m
k_m_wall = ${fparse 0.023 * (rho * 1.0 * 1.0 / mu)^0.8 * Sc_mol^(1/3) * D_mol / 1.0}
# Element-resolved corrosion of the Ni-based structural alloy (Hastelloy N class)
T_ref = 973.15        # K, reference temperature of the equilibria (core-average salt ~957 K)
c_eq_Cr = 4.0         # mol/m^3 (~50 wppm Cr in a 4284 kg/m^3 salt at U(III)/U ~ 1 %)
c_eq_Fe = 0.05        # mol/m^3 (FeF2 much less stable than CrF2)
c_eq_Ni = 1.0e-4      # mol/m^3 (NiF2 least stable: Ni barely dissolves)
f_Cr = 0.90           # seed fractions: selective Cr leaching of the Ni alloy (Cr 7 wt%, Fe 4 wt%)
f_Fe = 0.08
f_Ni = 0.02
dH_diss = 50.0        # kJ/mol dissolution enthalpy of the equilibria (assumed, all elements)
Ea_corr = 54.83553514725915   # kJ/mol exchange-current activation energy (corrosion database)
z_M = 2               # valence of Cr, Fe, Ni (database)
Faraday = 96485.33212
# Exchange currents, as seeded by the action (verbose output): the calibrated reduced-empirical rate
# for hastelloy_n / fluoride_fuel / msre_or_fuel_baseline at T_ref (3.905 um/y of alloy recession)
# times the seed fraction, converted to a current density and divided by the fresh-salt (c = 0) anodic
# bracket exp(alpha_a z f eta_eq) = sqrt(c_eq / c_ref) (alpha_a = alpha_c = 0.5, c_ref = 1 mol/m^3)
rate_cal = 3.905489   # um/y (printed by the action as planned_rate / seed fraction)
rho_alloy = 8890      # kg/m^3 Hastelloy N (database)
i_cal = ${fparse rate_cal * 1e-6 / 3.15576e7 * rho_alloy * z_M * Faraday}   # A/m^2 per kg/mol of metal
i0_Cr = ${fparse i_cal / 51.9961e-3 * f_Cr / sqrt(c_eq_Cr)}
i0_Fe = ${fparse i_cal / 55.845e-3 * f_Fe / sqrt(c_eq_Fe)}
i0_Ni = ${fparse i_cal / 58.6934e-3 * f_Ni / sqrt(c_eq_Ni)}
R_gas = 8.314462618

[Mesh]
  [fmg]
    type = FileMeshGenerator
    # converged flow of step 1 (msfr_flow.i), with the VTB Griffin power density and hx_top/hx_bot
    file = 'msfr_flow_restart.e'
    use_for_exodus_restart = true
  []
  coord_type = RZ
  rz_coord_axis = Y
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system TKE_system TKED_system energy_system '
                     'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys U_IV_sys U_III_sys '
                     'c_Cr_sys c_Fe_sys c_Ni_sys '
                     'I135_sys Xe135_sys Mo99_sys Ru103_sys Te132_sys'
  previous_nl_solution_required = true
  allow_initial_conditions_with_restart = true  # flow fields from the file, chemistry from ICs
[]

[FVInterpolationMethods]
  [upwind]
    type = FVAdvectedUpwind
  []
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method = 'upwind'
  advected_interp_method_name = 'upwind'
[]

[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
    block = ${fluid}
  []
[]

[Functions]
  [pump_profile]
    type = ConstantFunction
    value = ${pump_force}
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_from_file_var = vel_x
    block = ${fluid}
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_from_file_var = vel_y
    block = ${fluid}
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_from_file_var = pressure
    block = ${fluid}
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    initial_from_file_var = TKE
    block = ${fluid}
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    initial_from_file_var = TKED
    block = ${fluid}
  []
  [T_fluid]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
    initial_from_file_var = T_fluid
    block = ${fluid}
  []
[]

[AuxVariables]
  # wall-law fluxes [mol/m^2/s] and the local Cr equilibrium [mol/m^3], sampled on the walls
  [J_Cr_wall]
    type = MooseVariableFVReal
    block = ${fluid}
  []
  [J_Fe_wall]
    type = MooseVariableFVReal
    block = ${fluid}
  []
  [J_Ni_wall]
    type = MooseVariableFVReal
    block = ${fluid}
  []
  [c_eq_Cr_wall]
    type = MooseVariableFVReal
    block = ${fluid}
  []
  [wall_distance]
    type = MooseVariableFVReal
    initial_condition = 0.1
    block = ${fluid}
  []
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_from_file_var = mu_t
    block = ${fluid}
  []
  [power_density]
    type = MooseVariableFVReal
    initial_from_file_var = power_density
    block = ${fluid}
  []
[]

[AuxKernels]
  [J_Cr_aux]
    type = FunctorAux
    variable = J_Cr_wall
    functor = J_Cr
    execute_on = TIMESTEP_END
    block = ${fluid}
  []
  [J_Fe_aux]
    type = FunctorAux
    variable = J_Fe_wall
    functor = J_Fe
    execute_on = TIMESTEP_END
    block = ${fluid}
  []
  [J_Ni_aux]
    type = FunctorAux
    variable = J_Ni_wall
    functor = J_Ni
    execute_on = TIMESTEP_END
    block = ${fluid}
  []
  [Cr_eq_aux]
    type = FunctorAux
    variable = c_eq_Cr_wall
    functor = c_eq_Cr_T
    execute_on = TIMESTEP_END
    block = ${fluid}
  []
  [compute_wall_distance]
    type = WallDistanceAux
    variable = wall_distance
    walls = ${walls}
    execute_on = 'INITIAL'
    block = ${fluid}
  []
  [compute_mu_t]
    type = kEpsilonViscosityAux
    variable = mu_t
    C_mu = ${C_mu}
    tke = TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    bulk_wall_treatment = false
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    mu_t_ratio_max = 1e5
    execute_on = 'NONLINEAR'
    block = ${fluid}
  []
[]

[FunctorMaterials]
  [D_eff_mat]
    type = ParsedFunctorMaterial
    property_name = D_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${D_mol} + mut/(${rho}*${Sc_t})'
    block = ${fluid}
  []
  [dose_mat]
    type = ParsedFunctorMaterial
    property_name = dose_rate
    functor_names = 'power_density'
    functor_symbols = 'q'
    expression = '${dose_fraction} * q'
    block = ${fluid}
  []
  # Removal-rate maps [1/s]
  [bubbling_mat]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = bubbling_rate
    subdomain_to_prop_value = 'fuel 0 pump ${r_bubble} hx 0'
    block = ${fluid}
  []
  [hx_plating_mat]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = hx_plating_rate
    subdomain_to_prop_value = 'fuel 0 pump 0 hx ${r_hx_plate}'
    block = ${fluid}
  []
  [ctrl_mat]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = ctrl_rate
    subdomain_to_prop_value = 'fuel 0 pump ${k_ctrl} hx 0'
    block = ${fluid}
  []
  # net U(IV) -> U(III) reduction rate of the redox control [mol/m^3/s]
  [ctrl_reduction_mat]
    type = ParsedFunctorMaterial
    property_name = redox_ctrl_reduction
    functor_names = 'ctrl_rate U_III'
    functor_symbols = 'k u3'
    expression = 'k * (${U_III_0} - u3)'
    block = ${fluid}
  []
  [ctrl_UIII_source_mat]
    type = ParsedFunctorMaterial
    property_name = ctrl_UIII_source
    functor_names = 'ctrl_rate'
    functor_symbols = 'k'
    expression = 'k * ${U_III_0}'
    block = ${fluid}
  []
  [ctrl_UIV_source_mat]
    type = ParsedFunctorMaterial
    property_name = ctrl_UIV_source
    functor_names = 'ctrl_rate'
    functor_symbols = 'k'
    expression = 'k * ${U_IV_0}'
    block = ${fluid}
  []
  [fission_oxidation_mat]
    type = ParsedFunctorMaterial
    property_name = fission_oxidation
    functor_names = 'power_density'
    functor_symbols = 'q'
    expression = '${fission_eq_per_J} * q'
    block = ${fluid}
  []
  [Xe_strip_mat]
    type = ParsedFunctorMaterial
    property_name = Xe_strip_flux
    functor_names = 'bubbling_rate Xe135'
    functor_symbols = 'r N'
    expression = 'r * N'
    block = ${fluid}
  []
  [Mo_hx_plating_mat]
    type = ParsedFunctorMaterial
    property_name = Mo_hx_plating
    functor_names = 'hx_plating_rate Mo99'
    functor_symbols = 'rp N'
    expression = 'rp * N'
    block = ${fluid}
  []
  [Mo_wall_flux_mat]
    type = ParsedFunctorMaterial
    property_name = Mo_wall_flux
    functor_names = 'Mo99'
    functor_symbols = 'N'
    expression = '${k_m_wall} * N'
    block = ${fluid}
  []
  # Wall corrosion flux into the salt [mol/m^2/s] (+ dissolution, - plating): the BC law evaluated
  # with the wall-cell temperature and concentration,
  #   J_k = i0_k(T) / (z F sqrt(c_eq,k(T))) (c_eq,k(T) - c_k),  i0_k(T) = i0_k exp[Ea/R (1/T_ref - 1/T)]
  [J_Cr_mat]
    type = ParsedFunctorMaterial
    property_name = J_Cr
    functor_names = 'T_fluid c_Cr'
    functor_symbols = 'T c'
    expression = 'ceq := ${c_eq_Cr} * exp(-${fparse dH_diss * 1e3} / ${R_gas} * (1 / T - 1 / ${T_ref}));
                  ${i0_Cr} * exp(${fparse Ea_corr * 1e3} / ${R_gas} * (1 / ${T_ref} - 1 / T)) / (${z_M} * ${Faraday} * sqrt(ceq)) * (ceq - c)'
    block = ${fluid}
  []
  [J_Fe_mat]
    type = ParsedFunctorMaterial
    property_name = J_Fe
    functor_names = 'T_fluid c_Fe'
    functor_symbols = 'T c'
    expression = 'ceq := ${c_eq_Fe} * exp(-${fparse dH_diss * 1e3} / ${R_gas} * (1 / T - 1 / ${T_ref}));
                  ${i0_Fe} * exp(${fparse Ea_corr * 1e3} / ${R_gas} * (1 / ${T_ref} - 1 / T)) / (${z_M} * ${Faraday} * sqrt(ceq)) * (ceq - c)'
    block = ${fluid}
  []
  [J_Ni_mat]
    type = ParsedFunctorMaterial
    property_name = J_Ni
    functor_names = 'T_fluid c_Ni'
    functor_symbols = 'T c'
    expression = 'ceq := ${c_eq_Ni} * exp(-${fparse dH_diss * 1e3} / ${R_gas} * (1 / T - 1 / ${T_ref}));
                  ${i0_Ni} * exp(${fparse Ea_corr * 1e3} / ${R_gas} * (1 / ${T_ref} - 1 / T)) / (${z_M} * ${Faraday} * sqrt(ceq)) * (ceq - c)'
    block = ${fluid}
  []
  # 0D (well-mixed loop) wall balance: sum_walls w (c_eq - c) = 0 with w = i0(T) / sqrt(c_eq(T))
  [w_Cr_mat]
    type = ParsedFunctorMaterial
    property_name = w_Cr
    functor_names = 'T_fluid'
    functor_symbols = 'T'
    expression = 'ceq := ${c_eq_Cr} * exp(-${fparse dH_diss * 1e3} / ${R_gas} * (1 / T - 1 / ${T_ref}));
                  exp(${fparse Ea_corr * 1e3} / ${R_gas} * (1 / ${T_ref} - 1 / T)) / sqrt(ceq)'
    block = ${fluid}
  []
  [w_ceq_Cr_mat]
    type = ParsedFunctorMaterial
    property_name = w_ceq_Cr
    functor_names = 'T_fluid'
    functor_symbols = 'T'
    expression = 'ceq := ${c_eq_Cr} * exp(-${fparse dH_diss * 1e3} / ${R_gas} * (1 / T - 1 / ${T_ref}));
                  exp(${fparse Ea_corr * 1e3} / ${R_gas} * (1 / ${T_ref} - 1 / T)) * sqrt(ceq)'
    block = ${fluid}
  []
  # Cr equilibrium at the local temperature (for the plots)
  [Cr_eq_mat]
    type = ParsedFunctorMaterial
    property_name = c_eq_Cr_T
    functor_names = 'T_fluid'
    functor_symbols = 'T'
    expression = '${c_eq_Cr} * exp(-${fparse dH_diss * 1e3} / ${R_gas} * (1 / T - 1 / ${T_ref}))'
    block = ${fluid}
  []
  # Cr wall flux split into the dissolving (J > 0) and plating (J < 0) walls, and the HX duct walls
  [J_Cr_pos_mat]
    type = ParsedFunctorMaterial
    property_name = J_Cr_pos
    functor_names = 'J_Cr'
    functor_symbols = 'J'
    expression = 'max(J, 0)'
    block = ${fluid}
  []
  [J_Cr_neg_mat]
    type = ParsedFunctorMaterial
    property_name = J_Cr_neg
    functor_names = 'J_Cr'
    functor_symbols = 'J'
    expression = 'min(J, 0)'
    block = ${fluid}
  []
  [dissolving_mask_mat]
    type = ParsedFunctorMaterial
    property_name = dissolving_mask
    functor_names = 'J_Cr'
    functor_symbols = 'J'
    expression = 'if(J > 0, 1, 0)'
    block = ${fluid}
  []
  [J_Cr_hx_mat]
    type = ParsedFunctorMaterial
    property_name = J_Cr_hx
    functor_names = 'J_Cr'
    functor_symbols = 'J'
    expression = 'if(x > 1.93 & y < 0.03 & y > -0.68, J, 0)'
    block = ${fluid}
  []
  [k_eff_mat]
    type = ParsedFunctorMaterial
    property_name = k_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${k_salt} + mut*${cp}/${Pr_t}'
    block = ${fluid}
  []
  [mu_eff_mat]
    type = ParsedFunctorMaterial
    property_name = mu_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${mu} + mut'
    block = ${fluid}
  []
  [alpha_mat]
    type = GenericFunctorMaterial
    prop_names = 'alpha_b'
    prop_values = '${alpha_b}'
    block = ${fluid}
  []
  # Lagged speed for the linearised Forchheimer friction F |u| u of the heat exchanger
  [speed_mat]
    type = ParsedFunctorMaterial
    property_name = speed_lag
    functor_names = 'vel_x vel_y'
    functor_symbols = 'u v'
    expression = 'sqrt(u*u + v*v)'
    block = 'hx'
  []
  [forch_mat]
    type = GenericVectorFunctorMaterial
    prop_names = 'forch'
    prop_values = '${F_hx} ${F_hx} ${F_hx}'
    block = 'hx'
  []
  [hx_removal_mat]
    type = ParsedFunctorMaterial
    property_name = hx_removal
    functor_names = 'T_fluid'
    functor_symbols = 'T'
    expression = '${h_hx} * (T - ${T_HX})'
    block = 'hx'
  []
[]


[MoltenSaltRadiolysis]
  block = ${fluid}
  salt_type = fluoride
  metals = 'U'                        # U(IV)/U(III) redox buffer of the fuel salt
  temperature = T_fluid
  dose_rate = dose_rate
  # Charge-consistent radiolysis: every F atom is born with an electron. In the hot salt the pair
  # recombines (e_sol + F -> F-, or e_sol reduces U(IV) and the F atom / F2 re-oxidizes U(III)) within
  # ~1e-10 s, so the NET radiolytic oxidation is zero and no F2 survives (MSRE: F2 only from frozen
  # salt below ~100 C). The full e_sol / F set is too stiff for the segregated Picard iterations, so
  # its zero net effect is imposed by switching all three yields off (G_e = G_F = G_F2 = 0). The salt
  # oxidation that does exist is the fission oxidation (fission_oxidation kernels below).
  g_value_species = 'e_sol F_rad F2_diss'
  g_value_overrides = '0 0 0'
  rhie_chow_user_object = 'rc'
  diffusivity = D_eff
  time_derivative = false
  initial_condition_species = 'F_ion U_IV U_III'
  initial_condition_values = '${F_salt} ${U_IV_0} ${U_III_0}'
  verbose = true
[]

[CorrosionPlatingFlow]
  block = ${fluid}
  elements = 'Cr Fe Ni'
  reaction_boundary = ${walls}
  temperature = T_fluid
  reference_temperature = ${T_ref}
  material_class = hastelloy_n         # Ni-based structural alloy (Hastelloy N class, as in the VTB)
  salt_class = fluoride_fuel
  redox_class = msre_or_fuel_baseline  # U(III)/U ~ 1 % fuel salt
  reference_concentration = 1.0
  temperature_dependent_kinetics = true
  # Element-resolved thermodynamic wall law: each element relaxes toward its own equilibrium
  # c_eq,k(T) = c_eq,k exp[-dH/R (1/T - 1/T_ref)], so walls hotter than the loop equilibrium dissolve
  # and colder walls plate. The calibrated rate is split by selective Cr leaching (seed_fractions).
  equilibrium_concentrations = '${c_eq_Cr} ${c_eq_Fe} ${c_eq_Ni}'
  dissolution_enthalpies = '${dH_diss}'
  seed_fractions = '${f_Cr} ${f_Fe} ${f_Ni}'
  rhie_chow_user_object = 'rc'
  diffusivity = D_eff
  time_derivative = false
  initial_concentration = ${c_eq_Cr}    # start at the reference Cr equilibrium c_eq,Cr(T_ref) (see README)
  verbose = true
[]

[FissionProductInventory]
  block = ${fluid}
  nuclides = 'I135 Xe135 Mo99 Ru103 Te132'
  yield_set = U235_fast                # stand-in: no U-233 set in the database (MSFR start-up on 233U)
  power_density = power_density        # the Griffin fission power of the energy equation
  removal_groups = 'noble_gas noble_metal'
  removal_rates = 'bubbling_rate hx_plating_rate'
  removal_destinations = 'offgas plated'
  track_removed_inventories = false
  plating_boundaries = ${walls}
  plating_mass_transfer_coefficient = ${k_m_wall}
  rhie_chow_user_object = 'rc'
  diffusivity = D_eff
  time_derivative = false
  verbose = true
[]

[LinearFVKernels]
  # --- redox: fission oxidation (all fuel) and redox control (pump region) ---
  # Fission oxidation: 0.8 equivalents per fission convert U(III) to U(IV)
  [fission_oxidation_UIII]
    type = LinearFVSource
    variable = U_III
    source_density = power_density
    scaling_factor = ${fparse -fission_eq_per_J}
    block = ${fluid}
  []
  [fission_oxidation_UIV]
    type = LinearFVSource
    variable = U_IV
    source_density = power_density
    scaling_factor = ${fission_eq_per_J}
    block = ${fluid}
  []
  # Redox control: the processing unit returns the pump-region salt at the target redox state,
  #   dU_III/dt = k (U_III_0 - U_III),  dU_IV/dt = k (U_IV_0 - U_IV)  (pump region only)
  [control_UIII_sink]
    type = LinearFVReaction
    variable = U_III
    coeff = ctrl_rate
    block = ${fluid}
  []
  [control_UIII_feed]
    type = LinearFVSource
    variable = U_III
    source_density = ctrl_UIII_source
    block = ${fluid}
  []
  [control_UIV_sink]
    type = LinearFVReaction
    variable = U_IV
    coeff = ctrl_rate
    block = ${fluid}
  []
  [control_UIV_feed]
    type = LinearFVSource
    variable = U_IV
    source_density = ctrl_UIV_source
    block = ${fluid}
  []
  # --- momentum ---
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'x'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
    block = ${fluid}
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [u_rz]
    type = LinearFVRZViscousSource
    variable = vel_x
    mu = 'mu_eff'
    momentum_component = 'x'
    block = ${fluid}
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
    block = ${fluid}
  []
  [u_friction]
    type = LinearFVMomentumFriction
    variable = vel_x
    mu = 'speed_lag'
    Darcy_name = 'forch'
    momentum_component = 'x'
    block = 'hx'
  []
  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'y'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
    block = ${fluid}
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
    block = ${fluid}
  []
  [v_friction]
    type = LinearFVMomentumFriction
    variable = vel_y
    mu = 'speed_lag'
    Darcy_name = 'forch'
    momentum_component = 'y'
    block = 'hx'
  []
  [v_buoyancy]
    type = LinearFVMomentumBoussinesq
    variable = vel_y
    T_fluid = T_fluid
    gravity = '0 -9.81 0'
    rho = ${rho}
    alpha_name = 'alpha_b'
    ref_temperature = ${T_HX}
    momentum_component = 'y'
    block = ${fluid}
  []
  [pump]
    type = LinearFVSource
    variable = vel_y
    source_density = pump_profile
    block = 'pump'
  []
  [p_diffusion]
    type = LinearFVPressureCorrectionDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [HbyA_divergence]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = true
    block = ${fluid}
  []

  # --- k-epsilon ---
  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
    block = ${fluid}
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [TKE_source_sink]
    type = LinearFVTKESourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    block = ${fluid}
  []
  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = ${walls}
    block = ${fluid}
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    walls = ${walls}
    block = ${fluid}
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    walls = ${walls}
    block = ${fluid}
  []
  [TKED_source_sink]
    type = LinearFVTKEDSourceSink
    variable = TKED
    u = vel_x
    v = vel_y
    tke = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    block = ${fluid}
  []

  # --- energy ---
  [energy_advection]
    type = LinearFVEnergyAdvection
    variable = T_fluid
    advected_quantity = temperature
    cp = ${cp}
    block = ${fluid}
  []
  [energy_conduction]
    type = LinearFVDiffusion
    variable = T_fluid
    diffusion_coeff = 'k_eff'
    use_nonorthogonal_correction = false
    block = ${fluid}
  []
  [fission_heat]
    type = LinearFVSource
    variable = T_fluid
    source_density = power_density
    block = ${fluid}
  []
  [hx_sink]
    type = LinearFVReaction
    variable = T_fluid
    coeff = ${h_hx}
    block = 'hx'
  []
  [hx_secondary]
    type = LinearFVSource
    variable = T_fluid
    source_density = ${fparse h_hx * T_HX}
    block = 'hx'
  []
[]

[LinearFVBCs]
  [walls_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_x
    functor = 0.0
  []
  [walls_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_y
    functor = 0.0
  []
  [axis_u]
    type = LinearFVVelocitySymmetryBC
    boundary = fluid_symmetry
    variable = vel_x
    u = vel_x
    v = vel_y
    momentum_component = x
  []
  [axis_v]
    type = LinearFVVelocitySymmetryBC
    boundary = fluid_symmetry
    variable = vel_y
    u = vel_x
    v = vel_y
    momentum_component = y
  []
  [axis_p]
    type = LinearFVPressureSymmetryBC
    boundary = fluid_symmetry
    variable = pressure
    HbyA_flux = 'HbyA'
  []
  [walls_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = ${walls}
    variable = mu_t
    u = vel_x
    v = vel_y
    rho = ${rho}
    mu = ${mu}
    tke = TKE
    wall_treatment = ${wall_treatment}
  []
[]

[Postprocessors]
  [max_v]
    type = ElementExtremeValue
    variable = vel_x
    value_type = max
    block = ${fluid}
  []
  [max_vy]
    type = ElementExtremeValue
    variable = vel_y
    value_type = max
    block = ${fluid}
  []
  [min_vy]
    type = ElementExtremeValue
    variable = vel_y
    value_type = min
    block = ${fluid}
  []
  [T_max]
    type = ElementExtremeValue
    variable = T_fluid
    block = ${fluid}
  []
  [T_min]
    type = ElementExtremeValue
    variable = T_fluid
    value_type = min
    block = ${fluid}
  []
  [T_avg]
    type = ElementAverageValue
    variable = T_fluid
    block = ${fluid}
  []
  [p_max]
    type = ElementExtremeValue
    variable = pressure
    block = ${fluid}
  []
  [p_min]
    type = ElementExtremeValue
    variable = pressure
    value_type = min
    block = ${fluid}
  []
  [mu_t_max]
    type = ElementExtremeValue
    variable = mu_t
    block = ${fluid}
  []
  [flow_hx_top]
    type = VolumetricFlowRate
    boundary = 'hx_top'
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = 1
    rhie_chow_user_object = rc
    execute_on = 'TIMESTEP_END'
  []
  [mdot]
    type = ParsedPostprocessor
    expression = '-${rho} * flow_hx_top'
    pp_names = 'flow_hx_top'
  []
  [Tflow_hx_top]
    type = VolumetricFlowRate
    boundary = 'hx_top'
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = T_fluid
    rhie_chow_user_object = rc
    execute_on = 'TIMESTEP_END'
    outputs = none
  []
  [flow_hx_bot]
    type = VolumetricFlowRate
    boundary = 'hx_bot'
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = 1
    rhie_chow_user_object = rc
    execute_on = 'TIMESTEP_END'
    outputs = none
  []
  [Tflow_hx_bot]
    type = VolumetricFlowRate
    boundary = 'hx_bot'
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = T_fluid
    rhie_chow_user_object = rc
    execute_on = 'TIMESTEP_END'
    outputs = none
  []
  [T_hx_in]
    type = ParsedPostprocessor
    expression = 'Tflow_hx_top / flow_hx_top'
    pp_names = 'Tflow_hx_top flow_hx_top'
  []
  [T_hx_out]
    type = ParsedPostprocessor
    expression = 'Tflow_hx_bot / flow_hx_bot'
    pp_names = 'Tflow_hx_bot flow_hx_bot'
  []
  [fission_power]
    type = ElementIntegralVariablePostprocessor
    variable = power_density
    block = ${fluid}
  []
  [hx_removal_power]
    type = ElementIntegralFunctorPostprocessor
    functor = hx_removal
    block = 'hx'
  []
  [energy_balance_error]
    type = ParsedPostprocessor
    pp_names = 'fission_power hx_removal_power'
    expression = 'abs(fission_power - hx_removal_power) / fission_power'
  []
  [salt_volume]
    type = VolumePostprocessor
    block = ${fluid}
  []
  # --- bulk (mass-flow-weighted) temperatures: hot leg = HX inlet, cold leg = HX outlet ---
  [T_bulk_hot]
    type = ParsedPostprocessor
    expression = 'T_hx_in'
    pp_names = 'T_hx_in'
  []
  [T_bulk_cold]
    type = ParsedPostprocessor
    expression = 'T_hx_out'
    pp_names = 'T_hx_out'
  []
  [loop_transit_time]
    type = ParsedPostprocessor
    expression = '${rho} * salt_volume / abs(mdot)'
    pp_names = 'salt_volume mdot'
  []
  # --- redox: fission oxidation balanced by the redox control ---
  [U_III_avg]
    type = ElementAverageValue
    variable = U_III
    block = ${fluid}
  []
  [U_IV_avg]
    type = ElementAverageValue
    variable = U_IV
    block = ${fluid}
  []
  [U_total_avg]
    type = ParsedPostprocessor
    expression = 'U_III_avg + U_IV_avg'
    pp_names = 'U_III_avg U_IV_avg'
  []
  [U_III_fraction]
    type = ParsedPostprocessor
    expression = 'U_III_avg / (U_III_avg + U_IV_avg)'
    pp_names = 'U_III_avg U_IV_avg'
  []
  [U_III_min]
    type = ElementExtremeValue
    variable = U_III
    value_type = min
    block = ${fluid}
  []
  [U_III_max]
    type = ElementExtremeValue
    variable = U_III
    block = ${fluid}
  []
  [fission_oxidation_rate]
    # mol eq / s for the whole salt (0.8 eq per fission)
    type = ElementIntegralFunctorPostprocessor
    functor = fission_oxidation
    block = ${fluid}
  []
  [redox_control_rate]
    # mol / s of U(IV) reduced back to U(III) by the redox control in the pump region
    type = ElementIntegralFunctorPostprocessor
    functor = redox_ctrl_reduction
    block = ${fluid}
  []
  [redox_balance_error]
    type = ParsedPostprocessor
    expression = 'abs(redox_control_rate - fission_oxidation_rate) / fission_oxidation_rate'
    pp_names = 'redox_control_rate fission_oxidation_rate'
  []
  [reductant_per_day]
    # mol/day of reducing equivalents the fuel-processing unit must supply (e.g. /4 for Th metal)
    type = ParsedPostprocessor
    expression = 'redox_control_rate * 86400'
    pp_names = 'redox_control_rate'
  []
  [dose_power]
    # radiolytic dose deposited in the salt [W] (0.87 x the fission power)
    type = ElementIntegralFunctorPostprocessor
    functor = dose_rate
    block = ${fluid}
  []
  [gross_radiolysis_pairs]
    # e_sol / F pairs made per second at G = 0.2 /100 eV [mol/s]; they recombine (net G = 0)
    type = ParsedPostprocessor
    expression = 'dose_power * 0.2 / (100 * ${eV} * ${N_A})'
    pp_names = 'dose_power'
  []
  [surviving_fraction_equivalent]
    # fraction of the gross pairs that, if it escaped recombination as F2, would equal the fission oxidation
    type = ParsedPostprocessor
    expression = 'fission_oxidation_rate / gross_radiolysis_pairs'
    pp_names = 'fission_oxidation_rate gross_radiolysis_pairs'
  []
  [F2_max]
    type = ElementExtremeValue
    variable = F2_diss
    block = ${fluid}
  []
  # --- corrosion: element-resolved wall exchange ---
  [Cr_avg]
    type = ElementAverageValue
    variable = c_Cr
    block = ${fluid}
  []
  [Fe_avg]
    type = ElementAverageValue
    variable = c_Fe
    block = ${fluid}
  []
  [Ni_avg]
    type = ElementAverageValue
    variable = c_Ni
    block = ${fluid}
  []
  [Cr_max]
    type = ElementExtremeValue
    variable = c_Cr
    block = ${fluid}
  []
  [Cr_min]
    type = ElementExtremeValue
    variable = c_Cr
    value_type = min
    block = ${fluid}
  []
  [Cr_wppm]
    type = ParsedPostprocessor
    expression = 'Cr_avg * 51.9961e-3 / ${rho} * 1e6'
    pp_names = 'Cr_avg'
  []
  [Cr_eq_hot]
    # Cr equilibrium at the bulk hot-leg temperature
    type = ParsedPostprocessor
    expression = '${c_eq_Cr} * exp(-${fparse dH_diss * 1e3} / ${R_gas} * (1 / T_hx_in - 1 / ${T_ref}))'
    pp_names = 'T_hx_in'
  []
  [Cr_eq_cold]
    type = ParsedPostprocessor
    expression = '${c_eq_Cr} * exp(-${fparse dH_diss * 1e3} / ${R_gas} * (1 / T_hx_out - 1 / ${T_ref}))'
    pp_names = 'T_hx_out'
  []
  [Cr_wall_net]
    # net Cr release of all walls [mol/s]; ~0 at the closed-loop steady state
    type = SideIntegralFunctorPostprocessor
    functor = J_Cr
    boundary = ${walls}
  []
  [Cr_dissolved]
    # Cr released by the dissolving (hot) walls [mol/s]
    type = SideIntegralFunctorPostprocessor
    functor = J_Cr_pos
    boundary = ${walls}
  []
  [Cr_plated]
    # Cr deposited on the plating (cold) walls [mol/s] (< 0)
    type = SideIntegralFunctorPostprocessor
    functor = J_Cr_neg
    boundary = ${walls}
  []
  [Cr_plated_hx]
    # part of it deposited in the HX ducts [mol/s]
    type = SideIntegralFunctorPostprocessor
    functor = J_Cr_hx
    boundary = ${walls}
  []
  [Cr_wall_balance]
    # VERIFICATION: closed loop, steady state -> what dissolves hot plates cold
    type = ParsedPostprocessor
    expression = 'abs(Cr_wall_net) / Cr_dissolved'
    pp_names = 'Cr_wall_net Cr_dissolved'
  []
  [dissolving_area]
    type = SideIntegralFunctorPostprocessor
    functor = dissolving_mask
    boundary = ${walls}
  []
  [dissolving_area_fraction]
    type = ParsedPostprocessor
    expression = 'dissolving_area / wall_area'
    pp_names = 'dissolving_area wall_area'
  []
  [Cr_loss_mg_cm2_y]
    # area-averaged Cr loss of the dissolving walls [mg/cm^2/y]
    type = ParsedPostprocessor
    expression = 'Cr_dissolved / dissolving_area * 51.9961 * 3.15576e7 * 0.1'
    pp_names = 'Cr_dissolved dissolving_area'
  []
  [Cr_redistribution_kg_y]
    # Cr transferred from the hot to the cold walls [kg/y]
    type = ParsedPostprocessor
    expression = 'Cr_dissolved * 51.9961e-3 * 3.15576e7'
    pp_names = 'Cr_dissolved'
  []
  [w_Cr_int]
    type = SideIntegralFunctorPostprocessor
    functor = w_Cr
    boundary = ${walls}
    outputs = none
  []
  [w_ceq_Cr_int]
    type = SideIntegralFunctorPostprocessor
    functor = w_ceq_Cr
    boundary = ${walls}
    outputs = none
  []
  [Cr_0d]
    # 0D well-mixed-loop Cr level at which the wall exchange balances (VERIFICATION target)
    type = ParsedPostprocessor
    expression = 'w_ceq_Cr_int / w_Cr_int'
    pp_names = 'w_ceq_Cr_int w_Cr_int'
  []
  [Cr_vs_0d]
    type = ParsedPostprocessor
    expression = 'Cr_avg / Cr_0d'
    pp_names = 'Cr_avg Cr_0d'
  []
  [Fe_wall_net]
    type = SideIntegralFunctorPostprocessor
    functor = J_Fe
    boundary = ${walls}
  []
  [Ni_wall_net]
    type = SideIntegralFunctorPostprocessor
    functor = J_Ni
    boundary = ${walls}
  []
  # --- fission products ---
  [bubbling_integral]
    type = ElementIntegralFunctorPostprocessor
    functor = bubbling_rate
    block = ${fluid}
    outputs = none
  []
  [hx_plating_integral]
    type = ElementIntegralFunctorPostprocessor
    functor = hx_plating_rate
    block = ${fluid}
    outputs = none
  []
  [wall_area]
    type = AreaPostprocessor
    boundary = ${walls}
  []
  [Xe135_avg]
    type = ElementAverageValue
    variable = Xe135
    block = ${fluid}
  []
  [I135_avg]
    type = ElementAverageValue
    variable = I135
    block = ${fluid}
  []
  [Mo99_avg]
    type = ElementAverageValue
    variable = Mo99
    block = ${fluid}
  []
  [Ru103_avg]
    type = ElementAverageValue
    variable = Ru103
    block = ${fluid}
  []
  [Te132_avg]
    type = ElementAverageValue
    variable = Te132
    block = ${fluid}
  []
  [Xe135_max]
    type = ElementExtremeValue
    variable = Xe135
    block = ${fluid}
  []
  [Xe135_min]
    type = ElementExtremeValue
    variable = Xe135
    value_type = min
    block = ${fluid}
  []
  # production [mol/s]: cumulative mass-135 and mass-99 U235_fast yields, 200 MeV/fission
  [P135]
    type = ParsedPostprocessor
    expression = 'fission_power * 0.0660108 * 5.182134828e-14'
    pp_names = 'fission_power'
    outputs = none
  []
  [P_Mo99]
    type = ParsedPostprocessor
    expression = 'fission_power * 0.0594313 * 5.182134828e-14'
    pp_names = 'fission_power'
    outputs = none
  []
  # 0D (well-mixed loop) equilibria with the loop-integrated removal rates (VERIFICATION)
  [Xe135_vs_0d]
    type = ParsedPostprocessor
    expression = 'Xe135_avg / (P135 / salt_volume / (2.106574218e-05 + bubbling_integral / salt_volume))'
    pp_names = 'Xe135_avg P135 salt_volume bubbling_integral'
  []
  [Mo99_vs_0d]
    type = ParsedPostprocessor
    expression = 'Mo99_avg / (P_Mo99 / salt_volume / (2.918347331e-06 + (hx_plating_integral + ${k_m_wall} * wall_area) / salt_volume))'
    pp_names = 'Mo99_avg P_Mo99 salt_volume hx_plating_integral wall_area'
  []
  [Xe135_stripped]
    type = ElementIntegralFunctorPostprocessor
    functor = Xe_strip_flux
    block = ${fluid}
    outputs = none
  []
  [Xe135_stripped_fraction]
    type = ParsedPostprocessor
    expression = 'Xe135_stripped / P135'
    pp_names = 'Xe135_stripped P135'
  []
  [Mo99_hx_plated]
    type = ElementIntegralFunctorPostprocessor
    functor = Mo_hx_plating
    block = ${fluid}
    outputs = none
  []
  [Mo99_wall_plated]
    type = SideIntegralFunctorPostprocessor
    functor = Mo_wall_flux
    boundary = ${walls}
    outputs = none
  []
  [Mo99_hx_plated_fraction]
    type = ParsedPostprocessor
    expression = 'Mo99_hx_plated / P_Mo99'
    pp_names = 'Mo99_hx_plated P_Mo99'
  []
  [Mo99_wall_plated_fraction]
    type = ParsedPostprocessor
    expression = 'Mo99_wall_plated / P_Mo99'
    pp_names = 'Mo99_wall_plated P_Mo99'
  []
  [Mo99_decayed_fraction]
    type = ParsedPostprocessor
    expression = '2.918347331e-06 * Mo99_avg * salt_volume / P_Mo99'
    pp_names = 'Mo99_avg salt_volume P_Mo99'
  []
  [Mo99_balance_error]
    type = ParsedPostprocessor
    expression = 'abs(1 - Mo99_hx_plated_fraction - Mo99_wall_plated_fraction - Mo99_decayed_fraction)'
    pp_names = 'Mo99_hx_plated_fraction Mo99_wall_plated_fraction Mo99_decayed_fraction'
  []
[]

[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  energy_system = 'energy_system'
  passive_scalar_systems = 'e_sol_sys F_ion_sys F_rad_sys F2_diss_sys U_IV_sys U_III_sys '
                           'c_Cr_sys c_Fe_sys c_Ni_sys '
                           'I135_sys Xe135_sys Mo99_sys Ru103_sys Te132_sys'
  passive_scalar_equation_relaxation = '1 1 1 1 1 1 1 1 1 1 1 1 1 1'
  passive_scalar_absolute_tolerance = '2e-9 2e-9 2e-9 2e-9 2e-9 2e-9 2e-9 2e-9 2e-9 2e-9 2e-9 2e-9 2e-9 2e-9'
  passive_scalar_l_tol = 1e-12
  passive_scalar_l_abs_tol = 1e-14
  passive_scalar_l_max_its = 500
  # direct solve: the closed loop with a weak wall exchange is nearly singular for the dissolved metals
  passive_scalar_petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  passive_scalar_petsc_options_value = 'lu mumps'
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.3 0.3'
  turbulence_field_relaxation = '0.3 0.3'
  energy_equation_relaxation = 0.9
  # The flow restarted from msfr_flow.i is converged: loose flow tolerances stop the SIMPLE loop after
  # its first (pressure-corrected) iteration, and num_iterations then caps the passive-scalar loop.
  num_iterations = 3000
  pressure_absolute_tolerance = 5e-2
  momentum_absolute_tolerance = 5e-3
  turbulence_absolute_tolerance = '5e-2 5e-2'
  energy_absolute_tolerance = 5e-3
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'
  turbulence_l_max_its = 30
  pin_pressure = true
  pressure_pin_value = 1e5
  pressure_pin_point = '0.5 0.0 0.0'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'
  print_fields = false
  continue_on_max_its = true
[]

[VectorPostprocessors]
  # Wall-cell values on the inner (shield) and outer (reflector) walls; plot.py orders them by arc length
  [shield_profile]
    type = SideValueSampler
    variable = 'T_fluid c_Cr c_eq_Cr_wall J_Cr_wall J_Fe_wall J_Ni_wall Mo99'
    boundary = shield_wall
    sort_by = id
    execute_on = TIMESTEP_END
  []
  [reflector_profile]
    type = SideValueSampler
    variable = 'T_fluid c_Cr c_eq_Cr_wall J_Cr_wall J_Fe_wall J_Ni_wall Mo99'
    boundary = reflector_wall
    sort_by = id
    execute_on = TIMESTEP_END
  []
  # External loop (x = 2.13 m), top to bottom: hot leg, pump (bubbling), HX, cold leg
  [loop]
    type = LineValueSampler
    warn_discontinuous_face_values = false  # FV cell values; a point may sit on a face
    variable = 'T_fluid vel_y Xe135 I135 Mo99 U_III c_Cr'
    start_point = '2.13 1.0 0'
    end_point = '2.13 -1.0 0'
    num_points = 101
    sort_by = y
    execute_on = TIMESTEP_END
  []
  # Core mid-plane (y = 0), axis to shield wall
  [midplane]
    type = LineValueSampler
    warn_discontinuous_face_values = false  # FV cell values; a point may sit on a face
    variable = 'T_fluid vel_x vel_y Xe135 Mo99'
    start_point = '0.01 0 0'
    end_point = '1.05 0 0'
    num_points = 53
    sort_by = x
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  csv = true
[]
