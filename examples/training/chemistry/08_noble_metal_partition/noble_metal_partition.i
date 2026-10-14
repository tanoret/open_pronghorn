# =====================================================================================================
# Training case 08: where do the noble metals go? Plating, off-gas and decay on the walls.
#
# The noble and semi-noble fission products (Mo, Tc, Ru, Rh, Te, Nb) do not form stable fluorides in
# a reducing fuel salt. They exist as metal atoms/colloids that leave the salt quickly: they plate on
# metal surfaces (heat exchanger, piping) and graphite, and are entrained into the off-gas with the
# sparging bubbles. At the MSRE only a few percent of the Mo-99, Ru-103 and Te-132 produced were found
# in salt samples; most were on the surfaces or in the off-gas (Compere et al., ORNL-4865, 1975).
#
# The removal is modelled as two first-order pathways from the salt, plated (r_p) and off-gas (r_o):
#
#   dN/dt = Y F / N_A + feeds - (lambda + r_p + r_o) N,   d N_plated / dt = r_p N - lambda N_plated + ...
#
# and the removed inventories keep decaying. This matters for the source term and the decay heat:
#   * Mo-99 (66 h) -> Tc-99m (6 h) -> Tc-99: plated Mo-99 produces plated Tc-99m;
#   * Te-132 (3.2 d) -> I-132 (2.3 h): plated tellurium produces *iodine* on the heat-exchanger walls,
#     although iodine itself is salt-seeking (halogen group, not removed).
#
# Analytic check: at steady state each salt inventory without a tracked parent is
# N = Y F / (N_A (lambda + r_p + r_o)), and the removed fraction splits as r_p : r_o.
# =====================================================================================================

P0 = 2.0e7 # W/m^3
r_plate = 5.0e-4 # 1/s: k_m A / V, e.g. 1e-5 m/s on 50 m^2 of heat-exchanger per m^3 of salt
r_offgas = 1.0e-4 # 1/s: entrainment with the sparge gas
day = 86400

[Mesh]
  [cell]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  linear_sys_names = 'Mo99_sys Tc99m_sys Tc99_sys Ru103_sys Ru106_sys Te132_sys I132_sys '
                     'Mo99_offgas_sys Tc99m_offgas_sys Tc99_offgas_sys Ru103_offgas_sys '
                     'Ru106_offgas_sys Te132_offgas_sys I132_offgas_sys '
                     'Mo99_plated_sys Tc99m_plated_sys Tc99_plated_sys Ru103_plated_sys '
                     'Ru106_plated_sys Te132_plated_sys I132_plated_sys'
[]

[FissionProductInventory]
  nuclides = 'Mo99 Tc99m Tc99 Ru103 Ru106 Te132 I132'
  power_density = ${P0}
  # Two pathways for the noble-metal group (Mo, Tc, Ru, Te); iodine (halogen) stays in the salt but
  # is born in the removed inventories from the decay of removed Te-132
  removal_groups = 'noble_metal noble_metal'
  removal_rates = '${r_plate} ${r_offgas}'
  removal_destinations = 'plated offgas'
  # Print the plan: note I132_plated and I132_offgas, created because I-132 is a daughter of Te-132
  verbose = true
[]

[Postprocessors]
  [time_d]
    type = ParsedPostprocessor
    expression = 't / ${day}'
    use_t = true
  []
  # ---- Mo-99 ----
  [Mo99_salt]
    type = ElementAverageValue
    variable = Mo99
  []
  [Mo99_pl]
    type = ElementAverageValue
    variable = Mo99_plated
  []
  [Mo99_og]
    type = ElementAverageValue
    variable = Mo99_offgas
  []
  [Mo99_salt_fraction]
    type = ParsedPostprocessor
    expression = 'Mo99_salt / (Mo99_salt + Mo99_pl + Mo99_og)'
    pp_names = 'Mo99_salt Mo99_pl Mo99_og'
  []
  [Mo99_plated_fraction]
    type = ParsedPostprocessor
    expression = 'Mo99_pl / (Mo99_salt + Mo99_pl + Mo99_og)'
    pp_names = 'Mo99_salt Mo99_pl Mo99_og'
  []
  [Mo99_salt_analytic]
    type = ParsedPostprocessor
    expression = '0.0610873 * ${P0} / (200e6 * 1.602176634e-19) / 6.02214076e23 '
                 '/ (log(2) / 237513.6 + ${r_plate} + ${r_offgas})'
  []
  # ---- Ru-103 and Ru-106 ----
  [Ru103_salt]
    type = ElementAverageValue
    variable = Ru103
  []
  [Ru103_pl]
    type = ElementAverageValue
    variable = Ru103_plated
  []
  [Ru103_og]
    type = ElementAverageValue
    variable = Ru103_offgas
  []
  [Ru103_salt_fraction]
    type = ParsedPostprocessor
    expression = 'Ru103_salt / (Ru103_salt + Ru103_pl + Ru103_og)'
    pp_names = 'Ru103_salt Ru103_pl Ru103_og'
  []
  [Ru106_pl]
    type = ElementAverageValue
    variable = Ru106_plated
  []
  # ---- Te-132 / I-132 ----
  [Te132_salt]
    type = ElementAverageValue
    variable = Te132
  []
  [Te132_pl]
    type = ElementAverageValue
    variable = Te132_plated
  []
  [I132_salt]
    type = ElementAverageValue
    variable = I132
  []
  [I132_pl]
    type = ElementAverageValue
    variable = I132_plated
  []
  [Te132_salt_fraction]
    type = ParsedPostprocessor
    expression = 'Te132_salt / (Te132_salt + Te132_pl + Te132_og)'
    pp_names = 'Te132_salt Te132_pl Te132_og'
  []
  [Te132_og]
    type = ElementAverageValue
    variable = Te132_offgas
  []
  # Fraction of the I-132 activity located on the walls
  [I132_plated_fraction]
    type = ParsedPostprocessor
    expression = 'I132_pl / (I132_salt + I132_pl + I132_og)'
    pp_names = 'I132_salt I132_pl I132_og'
  []
  [I132_og]
    type = ElementAverageValue
    variable = I132_offgas
  []
  # Plated activity per unit salt volume [Bq/m^3] -- the heat-exchanger source term
  [plated_activity]
    type = ParsedPostprocessor
    expression = '6.02214076e23 * log(2) * (Mo99_pl / 237513.6 + Tc99m_pl / 21624.12 '
                 '+ Ru103_pl / 3390941.0 + Ru106_pl / 32123520.0 + Te132_pl / 276825.6 '
                 '+ I132_pl / 8262.0)'
    pp_names = 'Mo99_pl Tc99m_pl Ru103_pl Ru106_pl Te132_pl I132_pl'
  []
  [Tc99m_pl]
    type = ElementAverageValue
    variable = Tc99m_plated
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  dt = 3600
  end_time = ${fparse 30 * day}
  l_tol = 1e-14
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
