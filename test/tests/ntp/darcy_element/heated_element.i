# One Ludewig Design 1 element (12 x 5 cone mesh) in the Darcy-Forchheimer form, with the gas
# heated by 1 MW/m^3 in the fuel bed and no particle phase. The steady outlet mass flow must equal
# the imposed inlet flow, and the outlet enthalpy flow the inlet one plus the heat input, so the
# gas leaves at T with h(T) = h(300 K) + Q / mdot (300.083 K). Checks NTPDarcyAdvection,
# NTPPackedBedFunctorMaterial, NTPDuctFlowFunctorMaterial and NTPDarcyFluxIntegral together.
L = 0.7294
r1_in = 0.00801
r1_out = 0.02129
r2_in = 0.01101
r2_out = 0.02429
r3_in = 0.01976
r3_out = 0.02975
r4_in = 0.02176
r4_out = 0.03175
r5_in = 0.02376
r5_out = 0.03375
n_elements = 37
eps_cf = 0.32
d_cf = 8.0e-5
k_cf = 190.0
eps_bed = 0.40
d_bed = 5.0e-4
k_bed = 22.0
eps_hf = 0.28
d_hf = 1.2e-4
k_hf = 33.0
f_ann = 0.02
nz = 12
n_ch = 1
n_hf = 1
n_bed = 1
n_cf = 1
n_pl = 1
mdot = ${fparse 21.054294673192363 / n_elements}
p_out = 7.0e6
[Mesh]
  [strip]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1 1 1 1 1'
    ix = '${n_ch} ${n_hf} ${n_bed} ${n_cf} ${n_pl}'
    dy = '${L}'
    iy = '${nz}'
    subdomain_id = '0 1 2 3 4'
  []
  # Boundaries are cut on the reference strip, where the region limits are fixed:
  # inlet = plenum at z = 0, outlet = chamber at z = L, closed_ends = the rest of both ends.
  [inlet]
    type = ParsedGenerateSideset
    input = strip
    combinatorial_geometry = 'y < 1e-9 & x > 4'
    new_sideset_name = inlet
  []
  [outlet]
    type = ParsedGenerateSideset
    input = inlet
    combinatorial_geometry = 'y > ${L} - 1e-9 & x < 1'
    new_sideset_name = outlet
  []
  [closed_ends]
    type = ParsedGenerateSideset
    input = outlet
    combinatorial_geometry = '(y < 1e-9 & x < 4) | (y > ${L} - 1e-9 & x > 1)'
    new_sideset_name = closed_ends
  []
  [cone]
    type = ParsedNodeTransformGenerator
    input = closed_ends
    constant_names = 'L r1a r1b r2a r2b r3a r3b r4a r4b r5a r5b'
    constant_expressions = '${L} ${r1_in} ${r1_out} ${r2_in} ${r2_out} ${r3_in} ${r3_out}
                            ${r4_in} ${r4_out} ${r5_in} ${r5_out}'
    # f = y/L is the axial fraction; within region i the node keeps its fraction x - i.
    x_function = 'f := y / L;
                  r1 := r1a + (r1b - r1a) * f; r2 := r2a + (r2b - r2a) * f;
                  r3 := r3a + (r3b - r3a) * f; r4 := r4a + (r4b - r4a) * f;
                  r5 := r5a + (r5b - r5a) * f;
                  if(x < 1, x * r1,
                  if(x < 2, r1 + (x - 1) * (r2 - r1),
                  if(x < 3, r2 + (x - 2) * (r3 - r2),
                  if(x < 4, r3 + (x - 3) * (r4 - r3),
                            r4 + (x - 4) * (r5 - r4)))))'
  []
  [names]
    type = RenameBlockGenerator
    input = cone
    old_block = '0 1 2 3 4'
    new_block = 'chamber hot_frit fuel_bed cold_frit plenum'
  []
  [axis]
    type = RenameBoundaryGenerator
    input = names
    old_boundary = 'left right'
    new_boundary = 'axis outer_wall'
  []
  coord_type = RZ
  rz_coord_axis = Y
[]
[FluidProperties]
  [h2]
    type = NTPHydrogenFluidProperties
  []
[]
area_in = ${fparse pi * (r5_in^2 - r4_in^2)}

[Variables]
  [pressure]
    type = MooseVariableFVReal
  []
[]
[ICs]
  # plenum at the 300 K inlet pressure, falling linearly to p_out in the chamber
  [p_guess]
    type = FunctionIC
    variable = pressure
    function = '${p_out}'
  []
[]
[Variables]
  [T_fluid]
    type = MooseVariableFVReal
    initial_condition = 300
  []
[]
[FVKernels]
  [mass_flux]
    type = NTPDarcyMassFlux
    variable = pressure
    pressure = pressure
    rho_mobility = rho_mobility
  []
[]
[FVKernels]
  [gas_advection]
    type = NTPDarcyAdvection
    variable = T_fluid
    pressure = pressure
    rho_mobility = rho_mobility
    advected_quantity = dh
    inflow_boundaries = inlet
    boundary_value = 0
  []
  [gas_conduction]
    type = FVDiffusion
    variable = T_fluid
    coeff = k
  []
  # a uniform 1 MW/m^3 source in the bed only, gas side (no solid): checks the gas energy alone
  [source]
    type = FVBodyForce
    variable = T_fluid
    value = 1e6
    block = fuel_bed
  []
[]
[FVBCs]
  [T_in]
    type = FVDirichletBC
    variable = T_fluid
    boundary = inlet
    value = 300
  []
  [mass_in]
    type = FVFunctorNeumannBC
    variable = pressure
    boundary = inlet
    functor = ${fparse mdot / area_in}
  []
  [p_out]
    type = FVDirichletBC
    variable = pressure
    boundary = outlet
    value = ${p_out}
  []
[]
[FunctorMaterials]
  [gas]
    type = GeneralFunctorFluidProps
    fp = h2
    pressure = pressure
    T_fluid = T_fluid
    speed = superficial_speed
    porosity = 1
    characteristic_length = ${d_bed}
  []
  [Dh]
    type = ADParsedFunctorMaterial
    property_name = Dh
    expression = 'f := y / ${L};
                  r1 := ${r1_in} + (${r1_out} - ${r1_in}) * f;
                  r4 := ${r4_in} + (${r4_out} - ${r4_in}) * f;
                  r5 := ${r5_in} + (${r5_out} - ${r5_in}) * f;
                  if(x < r1 + 1e-9, 2 * r1, 2 * (r5 - r4))'
    block = 'chamber plenum'
  []
  [channels]
    type = NTPDuctFlowFunctorMaterial
    block = 'chamber plenum'
    hydraulic_diameter = Dh
    friction_factor = ${f_ann}
  []
  [cold_frit]
    type = NTPPackedBedFunctorMaterial
    block = cold_frit
    porosity = ${eps_cf}
    particle_diameter = ${d_cf}
    solid_conductivity = ${k_cf}
  []
  [fuel_bed]
    type = NTPPackedBedFunctorMaterial
    block = fuel_bed
    porosity = ${eps_bed}
    particle_diameter = ${d_bed}
    solid_conductivity = ${k_bed}
  []
  [hot_frit]
    type = NTPPackedBedFunctorMaterial
    block = hot_frit
    porosity = ${eps_hf}
    particle_diameter = ${d_hf}
    solid_conductivity = ${k_hf}
  []
[]
[FunctorMaterials/dh]
  type = ADParsedFunctorMaterial
  property_name = dh
  functor_names = 'h'
  expression = 'h - (26470.08016246304)'
[]
[FunctorMaterials/enth]
  type = INSFVEnthalpyFunctorMaterial
  fp = h2
  rho = rho
  pressure = pressure
  temperature = T_fluid
  assumed_constant_cp = false
[]
[Postprocessors]
  [H_out]
    type = NTPDarcyFluxIntegral
    variable = pressure
    boundary = outlet
    rho_mobility = rho_mobility
    carried_quantity = dh
  []
  [Q]
    type = ElementIntegralFunctorPostprocessor
    functor = 1e6
    block = fuel_bed
  []
  [energy_error]
    type = ParsedPostprocessor
    pp_names = 'H_out Q'
    expression = 'H_out / Q - 1'
  []
  [T_out]
    type = SideAverageValue
    variable = T_fluid
    boundary = outlet
  []
  [mdot_out]
    type = NTPDarcyFluxIntegral
    variable = pressure
    boundary = outlet
    rho_mobility = rho_mobility
  []
  [p_in]
    type = SideAverageValue
    variable = pressure
    boundary = inlet
  []
[]
[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_max_its = 40
[]
[Outputs]
  csv = true
  exodus = true
[]
