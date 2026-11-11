# Steady thermal-hydraulics of one particle-bed fuel element, axisymmetric (r, z).
#
# Flow: Darcy-Forchheimer, u = -M grad p (u superficial), in all five regions. In the frits and
# the bed M follows from the KTA packed-bed correlation; in the plenum and chamber from duct
# friction along the axis. Mass:  div(rho u) = -div(rho M grad p) = 0.
# Energy: local thermal non-equilibrium between the gas (T_fluid) and the particles (T_solid)
# of the cold frit, fuel bed and hot frit; the gas enthalpy is carried by the same Darcy fluxes.
#
# Include after common/element_mesh.i and common/hydrogen_nasa.i. The including input sets:
#   mdot                  mass flow through the mesh [kg/s]
#   element_fraction      fraction of the element the mesh covers (1 in RZ, 0.5 for the 3-D half)
#   T_in                  inlet temperature [K]
#   p_out                 chamber outlet pressure [Pa]
#   f_ann                 Darcy friction factor of the plenum and chamber [-]
#   eps_cf d_cf k_cf      cold frit: porosity, particle diameter [m], particle conductivity [W/(m K)]
#   eps_bed d_bed k_bed   fuel bed
#   eps_hf d_hf k_hf      hot frit
#   kz_bed kr_bed         fuel-bed axial / radial conductivity multipliers
#   power_density         function: fuel-bed power density q'''(r, z, t) [W/m^3]
#   cold_frit_multiplier  functor or number multiplying the cold-frit drag, chi_cf(z) [-]
#   T_solid_init          initial temperature of the particles [K]
#   h_in                  specific enthalpy of the gas at T_in [J/kg]
#
# The steady state is reached by pseudo-transient continuation. The time terms only damp the
# Newton iterations; the particle heat capacity is set for speed, not realism.

# The mesh may cover a fraction of the element (half of it in 3-D, with symmetry planes); mdot is
# the flow through that fraction. The inlet mass flux uses the area of the inlet of the mesh (in 3-D
# the faceted annulus is smaller than the circular one), so that mdot enters exactly.

# energy closures of the packed beds (THEORY.md 4.3): particle emissivity of the radiative
# conductivity (0: off) and the coefficient of the gas thermal dispersion (0: off)
bed_emissivity = 0
bed_dispersion = 0

[Variables]
  [pressure]
    type = MooseVariableFVReal
    initial_condition = ${p_out}
  []
  [T_fluid]
    type = MooseVariableFVReal
    initial_condition = ${T_in}
  []
  [T_solid]
    type = MooseVariableFVReal
    block = 'cold_frit fuel_bed hot_frit'
    initial_condition = ${T_solid_init}
  []
[]

[FVKernels]
  # ---- gas mass
  [mass_time]
    type = NTPGasMassTimeDerivative
    variable = pressure
    fp = h2
    T_fluid = T_fluid
    porosity = porosity
  []
  [mass_flux]
    type = NTPDarcyMassFlux
    variable = pressure
    pressure = pressure
    rho_mobility = rho_mobility
  []

  # ---- gas energy
  [gas_time]
    type = PINSFVEnergyTimeDerivative
    variable = T_fluid
    rho = rho
    cp = cp
    porosity = porosity
    is_solid = false
  []
  [gas_advection]
    type = NTPDarcyAdvection
    variable = T_fluid
    pressure = pressure
    rho_mobility = rho_mobility
    # the enthalpy rise h - h_in; the convective form makes the reference level irrelevant
    advected_quantity = dh
    inflow_boundaries = inlet
    boundary_value = 0
  []
  [gas_conduction]
    type = FVDiffusion
    variable = T_fluid
    coeff = eps_k
  []
  [gas_from_solid]
    type = PINSFVEnergyAmbientConvection
    variable = T_fluid
    T_fluid = T_fluid
    T_solid = T_solid
    is_solid = false
    h_solid_fluid = alpha
    block = 'cold_frit fuel_bed hot_frit'
  []

  # ---- particles
  [solid_time]
    type = PINSFVEnergyTimeDerivative
    variable = T_solid
    rho = rho_s
    cp = cp_s
    porosity = porosity
    is_solid = true
  []
  [solid_conduction]
    type = FVAnisotropicDiffusion
    variable = T_solid
    coeff = kappa_s
  []
  [fission_heat]
    type = FVBodyForce
    variable = T_solid
    function = ${power_density}
    block = fuel_bed
  []
  [solid_to_gas]
    type = PINSFVEnergyAmbientConvection
    variable = T_solid
    T_fluid = T_fluid
    T_solid = T_solid
    is_solid = true
    h_solid_fluid = alpha
  []
[]

[FVBCs]
  # an inflow mass flux: -rho M grad p . n = -mdot / A at the inlet
  [mass_in]
    type = FVFunctorNeumannBC
    variable = pressure
    boundary = inlet
    functor = inlet_mass_flux
  []
  [p_out]
    type = FVDirichletBC
    variable = pressure
    boundary = outlet
    value = ${p_out}
  []
  [T_in]
    type = FVDirichletBC
    variable = T_fluid
    boundary = inlet
    value = ${T_in}
  []
[]

# Cell superficial mass flux reconstructed from the face fluxes, for the heat-transfer correlations:
# an AD functor of the solution, so that its dependence enters the Newton Jacobian (an auxiliary
# variable would lag it by one residual evaluation and degrade Newton next to the flow excursion).
# G_r and G_z copy it for the output.
[AuxVariables]
  [G_r]
    type = MooseVariableFVReal
  []
  [G_z]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [G_r]
    type = FunctorAux
    variable = G_r
    functor = Gx_cell
    execute_on = 'initial timestep_end'
  []
  [G_z]
    type = FunctorAux
    variable = G_z
    functor = Gy_cell
    execute_on = 'initial timestep_end'
  []
[]

[Postprocessors]
  [inlet_area]
    type = AreaPostprocessor
    boundary = inlet
    execute_on = initial
    outputs = none
  []
[]

[FunctorMaterials]
  [inlet_flux]
    type = ADParsedFunctorMaterial
    property_name = inlet_mass_flux
    functor_names = inlet_area
    expression = '${mdot} / inlet_area'
  []
  [mass_flux]
    type = NTPDarcyMassFluxFunctorMaterial
    pressure = pressure
    rho_mobility = rho_mobility
    mass_flux_x_name = Gx_cell
    mass_flux_y_name = Gy_cell
    mass_flux_z_name = Gz_cell
    inflow_boundaries = inlet
    inflow_mass_fluxes = inlet_mass_flux
  []
  [porosity]
    type = ADPiecewiseByBlockFunctorMaterial
    prop_name = porosity
    subdomain_to_prop_value = 'chamber 1 plenum 1 cold_frit ${eps_cf} fuel_bed ${eps_bed}
                               hot_frit ${eps_hf}'
  []
  [gas]
    type = GeneralFunctorFluidProps
    fp = h2
    pressure = pressure
    T_fluid = T_fluid
    speed = superficial_speed
    porosity = porosity
    characteristic_length = ${d_bed}
  []
  [enthalpy]
    type = INSFVEnthalpyFunctorMaterial
    fp = h2
    rho = rho
    pressure = pressure
    temperature = T_fluid
    assumed_constant_cp = false
  []
  [enthalpy_rise]
    type = ADParsedFunctorMaterial
    property_name = dh
    functor_names = 'h'
    expression = 'h - ${h_in}'
  []
  # gas thermal dispersion of the packed beds (0 in the free channels)
  [k_disp]
    type = ADPiecewiseByBlockFunctorMaterial
    prop_name = k_disp
    subdomain_to_prop_value = 'chamber 0 plenum 0 cold_frit k_disp_cold_frit fuel_bed k_disp_fuel_bed
                               hot_frit k_disp_hot_frit'
  []
  [gas_conduction]
    type = ADParsedFunctorMaterial
    property_name = eps_k
    functor_names = 'porosity k k_disp'
    expression = 'porosity * k + k_disp'
  []

  # ---- free channels: plenum and chamber
  [Dh]
    type = ADParsedFunctorMaterial
    property_name = Dh
    expression = 'f := y / ${L};
                  r1 := ${r1_in} + (${r1_out} - ${r1_in}) * f;
                  r4 := ${r4_in} + (${r4_out} - ${r4_in}) * f;
                  r5 := ${r5_in} + (${r5_out} - ${r5_in}) * f;
                  r := sqrt(x * x + z * z);
                  if(r < r1 + 1e-9, 2 * r1, 2 * (r5 - r4))'
    block = 'chamber plenum'
  []
  [channels]
    type = NTPDuctFlowFunctorMaterial
    block = 'chamber plenum'
    hydraulic_diameter = Dh
    friction_factor = ${f_ann}
  []

  # ---- packed beds
  [cold_frit]
    type = NTPPackedBedFunctorMaterial
    block = cold_frit
    mass_flux_x = Gx_cell
    mass_flux_y = Gy_cell
    mass_flux_z = Gz_cell
    porosity = ${eps_cf}
    particle_diameter = ${d_cf}
    pressure_drop_model = kta
    heat_transfer_model = kta
    orifice_multiplier = ${cold_frit_multiplier}
    solid_conductivity = ${k_cf}
    kappa_name = kappa_s
    particle_emissivity = ${bed_emissivity}
    dispersion_coefficient = ${bed_dispersion}
    dispersion_name = k_disp_cold_frit
  []
  [fuel_bed]
    type = NTPPackedBedFunctorMaterial
    block = fuel_bed
    mass_flux_x = Gx_cell
    mass_flux_y = Gy_cell
    mass_flux_z = Gz_cell
    porosity = ${eps_bed}
    particle_diameter = ${d_bed}
    pressure_drop_model = kta
    heat_transfer_model = achenbach
    solid_conductivity = ${k_bed}
    axial_conductivity_multiplier = ${kz_bed}
    radial_conductivity_multiplier = ${kr_bed}
    kappa_name = kappa_s
    particle_emissivity = ${bed_emissivity}
    dispersion_coefficient = ${bed_dispersion}
    dispersion_name = k_disp_fuel_bed
  []
  [hot_frit]
    type = NTPPackedBedFunctorMaterial
    block = hot_frit
    mass_flux_x = Gx_cell
    mass_flux_y = Gy_cell
    mass_flux_z = Gz_cell
    porosity = ${eps_hf}
    particle_diameter = ${d_hf}
    pressure_drop_model = kta
    heat_transfer_model = kta
    solid_conductivity = ${k_hf}
    kappa_name = kappa_s
    particle_emissivity = ${bed_emissivity}
    dispersion_coefficient = ${bed_dispersion}
    dispersion_name = k_disp_hot_frit
  []
  # pseudo heat capacity of the particles: only shapes the path to the steady state
  [solid_inertia]
    type = ADGenericFunctorMaterial
    prop_names = 'rho_s cp_s'
    prop_values = '1000 50'
  []
[]
