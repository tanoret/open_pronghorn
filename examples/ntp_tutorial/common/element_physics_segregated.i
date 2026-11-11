# Thermal-hydraulics of one particle-bed fuel element with the segregated SIMPLE (steady) and
# PIMPLE (pseudo-transient) solvers of the linear finite-volume Navier-Stokes module. It solves
# the closures of the Newton Darcy-Forchheimer model (element_physics.i) with the full momentum
# balance of the free channels (THEORY.md, Section 6.2):
#
#  - momentum is solved for the superficial velocity in all five regions. The packed-bed drag of
#    the frits and the bed is a Darcy friction mu D u, D = chi (A + B|u|) / mu
#    (NTPPackedBedSegregatedFunctorMaterial, KTA), the cell speed |u| from the face mass fluxes.
#    The plenum and the chamber keep inertia and viscous stress, plus the friction of a duct of
#    their hydraulic diameter along the axis (NTPDuctFrictionSegregatedFunctorMaterial);
#  - the pressure-velocity coupling across the porosity jumps uses the series face coupling and
#    pressure gradient of NTPRhieChowMassFlux;
#  - the gas energy is solved for the specific enthalpy (exact for the variable cp of hydrogen),
#    the gas-to-particle exchange with partial elimination (NTPLinearFVEnthalpyExchange); the
#    particle temperature T_solid lives in its own linear system on the porous regions.
#
# Lagged quantities. The segregated pressure equation has no compressibility term, and in the
# packed beds the momentum diagonal is ten orders of magnitude larger than in the free channels:
# a density that follows the pressure correction amplifies it without bound, and one that follows
# the gas temperature (a factor of 100 between inlet and outlet) makes the iterations cycle. The
# flow properties are evaluated at a lagged state (p_rho, T_rho), and the turbulent duct friction at
# a lagged axial speed; both move towards the solved state by 'lagged_state_relaxation' at the end
# of every solve: once per fixed-point iteration with SIMPLE, once per time step with PIMPLE. The
# converged state is the one of the fully coupled equations.
#
# Include after the mesh and common/hydrogen_nasa.i, then include one of the dimension files
# element_segregated_rz.i or element_segregated_3d.i after this one. The including input sets the
# dimension lists first (the dimension files document them):
#   segregated_momentum_systems, segregated_velocities, segregated_inlet_velocity,
#   segregated_inlet_fraction, segregated_walls, segregated_wall_types,
#   segregated_wall_energy_types, segregated_wall_energy_functors, segregated_initial_velocity.
# Parameters the including input must also set (SI units):
#   mdot, T_in, p_out, f_ann, eps_cf d_cf k_cf, eps_bed d_bed k_bed, eps_hf d_hf k_hf, kz_bed, kr_bed
#   power_density          function of the fuel-bed power density q'''(x, t) [W/m^3]
#   cold_frit_multiplier   functor or number multiplying the cold-frit drag, chi_cf(z) [-]
#   T_solid_init           initial particle temperature [K]

area_in = ${fparse pi * (r5_in^2 - r4_in^2) * segregated_inlet_fraction}
lagged_state_relaxation = 0.3
# energy closures of the packed beds (THEORY.md 4.3): particle emissivity of the radiative
# conductivity (0: off) and the coefficient of the gas thermal dispersion (0: off)
bed_emissivity = 0
bed_dispersion = 0
# wall friction of the free channels as a body force (one cell across them); a model that resolves
# the channels (channel_turbulence.i) switches it off and sets f_ann = 0
channel_wall_friction = true
# uniform axial velocity added to the plug-flow start [m/s]: every cell then has an advection
# diagonal in the first momentum solve (needed on coarse azimuthal 3-D meshes)
segregated_start_floor = 5

[Problem]
  linear_sys_names = '${segregated_momentum_systems} pressure_system energy_system
                      solid_energy_system'
[]

[UserObjects]
  # Rhie-Chow face fluxes with the series coupling of the half cells. Declared under the Physics'
  # default name and before [Physics]: the objects of a task are created in input order, and the
  # Physics does not add a Rhie-Chow object when one exists on its blocks.
  [ins_rhie_chow_interpolator]
    type = NTPRhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = rho
    p_diffusion_kernel = flow_p_diffusion
  []
[]

[Physics]
  [NavierStokes]
    [FlowSegregated/flow]
      compressibility = weakly-compressible
      velocity_variable = '${segregated_velocities}'
      density = rho
      dynamic_viscosity = mu
      # '0 0 0' (at rest) or, along the axis y, the plug-flow estimate segregated_plug_velocity
      initial_velocity = '${segregated_initial_velocity}'
      initial_pressure = ${p_out}

      # The inlet velocity carries mdot at the face density the Rhie-Chow object uses, so the inlet
      # mass flow is mdot exactly.
      inlet_boundaries = inlet
      momentum_inlet_types = fixed-velocity
      momentum_inlet_functors = '${segregated_inlet_velocity}'
      outlet_boundaries = outlet
      momentum_outlet_types = fixed-pressure
      pressure_functors = ${p_out}
      wall_boundaries = '${segregated_walls}'
      momentum_wall_types = '${segregated_wall_types}'

      friction_blocks = 'cold_frit; fuel_bed; hot_frit; plenum chamber'
      friction_types = 'darcy; darcy; darcy; darcy'
      friction_coeffs = 'drag_cf; drag_bed; drag_hf; drag_duct'

      momentum_advection_interpolation = upwind
    []
    [FluidHeatTransferSegregated/gas]
      coupled_flow_physics = flow
      solve_for_enthalpy = true
      fp = h2
      fluid_temperature_variable = T_fluid
      thermal_conductivity = eps_k
      specific_heat = cp
      initial_temperature = ${T_in}
      energy_inlet_types = fixed-temperature
      energy_inlet_functors = ${T_in}
      energy_wall_types = '${segregated_wall_energy_types}'
      energy_wall_functors = '${segregated_wall_energy_functors}'
      energy_advection_interpolation = upwind
    []
  []
[]

[Functions]
  # plug-flow estimate of the cold axial velocity, a start for SIMPLE: the inlet speed in the
  # plenum falling linearly to zero at the closed end, the chamber speed rising linearly from zero
  # to the speed of the cold gas at the outlet, plus a uniform floor segregated_start_floor in
  # every cell. Started from rest, the first momentum solve sees only the gas viscosity in the
  # free channels, and in 3-D the iterations diverge; a cell left at zero (the beds, without the
  # floor) does the same on coarse azimuthal meshes (n_theta <= 5).
  [segregated_plug_velocity]
    type = ParsedFunction
    # cold gas (ideal gas at p_out, T_in) through the inlet annulus and the outlet of the chamber
    expression = 'f := y / ${L}; r := sqrt(x * x + z * z);
                  r1 := ${r1_in} + (${r1_out} - ${r1_in}) * f;
                  r4 := ${r4_in} + (${r4_out} - ${r4_in}) * f;
                  rho := ${p_out} * 2.01588e-3 / (8.3144598 * ${T_in});
                  u_in := ${mdot} / (pi * (${r5_in}^2 - ${r4_in}^2) * ${segregated_inlet_fraction} * rho);
                  u_out := ${mdot} / (pi * ${r1_out}^2 * ${segregated_inlet_fraction} * rho);
                  ${segregated_start_floor} + if(r > r4, u_in * (1 - f), if(r < r1, u_out * f, 0))'
  []
[]

[Variables]
  [T_solid]
    type = MooseLinearVariableFVReal
    solver_sys = solid_energy_system
    block = 'cold_frit fuel_bed hot_frit'
    initial_condition = ${T_solid_init}
  []
[]

[FVInterpolationMethods]
  # face conductivity of the particle conduction (THEORY.md 4.3)
  [particle_conductivity_interpolation]
    type = FVHarmonicAverage
  []
[]

[LinearFVKernels]
  # gas side of the particle-to-gas exchange, consistent with the enthalpy unknown; the particle
  # response of the cell is eliminated (solid_conductivity), which makes the gas-particle
  # iterations converge in tens instead of thousands of iterations
  [gas_from_particles]
    type = NTPLinearFVEnthalpyExchange
    variable = h
    alpha = alpha
    T_solid = T_solid
    T_fluid = T_fluid
    cp = cp
    solid_conductivity = kappa_s
    block = 'cold_frit fuel_bed hot_frit'
  []
  # particles. The face conductivity between two cells is their harmonic mean, as in the Newton
  # model (FVAnisotropicDiffusion). Evaluated directly on the face, the conductivity between the
  # cold frit and the bed (129 and 13 W/(m K)) would be the frit's, and the bed would heat the
  # cold frit's gas to twice the temperature of the Newton model.
  [particle_conduction]
    type = LinearFVAnisotropicDiffusion
    variable = T_solid
    diffusion_tensor = kappa_s
    coeff_interp_method = particle_conductivity_interpolation
  []
  [fission_heat]
    type = LinearFVSource
    variable = T_solid
    source_density = ${power_density}
    block = fuel_bed
  []
  [particles_to_gas]
    type = LinearFVVolumetricHeatTransfer
    variable = T_solid
    h_solid_fluid = alpha
    T_fluid = T_fluid
    T_solid = T_solid
    is_solid = true
  []
[]

[AuxVariables]
  # superficial mass flux reconstructed from the face fluxes, for the bed closures
  [G_x]
    type = MooseVariableFVReal
  []
  [G_y]
    type = MooseVariableFVReal
  []
  [G_z]
    type = MooseVariableFVReal
  []
  # lagged state of the flow properties and of the duct friction
  [p_rho]
    type = MooseVariableFVReal
    initial_condition = ${p_out}
  []
  [T_rho]
    type = MooseVariableFVReal
    initial_condition = ${T_in}
  []
  [speed_duct]
    type = MooseVariableFVReal
    block = 'plenum chamber'
    initial_condition = 0
  []
[]

[AuxKernels]
  [G_x]
    type = NTPRhieChowMassFluxAux
    variable = G_x
    rhie_chow_user_object = ins_rhie_chow_interpolator
    component = x
    execute_on = 'nonlinear timestep_end'
  []
  [G_y]
    type = NTPRhieChowMassFluxAux
    variable = G_y
    rhie_chow_user_object = ins_rhie_chow_interpolator
    component = y
    execute_on = 'nonlinear timestep_end'
  []
  [G_z]
    type = NTPRhieChowMassFluxAux
    variable = G_z
    rhie_chow_user_object = ins_rhie_chow_interpolator
    component = z
    execute_on = 'nonlinear timestep_end'
  []
  [p_rho]
    type = ParsedAux
    variable = p_rho
    coupled_variables = 'p_rho pressure'
    expression = 'p_rho + ${lagged_state_relaxation} * (pressure - p_rho)'
    execute_on = 'timestep_end'
  []
  [T_rho]
    type = ParsedAux
    variable = T_rho
    coupled_variables = T_rho
    functor_names = T_fluid
    functor_symbols = T
    expression = 'T_rho + ${lagged_state_relaxation} * (T - T_rho)'
    execute_on = 'timestep_end'
  []
  # axial speed of the duct friction (the element axis is y)
  [speed_duct]
    type = ParsedAux
    variable = speed_duct
    coupled_variables = 'speed_duct vel_y'
    expression = 'speed_duct + 0.5 * (abs(vel_y) - speed_duct)'
    execute_on = 'timestep_end'
    block = 'plenum chamber'
  []
[]

[FunctorMaterials]
  # flow properties at the lagged state
  [gas_flow]
    type = NTPGasPropertiesFunctorMaterial
    fp = h2
    pressure = p_rho
    temperature = T_rho
    specific_heat_name = cp_lagged
    conductivity_name = k_lagged
  []
  # energy properties at the current gas temperature
  [gas_energy]
    type = NTPGasPropertiesFunctorMaterial
    fp = h2
    pressure = p_rho
    temperature = T_fluid
    density_name = rho_current
    viscosity_name = mu_current
    speed_of_sound_name = speed_of_sound_current
  []
  [inlet_velocity]
    type = ParsedFunctorMaterial
    property_name = u_inlet
    functor_names = rho
    expression = '${mdot} / (${area_in} * rho)'
  []
  # gas conduction eps k, as in the Newton model
  [porosity]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = porosity
    subdomain_to_prop_value = 'chamber 1 plenum 1 cold_frit ${eps_cf} fuel_bed ${eps_bed}
                               hot_frit ${eps_hf}'
  []
  # gas thermal dispersion of the packed beds (0 in the free channels)
  [k_disp]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = k_disp
    subdomain_to_prop_value = 'chamber 0 plenum 0 cold_frit k_disp_cold_frit fuel_bed k_disp_fuel_bed
                               hot_frit k_disp_hot_frit'
  []
  [eps_k]
    type = ParsedFunctorMaterial
    property_name = eps_k
    functor_names = 'porosity k k_disp'
    expression = 'porosity * k + k_disp'
  []
  # ---- free channels: duct friction along the axis
  [Dh]
    type = ParsedFunctorMaterial
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
    type = NTPDuctFrictionSegregatedFunctorMaterial
    block = 'chamber plenum'
    hydraulic_diameter = Dh
    friction_factor = ${f_ann}
    laminar_term = ${channel_wall_friction}
    speed = speed_duct
    drag_name = drag_duct
  []
  # ---- packed beds
  [cold_frit]
    type = NTPPackedBedSegregatedFunctorMaterial
    block = cold_frit
    porosity = ${eps_cf}
    particle_diameter = ${d_cf}
    pressure_drop_model = kta
    heat_transfer_model = kta
    orifice_multiplier = ${cold_frit_multiplier}
    solid_conductivity = ${k_cf}
    mass_flux_x = G_x
    mass_flux_y = G_y
    mass_flux_z = G_z
    drag_name = drag_cf
    alpha_name = alpha_cf
    kappa_name = kappa_cf
    speed_name = speed_cf
    particle_emissivity = ${bed_emissivity}
    dispersion_coefficient = ${bed_dispersion}
    dispersion_name = k_disp_cold_frit
  []
  [fuel_bed]
    type = NTPPackedBedSegregatedFunctorMaterial
    block = fuel_bed
    porosity = ${eps_bed}
    particle_diameter = ${d_bed}
    pressure_drop_model = kta
    heat_transfer_model = achenbach
    solid_conductivity = ${k_bed}
    axial_conductivity_multiplier = ${kz_bed}
    radial_conductivity_multiplier = ${kr_bed}
    mass_flux_x = G_x
    mass_flux_y = G_y
    mass_flux_z = G_z
    drag_name = drag_bed
    alpha_name = alpha_bed
    kappa_name = kappa_bed
    speed_name = speed_bed
    particle_emissivity = ${bed_emissivity}
    dispersion_coefficient = ${bed_dispersion}
    dispersion_name = k_disp_fuel_bed
  []
  [hot_frit]
    type = NTPPackedBedSegregatedFunctorMaterial
    block = hot_frit
    porosity = ${eps_hf}
    particle_diameter = ${d_hf}
    pressure_drop_model = kta
    heat_transfer_model = kta
    solid_conductivity = ${k_hf}
    mass_flux_x = G_x
    mass_flux_y = G_y
    mass_flux_z = G_z
    drag_name = drag_hf
    alpha_name = alpha_hf
    kappa_name = kappa_hf
    speed_name = speed_hf
    particle_emissivity = ${bed_emissivity}
    dispersion_coefficient = ${bed_dispersion}
    dispersion_name = k_disp_hot_frit
  []
  [alpha]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = alpha
    subdomain_to_prop_value = 'cold_frit alpha_cf fuel_bed alpha_bed hot_frit alpha_hf'
  []
  [kappa_s]
    type = PiecewiseByBlockVectorFunctorMaterial
    prop_name = kappa_s
    subdomain_to_prop_value = 'cold_frit kappa_cf fuel_bed kappa_bed hot_frit kappa_hf'
  []
[]
