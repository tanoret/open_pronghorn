################################################################################
# GFMR training -- Stage 2: 8-group neutron diffusion eigenvalue, 2-D R-Z
#
# Continuous Galerkin (CFEM) multigroup diffusion on an R-Z model of the GFMR. The cross
# sections are the OpenMC (ENDF/B-VIII.0) flux-volume homogenized constants of the lattice and
# the graphite reflector (02_neutronics/openmc -> xs/gfmr_8g.json -> xs/gfmr_xs.i).
#
#   -div(D grad phi) + R phi = (1/k) F phi        (phi: 8-component array variable)
#
# Regions (R-Z, z = 0 at the bottom of the graphite):
#   core       0 <= r <= R_core, 0 <= z <= H  homogenized salt + graphite lattice; R_core is
#              the radius of the circle with the area of the 397 unit cells (379 + 18 plugs)
#   reflector  R_core <= r <= 1.75 m          graphite
#   rodded     the top 36 cm of an annulus at the radius of the inner ring-4 rod circle with the
#              area of the six ring-4 rods + 3 central rods smeared (see README; the 3-D rod worth
#              of OpenMC is the reference and the smeared Hf density is calibrated to it)
# Boundaries: Marshak vacuum on the top, bottom and outer surfaces.
#
# Outputs: k_eff, the power density (W/m^3, normalized to 150 MW), the fast flux (> 0.1 MeV,
# groups 1-2) for the graphite lifetime, and axial/radial profiles.
################################################################################

!include xs/gfmr_xs.i

H = 3.6
R_refl = 1.75
# 397 hexagonal unit cells of pitch 0.14 m: area 397 * sqrt(3)/2 * 0.14^2 = 6.739 m^2
R_core = ${fparse sqrt(397 * sqrt(3) / 2 * 0.14^2 / pi)}
power = 150e6
# fraction of the total heating deposited in the fuel salt (OpenMC coupled n-gamma heating)
f_salt = 0.92

nr_core = 30
nr_refl = 8
nz = 36

[Mesh]
  [gen]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${R_core} ${fparse R_refl - R_core}'
    ix = '${nr_core} ${nr_refl}'
    dy = '${H}'
    iy = '${nz}'
    subdomain_id = '1 2'
  []
  [names]
    type = RenameBlockGenerator
    input = gen
    old_block = '1 2'
    new_block = 'core reflector'
  []
  coord_type = RZ
  rz_coord_axis = Y
[]

[Variables]
  [phi]
    order = FIRST
    family = LAGRANGE
    components = ${num_groups}
    initial_condition = '1 1 1 1 1 1 1 1'
  []
[]

[Kernels]
  [diffusion]
    type = ArrayDiffusion
    variable = phi
    diffusion_coefficient = D
  []
  [removal]
    type = ArrayReaction
    variable = phi
    reaction_coefficient = R
  []
  # MOOSE solves A phi = lambda B phi with B assembled from the 'eigen' tag; ArrayReaction adds
  # +F phi, so B = -F and lambda = -1/k. F_eig = -F_total makes lambda = 1/k.
  [fission]
    type = ArrayReaction
    variable = phi
    reaction_coefficient = F_eig
    block = core
    extra_vector_tags = 'eigen'
  []
[]

[BCs]
  [vacuum]
    type = ArrayVacuumBC
    variable = phi
    boundary = 'top bottom right'
    alpha = ${VAC_ALPHA}
  []
[]

[Materials]
  [core]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_core}
    block = core
  []
  [core_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_core}
    block = core
  []
  [core_F]
    type = GenericConstant2DArray
    prop_name = F_eig
    prop_value = ${F_total_neg}
    block = core
  []
  [refl]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_refl}
    block = reflector
  []
  [refl_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_refl}
    block = reflector
  []
[]

[AuxVariables]
  # local thermal power deposited in the salt-graphite lattice, W/m^3 (homogenized)
  [power_density]
    order = CONSTANT
    family = MONOMIAL
    block = core
  []
  # fission neutron production rate, for the precursor source
  [nu_fission_rate]
    order = CONSTANT
    family = MONOMIAL
    block = core
  []
  # fast flux, E > 0.821 MeV (group 1) and 5.53 keV - 0.821 MeV (group 2), for the graphite lifetime
  [fast_flux]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [power_density]
    type = ArrayVarReductionAux
    variable = power_density
    array_variable = phi
    value_type = sum
    weights = ${KSF_core}
    block = core
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [nu_fission_rate]
    type = ArrayVarReductionAux
    variable = nu_fission_rate
    array_variable = phi
    value_type = sum
    weights = ${NSF_core}
    block = core
    execute_on = TIMESTEP_END
  []
  [fast_flux]
    type = ArrayVarReductionAux
    variable = fast_flux
    array_variable = phi
    value_type = sum
    weights = '1 1 0 0 0 0 0 0'
    execute_on = TIMESTEP_END
  []
[]

[Executioner]
  type = Eigenvalue
  solve_type = PJFNK
  # Scale the eigenvector so that the kappa-fission power integrates to 150 MW. The kappa
  # values include all recoverable energy (fragments, betas, gammas, neutrons), so this is the
  # thermal power of the core; OpenMC distributes 8 % of it to the graphite (see README).
  normalization = total_power
  normal_factor = ${power}
  free_power_iterations = 4
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Postprocessors]
  [total_power]
    type = ElementIntegralVariablePostprocessor
    variable = power_density
    block = core
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [lambda]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_real
    index = 0
    execute_on = TIMESTEP_END
  []
  [k_eff]
    type = ParsedPostprocessor
    expression = '1 / lambda'
    pp_names = 'lambda'
    execute_on = TIMESTEP_END
  []
[]

[VectorPostprocessors]
  [eigenvalues]
    type = Eigenvalues
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
