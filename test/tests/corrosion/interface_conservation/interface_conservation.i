# Two-block (salt + solid) corrosion interface: the Butler-Volmer interface kernel transfers
# chromium from the solid metal into the salt. With no-flux outer boundaries and no volumetric
# source, the total chromium inventory (salt + solid) must be conserved: every mole that leaves the
# solid appears in the salt. This verifies the Element/Neighbor flux balance of ButlerVolmerInterface.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 40
    xmax = 2.0
  []
  [solid_block]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    bottom_left = '1.0 0 0'
    top_right = '2.0 0 0'
  []
  [names]
    type = RenameBlockGenerator
    input = solid_block
    old_block = '0 1'
    new_block = 'salt solid'
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = names
    primary_block = salt
    paired_block = solid
    new_boundary = interface
  []
[]

[Variables]
  [c_Cr]
    block = salt
    initial_condition = 1.0
  []
  [cs_Cr]
    block = solid
    initial_condition = 1000.0
  []
[]

[Kernels]
  [c_time]
    type = ADTimeDerivative
    variable = c_Cr
    block = salt
  []
  [c_diffusion]
    type = CorrosionNernstPlanckFlux
    variable = c_Cr
    diffusivity = 1.0e-9
    valence = 2
    temperature = 923.15
    block = salt
  []
  [cs_time]
    type = ADTimeDerivative
    variable = cs_Cr
    block = solid
  []
  [cs_diffusion]
    type = CorrosionSolidDiffusion
    variable = cs_Cr
    diffusivity = 1.0e-11
    block = solid
  []
[]

[InterfaceKernels]
  [bv]
    type = ButlerVolmerInterface
    variable = c_Cr
    neighbor_var = cs_Cr
    boundary = interface
    valence = 2
    exchange_current_density = 50.0
    temperature = 923.15
    metal_potential = 0.05
    salt_potential = 0.0
    c_ref = 1.0
  []
[]

[Postprocessors]
  [cr_salt]
    type = ElementIntegralVariablePostprocessor
    variable = c_Cr
    block = salt
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_solid]
    type = ElementIntegralVariablePostprocessor
    variable = cs_Cr
    block = solid
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cr_total]
    type = LinearCombinationPostprocessor
    pp_names = 'cr_salt cr_solid'
    pp_coefs = '1.0 1.0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 2.0e4
  num_steps = 5
  automatic_scaling = true
  nl_abs_tol = 1.0e-12
[]

[Outputs]
  csv = true
[]
