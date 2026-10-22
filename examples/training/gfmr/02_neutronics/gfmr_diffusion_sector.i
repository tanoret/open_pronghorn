################################################################################
# GFMR training -- Stage 2b: 8-group heterogeneous neutron diffusion on the 3-D CFD sector mesh
#
# Same physics as gfmr_diffusion_rz.i, but on the 60-degree sector mesh of the CFD
# (../01_mesh/gfmr_sector_3d_cfd.e): every fuel channel, control-rod guide tube, absorber, helium gap
# and the graphite around them is resolved. Each mesh block carries the 8-group constants of the
# matching OpenMC cell region (cell-domain MGXS of the nominal case, xs/gfmr_xs.i), so the power
# density lands on the CFD cells without any mapping.
#
# Boundaries: Marshak vacuum on the top, bottom and the reflector side; reflective on the two
# symmetry planes (natural boundary condition of the diffusion operator).
#
# The Hf absorber is close to black for thermal neutrons, where diffusion theory is poor; its
# removal matrix is scaled (make_diffusion_input.py --rod-factor) so that the worth of the full
# insertion matches the OpenMC reference: a calibration, not a prediction (see README).
#
# Run (24 ranks, ~3.6e6 unknowns):
#   mpirun -n 24 open_pronghorn-opt -i gfmr_diffusion_sector.i
################################################################################

!include xs/gfmr_xs.i

power = 150e6
# 1/6 of the core is modelled
power_sector = ${fparse power / 6}

salt = 'salt_inner_tri salt_inner salt_outer_tri salt_outer salt_cr'
# mesh of the rod position: gfmr_sector_3d_cfd.e (nominal, 36 cm), gfmr_sector_3d_rods_in.e or _rods_out.e
mesh_file = '../01_mesh/gfmr_sector_3d_cfd.e'
# blocks holding Hf / He in that mesh (rods out: absorber = '' and helium also the core)
absorber = 'absorber_tri absorber'
helium = 'helium helium_tri'

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = ${mesh_file}
  []
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
  # eigen tag: B = -F so that MOOSE's lambda = 1/k (see gfmr_diffusion_rz.i)
  [fission]
    type = ArrayReaction
    variable = phi
    reaction_coefficient = F_eig
    block = ${salt}
    extra_vector_tags = 'eigen'
  []
[]

[BCs]
  [vacuum]
    type = ArrayVacuumBC
    variable = phi
    boundary = 'top bottom reflector_outer'
    alpha = ${VAC_ALPHA}
  []
[]

[Materials]
  # ---- fuel salt: the three channel types have their own constants
  [inner_salt]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_inner_salt}
    block = 'salt_inner_tri salt_inner'
  []
  [inner_salt_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_inner_salt}
    block = 'salt_inner_tri salt_inner'
  []
  [inner_salt_F]
    type = GenericConstant2DArray
    prop_name = F_eig
    prop_value = ${FN_inner_salt}
    block = 'salt_inner_tri salt_inner'
  []
  [outer_salt]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_outer_salt}
    block = 'salt_outer_tri salt_outer'
  []
  [outer_salt_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_outer_salt}
    block = 'salt_outer_tri salt_outer'
  []
  [outer_salt_F]
    type = GenericConstant2DArray
    prop_name = F_eig
    prop_value = ${FN_outer_salt}
    block = 'salt_outer_tri salt_outer'
  []
  [cr_salt]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_cr_salt}
    block = 'salt_cr'
  []
  [cr_salt_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_cr_salt}
    block = 'salt_cr'
  []
  [cr_salt_F]
    type = GenericConstant2DArray
    prop_name = F_eig
    prop_value = ${FN_cr_salt}
    block = 'salt_cr'
  []
  # ---- graphite: the mesh does not split the moderator by zone, the inner/outer/cr graphite
  # D differs by < 0.2 % and the removal by up to 2.4 % (cr graphite, fast groups); the outer-zone
  # set (88 % of the moderator) is used
  [graphite]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_outer_graphite}
    block = 'graphite'
  []
  [graphite_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_outer_graphite}
    block = 'graphite'
  []
  [plug]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_plug_graphite}
    block = 'graphite_plug_tri graphite_plug'
  []
  [plug_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_plug_graphite}
    block = 'graphite_plug_tri graphite_plug'
  []
  [reflector]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_reflector_lattice}
    block = 'reflector'
  []
  [reflector_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_reflector_lattice}
    block = 'reflector'
  []
  # ---- control-rod channel
  [guide]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_cr_guide}
    block = 'guide_tube'
  []
  [guide_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_cr_guide}
    block = 'guide_tube'
  []
  [helium]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_cr_gap}
    block = ${helium}
  []
  [helium_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_cr_gap}
    block = ${helium}
  []
  [absorber]
    type = GenericConstantArray
    prop_name = D
    prop_value = ${D_cr_rod}
    block = ${absorber}
  []
  # absorber removal matrix with the transport correction of the black Hf rod
  # (make_diffusion_input.py --rod-factor 0.602, calibrated by hand against the OpenMC total
  # worth; see README)
  [absorber_R]
    type = GenericConstant2DArray
    prop_name = R
    prop_value = ${REM_cr_rod}
    block = ${absorber}
  []
[]

[AuxVariables]
  [power_density]
    order = CONSTANT
    family = MONOMIAL
    block = ${salt}
  []
  [fast_flux]
    # E > 0.821 MeV + 5.53 keV - 0.821 MeV (groups 1-2); the graphite lifetime limit is set on
    # E > 0.1 MeV, so this is an upper bound (see README)
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [power_inner]
    type = ArrayVarReductionAux
    variable = power_density
    array_variable = phi
    value_type = sum
    weights = ${KSF_inner_salt}
    block = 'salt_inner_tri salt_inner'
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [power_outer]
    type = ArrayVarReductionAux
    variable = power_density
    array_variable = phi
    value_type = sum
    weights = ${KSF_outer_salt}
    block = 'salt_outer_tri salt_outer'
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [power_cr]
    type = ArrayVarReductionAux
    variable = power_density
    array_variable = phi
    value_type = sum
    weights = ${KSF_cr_salt}
    block = 'salt_cr'
    execute_on = 'LINEAR TIMESTEP_END'
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
  # kappa-fission power of the sector normalized to 150 MW / 6
  normalization = total_power
  normal_factor = ${power_sector}
  free_power_iterations = 6
  # Newton converges to |R| ~ 1e-8 in four steps (eigenvalue to 1e-9). PETSc's default
  # nonlinear step tolerance then stops the next step with SNES_DIVERGED_LINE_SEARCH-like
  # stagnation; accept convergence on the step size as well.
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-7
  nl_abs_step_tol = 1e-10
  # GMRES restart: the default (30) stalls on the 8-group block system at ~60 iterations
  l_tol = 1e-5
  l_max_its = 400
  petsc_options = '-ksp_gmres_modifiedgramschmidt'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 200'
[]

[Postprocessors]
  [total_power]
    type = ElementIntegralVariablePostprocessor
    variable = power_density
    block = ${salt}
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [power_inner]
    type = ElementIntegralVariablePostprocessor
    variable = power_density
    block = 'salt_inner_tri salt_inner'
    execute_on = TIMESTEP_END
  []
  [power_outer]
    type = ElementIntegralVariablePostprocessor
    variable = power_density
    block = 'salt_outer_tri salt_outer'
    execute_on = TIMESTEP_END
  []
  [power_cr]
    type = ElementIntegralVariablePostprocessor
    variable = power_density
    block = 'salt_cr'
    execute_on = TIMESTEP_END
  []
  [peak_power_density]
    type = ElementExtremeValue
    variable = power_density
    block = ${salt}
    execute_on = TIMESTEP_END
  []
  [peak_fast_flux_graphite]
    type = ElementExtremeValue
    variable = fast_flux
    block = 'graphite'
    execute_on = TIMESTEP_END
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
  # axial power along the central channel (cell (0,0), inner zone)
  [axial_center]
    type = LineValueSampler
    variable = power_density
    start_point = '0.005 0.003 0.045'
    end_point = '0.005 0.003 3.555'
    num_points = 40
    sort_by = z
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
