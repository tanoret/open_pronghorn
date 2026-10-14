# =====================================================================================================
# Training case 02: the chloride radical chain in irradiated LiCl-KCl (eutectic, 400 C).
#
# Gamma radiolysis of a chloride melt creates solvated electrons and chlorine atoms (G = 0.3 each):
#   Cl_rad + Cl-      -> Cl2m_rad        (fast, pseudo-first order in a 2e4 mol/m^3 chloride melt)
#   e_sol + Cl_rad    -> Cl-
#   2 Cl_rad          -> Cl2
#   2 Cl2m_rad        -> Cl3- + Cl-      (second order: sets the radical steady state)
#   Cl3-              -> Cl2 + Cl-       (k = 1 1/s)
#
# Every cell of this 1D mesh is an independent well-mixed (0D) sample: there is no transport, and the
# dose rate is a functor that steps by a decade from cell to cell, so a single run is a dose scan
# from 1e2 to 1e5 Gy/s (J/kg/s ~ J/m^3/s here per unit density; the model takes J/m^3/s).
#
# Analytic checks (S = G D / (100 eV N_A) is the radical source in mol/m^3/s):
#   Cl2m_rad -> sqrt(S / (2 k))    (scales as dose^0.5)
#   Cl3-     -> S / 2              (scales as dose^1)
#   e_sol    =  S t                (no electron scavenger: accumulates linearly)
#   Cl2      -> S t / 2            (accumulates: chlorine is the stable radiolysis product)
# =====================================================================================================

n_doses = 4

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = ${n_doses}
    xmax = ${n_doses}
  []
[]

[Functions]
  # Cell i (centroid x = i + 0.5) receives 10^(2+i) J/m^3/s
  [dose_rate]
    type = ParsedFunction
    expression = '10^(2 + floor(x))'
  []
[]

[Problem]
  linear_sys_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys'
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  temperature = 673.15
  dose_rate = dose_rate
  initial_condition_species = 'Cl_ion'
  initial_condition_values = '2.0e4'
[]

[VectorPostprocessors]
  # One row per cell (= per dose rate)
  [dose_scan]
    type = ElementValueSampler
    variable = 'e_sol Cl_rad Cl2m_rad Cl3_ion Cl2_diss'
    sort_by = x
    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]
  # Time histories in the highest-dose cell (1e5 J/m^3/s)
  [c_e_sol]
    type = PointValue
    variable = e_sol
    point = '${fparse n_doses - 0.5} 0 0'
  []
  [c_Cl_rad]
    type = PointValue
    variable = Cl_rad
    point = '${fparse n_doses - 0.5} 0 0'
  []
  [c_Cl2m_rad]
    type = PointValue
    variable = Cl2m_rad
    point = '${fparse n_doses - 0.5} 0 0'
  []
  [c_Cl3_ion]
    type = PointValue
    variable = Cl3_ion
    point = '${fparse n_doses - 0.5} 0 0'
  []
  [c_Cl2_diss]
    type = PointValue
    variable = Cl2_diss
    point = '${fparse n_doses - 0.5} 0 0'
  []
[]

[Convergence]
  # The species systems are coupled through lagged reaction terms: sweep them a fixed number of
  # times per step (Gauss-Seidel over the systems in 'system_names' order)
  [fp]
    type = IterationCountConvergence
    max_iterations = 20
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys Cl_ion_sys Cl_rad_sys Cl2m_rad_sys Cl3_ion_sys Cl2_diss_sys'
  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp
  scheme = implicit-euler
  start_time = 0
  end_time = 1e3
  # Radical lifetimes are microseconds: start tiny and let the step grow geometrically
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-9
    growth_factor = 1.5
  []
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  [out]
    type = CSV
    # Write the dose-scan table only once, at the end of the run
    execute_vector_postprocessors_on = FINAL
  []
[]
