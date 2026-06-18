//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* OpenPronghorn is powered by the MOOSE Framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Copyright 2025, Battelle Energy Alliance, LLC
//* ALL RIGHTS RESERVED
//*

#pragma once

#include "Action.h"

#include <string>
#include <vector>

/**
 * Sets up a molten salt corrosion and plating problem from the corrosion database.
 *
 * For each tracked element the action creates the salt-side cation variable and (where the solid is
 * modeled) the solid-side metal variable, the Nernst-Planck transport (diffusion plus
 * electromigration) of the cations, the optional current-continuity equation for the electric
 * potential (the electrophoresis driving force), and the Butler-Volmer electrode kinetics. The
 * kinetics are parameterized from the validated effective-correlation port (MoltenSaltCorrosionModel)
 * using the case feature selectors, so a single interface or boundary reproduces the calibrated
 * dissolution rate.
 *
 * Three topologies are supported:
 *   - two_block: salt and solid are separate mesh blocks; Butler-Volmer is an interface kernel.
 *   - salt_only: only the salt is modeled; Butler-Volmer is a boundary condition (metal is external).
 *   - solid_only: only the solid is modeled; Butler-Volmer is a boundary condition (salt is external).
 */
class CorrosionPlatingAction : public Action
{
public:
  static InputParameters validParams();

  CorrosionPlatingAction(const InputParameters & parameters);

  virtual void act() override;

protected:
  /// Topology of the problem.
  enum class Topology
  {
    TwoBlock,
    SaltOnly,
    SolidOnly
  };

  /// Per-element transport and kinetics plan, assembled from the database and the correlation port.
  struct ElementPlan
  {
    std::string name;        ///< Element symbol, e.g. "Cr".
    std::string salt_var;    ///< Salt-side cation variable name, e.g. "c_Cr".
    std::string solid_var;   ///< Solid-side metal variable name, e.g. "cs_Cr".
    Real valence = 2.0;      ///< Charge number z.
    Real molar_mass = 55.0;  ///< Molar mass [g/mol].
    Real diffusivity = 1.0e-9;      ///< Salt diffusivity [m^2/s].
    Real solid_diffusivity = 1.0e-15; ///< Solid-state diffusivity [m^2/s].
    Real i0 = 0.0;           ///< Exchange current density [A/m^2].
    Real alpha_a = 0.5;      ///< Anodic transfer coefficient.
    Real alpha_c = 0.5;      ///< Cathodic transfer coefficient.
    Real E0 = 0.0;           ///< Standard electrode potential [V].
    Real c_ref = 1.0;        ///< Nernst reference concentration [mol/m^3].
    Real initial_salt = 1.0;       ///< Initial salt concentration [mol/m^3].
    Real initial_solid = 1.0e4;    ///< Initial solid concentration [mol/m^3].
    bool is_reference = false;     ///< Eliminated by electroneutrality (no salt transport variable).
  };

  /// Assemble the element plan from the database and case features (called from the constructor).
  void buildPlan();

  /// Apply the configured block restriction for the salt or solid phase to a parameter set.
  void assignBlocks(InputParameters & params, bool solid) const;

  // Task handlers.
  void addVariables();
  void addAuxVariables();
  void addInitialConditions();
  void addKernels();
  void addInterfaceKernels();
  void addBoundaryConditions();
  void addAuxKernels();

  /// True when the salt phase is modeled (two_block or salt_only).
  bool modelSalt() const { return _topology != Topology::SolidOnly; }
  /// True when the solid phase is modeled (two_block or solid_only).
  bool modelSolid() const { return _topology != Topology::SaltOnly; }

  /// Name of the salt electric potential variable.
  std::string saltPotential() const;
  /// Name of the solid electric potential variable.
  std::string solidPotential() const;

  const Topology _topology;
  const MooseFunctorName _temperature;
  /// Scalar temperature used for the kinetics seed and the Butler-Volmer exponent [K].
  const Real _reference_temperature;
  const bool _solve_potential;
  const bool _supporting_electrolyte;
  const bool _transient;

  /// Tracked elements.
  std::vector<ElementPlan> _elements;
  /// Index of the recession-controlling element in _elements.
  unsigned int _recession_index = 0;
};
