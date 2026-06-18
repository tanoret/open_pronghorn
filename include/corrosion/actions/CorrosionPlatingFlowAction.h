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
 * Sets up salt-side molten salt corrosion of a wall on the linear finite-volume segregated basis, so
 * that the dissolved corrosion products are transported as passive scalars by the same flow solve
 * that carries the radiolysis chemistry and energy. This is the corrosion counterpart of
 * MoltenSaltRadiolysisAction, and it is the natural way to couple corrosion into a flowing MSR
 * (Navier-Stokes + k-epsilon + energy + radiolysis) solved with SIMPLE/PIMPLE.
 *
 * For each tracked element the action creates a linear finite-volume cation concentration, its time
 * derivative, advection by the solved flow (Rhie-Chow) or a prescribed velocity, and molecular plus
 * turbulent diffusion, and applies the Butler-Volmer electrode reaction at the wall through
 * CorrosionLinearFVButlerVolmerBC. The exchange current is seeded from the validated effective
 * corrosion correlation (MoltenSaltCorrosionModel) and is Arrhenius in the local temperature, so the
 * wall corrodes faster where the salt is hotter. The full electromigration / electric-potential solve
 * is provided by the finite-element CorrosionPlating action; here the bulk transport is advective and
 * diffusive, as in the radiolysis framework.
 */
class CorrosionPlatingFlowAction : public Action
{
public:
  static InputParameters validParams();

  CorrosionPlatingFlowAction(const InputParameters & parameters);

  virtual void act() override;

protected:
  /// Per-element transport and kinetics plan.
  struct ElementPlan
  {
    std::string name;
    std::string salt_var;
    Real valence = 2.0;
    Real molar_mass = 55.0;
    Real i0 = 0.0;
    Real alpha_a = 0.5;
    Real alpha_c = 0.5;
    Real E0 = 0.0;
    Real c_ref = 1.0;
    Real initial_salt = 0.0;
    /// Whether this action owns the variable (creates it and its transport) or releases the
    /// corrosion product into an existing variable (e.g. the radiolysis-tracked cation).
    bool owns_variable = true;
  };

  void buildPlan();
  void maybeAssignBlocks(InputParameters & params) const;
  std::string systemName(const std::string & var) const { return var + "_sys"; }

  void addVariables();
  void addInitialConditions();
  void addKernels();
  void addBoundaryConditions();

  const MooseFunctorName _temperature;
  const Real _reference_temperature;
  const bool _temperature_dependent;
  const bool _transient;

  std::vector<ElementPlan> _elements;
};
