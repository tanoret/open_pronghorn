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

#include "ADInterfaceKernel.h"

/**
 * Butler-Volmer charge flux coupling the interfacial current into the potential (current-continuity)
 * equations on the two blocks.
 *
 * The primary (this-side) variable is the salt-phase potential phi_salt; the neighbor variable is
 * the metal-phase potential phi_solid. The interfacial Butler-Volmer current i_BV enters the salt
 * charge balance as a source (-test i_BV) and leaves the metal charge balance as a sink
 * (+test_neighbor i_BV), conserving charge across the double layer. Add one instance per tracked
 * element; the currents are additive.
 */
class ButlerVolmerPotentialInterface : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

  ButlerVolmerPotentialInterface(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;

  /// Standard electrode potential E0 [V].
  const Real _E0;
  /// Charge number z of the cation.
  const Real _valence;
  /// Exchange current density i0 [A/m^2].
  const Real _i0;
  /// Anodic charge-transfer coefficient.
  const Real _alpha_a;
  /// Cathodic charge-transfer coefficient.
  const Real _alpha_c;
  /// Temperature [K].
  const Real _temperature;
  /// Reference concentration for the Nernst term [mol/m^3].
  const Real _c_ref;
  /// Salt-side concentration of the cation on this block [mol/m^3].
  const ADVariableValue & _concentration;
};
