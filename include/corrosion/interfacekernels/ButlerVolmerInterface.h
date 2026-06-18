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
 * Butler-Volmer mass flux across a metal/salt interface modeled as two mesh blocks.
 *
 * The primary (this-side) variable is the salt-side cation concentration c [mol/m^3]; the neighbor
 * (other-side) variable is the solid-side metal concentration c_s. The interfacial molar flux is
 *   J = i_BV / (z F),   i_BV = Butler-Volmer current density,
 * with i_BV anodic (positive) when the metal dissolves. The flux enters the salt as a source
 * (-test J on the Element side) and leaves the solid as a sink (+test_neighbor J on the Neighbor
 * side), so mass is conserved across the interface and a single object handles both dissolution and
 * plating.
 *
 * The metal and salt phase potentials drive the overpotential. They may be supplied as coupled
 * variables (phi_solid on the neighbor block, phi_salt on this block; fully coupled, AD-exact) or as
 * fixed constants (metal_potential, salt_potential) for a supporting-electrolyte problem in which
 * the potential field is not solved.
 */
class ButlerVolmerInterface : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

  ButlerVolmerInterface(const InputParameters & parameters);

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

  /// Whether the salt-phase potential is coupled (solved) on this side.
  const bool _has_salt_potential;
  /// Salt-phase potential on this side [V] (when coupled).
  const ADVariableValue * const _phi_salt;
  /// Whether the metal-phase potential is coupled (solved) on the neighbor side.
  const bool _has_metal_potential;
  /// Metal-phase potential on the neighbor side [V] (when coupled).
  const ADVariableValue * const _phi_solid;
  /// Fixed salt-phase potential when not coupled [V].
  const Real _salt_potential;
  /// Fixed metal-phase potential when not coupled [V].
  const Real _metal_potential;
};
