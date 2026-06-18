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

#include "ADIntegratedBC.h"

/**
 * Butler-Volmer electrode reaction applied as a boundary condition when only one phase (salt or
 * solid) is modeled and the metal/salt interface is an external boundary.
 *
 * The same kinetics as ButlerVolmerInterface are used; the absent counter-phase is supplied through
 * functors. The boundary can act on a concentration variable (flux_type = species, molar flux
 * J = i_BV/(z F)) or on a potential variable (flux_type = charge, current i_BV). When the salt is
 * modeled the reacting concentration is the boundary variable itself (or a coupled concentration for
 * a charge boundary); when the solid is modeled the salt concentration is a fixed functor.
 */
class ButlerVolmerBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ButlerVolmerBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// True when the boundary injects charge (i_BV) onto a potential variable, false for molar flux.
  const bool _charge_flux;
  /// True when the modeled phase is the solid metal, false when it is the salt.
  const bool _metal_domain;

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

  /// Whether the modeled-phase potential is a coupled (solved) variable.
  const bool _has_potential;
  /// Modeled-phase electric potential [V] (when coupled).
  const ADVariableValue * const _potential;
  /// Potential of the absent counter-phase [V].
  const Moose::Functor<Real> & _applied_potential;

  /// Whether a salt-ion concentration is coupled (for a charge boundary on the salt).
  const bool _has_concentration;
  /// Coupled salt-ion concentration [mol/m^3].
  const ADVariableValue * const _concentration;
  /// Fixed salt-ion concentration of the absent salt phase [mol/m^3] (solid-domain case).
  const Moose::Functor<Real> * const _counter_concentration;
};
