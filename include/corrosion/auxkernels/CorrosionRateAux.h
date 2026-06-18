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

#include "AuxKernel.h"

/**
 * Diagnostic Butler-Volmer output for an electrode boundary: the interfacial current density
 * [A/m^2], the equivalent penetration rate [um/y], or the time-integrated metal recession / plating
 * thickness [um]. The same kinetics as the Butler-Volmer objects are evaluated (non-AD), so the
 * value tracks the corrosion or plating the interface is driving and provides a directly comparable
 * quantity for validation against the reference correlation.
 */
class CorrosionRateAux : public AuxKernel
{
public:
  static InputParameters validParams();

  CorrosionRateAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Output quantity.
  enum class Mode
  {
    Current,
    PenetrationRate,
    Recession
  };
  const Mode _mode;

  /// Standard electrode potential E0 [V].
  const Real _E0;
  /// Charge number z of the cation.
  const Real _valence;
  /// Molar mass of the metal [g/mol].
  const Real _molar_mass;
  /// Metal density [g/cm^3].
  const Real _density;
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

  /// Salt-side concentration of the cation [mol/m^3].
  const VariableValue & _concentration;
  /// Salt-phase potential [V] (coupled) or constant fallback.
  const bool _has_salt_potential;
  const VariableValue & _salt_potential_var;
  const Real _salt_potential;
  /// Metal-phase potential [V] (coupled) or constant fallback.
  const bool _has_metal_potential;
  const VariableValue & _metal_potential_var;
  const Real _metal_potential;

  /// Old value of the recession variable (for time accumulation).
  const VariableValue & _u_old;
};
