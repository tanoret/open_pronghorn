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
 *
 * The salt concentration and either phase potential may be supplied as a coupled variable or a
 * functor. This lets a solid-only diagnostic consume the exact external-salt functors used by its
 * ButlerVolmerBC while retaining the original coupled-variable API.
 */
class CorrosionRateAux : public AuxKernel
{
public:
  static InputParameters validParams();

  CorrosionRateAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Evaluate a configured functor at the current node, boundary quadrature point, or element QP.
  Real functorValue(const Moose::Functor<Real> & functor) const;

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

  /// Salt-side concentration of the cation [mol/m^3], as one exclusive variable/functor source.
  const bool _has_concentration;
  const VariableValue & _concentration;
  const Moose::Functor<Real> * const _concentration_functor;
  /// Salt-phase potential [V], as one exclusive coupled/functor/fixed source.
  const bool _has_salt_potential;
  const VariableValue & _salt_potential_var;
  const Moose::Functor<Real> * const _salt_potential_functor;
  const Real _salt_potential;
  /// Metal-phase potential [V], as one exclusive coupled/functor/fixed source.
  const bool _has_metal_potential;
  const VariableValue & _metal_potential_var;
  const Moose::Functor<Real> * const _metal_potential_functor;
  const Real _metal_potential;

  /// Old value of the recession variable (for time accumulation).
  const VariableValue & _u_old;
};

