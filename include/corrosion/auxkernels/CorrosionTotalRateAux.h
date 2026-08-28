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

#include <vector>

/**
 * Sum species-resolved Butler-Volmer currents after converting each one to its own equivalent
 * penetration rate, or integrate that total rate as a recession thickness.
 *
 * This is intentionally separate from CorrosionRateAux: the latter retains its legacy
 * one-controlling-element API and behavior, while this object provides the physically meaningful
 * total-front diagnostic for a source-resolved multi-element kinetics model.
 */
class CorrosionTotalRateAux : public AuxKernel
{
public:
  static InputParameters validParams();

  CorrosionTotalRateAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Evaluate a configured functor at the current node, boundary quadrature point, or element QP.
  Real functorValue(const Moose::Functor<Real> & functor) const;

  enum class Mode
  {
    PenetrationRate,
    Recession
  };
  const Mode _mode;

  /// Per-species electrochemical and Faradaic-conversion properties.
  const std::vector<Real> _valences;
  const std::vector<Real> _molar_masses;
  const std::vector<Real> _exchange_current_densities;
  const std::vector<Real> _E0_values;
  const std::vector<Real> _alpha_a_values;
  const std::vector<Real> _alpha_c_values;
  const std::vector<Real> _c_ref_values;
  const Real _density;
  const Real _temperature;

  /// Salt concentrations, either one coupled variable per species or one shared external functor.
  std::vector<const VariableValue *> _concentrations;
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

  /// Old recession value for time integration.
  const VariableValue & _u_old;
};
