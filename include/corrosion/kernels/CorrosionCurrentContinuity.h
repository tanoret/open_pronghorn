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

#include "ADKernelGrad.h"

#include <vector>

/**
 * Current-continuity (charge-conservation) equation for the salt electric potential under the
 * dilute-solution / MacInnes closure.
 *
 * The ionic current density is
 *   i = -F sum_k z_k D_k grad c_k - kappa grad(phi),   kappa = (F^2 / R T) sum_k z_k^2 D_k c_k,
 * and div(i) = 0 yields the weak (grad-test) residual
 *   (grad test, kappa grad phi) + (grad test, F sum_k z_k D_k grad c_k).
 * Electroneutrality is enforced separately by eliminating one ion as the reference species.
 */
class CorrosionCurrentContinuity : public ADKernelGrad
{
public:
  static InputParameters validParams();

  CorrosionCurrentContinuity(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  /// Number of coupled salt-ion concentrations.
  const unsigned int _n_ions;
  /// Temperature functor [K].
  const Moose::Functor<Real> & _temperature;
  /// Charge number z_k of each coupled ion.
  const std::vector<Real> _valences;
  /// Diffusivity D_k of each coupled ion [m^2/s].
  const std::vector<Real> _diffusivities;
  /// Concentration of each coupled ion [mol/m^3].
  std::vector<const ADVariableValue *> _c;
  /// Concentration gradient of each coupled ion [mol/m^4].
  std::vector<const ADVariableGradient *> _grad_c;
};
