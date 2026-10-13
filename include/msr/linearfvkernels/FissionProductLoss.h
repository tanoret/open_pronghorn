//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVElementalKernel.h"

/**
 * Implicit first-order loss of a fission-product nuclide,
 *
 *   L = (lambda + sigma_c * phi + sum_k r_k) * N,
 *
 * combining radioactive decay (lambda), neutron capture (burnout, sigma_c * phi) and any number of
 * first-order removal pathways r_k (off-gas sparging, plating, chemical processing).
 */
class FissionProductLoss : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  FissionProductLoss(const InputParameters & params);

  virtual Real computeMatrixContribution() override;
  virtual Real computeRightHandSideContribution() override;

protected:
  /// Decay constant [1/s]
  const Real _lambda;
  /// Capture cross section [m^2]
  const Real _sigma;
  /// One-group neutron flux [1/(m^2 s)]
  const Moose::Functor<Real> & _flux;
  /// First-order removal rates [1/s]
  std::vector<const Moose::Functor<Real> *> _removal;
};
