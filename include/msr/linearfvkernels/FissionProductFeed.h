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
 * Production of a nuclide from a parent inventory,
 *
 *   S = (lambda_p * b + sigma_p * phi + r) * N_p,
 *
 * where lambda_p * b is the parent's partial decay constant to this nuclide, sigma_p * phi its
 * capture rate (when this nuclide is the capture product) and r a first-order transfer rate
 * (e.g. the rate at which a salt nuclide is moved to a plated or off-gas inventory). The parent
 * is evaluated at its latest iterate, so ordering the parent's solver system before the
 * daughter's makes the coupling implicit within a single fixed-point sweep.
 */
class FissionProductFeed : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  FissionProductFeed(const InputParameters & params);

  virtual Real computeMatrixContribution() override;
  virtual Real computeRightHandSideContribution() override;

protected:
  /// Parent inventory [mol/m^3]
  const Moose::Functor<Real> & _parent;
  /// Partial decay constant lambda_p * b [1/s]
  const Real _decay_coefficient;
  /// Parent capture cross section [m^2]
  const Real _sigma;
  /// One-group neutron flux [1/(m^2 s)]
  const Moose::Functor<Real> & _flux;
  /// First-order transfer rate [1/s]
  const Moose::Functor<Real> & _transfer_rate;
};
