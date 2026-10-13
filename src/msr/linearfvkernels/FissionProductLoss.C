//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FissionProductLoss.h"
#include "FissionProductData.h"

registerMooseObject("OpenPronghornApp", FissionProductLoss);

InputParameters
FissionProductLoss::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription(
      "Implicit loss of a fission-product nuclide by decay, neutron capture and first-order "
      "removal: (lambda + sigma_c phi + sum_k r_k) N.");
  params.addRangeCheckedParam<Real>(
      "decay_constant", 0.0, "decay_constant>=0", "Decay constant [1/s].");
  params.addRangeCheckedParam<Real>(
      "capture_cross_section", 0.0, "capture_cross_section>=0", "Capture cross section [barn].");
  params.addParam<MooseFunctorName>(
      "neutron_flux", "0", "One-group neutron flux [1/(m^2 s)] driving the capture loss.");
  params.addParam<std::vector<MooseFunctorName>>(
      "removal_rates", {}, "First-order removal rates [1/s] (off-gas, plating, processing).");
  return params;
}

FissionProductLoss::FissionProductLoss(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _lambda(getParam<Real>("decay_constant")),
    _sigma(getParam<Real>("capture_cross_section") * FissionProducts::barn),
    _flux(getFunctor<Real>("neutron_flux"))
{
  for (const auto & name : getParam<std::vector<MooseFunctorName>>("removal_rates"))
    _removal.push_back(&getFunctor<Real>(name));
}

Real
FissionProductLoss::computeMatrixContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const auto state = determineState();
  Real coeff = _lambda + _sigma * _flux(elem_arg, state);
  for (const auto * r : _removal)
    coeff += (*r)(elem_arg, state);
  return coeff * _current_elem_volume;
}

Real
FissionProductLoss::computeRightHandSideContribution()
{
  return 0.0;
}
