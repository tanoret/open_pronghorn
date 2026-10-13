//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FissionProductFeed.h"
#include "FissionProductData.h"

registerMooseObject("OpenPronghornApp", FissionProductFeed);

InputParameters
FissionProductFeed::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription(
      "Production of a nuclide from a parent inventory by decay, neutron "
      "capture or first-order transfer: (lambda_p b + sigma_p phi + r) N_p.");
  params.addRequiredParam<MooseFunctorName>("parent", "The parent inventory [mol/m^3].");
  params.addRangeCheckedParam<Real>("decay_coefficient",
                                    0.0,
                                    "decay_coefficient>=0",
                                    "Partial decay constant lambda_p * branching [1/s].");
  params.addRangeCheckedParam<Real>("capture_cross_section",
                                    0.0,
                                    "capture_cross_section>=0",
                                    "Parent capture cross section [barn] when this nuclide is "
                                    "the capture product.");
  params.addParam<MooseFunctorName>("neutron_flux", "0", "One-group neutron flux [1/(m^2 s)].");
  params.addParam<MooseFunctorName>(
      "transfer_rate", "0", "First-order transfer rate from the parent inventory [1/s].");
  return params;
}

FissionProductFeed::FissionProductFeed(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _parent(getFunctor<Real>("parent")),
    _decay_coefficient(getParam<Real>("decay_coefficient")),
    _sigma(getParam<Real>("capture_cross_section") * FissionProducts::barn),
    _flux(getFunctor<Real>("neutron_flux")),
    _transfer_rate(getFunctor<Real>("transfer_rate"))
{
}

Real
FissionProductFeed::computeMatrixContribution()
{
  return 0.0;
}

Real
FissionProductFeed::computeRightHandSideContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const auto state = determineState();
  const Real rate =
      _decay_coefficient + _sigma * _flux(elem_arg, state) + _transfer_rate(elem_arg, state);
  return rate * _parent(elem_arg, state) * _current_elem_volume;
}
