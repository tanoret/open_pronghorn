//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MSRReaction.h"

#include <cmath>

registerMooseObject("OpenPronghornApp", MSRReaction);

InputParameters
MSRReaction::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription("Adds the mass-action source or sink of a single molten salt "
                             "radiolysis reaction to one species in a linear finite-volume system.");

  params.addRequiredParam<MooseFunctorName>("temperature", "The temperature functor [K].");

  MooseEnum mode("consumption production");
  params.addRequiredParam<MooseEnum>(
      "mode",
      mode,
      "Whether the species is consumed (reactant, implicit sink) or produced (product, source).");

  params.addRequiredParam<Real>("stoichiometric_coefficient",
                                "The net stoichiometric coefficient of the species in this "
                                "reaction (always positive).");
  params.addParam<Real>(
      "self_order",
      1.0,
      "The kinetic order of the species itself in this reaction (used in consumption mode).");

  params.addParam<MooseFunctorName>("partner1", "First partner reactant species concentration.");
  params.addParam<Real>("partner1_order", 1.0, "Kinetic order of the first partner reactant.");
  params.addParam<MooseFunctorName>("partner2", "Second partner reactant species concentration.");
  params.addParam<Real>("partner2_order", 1.0, "Kinetic order of the second partner reactant.");

  // Arrhenius rate-constant parameters (mirror the molten salt radiolysis database). A value of
  // zero marks an unsupplied parameter.
  params.addParam<Real>("A", 0.0, "Arrhenius pre-exponential factor [SI].");
  params.addParam<Real>("k_ref", 0.0, "Reference rate constant [SI].");
  params.addParam<Real>("T_ref", 0.0, "Reference temperature for k_ref [K].");
  params.addParam<Real>("Ea", 0.0, "Activation energy [J/mol].");

  return params;
}

MSRReaction::MSRReaction(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _temperature(getFunctor<Real>("temperature")),
    _mode(getParam<MooseEnum>("mode") == "consumption" ? Mode::Consumption : Mode::Production),
    _stoichiometric_coefficient(getParam<Real>("stoichiometric_coefficient")),
    _self_order(getParam<Real>("self_order")),
    _partner1(isParamValid("partner1") ? &getFunctor<Real>("partner1") : nullptr),
    _partner1_order(getParam<Real>("partner1_order")),
    _partner2(isParamValid("partner2") ? &getFunctor<Real>("partner2") : nullptr),
    _partner2_order(getParam<Real>("partner2_order"))
{
  _kinetics.A = getParam<Real>("A");
  _kinetics.k_ref = getParam<Real>("k_ref");
  _kinetics.T_ref = getParam<Real>("T_ref");
  _kinetics.Ea = getParam<Real>("Ea");
}

Real
MSRReaction::partnerProduct(const Moose::ElemArg & elem_arg, const Moose::StateArg & state) const
{
  Real product = 1.0;
  if (_partner1)
    product *= std::pow(std::max((*_partner1)(elem_arg, state), 0.0), _partner1_order);
  if (_partner2)
    product *= std::pow(std::max((*_partner2)(elem_arg, state), 0.0), _partner2_order);
  return product;
}

Real
MSRReaction::computeMatrixContribution()
{
  // Only consumption is treated implicitly on the diagonal.
  if (_mode != Mode::Consumption)
    return 0.0;

  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());

  const Real k = _kinetics.rate(_temperature(elem_arg, state));
  const Real self = std::max(_var.getElemValue(*_current_elem_info, state), 0.0);

  // coeff * [self] reproduces stoich * k * [self]^self_order * partner_product
  const Real self_factor = (_self_order > 1.0) ? std::pow(self, _self_order - 1.0) : 1.0;
  const Real coeff = _stoichiometric_coefficient * k * partnerProduct(elem_arg, state) * self_factor;

  return coeff * _current_elem_volume;
}

Real
MSRReaction::computeRightHandSideContribution()
{
  // Only production is treated explicitly on the right-hand side.
  if (_mode != Mode::Production)
    return 0.0;

  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());

  const Real k = _kinetics.rate(_temperature(elem_arg, state));
  const Real source = _stoichiometric_coefficient * k * partnerProduct(elem_arg, state);

  return source * _current_elem_volume;
}
