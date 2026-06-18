//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MSRGasExchange.h"
#include "MoltenSaltRadiolysisData.h"

registerMooseObject("OpenPronghornApp", MSRGasExchange);

InputParameters
MSRGasExchange::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription("Adds the gas-liquid exchange of a dissolved diatomic between the melt "
                             "and a well-mixed headspace to a linear finite-volume system.");

  MooseEnum mode("liquid gas");
  params.addRequiredParam<MooseEnum>(
      "mode", mode, "Whether this instance assembles the liquid or the gas-phase equation.");

  params.addRequiredParam<MooseFunctorName>("temperature", "The temperature functor [K].");
  params.addRequiredParam<MooseFunctorName>(
      "partner",
      "The partner concentration: the gas-phase concentration in liquid mode, the dissolved "
      "concentration in gas mode.");

  params.addRequiredParam<Real>("kLa", "Overall mass-transfer coefficient kLa [1/s].");
  params.addRequiredParam<Real>("kH", "Henry coefficient kH [mol/(m^3 Pa)].");
  params.addParam<Real>(
      "volume_ratio", 1.0, "Ratio of liquid to headspace volume V_liq / V_gas (gas mode).");

  return params;
}

MSRGasExchange::MSRGasExchange(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _mode(getParam<MooseEnum>("mode") == "liquid" ? Mode::Liquid : Mode::Gas),
    _temperature(getFunctor<Real>("temperature")),
    _kLa(getParam<Real>("kLa")),
    _kH(getParam<Real>("kH")),
    _volume_ratio(getParam<Real>("volume_ratio")),
    _partner(getFunctor<Real>("partner"))
{
}

Real
MSRGasExchange::computeMatrixContribution()
{
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());

  // Liquid equation: dC_liq/dt + kLa * C_liq = kLa * kH * R T * C_gas
  if (_mode == Mode::Liquid)
    return _kLa * _current_elem_volume;

  // Gas equation: dC_gas/dt + (V_liq/V_gas) kLa kH R T * C_gas = (V_liq/V_gas) kLa * C_liq
  const Real T = _temperature(elem_arg, state);
  return _volume_ratio * _kLa * _kH * MSR::R_gas * T * _current_elem_volume;
}

Real
MSRGasExchange::computeRightHandSideContribution()
{
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());

  const Real T = _temperature(elem_arg, state);
  const Real partner = std::max(_partner(elem_arg, state), 0.0);

  // Liquid equation: explicit equilibrium-restoring source from the lagged gas concentration
  if (_mode == Mode::Liquid)
    return _kLa * _kH * MSR::R_gas * T * partner * _current_elem_volume;

  // Gas equation: explicit gain from the lagged dissolved concentration
  return _volume_ratio * _kLa * partner * _current_elem_volume;
}
