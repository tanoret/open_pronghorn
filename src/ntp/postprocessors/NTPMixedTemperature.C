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

#include "NTPMixedTemperature.h"

#include <cmath>
#include "SinglePhaseFluidProperties.h"

registerMooseObject("OpenPronghornApp", NTPMixedTemperature);

InputParameters
NTPMixedTemperature::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription("Mixed-mean temperature T(p, H/mdot) of a stream from its enthalpy "
                             "and mass flow rates.");
  params.addRequiredParam<UserObjectName>("fp", "Fluid properties");
  params.addRequiredParam<PostprocessorName>("enthalpy_flow", "Enthalpy flow rate H [W]");
  params.addRequiredParam<PostprocessorName>("mass_flow", "Mass flow rate [kg/s]");
  params.addParam<PostprocessorName>("pressure", 1e6, "Pressure for the inversion [Pa]");
  return params;
}

NTPMixedTemperature::NTPMixedTemperature(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _fp(getUserObject<SinglePhaseFluidProperties>("fp")),
    _H(getPostprocessorValue("enthalpy_flow")),
    _mdot(getPostprocessorValue("mass_flow")),
    _p(getPostprocessorValue("pressure"))
{
}

Real
NTPMixedTemperature::getValue() const
{
  // no defined mixed temperature without a through-flow (e.g. at the initial state)
  if (std::abs(_mdot) < 1e-12 * std::max(std::abs(_H), 1.0))
    return 0;
  return _fp.T_from_p_h(_p, _H / _mdot);
}
