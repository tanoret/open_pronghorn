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

#include "NTPGasMassTimeDerivative.h"
#include "SinglePhaseFluidProperties.h"
#include "NS.h"

registerMooseObject("OpenPronghornApp", NTPGasMassTimeDerivative);

InputParameters
NTPGasMassTimeDerivative::validParams()
{
  InputParameters params = FVTimeKernel::validParams();
  params.addClassDescription("Gas mass storage eps d(rho(p, T))/dt on the pressure equation.");
  params.addRequiredParam<UserObjectName>(NS::fluid, "Fluid properties");
  params.addRequiredParam<MooseFunctorName>(NS::T_fluid, "Gas temperature");
  params.addRequiredParam<MooseFunctorName>(NS::porosity, "Porosity");
  return params;
}

NTPGasMassTimeDerivative::NTPGasMassTimeDerivative(const InputParameters & parameters)
  : FVTimeKernel(parameters),
    _fp(getUserObject<SinglePhaseFluidProperties>(NS::fluid)),
    _T(getFunctor<ADReal>(NS::T_fluid)),
    _eps(getFunctor<ADReal>(NS::porosity))
{
}

ADReal
NTPGasMassTimeDerivative::computeQpResidual()
{
  const auto elem = makeElemArg(_current_elem);
  const auto state = determineState();
  ADReal rho, drho_dp, drho_dT;
  _fp.rho_from_p_T(_var(elem, state), _T(elem, state), rho, drho_dp, drho_dT);
  return _eps(elem, state) * (drho_dp * _u_dot[_qp] + drho_dT * _T.dot(elem, state));
}
