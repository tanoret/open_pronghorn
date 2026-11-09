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

#include "NTPGasPropertiesFunctorMaterial.h"
#include "SinglePhaseFluidProperties.h"
#include "NS.h"

registerMooseObject("OpenPronghornApp", NTPGasPropertiesFunctorMaterial);
registerMooseObject("OpenPronghornApp", ADNTPGasPropertiesFunctorMaterial);

template <bool is_ad>
InputParameters
NTPGasPropertiesFunctorMaterialTempl<is_ad>::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription("Gas density, viscosity, specific heat and thermal conductivity from "
                             "fluid properties at given pressure and temperature functors, under "
                             "chosen names.");
  params.addRequiredParam<UserObjectName>(NS::fluid, "Fluid properties");
  params.addRequiredParam<MooseFunctorName>(NS::pressure, "Pressure functor");
  params.addRequiredParam<MooseFunctorName>("temperature", "Temperature functor");
  params.addParam<MooseFunctorName>("density_name", NS::density, "Name of the density");
  params.addParam<MooseFunctorName>("viscosity_name", NS::mu, "Name of the dynamic viscosity");
  params.addParam<MooseFunctorName>(
      "specific_heat_name", NS::cp, "Name of the isobaric specific heat");
  params.addParam<MooseFunctorName>("conductivity_name", NS::k, "Name of the thermal conductivity");
  return params;
}

template <bool is_ad>
NTPGasPropertiesFunctorMaterialTempl<is_ad>::NTPGasPropertiesFunctorMaterialTempl(
    const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _fp(getUserObject<SinglePhaseFluidProperties>(NS::fluid)),
    _p(getFunctor<GenericReal<is_ad>>(NS::pressure)),
    _T(getFunctor<GenericReal<is_ad>>("temperature"))
{
  using T = GenericReal<is_ad>;
  addFunctorProperty<T>(getParam<MooseFunctorName>("density_name"),
                        [this](const auto & r, const auto & t) -> T
                        { return _fp.rho_from_p_T(_p(r, t), _T(r, t)); });
  addFunctorProperty<T>(getParam<MooseFunctorName>("viscosity_name"),
                        [this](const auto & r, const auto & t) -> T
                        { return _fp.mu_from_p_T(_p(r, t), _T(r, t)); });
  addFunctorProperty<T>(getParam<MooseFunctorName>("specific_heat_name"),
                        [this](const auto & r, const auto & t) -> T
                        { return _fp.cp_from_p_T(_p(r, t), _T(r, t)); });
  addFunctorProperty<T>(getParam<MooseFunctorName>("conductivity_name"),
                        [this](const auto & r, const auto & t) -> T
                        { return _fp.k_from_p_T(_p(r, t), _T(r, t)); });
}

template class NTPGasPropertiesFunctorMaterialTempl<false>;
template class NTPGasPropertiesFunctorMaterialTempl<true>;
