//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Copyright 2025, Battelle Energy Alliance, LLC
//* ALL RIGHTS RESERVED
//*

#include "NTPDuctFrictionSegregatedFunctorMaterial.h"
#include "NS.h"

registerMooseObject("OpenPronghornApp", NTPDuctFrictionSegregatedFunctorMaterial);

InputParameters
NTPDuctFrictionSegregatedFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription("Axial duct wall friction of a free channel of a particle-bed element "
                             "as a Darcy coefficient for the segregated solvers.");
  params.addParam<MooseFunctorName>(NS::density, NS::density, "Fluid density");
  params.addParam<MooseFunctorName>(NS::mu, NS::mu, "Fluid dynamic viscosity");
  params.addParam<MooseFunctorName>("v", "vel_y", "Axial velocity");
  params.addParam<MooseFunctorName>(
      "speed", "Axial speed |u_z| of the turbulent term, instead of the one of the cell velocity");
  params.addRequiredParam<MooseFunctorName>("hydraulic_diameter", "Hydraulic diameter [m]");
  params.addRequiredRangeCheckedParam<Real>(
      "friction_factor", "friction_factor >= 0", "Darcy-Weisbach friction factor f [-]");
  params.addParam<bool>("laminar_term", true, "Include the Hagen-Poiseuille term 32 mu / D_h^2");
  params.addParam<MooseFunctorName>(
      "drag_name", "duct_drag", "Name of the Darcy coefficient vector");
  return params;
}

NTPDuctFrictionSegregatedFunctorMaterial::NTPDuctFrictionSegregatedFunctorMaterial(
    const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _rho(getFunctor<Real>(NS::density)),
    _mu(getFunctor<Real>(NS::mu)),
    _v(getFunctor<Real>("v")),
    _speed(isParamValid("speed") ? &getFunctor<Real>("speed") : nullptr),
    _Dh(getFunctor<Real>("hydraulic_diameter")),
    _f(getParam<Real>("friction_factor")),
    _laminar(getParam<bool>("laminar_term"))
{
  addFunctorProperty<RealVectorValue>(getParam<MooseFunctorName>("drag_name"),
                                      [this](const auto & r, const auto & t) -> RealVectorValue
                                      {
                                        const Real mu = _mu(r, t);
                                        const Real Dh = _Dh(r, t);
                                        const Real speed =
                                            _speed ? (*_speed)(r, t) : std::abs(_v(r, t));
                                        const Real a = _laminar ? 32 * mu / (Dh * Dh) : 0;
                                        const Real b = _f * _rho(r, t) / (2 * Dh);
                                        return RealVectorValue(0, (a + b * speed) / mu, 0);
                                      });
}
