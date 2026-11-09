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

#include "NTPDuctFlowFunctorMaterial.h"
#include "NS.h"

registerMooseObject("OpenPronghornApp", NTPDuctFlowFunctorMaterial);

InputParameters
NTPDuctFlowFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "Darcy-type mobility of a free channel: Hagen-Poiseuille plus Darcy-Weisbach friction "
      "along the axis, a large mobility across it.");
  params.addParam<MooseFunctorName>(NS::density, NS::density, "Fluid density");
  params.addParam<MooseFunctorName>(NS::mu, NS::mu, "Fluid dynamic viscosity");
  params.addParam<MooseFunctorName>(NS::pressure, NS::pressure, "Pressure");
  params.addRequiredParam<MooseFunctorName>("hydraulic_diameter", "Hydraulic diameter [m]");
  params.addRequiredRangeCheckedParam<Real>(
      "friction_factor", "friction_factor >= 0", "Darcy friction factor [-]");
  params.addRangeCheckedParam<Real>("radial_mobility_ratio",
                                    100,
                                    "radial_mobility_ratio > 0",
                                    "Ratio of the radial to the axial mobility [-]");
  params.addParam<MooseFunctorName>(
      "rho_mobility_name", "rho_mobility", "Name of the mass mobility vector rho M [s]");
  params.addParam<MooseFunctorName>(
      "speed_name", "superficial_speed", "Name of the axial speed |u_z| [m/s]");
  params.addParam<MooseFunctorName>(
      "velocity_name", "superficial_velocity", "Name of the velocity vector u = -M grad p [m/s]");
  return params;
}

NTPDuctFlowFunctorMaterial::NTPDuctFlowFunctorMaterial(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _rho(getFunctor<ADReal>(NS::density)),
    _mu(getFunctor<ADReal>(NS::mu)),
    _p(getFunctor<ADReal>(NS::pressure)),
    _Dh(getFunctor<ADReal>("hydraulic_diameter")),
    _f(getParam<Real>("friction_factor")),
    _radial_ratio(getParam<Real>("radial_mobility_ratio"))
{
  const auto axial_mobility = [this](const auto & r, const auto & t) -> ADReal
  {
    using std::abs;
    using std::sqrt;
    const auto Dh = _Dh(r, t);
    const auto a = 32 * _mu(r, t) / (Dh * Dh);
    const auto b = _f * _rho(r, t) / (2 * Dh);
    const auto dpdz = _p.gradient(r, t)(1);
    const auto G = sqrt(dpdz * dpdz + 1e-6);
    return 2 / (a + sqrt(a * a + 4 * b * G));
  };

  addFunctorProperty<ADRealVectorValue>(
      getParam<MooseFunctorName>("rho_mobility_name"),
      [this, axial_mobility](const auto & r, const auto & t) -> ADRealVectorValue
      {
        const auto rhoM = _rho(r, t) * axial_mobility(r, t);
        // y is the channel axis; x (and z in 3-D) are transverse
        return ADRealVectorValue(_radial_ratio * rhoM, rhoM, _radial_ratio * rhoM);
      });

  addFunctorProperty<ADRealVectorValue>(
      getParam<MooseFunctorName>("velocity_name"),
      [this, axial_mobility](const auto & r, const auto & t) -> ADRealVectorValue
      {
        const auto Mz = axial_mobility(r, t);
        const auto grad_p = _p.gradient(r, t);
        return ADRealVectorValue(
            -_radial_ratio * Mz * grad_p(0), -Mz * grad_p(1), -_radial_ratio * Mz * grad_p(2));
      });

  addFunctorProperty<ADReal>(getParam<MooseFunctorName>("speed_name"),
                             [this, axial_mobility](const auto & r, const auto & t) -> ADReal
                             {
                               using std::abs;
                               return axial_mobility(r, t) * abs(_p.gradient(r, t)(1));
                             });
}
