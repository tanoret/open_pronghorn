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

#include "NTPPackedBedSegregatedFunctorMaterial.h"
#include "NTPPackedBedCorrelations.h"
#include "NS.h"

registerMooseObject("OpenPronghornApp", NTPPackedBedSegregatedFunctorMaterial);

InputParameters
NTPPackedBedSegregatedFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription("Packed-bed drag, interphase heat-transfer and particle-conduction "
                             "closures of a particle-bed region, for the segregated "
                             "(SIMPLE/PIMPLE) solvers.");
  params.addRequiredRangeCheckedParam<Real>(
      "porosity", "porosity > 0 & porosity < 1", "Bed void fraction [-]");
  params.addRequiredRangeCheckedParam<Real>(
      "particle_diameter", "particle_diameter > 0", "Particle diameter [m]");
  params.addParam<MooseEnum>(
      "pressure_drop_model", MooseEnum("kta ergun", "kta"), "Packed-bed pressure-drop closure");
  params.addParam<MooseEnum>("heat_transfer_model",
                             MooseEnum("achenbach kta", "achenbach"),
                             "Particle-to-fluid Nusselt-number correlation");
  params.addParam<Real>("pressure_drop_multiplier", 1, "Multiplier on the pressure gradient");
  params.addParam<MooseFunctorName>(
      "orifice_multiplier", 1, "Spatially varying multiplier on the drag (orifice shaping) [-]");
  params.addParam<Real>(
      "heat_transfer_multiplier", 1, "Multiplier on the heat-transfer coefficient");
  params.addParam<Real>("surface_area_multiplier", 1, "Multiplier on the specific surface area");
  params.addRequiredParam<Real>("solid_conductivity", "Particle material conductivity [W/(m K)]");
  params.addParam<Real>("axial_conductivity_multiplier", 1, "Multiplier on the axial conductivity");
  params.addParam<Real>(
      "radial_conductivity_multiplier", 1, "Multiplier on the radial conductivity");
  params.addParam<Real>(
      "minimum_reynolds", 1e-3, "Floor on the Reynolds numbers used in the correlations");

  params.addParam<MooseFunctorName>(NS::density, NS::density, "Fluid density");
  params.addParam<MooseFunctorName>(NS::mu, NS::mu, "Fluid dynamic viscosity");
  params.addParam<MooseFunctorName>(NS::k, NS::k, "Fluid thermal conductivity");
  params.addParam<MooseFunctorName>(NS::cp, NS::cp, "Fluid isobaric specific heat");
  params.addParam<MooseFunctorName>("u", "vel_x", "Superficial radial velocity");
  params.addParam<MooseFunctorName>("v", "vel_y", "Superficial axial velocity");
  params.addParam<MooseFunctorName>(
      "mass_flux_x",
      "Radial superficial mass flux, e.g. from NTPRhieChowMassFluxAux; with 'mass_flux_y', the "
      "superficial speed is |G| / rho instead of the cell velocity");
  params.addParam<MooseFunctorName>("mass_flux_y", "Axial superficial mass flux");
  params.addParam<MooseFunctorName>("mass_flux_z", "Third superficial mass flux component (3-D)");

  params.addParam<MooseFunctorName>(
      "drag_name", "darcy_drag", "Name of the Darcy coefficient vector");
  params.addParam<MooseFunctorName>(
      "alpha_name", "alpha", "Name of the interphase coefficient h a_s [W/(m^3 K)]");
  params.addParam<MooseFunctorName>(
      "kappa_name", "kappa_solid", "Name of the particle conductivity vector");
  params.addParam<MooseFunctorName>(
      "speed_name", "superficial_speed", "Name of the superficial speed |u| [m/s]");
  return params;
}

NTPPackedBedSegregatedFunctorMaterial::NTPPackedBedSegregatedFunctorMaterial(
    const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _eps(getParam<Real>("porosity")),
    _d(getParam<Real>("particle_diameter")),
    _pressure_drop_model(getParam<MooseEnum>("pressure_drop_model")),
    _heat_transfer_model(getParam<MooseEnum>("heat_transfer_model")),
    _dp_mult(getParam<Real>("pressure_drop_multiplier")),
    _htc_mult(getParam<Real>("heat_transfer_multiplier")),
    _area_mult(getParam<Real>("surface_area_multiplier")),
    _k_solid(getParam<Real>("solid_conductivity")),
    _k_axial_mult(getParam<Real>("axial_conductivity_multiplier")),
    _k_radial_mult(getParam<Real>("radial_conductivity_multiplier")),
    _Re_floor(getParam<Real>("minimum_reynolds")),
    _rho(getFunctor<Real>(NS::density)),
    _mu(getFunctor<Real>(NS::mu)),
    _k_fluid(getFunctor<Real>(NS::k)),
    _cp(getFunctor<Real>(NS::cp)),
    _u(getFunctor<Real>("u")),
    _v(getFunctor<Real>("v")),
    _orifice(getFunctor<Real>("orifice_multiplier")),
    _Gx(isParamValid("mass_flux_x") ? &getFunctor<Real>("mass_flux_x") : nullptr),
    _Gy(isParamValid("mass_flux_y") ? &getFunctor<Real>("mass_flux_y") : nullptr),
    _Gz(isParamValid("mass_flux_z") ? &getFunctor<Real>("mass_flux_z") : nullptr)
{
  if (bool(_Gx) != bool(_Gy))
    paramError("mass_flux_x", "'mass_flux_x' and 'mass_flux_y' go together");
  using NTPPackedBed::HeatTransfer;
  using NTPPackedBed::PressureDrop;
  const auto dp_model = _pressure_drop_model == "kta" ? PressureDrop::KTA : PressureDrop::ERGUN;
  const auto ht_model =
      _heat_transfer_model == "achenbach" ? HeatTransfer::ACHENBACH : HeatTransfer::KTA;

  addFunctorProperty<Real>(getParam<MooseFunctorName>("speed_name"),
                           [this](const auto & r, const auto & t) -> Real { return speed(r, t); });

  addFunctorProperty<RealVectorValue>(
      getParam<MooseFunctorName>("drag_name"),
      [this, dp_model](const auto & r, const auto & t) -> RealVectorValue
      {
        const Real rho = _rho(r, t);
        const Real mu = _mu(r, t);
        const Real U = speed(r, t);
        const Real A = NTPPackedBed::viscousCoefficient(dp_model, mu, _eps, _d);
        const Real B = NTPPackedBed::inertialCoefficient(dp_model, rho, mu, U, _eps, _d, _Re_floor);
        const Real D = _dp_mult * _orifice(r, t) * (A + B * U) / mu;
        return RealVectorValue(D, D, D);
      });

  const Real a_s = _area_mult * NTPPackedBed::specificSurface(_eps, _d);
  addFunctorProperty<Real>(getParam<MooseFunctorName>("alpha_name"),
                           [this, ht_model, a_s](const auto & r, const auto & t) -> Real
                           {
                             const Real Nu = NTPPackedBed::nusselt(ht_model,
                                                                   _rho(r, t),
                                                                   _mu(r, t),
                                                                   _k_fluid(r, t),
                                                                   _cp(r, t),
                                                                   speed(r, t),
                                                                   _eps,
                                                                   _d,
                                                                   _Re_floor);
                             return _htc_mult * Nu * _k_fluid(r, t) / _d * a_s;
                           });

  addFunctorProperty<RealVectorValue>(
      getParam<MooseFunctorName>("kappa_name"),
      [this](const auto & r, const auto & t) -> RealVectorValue
      {
        const Real k = (1 - _eps) * _k_solid + _eps * _k_fluid(r, t);
        return RealVectorValue(_k_radial_mult * k, _k_axial_mult * k, _k_radial_mult * k);
      });
}

template <typename Space, typename Time>
Real
NTPPackedBedSegregatedFunctorMaterial::speed(const Space & r, const Time & t) const
{
  if (_Gx)
  {
    const Real Gx = (*_Gx)(r, t);
    const Real Gy = (*_Gy)(r, t);
    const Real Gz = _Gz ? (*_Gz)(r, t) : 0;
    return std::sqrt(Gx * Gx + Gy * Gy + Gz * Gz) / _rho(r, t);
  }
  const Real u = _u(r, t);
  const Real v = _v(r, t);
  return std::sqrt(u * u + v * v);
}
