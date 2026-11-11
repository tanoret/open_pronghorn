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

#include "NTPPackedBedFunctorMaterial.h"
#include "NS.h"
#include "NTPPackedBedCorrelations.h"

registerMooseObject("OpenPronghornApp", NTPPackedBedFunctorMaterial);

InputParameters
NTPPackedBedFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription("Darcy-Forchheimer mobility, interphase heat-transfer and "
                             "solid-conduction closures of one packed-bed region of a "
                             "particle-bed fuel element.");
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
      "orifice_multiplier",
      1,
      "Spatially varying multiplier on the pressure gradient, for orifice shaping [-]");
  params.addParam<Real>(
      "heat_transfer_multiplier", 1, "Multiplier on the heat-transfer coefficient");
  params.addParam<Real>("surface_area_multiplier", 1, "Multiplier on the specific surface area");
  params.addRequiredParam<Real>("solid_conductivity", "Particle material conductivity [W/(m K)]");
  params.addParam<Real>("axial_conductivity_multiplier", 1, "Multiplier on the axial conductivity");
  params.addParam<Real>(
      "radial_conductivity_multiplier", 1, "Multiplier on the radial conductivity");
  params.addRangeCheckedParam<Real>(
      "particle_emissivity",
      0,
      "particle_emissivity >= 0 & particle_emissivity < 2",
      "Emissivity of the particles; > 0 adds the radiative conductivity of the bed "
      "4 sigma e / (2 - e) d T_s^3 to the particle conductivity");
  params.addParam<MooseFunctorName>(
      "solid_temperature", "T_solid", "Particle temperature (for the radiative conductivity)");
  params.addRangeCheckedParam<Real>(
      "dispersion_coefficient",
      0,
      "dispersion_coefficient >= 0",
      "Coefficient C of the gas thermal dispersion conductivity C Pe k = C |G| cp d (about 0.1 "
      "transverse, 0.5 longitudinal); 0 switches it off");
  params.addParam<MooseFunctorName>(
      "dispersion_name", "k_dispersion", "Name of the gas thermal dispersion conductivity");
  params.addParam<Real>(
      "minimum_reynolds", 1e-3, "Floor on the Reynolds numbers used in the correlations");

  params.addParam<MooseFunctorName>(NS::density, NS::density, "Fluid density");
  params.addParam<MooseFunctorName>(NS::mu, NS::mu, "Fluid dynamic viscosity");
  params.addParam<MooseFunctorName>(NS::k, NS::k, "Fluid thermal conductivity");
  params.addParam<MooseFunctorName>(NS::cp, NS::cp, "Fluid isobaric specific heat");
  params.addParam<MooseFunctorName>(NS::pressure, NS::pressure, "Pressure");
  params.addParam<MooseFunctorName>(
      "mass_flux_x",
      "Radial superficial mass flux rho u_r [kg/(m^2 s)], e.g. from "
      "NTPDarcyMassFluxFunctorMaterial; with "
      "mass_flux_y it sets the speed of the heat-transfer correlation (recommended)");
  params.addParam<MooseFunctorName>("mass_flux_y", "Axial superficial mass flux rho u_z");
  params.addParam<MooseFunctorName>("mass_flux_z", "Third superficial mass flux component (3-D)");

  params.addParam<MooseFunctorName>(
      "rho_mobility_name", "rho_mobility", "Name of the mass mobility vector rho M [s]");
  params.addParam<MooseFunctorName>(
      "speed_name", "superficial_speed", "Name of the superficial speed |u| [m/s]");
  params.addParam<MooseFunctorName>("velocity_name",
                                    "superficial_velocity",
                                    "Name of the superficial velocity vector u = -M grad p [m/s]");
  params.addParam<MooseFunctorName>(
      "alpha_name", "alpha", "Name of the interphase coefficient h a_s [W/(m^3 K)]");
  params.addParam<MooseFunctorName>(
      "kappa_name", "kappa_solid", "Name of the solid effective conductivity vector");
  params.addParam<MooseFunctorName>("nusselt_name", "Nu", "Name of the Nusselt number");
  return params;
}

NTPPackedBedFunctorMaterial::NTPPackedBedFunctorMaterial(const InputParameters & parameters)
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
    _emissivity(getParam<Real>("particle_emissivity")),
    _T_solid(_emissivity > 0 ? &getFunctor<ADReal>("solid_temperature") : nullptr),
    _dispersion(getParam<Real>("dispersion_coefficient")),
    _Re_floor(getParam<Real>("minimum_reynolds")),
    _rho(getFunctor<ADReal>(NS::density)),
    _mu(getFunctor<ADReal>(NS::mu)),
    _k_fluid(getFunctor<ADReal>(NS::k)),
    _cp(getFunctor<ADReal>(NS::cp)),
    _p(getFunctor<ADReal>(NS::pressure)),
    _orifice(getFunctor<ADReal>("orifice_multiplier")),
    _G_x(isParamValid("mass_flux_x") ? &getFunctor<ADReal>("mass_flux_x") : nullptr),
    _G_y(isParamValid("mass_flux_y") ? &getFunctor<ADReal>("mass_flux_y") : nullptr),
    _G_z(isParamValid("mass_flux_z") ? &getFunctor<ADReal>("mass_flux_z") : nullptr)
{
  if (bool(_G_x) != bool(_G_y))
    paramError("mass_flux_x", "'mass_flux_x' and 'mass_flux_y' go together");
  addFunctorProperty<ADRealVectorValue>(getParam<MooseFunctorName>("rho_mobility_name"),
                                        [this](const auto & r, const auto & t) -> ADRealVectorValue
                                        {
                                          ADReal M, speed;
                                          mobility(r, t, M, speed);
                                          const auto rhoM = _rho(r, t) * M;
                                          return ADRealVectorValue(rhoM, rhoM, rhoM);
                                        });

  addFunctorProperty<ADRealVectorValue>(getParam<MooseFunctorName>("velocity_name"),
                                        [this](const auto & r, const auto & t) -> ADRealVectorValue
                                        {
                                          ADReal M, speed;
                                          mobility(r, t, M, speed);
                                          return -M * _p.gradient(r, t);
                                        });

  // Speed of the heat-transfer correlation: from the reconstructed cell mass flux when given (exact
  // also in porous layers one cell thick), else from the cell pressure gradient
  const auto & speed =
      addFunctorProperty<ADReal>(getParam<MooseFunctorName>("speed_name"),
                                 [this](const auto & r, const auto & t) -> ADReal
                                 {
                                   if (_G_x)
                                   {
                                     using std::sqrt;
                                     const auto Gx = (*_G_x)(r, t);
                                     const auto Gy = (*_G_y)(r, t);
                                     const ADReal Gz = _G_z ? (*_G_z)(r, t) : ADReal(0);
                                     return sqrt(Gx * Gx + Gy * Gy + Gz * Gz + 1e-20) / _rho(r, t);
                                   }
                                   ADReal M, speed;
                                   mobility(r, t, M, speed);
                                   return speed;
                                 });

  const auto & Nu = addFunctorProperty<ADReal>(
      getParam<MooseFunctorName>("nusselt_name"),
      [this, &speed](const auto & r, const auto & t) -> ADReal
      {
        return NTPPackedBed::nusselt(_heat_transfer_model == "achenbach"
                                         ? NTPPackedBed::HeatTransfer::ACHENBACH
                                         : NTPPackedBed::HeatTransfer::KTA,
                                     _rho(r, t),
                                     _mu(r, t),
                                     _k_fluid(r, t),
                                     _cp(r, t),
                                     speed(r, t),
                                     _eps,
                                     _d,
                                     _Re_floor);
      });

  const Real a_s = _area_mult * NTPPackedBed::specificSurface(_eps, _d);
  addFunctorProperty<ADReal>(getParam<MooseFunctorName>("alpha_name"),
                             [this, &Nu, a_s](const auto & r, const auto & t) -> ADReal
                             { return _htc_mult * Nu(r, t) * _k_fluid(r, t) / _d * a_s; });

  addFunctorProperty<ADRealVectorValue>(
      getParam<MooseFunctorName>("kappa_name"),
      [this](const auto & r, const auto & t) -> ADRealVectorValue
      {
        ADReal k = (1 - _eps) * _k_solid + _eps * _k_fluid(r, t);
        if (_T_solid)
          k += NTPPackedBed::radiativeConductivity((*_T_solid)(r, t), _emissivity, _d);
        return ADRealVectorValue(_k_radial_mult * k, _k_axial_mult * k, _k_radial_mult * k);
      });

  // gas thermal dispersion, C |G| cp d
  addFunctorProperty<ADReal>(getParam<MooseFunctorName>("dispersion_name"),
                             [this, &speed](const auto & r, const auto & t) -> ADReal
                             {
                               if (_dispersion == 0)
                                 return ADReal(0);
                               return NTPPackedBed::dispersionConductivity(
                                   speed(r, t) * _rho(r, t), _cp(r, t), _d, _dispersion);
                             });
}

template <typename Space, typename Time>
void
NTPPackedBedFunctorMaterial::mobility(const Space & r,
                                      const Time & t,
                                      ADReal & M,
                                      ADReal & speed) const
{
  using std::sqrt;

  const Real eps = _eps;
  const Real d = _d;
  const auto rho = _rho(r, t);
  const auto mu = _mu(r, t);
  const auto mult = _dp_mult * _orifice(r, t);
  const auto grad_p = _p.gradient(r, t);
  // the small offset keeps the derivative of |grad p| finite where the gradient vanishes
  const auto G = sqrt(grad_p * grad_p + 1e-6);

  using NTPPackedBed::PressureDrop;
  const auto model = _pressure_drop_model == "kta" ? PressureDrop::KTA : PressureDrop::ERGUN;
  const ADReal A = mult * NTPPackedBed::viscousCoefficient(model, mu, eps, d);
  // |u| = M G and -grad p = (A + B |u|) u give B G M^2 + A M - 1 = 0
  const auto solve = [&](const ADReal & B) { return 2 / (A + sqrt(A * A + 4 * B * G)); };

  // Ergun is explicit; the weak Re_m^-0.1 dependence of KTA is a fixed point started from Ergun
  const ADReal zero_speed = 0;
  M = solve(mult * NTPPackedBed::inertialCoefficient(
                       PressureDrop::ERGUN, rho, mu, zero_speed, eps, d, _Re_floor));
  if (model == PressureDrop::KTA)
    for (unsigned int it = 0; it < 12; ++it)
    {
      const ADReal u = M * G;
      M = solve(mult * NTPPackedBed::inertialCoefficient(model, rho, mu, u, eps, d, _Re_floor));
    }
  speed = M * G;
}
