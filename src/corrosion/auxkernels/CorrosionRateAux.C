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

#include "CorrosionRateAux.h"
#include "ButlerVolmerKinetics.h"
#include "CorrosionChemistry.h"

registerMooseObject("OpenPronghornApp", CorrosionRateAux);

InputParameters
CorrosionRateAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Butler-Volmer interfacial current density [A/m^2], penetration rate "
                             "[um/y], or accumulated recession/plating thickness [um].");

  MooseEnum mode("current penetration_rate recession", "penetration_rate");
  params.addParam<MooseEnum>("mode", mode, "The diagnostic quantity to output.");

  params.addRequiredCoupledVar("concentration", "Salt-side cation concentration [mol/m^3].");
  params.addRequiredParam<Real>("valence", "Charge number z of the cation.");
  params.addRequiredParam<Real>("molar_mass", "Molar mass of the metal [g/mol].");
  params.addRequiredParam<Real>("density", "Metal density [g/cm^3].");
  params.addRequiredParam<Real>("exchange_current_density", "Exchange current density i0 [A/m^2].");
  params.addParam<Real>("E0", 0.0, "Standard electrode potential [V].");
  params.addParam<Real>("alpha_a", 0.5, "Anodic charge-transfer coefficient.");
  params.addParam<Real>("alpha_c", 0.5, "Cathodic charge-transfer coefficient.");
  params.addRequiredParam<Real>("temperature", "Temperature [K].");
  params.addParam<Real>("c_ref", 1.0, "Reference concentration for the Nernst term [mol/m^3].");

  params.addCoupledVar("salt_potential", "Salt-phase electric potential [V].");
  params.addCoupledVar("metal_potential", "Metal-phase electric potential [V].");
  params.addParam<Real>(
      "salt_potential_value", 0.0, "Fixed salt-phase potential [V] when not coupled.");
  params.addParam<Real>(
      "metal_potential_value", 0.0, "Fixed metal-phase potential [V] when not coupled.");
  return params;
}

CorrosionRateAux::CorrosionRateAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _mode(getParam<MooseEnum>("mode") == "current"
              ? Mode::Current
              : (getParam<MooseEnum>("mode") == "recession" ? Mode::Recession
                                                            : Mode::PenetrationRate)),
    _E0(getParam<Real>("E0")),
    _valence(getParam<Real>("valence")),
    _molar_mass(getParam<Real>("molar_mass")),
    _density(getParam<Real>("density")),
    _i0(getParam<Real>("exchange_current_density")),
    _alpha_a(getParam<Real>("alpha_a")),
    _alpha_c(getParam<Real>("alpha_c")),
    _temperature(getParam<Real>("temperature")),
    _c_ref(getParam<Real>("c_ref")),
    _concentration(coupledValue("concentration")),
    _has_salt_potential(isCoupled("salt_potential")),
    _salt_potential_var(_has_salt_potential ? coupledValue("salt_potential") : _zero),
    _salt_potential(getParam<Real>("salt_potential_value")),
    _has_metal_potential(isCoupled("metal_potential")),
    _metal_potential_var(_has_metal_potential ? coupledValue("metal_potential") : _zero),
    _metal_potential(getParam<Real>("metal_potential_value")),
    _u_old(uOld())
{
  if (_mode == Mode::Recession && !isNodal())
    paramError("mode", "Recession accumulation requires a nodal auxiliary variable.");
}

Real
CorrosionRateAux::computeValue()
{
  const Real phi_salt = _has_salt_potential ? _salt_potential_var[_qp] : _salt_potential;
  const Real phi_metal = _has_metal_potential ? _metal_potential_var[_qp] : _metal_potential;
  const Real c_ion = _concentration[_qp];

  const Real i_bv = Corrosion::bvCurrent(
      phi_metal, phi_salt, c_ion, _E0, _valence, _i0, _alpha_a, _alpha_c, _temperature, _c_ref);

  if (_mode == Mode::Current)
    return i_bv;

  // Convert the current density [A/m^2] to a penetration rate [um/y] (the conversion expects
  // A/cm^2, so divide by 1e4).
  const Real rate_um_y =
      Corrosion::corrosionCurrentToUmY(i_bv * 1.0e-4, _valence, _molar_mass, _density);
  if (_mode == Mode::PenetrationRate)
    return rate_um_y;

  // Recession / plating thickness [um], accumulated over time.
  const Real dt_years = _dt / Corrosion::seconds_per_year;
  return _u_old[_qp] + rate_um_y * dt_years;
}
