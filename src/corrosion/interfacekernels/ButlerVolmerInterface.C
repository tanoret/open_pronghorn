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

#include "ButlerVolmerInterface.h"
#include "ButlerVolmerKinetics.h"
#include "CorrosionChemistry.h"

registerMooseObject("OpenPronghornApp", ButlerVolmerInterface);

InputParameters
ButlerVolmerInterface::validParams()
{
  InputParameters params = ADInterfaceKernel::validParams();
  params.addClassDescription("Butler-Volmer mass flux of a metal cation across a metal/salt "
                             "interface modeled as two mesh blocks.");

  params.addRequiredParam<Real>("valence", "Charge number z of the cation.");
  params.addRequiredParam<Real>("exchange_current_density", "Exchange current density i0 [A/m^2].");
  params.addParam<Real>("E0", 0.0, "Standard electrode potential [V].");
  params.addParam<Real>("alpha_a", 0.5, "Anodic charge-transfer coefficient.");
  params.addParam<Real>("alpha_c", 0.5, "Cathodic charge-transfer coefficient.");
  params.addRequiredParam<Real>("temperature", "Temperature [K].");
  params.addParam<Real>("c_ref", 1.0, "Reference concentration for the Nernst term [mol/m^3].");

  params.addCoupledVar("phi_salt", "Salt-phase electric potential on this block [V].");
  params.addCoupledVar("phi_solid", "Metal-phase electric potential on the neighbor block [V].");
  params.addParam<Real>(
      "salt_potential", 0.0, "Fixed salt-phase potential [V] used when phi_salt is not coupled.");
  params.addParam<Real>(
      "metal_potential", 0.0, "Fixed metal-phase potential [V] used when phi_solid is not coupled.");
  return params;
}

ButlerVolmerInterface::ButlerVolmerInterface(const InputParameters & parameters)
  : ADInterfaceKernel(parameters),
    _E0(getParam<Real>("E0")),
    _valence(getParam<Real>("valence")),
    _i0(getParam<Real>("exchange_current_density")),
    _alpha_a(getParam<Real>("alpha_a")),
    _alpha_c(getParam<Real>("alpha_c")),
    _temperature(getParam<Real>("temperature")),
    _c_ref(getParam<Real>("c_ref")),
    _has_salt_potential(isCoupled("phi_salt")),
    _phi_salt(_has_salt_potential ? &adCoupledValue("phi_salt") : nullptr),
    _has_metal_potential(isCoupled("phi_solid")),
    _phi_solid(_has_metal_potential ? &adCoupledNeighborValue("phi_solid") : nullptr),
    _salt_potential(getParam<Real>("salt_potential")),
    _metal_potential(getParam<Real>("metal_potential"))
{
}

ADReal
ButlerVolmerInterface::computeQpResidual(Moose::DGResidualType type)
{
  const ADReal phi_salt = _has_salt_potential ? (*_phi_salt)[_qp] : ADReal(_salt_potential);
  const ADReal phi_metal = _has_metal_potential ? (*_phi_solid)[_qp] : ADReal(_metal_potential);
  const ADReal c_ion = _u[_qp];

  const ADReal i_bv = Corrosion::bvCurrent(
      phi_metal, phi_salt, c_ion, _E0, _valence, _i0, _alpha_a, _alpha_c, _temperature, _c_ref);
  // Faradaic molar flux of the cation into the salt [mol/m^2/s].
  const ADReal flux = i_bv / (_valence * Corrosion::faraday);

  switch (type)
  {
    case Moose::Element:
      // Salt side gains the cation (influx source).
      return -_test[_i][_qp] * flux;

    case Moose::Neighbor:
      // Solid side loses metal (sink).
      return _test_neighbor[_i][_qp] * flux;
  }
  return 0.0;
}
