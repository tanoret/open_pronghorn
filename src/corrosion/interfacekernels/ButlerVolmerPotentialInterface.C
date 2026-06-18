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

#include "ButlerVolmerPotentialInterface.h"
#include "ButlerVolmerKinetics.h"
#include "CorrosionChemistry.h"

registerMooseObject("OpenPronghornApp", ButlerVolmerPotentialInterface);

InputParameters
ButlerVolmerPotentialInterface::validParams()
{
  InputParameters params = ADInterfaceKernel::validParams();
  params.addClassDescription("Butler-Volmer charge flux coupling the interfacial current into the "
                             "potential equations of the salt and metal blocks.");

  params.addRequiredCoupledVar("concentration",
                               "Salt-side concentration of the reacting cation on this block "
                               "[mol/m^3].");
  params.addRequiredParam<Real>("valence", "Charge number z of the cation.");
  params.addRequiredParam<Real>("exchange_current_density", "Exchange current density i0 [A/m^2].");
  params.addParam<Real>("E0", 0.0, "Standard electrode potential [V].");
  params.addParam<Real>("alpha_a", 0.5, "Anodic charge-transfer coefficient.");
  params.addParam<Real>("alpha_c", 0.5, "Cathodic charge-transfer coefficient.");
  params.addRequiredParam<Real>("temperature", "Temperature [K].");
  params.addParam<Real>("c_ref", 1.0, "Reference concentration for the Nernst term [mol/m^3].");
  return params;
}

ButlerVolmerPotentialInterface::ButlerVolmerPotentialInterface(const InputParameters & parameters)
  : ADInterfaceKernel(parameters),
    _E0(getParam<Real>("E0")),
    _valence(getParam<Real>("valence")),
    _i0(getParam<Real>("exchange_current_density")),
    _alpha_a(getParam<Real>("alpha_a")),
    _alpha_c(getParam<Real>("alpha_c")),
    _temperature(getParam<Real>("temperature")),
    _c_ref(getParam<Real>("c_ref")),
    _concentration(adCoupledValue("concentration"))
{
}

ADReal
ButlerVolmerPotentialInterface::computeQpResidual(Moose::DGResidualType type)
{
  // This-side variable is the salt potential; neighbor variable is the metal potential.
  const ADReal phi_salt = _u[_qp];
  const ADReal phi_metal = _neighbor_value[_qp];
  const ADReal c_ion = _concentration[_qp];

  const ADReal i_bv = Corrosion::bvCurrent(
      phi_metal, phi_salt, c_ion, _E0, _valence, _i0, _alpha_a, _alpha_c, _temperature, _c_ref);

  switch (type)
  {
    case Moose::Element:
      // Current enters the salt phase.
      return -_test[_i][_qp] * i_bv;

    case Moose::Neighbor:
      // Current leaves the metal phase.
      return _test_neighbor[_i][_qp] * i_bv;
  }
  return 0.0;
}
