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

#include "ButlerVolmerBC.h"
#include "ButlerVolmerKinetics.h"
#include "CorrosionChemistry.h"

registerMooseObject("OpenPronghornApp", ButlerVolmerBC);

InputParameters
ButlerVolmerBC::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("Butler-Volmer electrode reaction as a boundary condition for a "
                             "single-domain (salt-only or solid-only) corrosion/plating problem.");

  MooseEnum flux_type("species charge", "species");
  params.addParam<MooseEnum>(
      "flux_type",
      flux_type,
      "Whether the boundary adds the molar flux i_BV/(zF) to a concentration variable ('species') "
      "or the current i_BV to a potential variable ('charge').");
  params.addParam<bool>("metal_domain",
                        false,
                        "True when the modeled phase is the solid metal, false when it is the salt.");

  params.addRequiredParam<Real>("valence", "Charge number z of the cation.");
  params.addRequiredParam<Real>("exchange_current_density", "Exchange current density i0 [A/m^2].");
  params.addParam<Real>("E0", 0.0, "Standard electrode potential [V].");
  params.addParam<Real>("alpha_a", 0.5, "Anodic charge-transfer coefficient.");
  params.addParam<Real>("alpha_c", 0.5, "Cathodic charge-transfer coefficient.");
  params.addRequiredParam<Real>("temperature", "Temperature [K].");
  params.addParam<Real>("c_ref", 1.0, "Reference concentration for the Nernst term [mol/m^3].");

  params.addCoupledVar("potential", "The modeled-phase electric potential [V], if solved.");
  params.addParam<MooseFunctorName>(
      "applied_potential", "0", "Potential of the absent counter-phase [V].");
  params.addCoupledVar("concentration",
                       "Salt-ion concentration [mol/m^3] for a charge boundary on the salt.");
  params.addParam<MooseFunctorName>(
      "counter_concentration",
      "Fixed salt-ion concentration [mol/m^3] of the absent salt phase (solid-domain case).");
  return params;
}

ButlerVolmerBC::ButlerVolmerBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _charge_flux(getParam<MooseEnum>("flux_type") == "charge"),
    _metal_domain(getParam<bool>("metal_domain")),
    _E0(getParam<Real>("E0")),
    _valence(getParam<Real>("valence")),
    _i0(getParam<Real>("exchange_current_density")),
    _alpha_a(getParam<Real>("alpha_a")),
    _alpha_c(getParam<Real>("alpha_c")),
    _temperature(getParam<Real>("temperature")),
    _c_ref(getParam<Real>("c_ref")),
    _has_potential(isCoupled("potential")),
    _potential(_has_potential ? &adCoupledValue("potential") : nullptr),
    _applied_potential(getFunctor<Real>("applied_potential")),
    _has_concentration(isCoupled("concentration")),
    _concentration(_has_concentration ? &adCoupledValue("concentration") : nullptr),
    _counter_concentration(isParamValid("counter_concentration")
                               ? &getFunctor<Real>("counter_concentration")
                               : nullptr)
{
  // A species flux on the salt reads the concentration from the boundary variable itself; any other
  // combination needs an explicit salt-ion concentration source.
  const bool concentration_is_variable = !_charge_flux && !_metal_domain;
  if (!concentration_is_variable && !_has_concentration && !_counter_concentration)
    paramError("concentration",
               "Couple 'concentration' or supply 'counter_concentration' so the Butler-Volmer "
               "overpotential has a salt-ion concentration to use.");
}

ADReal
ButlerVolmerBC::computeQpResidual()
{
  const Moose::ElemQpArg qp_arg = {_current_elem, _qp, _qrule, _q_point[_qp]};
  const auto state = Moose::currentState();

  // Salt-ion concentration driving the Nernst term.
  ADReal c_ion;
  if (!_charge_flux && !_metal_domain)
    c_ion = _u[_qp];
  else if (_has_concentration)
    c_ion = (*_concentration)[_qp];
  else
    c_ion = ADReal((*_counter_concentration)(qp_arg, state));

  // Phase potentials: the counter phase is the applied functor. The modeled-phase potential is the
  // boundary variable itself for a charge boundary, or the optionally coupled 'potential' otherwise.
  const ADReal applied = ADReal(_applied_potential(qp_arg, state));
  ADReal modeled;
  if (_charge_flux)
    modeled = _u[_qp];
  else
    modeled = _has_potential ? (*_potential)[_qp] : ADReal(0.0);
  const ADReal phi_metal = _metal_domain ? modeled : applied;
  const ADReal phi_salt = _metal_domain ? applied : modeled;

  const ADReal i_bv = Corrosion::bvCurrent(
      phi_metal, phi_salt, c_ion, _E0, _valence, _i0, _alpha_a, _alpha_c, _temperature, _c_ref);
  const ADReal flux = _charge_flux ? i_bv : i_bv / (_valence * Corrosion::faraday);

  // Salt domain gains the cation / current (source, -test); solid domain loses metal (sink, +test).
  const Real sign = _metal_domain ? 1.0 : -1.0;
  return sign * _test[_i][_qp] * flux;
}
