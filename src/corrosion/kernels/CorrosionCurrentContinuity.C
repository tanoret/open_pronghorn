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

#include "CorrosionCurrentContinuity.h"
#include "CorrosionChemistry.h"

registerMooseObject("OpenPronghornApp", CorrosionCurrentContinuity);

InputParameters
CorrosionCurrentContinuity::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addClassDescription("Current-continuity equation div(i) = 0 for the salt electric "
                             "potential, with the ionic conductivity built from the coupled "
                             "concentrations (dilute-solution closure).");

  params.addRequiredCoupledVar("concentrations", "The coupled salt-ion concentrations c_k [mol/m^3].");
  params.addRequiredParam<std::vector<Real>>("valences",
                                             "Charge number z_k of each coupled ion (same order as "
                                             "'concentrations').");
  params.addRequiredParam<std::vector<Real>>("diffusivities",
                                             "Diffusivity D_k of each coupled ion [m^2/s] (same "
                                             "order as 'concentrations').");
  params.addRequiredParam<MooseFunctorName>("temperature", "Temperature functor [K].");
  return params;
}

CorrosionCurrentContinuity::CorrosionCurrentContinuity(const InputParameters & parameters)
  : ADKernelGrad(parameters),
    _n_ions(coupledComponents("concentrations")),
    _temperature(getFunctor<Real>("temperature")),
    _valences(getParam<std::vector<Real>>("valences")),
    _diffusivities(getParam<std::vector<Real>>("diffusivities"))
{
  if (_valences.size() != _n_ions)
    paramError("valences", "Must list one charge number per coupled concentration.");
  if (_diffusivities.size() != _n_ions)
    paramError("diffusivities", "Must list one diffusivity per coupled concentration.");

  for (unsigned int k = 0; k < _n_ions; ++k)
  {
    _c.push_back(&adCoupledValue("concentrations", k));
    _grad_c.push_back(&adCoupledGradient("concentrations", k));
  }
}

ADRealVectorValue
CorrosionCurrentContinuity::precomputeQpResidual()
{
  const Moose::ElemQpArg qp_arg = {_current_elem, _qp, _qrule, _q_point[_qp]};
  const Real temperature = _temperature(qp_arg, Moose::currentState());
  const Real f2_rt = Corrosion::faraday * Corrosion::faraday / (Corrosion::R_gas * temperature);

  // Ionic conductivity kappa = (F^2 / R T) sum_k z_k^2 D_k c_k (kept positive for ellipticity).
  ADReal kappa = 0.0;
  // Diffusion-current source F sum_k z_k D_k grad c_k.
  ADRealVectorValue diffusion_current(0.0, 0.0, 0.0);
  for (unsigned int k = 0; k < _n_ions; ++k)
  {
    const Real z = _valences[k];
    const Real d = _diffusivities[k];
    const ADReal c = (*_c[k])[_qp];
    kappa += f2_rt * z * z * d * (c > 0.0 ? c : ADReal(0.0));
    diffusion_current += (Corrosion::faraday * z * d) * (*_grad_c[k])[_qp];
  }

  // Residual of div(i) = 0: (grad test, kappa grad phi + F sum_k z_k D_k grad c_k).
  return kappa * _grad_u[_qp] + diffusion_current;
}
