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

#include "CorrosionNernstPlanckFlux.h"
#include "CorrosionChemistry.h"

registerMooseObject("OpenPronghornApp", CorrosionNernstPlanckFlux);

InputParameters
CorrosionNernstPlanckFlux::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addClassDescription("Nernst-Planck diffusion and electromigration flux of a dissolved "
                             "cation in the molten salt.");

  params.addRequiredParam<Real>("diffusivity", "Molecular diffusivity of the cation [m^2/s].");
  params.addRequiredParam<Real>("valence", "Charge number z of the cation.");
  params.addRequiredParam<MooseFunctorName>("temperature", "Temperature functor [K].");
  params.addCoupledVar(
      "potential",
      "The electric potential [V] driving migration. When omitted only diffusion is added.");
  return params;
}

CorrosionNernstPlanckFlux::CorrosionNernstPlanckFlux(const InputParameters & parameters)
  : ADKernelGrad(parameters),
    _diffusivity(getParam<Real>("diffusivity")),
    _valence(getParam<Real>("valence")),
    _temperature(getFunctor<Real>("temperature")),
    _has_potential(isCoupled("potential")),
    _grad_potential(_has_potential ? adCoupledGradient("potential") : _ad_grad_zero)
{
}

ADRealVectorValue
CorrosionNernstPlanckFlux::precomputeQpResidual()
{
  const Moose::ElemQpArg qp_arg = {_current_elem, _qp, _qrule, _q_point[_qp]};
  const Real temperature = _temperature(qp_arg, Moose::currentState());

  // Diffusion flux contribution (grad-test, D grad c).
  ADRealVectorValue flux = _diffusivity * _grad_u[_qp];

  // Migration flux contribution (grad-test, (z F / R T) D c grad phi).
  if (_has_potential)
  {
    const Real z_f_rt = _valence * Corrosion::faraday / (Corrosion::R_gas * temperature);
    flux += (z_f_rt * _diffusivity * _u[_qp]) * _grad_potential[_qp];
  }

  return flux;
}
