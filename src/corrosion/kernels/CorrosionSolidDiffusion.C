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

#include "CorrosionSolidDiffusion.h"

registerMooseObject("OpenPronghornApp", CorrosionSolidDiffusion);

InputParameters
CorrosionSolidDiffusion::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addClassDescription("Solid-state diffusion of a metal species in the alloy.");
  params.addRequiredParam<Real>("diffusivity", "Solid-state diffusivity [m^2/s].");
  return params;
}

CorrosionSolidDiffusion::CorrosionSolidDiffusion(const InputParameters & parameters)
  : ADKernelGrad(parameters), _diffusivity(getParam<Real>("diffusivity"))
{
}

ADRealVectorValue
CorrosionSolidDiffusion::precomputeQpResidual()
{
  return _diffusivity * _grad_u[_qp];
}
