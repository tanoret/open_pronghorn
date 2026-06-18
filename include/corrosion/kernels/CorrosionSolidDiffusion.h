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

#pragma once

#include "ADKernelGrad.h"

/**
 * Solid-state diffusion of a metal species in the alloy, e.g. chromium depletion feeding the
 * interface dissolution. The residual is the standard (grad test, D_s grad c_s); the solid-state
 * diffusivity is supplied as a constant (the CorrosionPlatingAction sets it from the Arrhenius
 * correlation D_s(T) at the operating temperature).
 */
class CorrosionSolidDiffusion : public ADKernelGrad
{
public:
  static InputParameters validParams();

  CorrosionSolidDiffusion(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  /// Solid-state diffusivity [m^2/s].
  const Real _diffusivity;
};
