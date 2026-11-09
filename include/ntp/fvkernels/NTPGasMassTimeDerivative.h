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

#include "FVTimeKernel.h"

/**
 * Mass storage of a gas in a porous medium, written on the pressure equation:
 *
 *   d(eps rho)/dt = eps (drho/dp dp/dt + drho/dT dT/dt),  rho = rho(p, T) from the fluid
 * properties.
 *
 * Used with the Darcy-Forchheimer mass equation of the particle-bed element. In steady
 * calculations the term only damps the pseudo-transient iterations, but it does so with the
 * physical compressibility, so pressure and mass flow relax together.
 */
class NTPGasMassTimeDerivative : public FVTimeKernel
{
public:
  static InputParameters validParams();

  NTPGasMassTimeDerivative(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  const SinglePhaseFluidProperties & _fp;
  const Moose::Functor<ADReal> & _T;
  const Moose::Functor<ADReal> & _eps;
};
