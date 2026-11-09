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

#include "NTPDarcyFluxInterface.h"

/**
 * Gas mass conservation div(m) = 0 on the pressure variable, with the face mass flux m of
 * NTPDarcyFluxInterface (Darcy-Forchheimer mobility, optional momentum flux in free channels).
 * Without a momentum flux it is identical to FVAnisotropicDiffusion with coefficient rho M.
 */
class NTPDarcyMassFlux : public NTPDarcyFluxInterface
{
public:
  static InputParameters validParams();

  NTPDarcyMassFlux(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;
};
