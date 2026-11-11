//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Copyright 2025, Battelle Energy Alliance, LLC
//* ALL RIGHTS RESERVED
//*

#pragma once

#include "GeneralUserObject.h"
#include "NTPChannelGradientFunctor.h"

/**
 * Registers, for each velocity component, the functor '<component>_channel'
 * (NTPChannelGradientFunctor): the velocity with a cell gradient that ignores the porous
 * neighbours of the free-channel blocks. Closures of the free channels that need the strain rate
 * (k-epsilon production, the viscous stress) take these functors instead of the velocity.
 */
class NTPChannelVelocity : public GeneralUserObject
{
public:
  static InputParameters validParams();

  NTPChannelVelocity(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual void finalize() override {}

protected:
  std::vector<std::unique_ptr<NTPChannelGradientFunctor>> _functors;
};
