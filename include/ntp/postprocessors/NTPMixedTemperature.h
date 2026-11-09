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

#include "GeneralPostprocessor.h"

class SinglePhaseFluidProperties;

/**
 * Mixed-mean (bulk) temperature of a stream from its enthalpy and mass flow rates:
 * T = T(p, h = H / mdot), inverted with the fluid properties. The pressure only matters for a
 * fluid whose enthalpy depends on it.
 */
class NTPMixedTemperature : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  NTPMixedTemperature(const InputParameters & parameters);

  void initialize() override {}
  void execute() override {}
  Real getValue() const override;

protected:
  const SinglePhaseFluidProperties & _fp;
  /// Enthalpy flow rate [W]
  const PostprocessorValue & _H;
  /// Mass flow rate [kg/s]
  const PostprocessorValue & _mdot;
  /// Pressure at which h is inverted [Pa]
  const PostprocessorValue & _p;
};
