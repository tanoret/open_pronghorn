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

#include "FunctorMaterial.h"

class SinglePhaseFluidProperties;

/**
 * Gas density, viscosity, isobaric specific heat and thermal conductivity as functors of a pressure
 * and a temperature functor, each under a chosen name, and nothing else. Two instances can
 * therefore coexist on the same blocks, e.g. the flow properties at a lagged state and the energy
 * properties at the current gas temperature of a segregated solve. Templated on AD:
 * NTPGasPropertiesFunctorMaterial gives Real functors (linear finite volume),
 * ADNTPGasPropertiesFunctorMaterial AD ones.
 */
template <bool is_ad>
class NTPGasPropertiesFunctorMaterialTempl : public FunctorMaterial
{
public:
  static InputParameters validParams();

  NTPGasPropertiesFunctorMaterialTempl(const InputParameters & parameters);

protected:
  const SinglePhaseFluidProperties & _fp;
  const Moose::Functor<GenericReal<is_ad>> & _p;
  const Moose::Functor<GenericReal<is_ad>> & _T;
};

typedef NTPGasPropertiesFunctorMaterialTempl<false> NTPGasPropertiesFunctorMaterial;
typedef NTPGasPropertiesFunctorMaterialTempl<true> ADNTPGasPropertiesFunctorMaterial;
