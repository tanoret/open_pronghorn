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

/**
 * Mobility of a free channel (the plenum or the chamber of a particle-bed element) for a
 * Darcy-type description of the flow, u = -M grad p, with y the channel axis (RZ, or 3-D with x
 * and z transverse).
 *
 * Along the axis the friction of a duct of hydraulic diameter D_h is
 *   -dp/dz = (32 mu / D_h^2 + f rho |u| / (2 D_h)) u,
 * the laminar Hagen-Poiseuille term plus the Darcy-Weisbach term with friction factor f, so
 *   M_z = 2 / (a + sqrt(a^2 + 4 b |dp/dz|)),  a = 32 mu / D_h^2,  b = f rho / (2 D_h).
 * Inertia (acceleration) is neglected, as in a reduced-order channel model. Across the channel
 * the mobility is 'radial_mobility_ratio' times M_z: a free channel offers almost no resistance
 * to the radial redistribution that feeds or drains the porous layers.
 */
class NTPDuctFlowFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  NTPDuctFlowFunctorMaterial(const InputParameters & parameters);

protected:
  const Moose::Functor<ADReal> & _rho;
  const Moose::Functor<ADReal> & _mu;
  const Moose::Functor<ADReal> & _p;
  const Moose::Functor<ADReal> & _Dh;
  const Real _f;
  const Real _radial_ratio;
};
