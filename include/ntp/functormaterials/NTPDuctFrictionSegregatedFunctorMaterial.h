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

#include "FunctorMaterial.h"

/**
 * Wall friction of a free channel (the plenum or the chamber of a particle-bed element) for the
 * segregated (SIMPLE / PIMPLE) solvers, as the Darcy coefficient of LinearFVMomentumFriction
 * (force per unit volume mu D u).
 *
 * The channels are resolved by a few cells across their gap, too few to resolve the turbulent
 * wall layer, so the wall shear is modelled as the friction of a duct of hydraulic diameter D_h,
 * the same law as NTPDuctFlowFunctorMaterial of the Newton path:
 *
 *   -dp/dz|_friction = (32 mu / D_h^2 + f rho |u_z| / (2 D_h)) u_z,
 *
 * so D_z = (32 mu / D_h^2 + f rho |u_z| / (2 D_h)) / mu, and D_r = 0: the friction acts on the flow
 * along the channel axis (y) only, at the axial speed |u_z|. The radial velocity of a channel cell
 * is the cross flow that feeds or drains the porous layer next to it, not duct flow; near the
 * closed end of the plenum it exceeds the axial velocity. |u_z| is lagged by one iteration, from
 * the cell velocity, or is the optional 'speed' functor: an under-relaxed speed, updated once per
 * outer iteration, keeps the Picard iteration on the turbulent term (whose coefficient grows with
 * the solution) from oscillating.
 */
class NTPDuctFrictionSegregatedFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  NTPDuctFrictionSegregatedFunctorMaterial(const InputParameters & parameters);

protected:
  const Moose::Functor<Real> & _rho;
  const Moose::Functor<Real> & _mu;
  const Moose::Functor<Real> & _v;
  const Moose::Functor<Real> * const _speed;
  const Moose::Functor<Real> & _Dh;
  const Real _f;
  const bool _laminar;
};
