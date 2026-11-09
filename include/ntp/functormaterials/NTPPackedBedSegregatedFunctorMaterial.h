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
 * Closures of one packed-bed region of a particle-bed fuel element for the segregated
 * (SIMPLE / PIMPLE, linear finite volume) solvers, in which the superficial velocity u is solved.
 * All outputs are Real functors, the type the linear finite-volume objects read.
 *
 *  - 'drag_name' (vector): Darcy coefficient D of LinearFVMomentumFriction, whose force per unit
 *    volume is mu D u. With the packed-bed law -grad p = (A + B |u|) u (NTPPackedBedCorrelations),
 *    D = chi (A + B |u|) / mu, isotropic, chi being 'pressure_drop_multiplier' times the functor
 *    'orifice_multiplier'. B depends on |u|, which the segregated iterations lag by one iteration.
 *  - 'alpha_name': interphase coefficient h a_s [W/(m^3 K)], Nu from Achenbach or KTA.
 *  - 'kappa_name' (vector): particle conductivity [(1 - eps) k_s + eps k_f] (chi_r, chi_z, chi_r),
 *    y being the element axis (in 3-D, x and z are both transverse).
 *
 * The superficial speed |u| is the cell velocity (u, v only), or |G| / rho with the optional
 * mass-flux inputs (NTPRhieChowMassFluxAux, three components in 3-D), which are exact in porous
 * layers one cell thick and are needed in 3-D.
 *
 * The correlations are the ones of NTPPackedBedFunctorMaterial (the Newton Darcy-Forchheimer
 * path), so the two solution methods solve the same closure equations.
 */
class NTPPackedBedSegregatedFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  NTPPackedBedSegregatedFunctorMaterial(const InputParameters & parameters);

protected:
  /// Superficial speed from the velocity components
  template <typename Space, typename Time>
  Real speed(const Space & r, const Time & t) const;

  const Real _eps;
  const Real _d;
  const MooseEnum _pressure_drop_model;
  const MooseEnum _heat_transfer_model;
  const Real _dp_mult;
  const Real _htc_mult;
  const Real _area_mult;
  const Real _k_solid;
  const Real _k_axial_mult;
  const Real _k_radial_mult;
  const Real _Re_floor;

  const Moose::Functor<Real> & _rho;
  const Moose::Functor<Real> & _mu;
  const Moose::Functor<Real> & _k_fluid;
  const Moose::Functor<Real> & _cp;
  const Moose::Functor<Real> & _u;
  const Moose::Functor<Real> & _v;
  const Moose::Functor<Real> & _orifice;
  /// Optional superficial mass flux components setting the speed
  const Moose::Functor<Real> * const _Gx;
  const Moose::Functor<Real> * const _Gy;
  const Moose::Functor<Real> * const _Gz;
};
