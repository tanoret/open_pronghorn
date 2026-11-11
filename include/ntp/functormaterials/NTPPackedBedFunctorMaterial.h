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
 * Closures of one packed-bed region (frit or fuel bed) of a particle-bed fuel element, for a
 * Darcy-Forchheimer description of the flow, u = -M grad p, with u the superficial velocity.
 *
 * Pressure drop. The packed-bed correlations are written as  -grad p = (A + B |u|) u  with
 *   KTA 3102.3:  A = 160 mu (1-eps)^2 / (eps^3 d^2),  B = 3 rho (1-eps) / (eps^3 d) Re_m^-0.1
 *                (that is Psi (1-eps)/(eps^3 d) rho |u| u / 2, Psi = 320/Re_m + 6/Re_m^0.1)
 *   Ergun:       A = 150 mu (1-eps)^2 / (eps^3 d^2),  B = 1.75 rho (1-eps) / (eps^3 d)
 * with Re_m = rho |u| d / (mu (1 - eps)). Both are scaled by 'pressure_drop_multiplier' and by the
 * functor 'orifice_multiplier' (orifice shaping that varies along the element). An orifice solved
 * together with the flow is applied by the flux kernels instead (NTPDarcyFluxInterface). Since
 * |u| = M |grad p|, the mobility solves  B |grad p| M^2 + A M - 1 = 0, so
 *   M = 2 / (A + sqrt(A^2 + 4 B |grad p|)),
 * explicit and smooth. The weak Re_m^-0.1 dependence of the KTA B is resolved by twelve fixed-point
 * iterations started from the Ergun value (converged to round-off; the contraction factor is 0.1).
 *
 * Heat transfer. alpha = h a_s [W/(m^3 K)], a_s = chi_a 6 (1-eps)/d, h = chi_h Nu k_f / d, Nu from
 * Achenbach (1995) or KTA 3102.2 with Re = rho |u| d / mu. The speed |u| should come from the cell
 * mass flux reconstructed from the face fluxes (mass_flux_x/y/z, NTPDarcyMassFluxFunctorMaterial):
 * the cell-pressure-gradient speed is wrong by a large factor in a layer one cell thick.
 *
 * Conduction. kappa_s = [(1 - eps) k_s + eps k_f] (chi_r, chi_z, chi_r), y being the element axis
 * (x the radius in RZ; x and z transverse in 3-D).
 *
 * The forms and multipliers are those of the reduced-order channel model the NTP tutorial
 * reproduces.
 */
class NTPPackedBedFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  NTPPackedBedFunctorMaterial(const InputParameters & parameters);

protected:
  /// Mobility M [m^3 s / kg] and superficial speed |u| from the local pressure gradient
  template <typename Space, typename Time>
  void mobility(const Space & r, const Time & t, ADReal & M, ADReal & speed) const;

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
  /// Particle emissivity (0: no radiative conductivity) and the particle temperature
  const Real _emissivity;
  const Moose::Functor<ADReal> * const _T_solid;
  /// Coefficient of the gas thermal dispersion conductivity (0: none)
  const Real _dispersion;
  /// Smallest Reynolds number used in the correlations
  const Real _Re_floor;

  const Moose::Functor<ADReal> & _rho;
  const Moose::Functor<ADReal> & _mu;
  const Moose::Functor<ADReal> & _k_fluid;
  const Moose::Functor<ADReal> & _cp;
  const Moose::Functor<ADReal> & _p;
  /// Spatially varying drag multiplier (orifice shaping)
  const Moose::Functor<ADReal> & _orifice;
  /// Optional reconstructed cell mass flux for the heat-transfer speed
  const Moose::Functor<ADReal> * const _G_x;
  const Moose::Functor<ADReal> * const _G_y;
  const Moose::Functor<ADReal> * const _G_z;
};
