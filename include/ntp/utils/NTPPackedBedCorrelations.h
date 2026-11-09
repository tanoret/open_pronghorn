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

#include "MooseTypes.h"

#include <cmath>

/**
 * Packed-bed correlations of the particle-bed fuel element, shared by the Newton (Darcy-
 * Forchheimer) and the segregated (SIMPLE/PIMPLE) closure materials. Templated on the number type
 * so that both Real and ADReal evaluations use the same formulas.
 *
 * Pressure drop, superficial velocity u, porosity eps, particle diameter d:
 *   -grad p = (A + B |u|) u
 *   KTA 3102.3:  A = 160 mu (1-eps)^2 / (eps^3 d^2),  B = 3 rho (1-eps) / (eps^3 d) Re_m^-0.1
 *   Ergun:       A = 150 mu (1-eps)^2 / (eps^3 d^2),  B = 1.75 rho (1-eps) / (eps^3 d)
 *   Re_m = rho |u| d / (mu (1 - eps))
 * (KTA: Psi = 320/Re_m + 6/Re_m^0.1 multiplying (1-eps)/(eps^3 d) rho |u| u / 2.)
 *
 * Particle-to-gas heat transfer, Re = rho |u| d / mu, Pr = cp mu / k:
 *   Achenbach (1995):  Nu = [(1.18 Re^0.58)^4 + (0.23 (Re/(1-eps))^0.75)^4]^(1/4)
 *   KTA 3102.2:        Nu = 1.27 Pr^(1/3) Re^0.36 / eps^1.18 + 0.033 Pr^(1/2) Re^0.86 / eps^1.07
 */
namespace NTPPackedBed
{
enum class PressureDrop
{
  KTA,
  ERGUN
};
enum class HeatTransfer
{
  ACHENBACH,
  KTA
};

/// Viscous coefficient A of -grad p = (A + B|u|) u [Pa s/m^2]
template <typename T>
T
viscousCoefficient(PressureDrop model, const T & mu, Real eps, Real d)
{
  return (model == PressureDrop::KTA ? 160.0 : 150.0) * mu * (1 - eps) * (1 - eps) /
         (eps * eps * eps * d * d);
}

/// Inertial coefficient B of -grad p = (A + B|u|) u [kg/m^4], at superficial speed 'speed'
template <typename T>
T
inertialCoefficient(PressureDrop model,
                    const T & rho,
                    const T & mu,
                    const T & speed,
                    Real eps,
                    Real d,
                    Real Re_floor)
{
  using std::max;
  using std::pow;
  const Real geom = (1 - eps) / (eps * eps * eps * d);
  if (model == PressureDrop::ERGUN)
    return 1.75 * rho * geom;
  const T Re_m = max(T(rho * speed * d / (mu * (1 - eps))), T(Re_floor));
  return 3.0 * rho * geom * pow(Re_m, -0.1);
}

/// Particle Nusselt number at superficial speed 'speed'
template <typename T>
T
nusselt(HeatTransfer model,
        const T & rho,
        const T & mu,
        const T & k,
        const T & cp,
        const T & speed,
        Real eps,
        Real d,
        Real Re_floor)
{
  using std::max;
  using std::pow;
  const T Re = max(T(rho * speed * d / mu), T(Re_floor));
  if (model == HeatTransfer::ACHENBACH)
    return pow(pow(1.18 * pow(Re, 0.58), 4) + pow(0.23 * pow(Re / (1 - eps), 0.75), 4), 0.25);
  const T Pr = cp * mu / k;
  return 1.27 * pow(Pr, 1.0 / 3.0) * pow(Re, 0.36) / std::pow(eps, 1.18) +
         0.033 * pow(Pr, 0.5) * pow(Re, 0.86) / std::pow(eps, 1.07);
}

/// Specific particle surface area a_s = 6 (1 - eps) / d [1/m]
inline Real
specificSurface(Real eps, Real d)
{
  return 6 * (1 - eps) / d;
}
}
