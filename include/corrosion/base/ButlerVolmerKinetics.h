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
#include "CorrosionChemistry.h"

#include <cmath>

namespace Corrosion
{

/**
 * Shared Butler-Volmer electrode kinetics, used identically by the two-block interface kernel and
 * the single-domain boundary condition so that both topologies solve the same physics.
 *
 * The exponent is clipped to +/-60 (matching the reference model.py::exp_clip) to keep the current
 * finite and Newton stable when the overpotential is large.
 */

/// exp(arg) with the exponent limited to [-60, 60]; the AD derivative is preserved inside the range.
template <typename T>
T
clampedExp(const T & arg)
{
  using std::exp;
  const Real raw = MetaPhysicL::raw_value(arg);
  if (raw > 60.0)
    return exp(T(60.0));
  if (raw < -60.0)
    return exp(T(-60.0));
  return exp(arg);
}

/**
 * Butler-Volmer current density for the reaction M(s) <-> M^{z+} + z e-, in concentration-explicit
 * form.
 *
 *   eta    = (phi_metal - phi_salt) - E0                              (activation overpotential)
 *   i_BV   = i0 [ exp(alpha_a z F eta / R T)
 *                 - (c_ion / c_ref) exp(-alpha_c z F eta / R T) ]
 *
 * The anodic (dissolution) branch carries the solid metal activity, taken constant (= 1), so the
 * corrosion rate is set by the overpotential, not by how much metal has already dissolved. The
 * cathodic (plating) branch carries the dissolved-ion activity c_ion / c_ref, so plating accelerates
 * as the salt loads and the reaction reaches equilibrium (i_BV = 0) when
 * c_ion = c_ref exp[(alpha_a + alpha_c) z F eta / R T]. A positive (anodic) current dissolves metal
 * into the salt; a negative (cathodic) current plates it out, so one expression covers both
 * corrosion and deposition. This is the form used for metal electrodes, and it keeps the dissolution
 * rate bounded for fresh (nearly clean) salt.
 *
 * @param phi_metal Electric potential of the metal phase [V].
 * @param phi_salt  Electric potential of the salt phase [V].
 * @param c_ion     Salt-side concentration of M^{z+} [mol/m^3].
 * @param E0        Standard electrode potential [V].
 * @param valence   Charge number z of the cation.
 * @param i0        Exchange current density [A/m^2].
 * @param alpha_a   Anodic charge-transfer coefficient.
 * @param alpha_c   Cathodic charge-transfer coefficient.
 * @param temperature Absolute temperature [K].
 * @param c_ref     Reference concentration at which i0 is defined [mol/m^3].
 * @return Butler-Volmer current density [A/m^2].
 */
template <typename T>
T
bvCurrent(const T & phi_metal,
          const T & phi_salt,
          const T & c_ion,
          Real E0,
          Real valence,
          Real i0,
          Real alpha_a,
          Real alpha_c,
          Real temperature,
          Real c_ref)
{
  const Real f = faraday / (R_gas * temperature);
  const T eta = (phi_metal - phi_salt) - E0;
  // Dissolved-ion activity for the cathodic (plating) branch, floored at zero so a Newton overshoot
  // to a slightly negative concentration cannot create a spurious deposition current.
  const T activity = ((c_ion > T(0.0)) ? c_ion : T(0.0)) / c_ref;
  const T anodic = clampedExp(static_cast<T>(alpha_a * valence * f * eta));
  const T cathodic = activity * clampedExp(static_cast<T>(-alpha_c * valence * f * eta));
  return i0 * (anodic - cathodic);
}

} // namespace Corrosion
