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

#include "CorrosionChemistry.h"

#include <algorithm>
#include <cmath>

namespace Corrosion
{

Real
expClip(Real x, Real lo, Real hi)
{
  return std::exp(std::min(std::max(x, lo), hi));
}

Real
safeLog(Real x, Real floor)
{
  return std::log(std::max(x, floor));
}

Real
mgCm2ToUm(Real value_mg_cm2, Real density_g_cm3)
{
  // 1 mg/cm^2 over a layer of density rho [g/cm^3] is a thickness of (value / (rho * 0.1)) um.
  return value_mg_cm2 / (density_g_cm3 * 0.1);
}

Real
umToMgCm2(Real value_um, Real density_g_cm3)
{
  return value_um * density_g_cm3 * 0.1;
}

Real
corrosionCurrentToUmY(Real current_a_cm2, Real valence, Real molar_mass_g_mol, Real density_g_cm3)
{
  const Real cm_per_s = current_a_cm2 * molar_mass_g_mol / (valence * faraday * density_g_cm3);
  return cm_per_s * 1.0e4 * seconds_per_year;
}

Real
umYToCorrosionCurrent(Real rate_um_y, Real valence, Real molar_mass_g_mol, Real density_g_cm3)
{
  const Real cm_per_s = rate_um_y / (1.0e4 * seconds_per_year);
  return cm_per_s * valence * faraday * density_g_cm3 / molar_mass_g_mol;
}

} // namespace Corrosion
