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

#include "ButlerVolmerKinetics.h"

#include <cmath>
#include <sstream>
#include <string>

namespace Corrosion
{
namespace ActionKinetics
{

/**
 * Butler-Volmer bracket used to seed an exchange current at c_ion = c_ref.
 *
 * Calling the production bvCurrent implementation with unit exchange current and unit activity
 * keeps the Action seed exactly aligned with the generated boundary objects, including the
 * metal-minus-salt potential convention, E0 subtraction, and exponent clipping.
 */
inline Real
butlerVolmerSeedBracket(Real metal_minus_salt_potential,
                        Real E0,
                        Real valence,
                        Real alpha_a,
                        Real alpha_c,
                        Real temperature)
{
  return bvCurrent(metal_minus_salt_potential,
                   0.0,
                   1.0,
                   E0,
                   valence,
                   1.0,
                   alpha_a,
                   alpha_c,
                   temperature,
                   1.0);
}

/**
 * Absolute metal potential corresponding to a reference metal-minus-salt difference.
 *
 * When the salt potential is solved, its constant mode is pinned to an absolute value. The default
 * metal potential must carry the same offset so their difference remains the seed overpotential.
 */
inline Real
defaultMetalPotential(Real metal_minus_salt_potential,
                      bool solve_salt_potential,
                      Real pinned_salt_potential)
{
  return metal_minus_salt_potential + (solve_salt_potential ? pinned_salt_potential : 0.0);
}

enum class ReferencePotentialStatus
{
  Valid,
  NonfiniteOverpotential,
  NonfiniteSaltPin,
  NonfiniteDefaultMetalPotential
};

/**
 * Classify reference-potential inputs before an Action creates numeric functors.
 *
 * A finite pinned salt potential is required only when the salt potential is solved. The final
 * addition is checked separately because two finite inputs can overflow.
 */
inline ReferencePotentialStatus
referencePotentialStatus(Real metal_minus_salt_potential,
                         bool solve_salt_potential,
                         Real pinned_salt_potential)
{
  if (!std::isfinite(metal_minus_salt_potential))
    return ReferencePotentialStatus::NonfiniteOverpotential;
  if (solve_salt_potential && !std::isfinite(pinned_salt_potential))
    return ReferencePotentialStatus::NonfiniteSaltPin;
  if (solve_salt_potential &&
      !std::isfinite(
          defaultMetalPotential(metal_minus_salt_potential, true, pinned_salt_potential)))
    return ReferencePotentialStatus::NonfiniteDefaultMetalPotential;
  return ReferencePotentialStatus::Valid;
}

/// Format a Real without the six-digit precision loss of std::to_string when it names a functor.
inline std::string
realFunctorName(Real value)
{
  std::ostringstream oss;
  oss.precision(17);
  oss << value;
  return oss.str();
}

} // namespace ActionKinetics
} // namespace Corrosion
