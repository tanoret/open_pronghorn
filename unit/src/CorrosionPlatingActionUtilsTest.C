//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "CorrosionPlatingActionUtils.h"

#include "gtest/gtest.h"

#include <cmath>
#include <initializer_list>
#include <string>

using namespace Corrosion;

TEST(CorrosionPlatingActionUtils, NonzeroE0SeedReproducesProductionButlerVolmerCurrent)
{
  const Real potential_difference = 0.18765432101234567;
  const Real E0 = 0.03125;
  const Real valence = 3.0;
  const Real alpha_a = 0.37;
  const Real alpha_c = 0.63;
  const Real temperature = 948.25;
  const Real c_ref = 42.0;
  const Real reference_current = 3.75;

  const Real effective_eta = potential_difference - E0;
  const Real f = faraday / (R_gas * temperature);
  const Real expected_bracket = expClip(alpha_a * valence * f * effective_eta) -
                                expClip(-alpha_c * valence * f * effective_eta);
  const Real bracket = ActionKinetics::butlerVolmerSeedBracket(
      potential_difference, E0, valence, alpha_a, alpha_c, temperature);
  EXPECT_DOUBLE_EQ(bracket, expected_bracket);

  // Both Actions divide their planned reference current by this shared bracket. At c=c_ref, the
  // production boundary law must therefore recover the planned current, including a nonzero E0.
  // Exercise the same 17-digit numeric-functor serialization used by the linear-FV Action.
  const Real i0 = reference_current / bracket;
  const Real runtime_potential =
      std::stod(ActionKinetics::realFunctorName(potential_difference));
  EXPECT_DOUBLE_EQ(runtime_potential, potential_difference);
  EXPECT_NEAR(bvCurrent(runtime_potential,
                        0.0,
                        c_ref,
                        E0,
                        valence,
                        i0,
                        alpha_a,
                        alpha_c,
                        temperature,
                        c_ref),
              reference_current,
              1.0e-13);

  const Real incorrect_unshifted = expClip(alpha_a * valence * f * potential_difference) -
                                   expClip(-alpha_c * valence * f * potential_difference);
  EXPECT_NE(bracket, incorrect_unshifted);
}

TEST(CorrosionPlatingActionUtils, ZeroE0PreservesLegacySeedBracket)
{
  const Real potential_difference = 0.1;
  const Real valence = 2.0;
  const Real alpha_a = 0.5;
  const Real alpha_c = 0.5;
  const Real temperature = 923.15;
  const Real f = faraday / (R_gas * temperature);
  const Real exponent = alpha_a * valence * f * potential_difference;
  const Real legacy_bracket = std::exp(exponent) - std::exp(-exponent);

  EXPECT_DOUBLE_EQ(ActionKinetics::butlerVolmerSeedBracket(
                       potential_difference, 0.0, valence, alpha_a, alpha_c, temperature),
                   legacy_bracket);
}

TEST(CorrosionPlatingActionUtils, RealFunctorNamesRoundTripWithoutPrecisionLoss)
{
  for (const Real value : {0.12345678901234566, -0.12345678901234566, 1.0000000000000002})
    EXPECT_DOUBLE_EQ(std::stod(ActionKinetics::realFunctorName(value)), value);

  EXPECT_NE(std::stod(std::to_string(0.12345678901234566)), 0.12345678901234566);
}

TEST(CorrosionPlatingActionUtils, DefaultMetalPotentialCarriesSolvedSaltPin)
{
  const Real metal_minus_salt = 0.1;
  const Real salt_pin = 0.0375;

  EXPECT_DOUBLE_EQ(
      ActionKinetics::defaultMetalPotential(metal_minus_salt, true, salt_pin), 0.1375);
  EXPECT_DOUBLE_EQ(
      ActionKinetics::defaultMetalPotential(metal_minus_salt, false, salt_pin), metal_minus_salt);
}
