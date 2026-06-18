//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"
#include "MoltenSaltRadiolysisData.h"
#include "DataFileUtils.h"

#include <cmath>

// These tests mirror the analytic targets in
// MoltenSaltRadiolysis/msr_radiolysis/data/reference_validation.yaml and the rate-law definitions
// in core.py, verifying that the bundled JSON chemistry database reproduces the validated model.

namespace
{
// Load the bundled database (resolved via the registered application data path).
MSR::MoltenSaltRadiolysisDatabase
loadDatabase()
{
  return MSR::MoltenSaltRadiolysisDatabase(Moose::DataFileUtils::getPath("msr_database.json").path);
}

const MSR::ReactionData &
findReaction(const std::vector<MSR::ReactionData> & network, const std::string & name)
{
  for (const auto & rxn : network)
    if (rxn.name == name)
      return rxn;
  throw std::runtime_error("reaction not found: " + name);
}
}

TEST(MoltenSaltRadiolysisData, gValueConversion)
{
  // reference_validation.yaml: g_value_conversion
  const Real edot = 1.0e6;
  const Real g = 0.005;
  const Real expected = 0.005 * 1.0e6 / (100.0 * 1.602176634e-19) / 6.02214076e23;
  EXPECT_NEAR(MSR::gToSource(g, edot), expected, std::abs(expected) * 1e-12);

  // A zero G value yields no source.
  EXPECT_DOUBLE_EQ(MSR::gToSource(0.0, edot), 0.0);
}

TEST(MoltenSaltRadiolysisData, arrheniusKRefAtReferenceTemperature)
{
  // For a reaction parametrized with (k_ref, T_ref, Ea), k(T_ref) == k_ref since exp(0) = 1.
  const auto db = loadDatabase();
  const auto chloride = db.coreReactions("chloride");
  const auto & dispro = findReaction(chloride, "2 Cl2m_rad -> Cl3_ion + Cl_ion");
  EXPECT_NEAR(dispro.rate(673.15), 2.2e6, 2.2e6 * 1e-12);
}

TEST(MoltenSaltRadiolysisData, arrheniusKRefTemperatureDependence)
{
  // k(T) = k_ref * exp(Ea/R * (1/T_ref - 1/T)) away from the reference temperature.
  const auto db = loadDatabase();
  const auto chloride = db.coreReactions("chloride");
  const auto & r1 = findReaction(chloride, "Cl_rad + Cl_ion -> Cl2m_rad");
  const Real T = 773.15;
  const Real expected = 1.0e7 * std::exp(2.0e4 / MSR::R_gas * (1.0 / 673.15 - 1.0 / T));
  EXPECT_NEAR(r1.rate(T), expected, expected * 1e-12);
}

TEST(MoltenSaltRadiolysisData, arrheniusPreExponentialForm)
{
  // For a reaction parametrized with (A, Ea), k(T) = A * exp(-Ea/(R T)).
  const auto db = loadDatabase();
  const auto zn = db.metalReactions("Zn", "chloride");
  const auto & capture = findReaction(zn, "e_sol + Zn_II -> Zn_I");
  const Real T = 673.15;
  const Real expected = 2.4e10 * std::exp(-3.56e4 / (MSR::R_gas * T));
  EXPECT_NEAR(capture.rate(T), expected, expected * 1e-12);
  // The capture rate should be on the order of 1e7 m^3/mol/s at 400 C (literature scoping value).
  EXPECT_GT(capture.rate(T), 1.0e7);
  EXPECT_LT(capture.rate(T), 1.0e8);
}

TEST(MoltenSaltRadiolysisData, rateConstantWithoutActivationEnergyIsTemperatureIndependent)
{
  // Cr reactions reported at a single temperature carry only k_ref (no Ea) -> T-independent.
  const auto db = loadDatabase();
  const auto cr = db.metalReactions("Cr", "chloride");
  const auto & r9 = findReaction(cr, "Cl2m_rad + Cr_II -> Cr_III + 2 Cl_ion");
  EXPECT_DOUBLE_EQ(r9.rate(673.15), 7.2e6);
  EXPECT_DOUBLE_EQ(r9.rate(873.15), 7.2e6);
}

TEST(MoltenSaltRadiolysisData, speciesAndNetworkSizes)
{
  const auto db = loadDatabase();
  EXPECT_EQ(db.coreSpecies("chloride").size(), 6u);
  EXPECT_EQ(db.coreSpecies("fluoride").size(), 4u);
  EXPECT_EQ(db.coreReactions("chloride").size(), 5u);
  EXPECT_EQ(db.coreReactions("fluoride").size(), 2u);

  EXPECT_EQ(db.metalSpecies("Zn").size(), 2u);
  EXPECT_EQ(db.metalSpecies("Cr").size(), 3u);
  EXPECT_EQ(db.metalSpecies("U").size(), 2u);

  EXPECT_EQ(db.metalReactions("Zn", "chloride").size(), 2u);
  EXPECT_EQ(db.metalReactions("Cr", "chloride").size(), 5u);
  EXPECT_EQ(db.metalReactions("U", "fluoride").size(), 2u);

  // Case-insensitive metal lookup (MultiMooseEnum upper-cases its values).
  EXPECT_EQ(db.metalSpecies("ZN").size(), 2u);
}

TEST(MoltenSaltRadiolysisData, defaultGValuesAndHenry)
{
  const auto db = loadDatabase();
  const auto g_chloride = db.defaultGValues("gamma", "chloride");
  EXPECT_DOUBLE_EQ(g_chloride.at("e_sol"), 0.30);
  EXPECT_DOUBLE_EQ(g_chloride.at("Cl_rad"), 0.30);

  const auto g_fluoride = db.defaultGValues("gamma", "fluoride");
  EXPECT_DOUBLE_EQ(g_fluoride.at("F2_diss"), 0.005);

  EXPECT_DOUBLE_EQ(db.henryCoefficient("Cl2"), 1.0e-5);
  EXPECT_DOUBLE_EQ(db.henryCoefficient("F2"), 1.0e-6);
  EXPECT_EQ(db.gasLiquidSpecies("Cl2"), "Cl2_diss");
  EXPECT_EQ(db.gasPhaseSpecies("Cl2"), "Cl2_gas");
}
