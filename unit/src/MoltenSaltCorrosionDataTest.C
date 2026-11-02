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

#include "gtest/gtest.h"
#include "MoltenSaltCorrosionData.h"
#include "CorrosionChemistry.h"
#include "DataFileUtils.h"

#include <cmath>

// These tests verify that the bundled corrosion database reproduces the engineering tables and the
// calibrated parameter set of the validated reference model (msr_corrosion_bv/chemistry.py and
// msr_corrosion_plating_model/results/parameters.json), and that the Faradaic conversions are exact
// inverses.

namespace
{
Corrosion::MoltenSaltCorrosionDatabase
loadDatabase()
{
  return Corrosion::MoltenSaltCorrosionDatabase(
      Moose::DataFileUtils::getPath("corrosion_database.json").path);
}
}

TEST(MoltenSaltCorrosionData, elementValencesAndMolarMasses)
{
  const auto db = loadDatabase();

  // Valences match chemistry.py VALENCE.
  EXPECT_DOUBLE_EQ(db.element("Cr").valence, 2.0);
  EXPECT_DOUBLE_EQ(db.element("Fe").valence, 2.0);
  EXPECT_DOUBLE_EQ(db.element("Ni").valence, 2.0);
  EXPECT_DOUBLE_EQ(db.element("Mo").valence, 3.0);
  EXPECT_DOUBLE_EQ(db.element("Nb").valence, 4.0);
  EXPECT_DOUBLE_EQ(db.element("Tc").valence, 4.0);
  EXPECT_DOUBLE_EQ(db.element("Ru").valence, 3.0);
  EXPECT_DOUBLE_EQ(db.element("Ag").valence, 1.0);
  EXPECT_DOUBLE_EQ(db.element("Sb").valence, 3.0);
  EXPECT_DOUBLE_EQ(db.element("Te").valence, 2.0);

  // Molar masses match chemistry.py MOLAR_MASS_G_MOL.
  EXPECT_DOUBLE_EQ(db.element("Cr").molar_mass_g_mol, 51.9961);
  EXPECT_DOUBLE_EQ(db.element("Ni").molar_mass_g_mol, 58.6934);
  EXPECT_DOUBLE_EQ(db.element("Ag").molar_mass_g_mol, 107.8682);

  // Case-insensitive lookup (MultiMooseEnum upper-cases its values).
  EXPECT_DOUBLE_EQ(db.element("cr").valence, 2.0);

  // An unknown element falls back to the generic-metal entry (z = 2, M = 55).
  EXPECT_FALSE(db.hasElement("Xx"));
  EXPECT_DOUBLE_EQ(db.element("Xx").valence, 2.0);
  EXPECT_DOUBLE_EQ(db.element("Xx").molar_mass_g_mol, 55.0);
}

TEST(MoltenSaltCorrosionData, densityAndChromiumFractionTables)
{
  const auto db = loadDatabase();

  // Densities and chromium weight fractions match chemistry.py.
  EXPECT_DOUBLE_EQ(db.density("hastelloy_n"), 8.89);
  EXPECT_DOUBLE_EQ(db.density("stainless_304"), 8.00);
  EXPECT_DOUBLE_EQ(db.density("graphite"), 1.80);
  EXPECT_DOUBLE_EQ(db.crWeightFraction("hastelloy_n"), 0.07);
  EXPECT_DOUBLE_EQ(db.crWeightFraction("stainless_304"), 0.19);
  EXPECT_DOUBLE_EQ(db.crWeightFraction("in625"), 0.215);

  // Unknown classes fall back to generic_metal (rho = 8.30, Cr fraction = 0.12).
  EXPECT_DOUBLE_EQ(db.density("unobtanium"), 8.30);
  EXPECT_DOUBLE_EQ(db.crWeightFraction("unobtanium"), 0.12);
}

TEST(MoltenSaltCorrosionData, solidDiffusivityProperties)
{
  const auto db = loadDatabase();

  // Hastelloy N has an explicit chromium diffusivity correlation from DeVan.
  EXPECT_TRUE(db.hasSolidDiffusivity("hastelloy_n", "Cr"));

  // Material and element lookup are case-insensitive.
  EXPECT_TRUE(db.hasSolidDiffusivity("HASTELLOY_N", "cr"));

  const auto props = db.solidDiffusivity("hastelloy_n", "Cr");

  EXPECT_DOUBLE_EQ(props.D0_cm2_s, 1.382126667841e-7);
  EXPECT_DOUBLE_EQ(props.Q_kJ_mol, 120.620464930);
  EXPECT_DOUBLE_EQ(props.measurement_temperature_min_K, 964.15);
  EXPECT_DOUBLE_EQ(props.measurement_temperature_max_K, 1143.15);
  EXPECT_EQ(props.source, "DeVan 1960, Table XVI, Loop 1248");
  EXPECT_EQ(props.status, "direct_tracer_measurement");

  // SUS316 has a direct chromium volume/lattice diffusivity correlation from Mizouchi et al.
  EXPECT_TRUE(db.hasSolidDiffusivity("stainless_316", "Cr"));

  const auto stainless_316 = db.solidDiffusivity("stainless_316", "Cr");
  EXPECT_DOUBLE_EQ(stainless_316.D0_cm2_s, 1.13e-3);
  EXPECT_DOUBLE_EQ(stainless_316.Q_kJ_mol, 234.0);
  EXPECT_DOUBLE_EQ(stainless_316.measurement_temperature_min_K, 888.0);
  EXPECT_DOUBLE_EQ(stainless_316.measurement_temperature_max_K, 1173.0);
  EXPECT_EQ(stainless_316.source,
            "Mizouchi et al. 2004, Sec. 3.1 and Fig. 2, 51Cr volume/lattice diffusion in SUS316");
  EXPECT_EQ(stainless_316.status, "direct_tracer_measurement");

  // Related stainless grades use the explicitly configured SUS316 engineering fallback.
  EXPECT_TRUE(db.hasSolidDiffusivity("stainless_316h", "Cr"));
  EXPECT_TRUE(db.hasSolidDiffusivity("stainless_316l", "Cr"));
  EXPECT_TRUE(db.hasSolidDiffusivity("stainless_304", "Cr"));
  EXPECT_TRUE(db.hasSolidDiffusivity("stainless_304l", "Cr"));

  const auto stainless_316h = db.solidDiffusivity("stainless_316h", "Cr");
  EXPECT_DOUBLE_EQ(stainless_316h.D0_cm2_s, stainless_316.D0_cm2_s);
  EXPECT_DOUBLE_EQ(stainless_316h.Q_kJ_mol, stainless_316.Q_kJ_mol);
  EXPECT_EQ(stainless_316h.source, stainless_316.source);
  EXPECT_EQ(stainless_316h.status, "direct_tracer_measurement");

  // Fallback lookup remains case-insensitive.
  EXPECT_TRUE(db.hasSolidDiffusivity("STAINLESS_316H", "cr"));

  // An unrelated material still has no alloy-specific or configured fallback entry.
  EXPECT_FALSE(db.hasSolidDiffusivity("unobtanium", "Cr"));

// Generic metal retains the old correlation as an explicit legacy fallback.
  EXPECT_TRUE(db.hasSolidDiffusivity("generic_metal", "Cr"));

  const auto legacy = db.solidDiffusivity("generic_metal", "Cr");
  EXPECT_DOUBLE_EQ(legacy.D0_cm2_s, 1.5195890862355468);
  EXPECT_DOUBLE_EQ(legacy.Q_kJ_mol, 240.0);
  EXPECT_EQ(legacy.source, "Legacy 61-parameter corrosion calibration");
  EXPECT_EQ(legacy.status, "legacy_fallback");
}

TEST(MoltenSaltCorrosionData, calibratedParameterLookups)
{
  const auto db = loadDatabase();

  // A spot-check of the fitted parameters from results/parameters.json and validation corrections.
  EXPECT_DOUBLE_EQ(db.parameter("log_rate0_um_y"), 0.26537174220781806);
  EXPECT_DOUBLE_EQ(db.parameter("redox_oxidizing_fef2"), 1.5026093263788216);
  EXPECT_DOUBLE_EQ(db.parameter("mat_stainless_304"), 1.5121011358409628);
  EXPECT_EQ(db.calibratedParameters().size(), 60u);
}

TEST(MoltenSaltCorrosionData, faradaicConversionsRoundTrip)
{
  using namespace Corrosion;
  // Chromium in Hastelloy N: z = 2, M = 51.9961 g/mol, rho = 8.89 g/cm^3.
  const Real z = 2.0, M = 51.9961, rho = 8.89;

  // The current<->rate conversions are exact inverses.
  const Real rate = 25.0; // um/y
  const Real current = umYToCorrosionCurrent(rate, z, M, rho);
  EXPECT_NEAR(corrosionCurrentToUmY(current, z, M, rho), rate, rate * 1.0e-12);

  // Closed-form check against chemistry.py: cm_per_s = rate / (1e4 * SEC_PER_YEAR).
  const Real cm_per_s = rate / (1.0e4 * seconds_per_year);
  const Real expected_current = cm_per_s * z * faraday * rho / M;
  EXPECT_NEAR(current, expected_current, std::abs(expected_current) * 1.0e-12);

  // Mass<->thickness conversion: 1 um of Hastelloy N is rho * 0.1 mg/cm^2.
  EXPECT_NEAR(umToMgCm2(1.0, rho), rho * 0.1, 1.0e-12);
  EXPECT_NEAR(mgCm2ToUm(umToMgCm2(3.0, rho), rho), 3.0, 1.0e-12);
}
