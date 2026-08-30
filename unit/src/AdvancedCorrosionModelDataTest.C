//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "AdvancedCorrosionModelData.h"
#include "MoltenSaltCorrosionData.h"

#include "DataFileUtils.h"
#include "gtest/gtest.h"

#include <array>
#include <fstream>
#include <stdexcept>
#include <string>

namespace
{
std::string
fixture(const std::string & name)
{
  const std::string source_path(__FILE__);
  const auto separator = source_path.find_last_of("/\\");
  const std::string source_relative =
      separator == std::string::npos
          ? ""
          : source_path.substr(0, separator) + "/../data/advanced/" + name;
  const std::array<std::string, 5> candidates = {source_relative,
                                                 "unit/data/advanced/" + name,
                                                 "../unit/data/advanced/" + name,
                                                 "data/advanced/" + name,
                                                 name};
  for (const auto & candidate : candidates)
  {
    if (!candidate.empty() && std::ifstream(candidate).good())
      return candidate;
  }
  throw std::runtime_error("Unable to locate advanced-corrosion fixture '" + name + "'.");
}
} // namespace

TEST(AdvancedCorrosionModelData, LoadsCorrectedCalibratedArtifactAndProvenance)
{
  const Corrosion::AdvancedCorrosionModelDatabase database(
      Moose::DataFileUtils::getPath("advanced_corrosion_models.json").path);
  EXPECT_EQ(database.schemaVersion(), "1.0");
  EXPECT_EQ(database.calibrationId(), "pr-corrosion-advanced-explicit-context-2026-08-28-v5");
  EXPECT_EQ(database.modelRevision(), "mstdb-nst-v2-fe-identity_dridn-v2-explicit-geometry");
  EXPECT_EQ(database.modelRevision(),
            Corrosion::AdvancedCorrosionModelDatabase::supportedModelRevision());
  EXPECT_EQ(database.expectedMSTDBVersion(), "3.1");
  EXPECT_EQ(database.expectedFluorideSHA256(),
            "a7ca53e6061c3aa16d1932e35f908fa61b328b88625e8361e63881bfc8592000");
  EXPECT_EQ(database.expectedChlorideSHA256(),
            "460e8d7148cc76c597064e167b4afa98ea762e59a666f9c42962f3d2d013e5c1");
  EXPECT_EQ(database.baseModelSourceSHA256(),
            "bfad8e09ada73900b800c4df262ae70192296e6c7b95cf5fc20943272637d0e6");
  EXPECT_EQ(database.baseModelParameters().size(), 14u);
  EXPECT_EQ(database.baseModelDensities().size(), 12u);
  EXPECT_EQ(database.baseModelElementProperties().size(), 3u);
  EXPECT_DOUBLE_EQ(database.baseModelParameter("Ea_corr_kJ_mol"), 54.83553514725915);
  EXPECT_DOUBLE_EQ(database.baseModelParameter("redox_oxidizing_fef2"), 1.5026092214000561);
  EXPECT_DOUBLE_EQ(database.baseModelDensity("stainless_316h"), 8.0);
  EXPECT_DOUBLE_EQ(database.baseModelElementProperty("Cr", "molar_mass_g_mol"), 51.9961);
  EXPECT_EQ(database.thermochemicalParameters().size(), 28u);
  EXPECT_EQ(database.dynamicParameters().size(), 26u);

  // Corrected-data anchors: these distinguish the 2026-08-28 refit from the stale extension.
  EXPECT_DOUBLE_EQ(database.thermochemicalParameter("log_front_rate0_um_y"), 0.0036741915643164766);
  EXPECT_DOUBLE_EQ(database.thermochemicalParameter("log_gamma_cr_flinak"), -8.056846922178943);
  EXPECT_DOUBLE_EQ(database.dynamicParameter("log_rate_scale"), 0.055365471992740514);
  EXPECT_DOUBLE_EQ(database.dynamicParameter("log_surface_reservoir_um"), 3.321684647830804);

  const Corrosion::MoltenSaltCorrosionDatabase base(
      Moose::DataFileUtils::getPath("corrosion_database.json").path);
  database.validateBaseModel(base);
}

TEST(AdvancedCorrosionModelData, RejectsLegacyModelLawBeforeParametersAreUsed)
{
  EXPECT_THROW(Corrosion::AdvancedCorrosionModelDatabase(
                   fixture("advanced_corrosion_models_legacy_revision.json")),
               std::exception);
}

TEST(AdvancedCorrosionModelData, RejectsMutatedOrImplicitBaseSemantics)
{
  const Corrosion::AdvancedCorrosionModelDatabase database(
      Moose::DataFileUtils::getPath("advanced_corrosion_models.json").path);

  const Corrosion::MoltenSaltCorrosionDatabase mutated_parameter(
      fixture("corrosion_database_mstdb_mutated.json"));
  EXPECT_THROW(database.validateBaseModel(mutated_parameter), std::exception);

  const Corrosion::MoltenSaltCorrosionDatabase mutated_element(
      fixture("corrosion_database_mutated_element.json"));
  EXPECT_THROW(database.validateBaseModel(mutated_element), std::exception);

  // element() has a generic_metal fallback.  Provenance requires explicit Cr/Fe/Ni records even
  // when generic_metal has been made numerically identical to the missing record.
  const Corrosion::MoltenSaltCorrosionDatabase implicit_chromium(
      fixture("corrosion_database_missing_explicit_cr.json"));
  EXPECT_THROW(database.validateBaseModel(implicit_chromium), std::exception);
}
