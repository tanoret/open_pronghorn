//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MSTDBTCStandardStateCorrosionModel.h"

#include "DataFileUtils.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace
{
constexpr const char * fluoride_fixture_sha256 =
    "e713af15f7172861bd35feef50352202838730bbe23696e3b1421d1ee0aa6fe6";
constexpr const char * chloride_fixture_sha256 =
    "4922c0b0ca3fc7655bf644d7b5a7ebbe6d6cb6c43a253a16ca4b19cb457c8bcb";

std::string
fixture(const std::string & name)
{
  const std::string source_path(__FILE__);
  const auto separator = source_path.find_last_of("/\\");
  const std::string source_relative =
      separator == std::string::npos
          ? ""
          : source_path.substr(0, separator) + "/../data/mstdb/" + name;
  const std::array<std::string, 4> candidates = {
      source_relative, "unit/data/mstdb/" + name, "../unit/data/mstdb/" + name, name};
  for (const auto & candidate : candidates)
  {
    if (!candidate.empty() && std::ifstream(candidate).good())
      return candidate;
  }
  throw std::runtime_error("Unable to locate synthetic corrosion fixture '" + name + "'.");
}

std::string
advancedFixture(const std::string & name)
{
  const std::string source_path(__FILE__);
  const auto separator = source_path.find_last_of("/\\");
  const std::string source_relative =
      separator == std::string::npos
          ? ""
          : source_path.substr(0, separator) + "/../data/advanced/" + name;
  const std::array<std::string, 4> candidates = {
      source_relative, "unit/data/advanced/" + name, "../unit/data/advanced/" + name, name};
  for (const auto & candidate : candidates)
  {
    if (!candidate.empty() && std::ifstream(candidate).good())
      return candidate;
  }
  throw std::runtime_error("Unable to locate advanced corrosion fixture '" + name + "'.");
}

class TestAdvancedCorrosionModelDatabase : public Corrosion::AdvancedCorrosionModelDatabase
{
public:
  explicit TestAdvancedCorrosionModelDatabase(const std::string & filename)
    : AdvancedCorrosionModelDatabase(filename)
  {
  }

  void bindMSTDB(const std::string & version,
                 const std::string & fluoride_sha256,
                 const std::string & chloride_sha256)
  {
    _expected_mstdb_version = version;
    _expected_fluoride_sha256 = fluoride_sha256;
    _expected_chloride_sha256 = chloride_sha256;
  }
};

struct ModelFixture
{
  ModelFixture()
    : base(Moose::DataFileUtils::getPath("corrosion_database.json").path),
      advanced(Moose::DataFileUtils::getPath("advanced_corrosion_models.json").path),
      pair(fixture("MSTDB-TC_V0.2_Fluorides_No_Func.dat"),
           fixture("MSTDB-TC_V0.2_Chlorides_No_Func.dat"),
           "0.2"),
      model(base, bindAndReturnAdvanced(), pair)
  {
  }

  TestAdvancedCorrosionModelDatabase & bindAndReturnAdvanced()
  {
    advanced.bindMSTDB("0.2", fluoride_fixture_sha256, chloride_fixture_sha256);
    return advanced;
  }

  Corrosion::MoltenSaltCorrosionDatabase base;
  TestAdvancedCorrosionModelDatabase advanced;
  Corrosion::MSTDBTCPair pair;
  Corrosion::MSTDBTCStandardStateCorrosionModel model;
};
} // namespace

TEST(MSTDBTCStandardStateCorrosionModel, CorrectedFeIdentityAffinity)
{
  ModelFixture f;
  Corrosion::MSTDBTCCorrosionFeatures input;
  input.material_class = "stainless_316h";
  input.salt_class = "flinak";
  input.redox_class = "purified_baseline";
  input.hot_temperature_K = 900.0;
  input.exposure_s = Corrosion::seconds_per_year;
  EXPECT_DOUBLE_EQ(f.model.reactionLogKOverQ(
                       Corrosion::MSTDBTCStandardStateCorrosionModel::Fe, input, 900.0),
                   0.0);

  input.salt_class = "flibe";
  EXPECT_DOUBLE_EQ(f.model.reactionLogKOverQ(
                       Corrosion::MSTDBTCStandardStateCorrosionModel::Fe, input, 900.0),
                   0.0);
  input.salt_class = "chloride";
  EXPECT_DOUBLE_EQ(f.model.reactionLogKOverQ(
                       Corrosion::MSTDBTCStandardStateCorrosionModel::Fe, input, 900.0),
                   0.0);
  input.redox_class = "oxidizing_fef2";
  EXPECT_DOUBLE_EQ(f.model.reactionLogKOverQ(
                       Corrosion::MSTDBTCStandardStateCorrosionModel::Fe, input, 900.0, "", 500.0),
                   0.0);
}

TEST(MSTDBTCStandardStateCorrosionModel, FractionsCaptureAndInventoryArePhysical)
{
  ModelFixture f;
  Corrosion::MSTDBTCCorrosionFeatures input;
  input.material_class = "stainless_316h";
  input.salt_class = "flinak";
  input.redox_class = "oxidizing_fef2";
  input.hot_temperature_K = 900.0;
  input.cold_temperature_K = 800.0;
  input.exposure_s = 2.0 * Corrosion::seconds_per_year;
  input.flow_factor = 1.2;
  input.area_to_salt_mass_cm2_g = 0.25;
  input.inventory_coupling_factor = 0.5;

  const auto output = f.model.evaluate(input);
  EXPECT_NEAR(std::accumulate(output.source_fraction.begin(), output.source_fraction.end(), 0.0),
              1.0,
              1.0e-14);
  EXPECT_NEAR(std::accumulate(output.deposit_fraction.begin(), output.deposit_fraction.end(), 0.0),
              1.0,
              1.0e-14);
  EXPECT_GT(output.front_rate_um_y, 0.0);
  EXPECT_GT(output.front_depth_um, 0.0);
  EXPECT_GE(output.mass_loss_mg_cm2, 0.0);
  EXPECT_GE(output.mass_gain_mg_cm2, 0.0);
  EXPECT_TRUE(std::isfinite(output.igc_depth_um));
  for (unsigned int i = 0; i < 3; ++i)
  {
    EXPECT_GE(output.source_fraction[i], 0.0);
    EXPECT_GE(output.cold_capture_fraction[i], 0.0);
    EXPECT_LE(output.cold_capture_fraction[i], 1.0);
    EXPECT_GE(output.dissolved_inventory_ppm[i], 0.0);
  }

  input.cold_temperature_K = input.hot_temperature_K;
  const auto isothermal = f.model.evaluate(input);
  for (const auto value : isothermal.cold_capture_fraction)
    EXPECT_DOUBLE_EQ(value, 0.0);
  for (const auto value : isothermal.deposit_fraction)
    EXPECT_DOUBLE_EQ(value, 0.0);
  EXPECT_DOUBLE_EQ(isothermal.mass_gain_mg_cm2, 0.0);
}

TEST(MSTDBTCStandardStateCorrosionModel, ZeroExposureHasZeroExtentsAndFiniteDiagnostics)
{
  ModelFixture f;
  Corrosion::MSTDBTCCorrosionFeatures input;
  input.material_class = "stainless_316h";
  input.salt_class = "flinak";
  input.redox_class = "purified_baseline";
  input.hot_temperature_K = 900.0;
  input.cold_temperature_K = 800.0;
  input.exposure_s = 0.0;
  input.flow_factor = 1.0e-4;
  input.area_to_salt_mass_cm2_g = 0.25;
  input.inventory_coupling_factor = 0.5;

  const auto output = f.model.evaluate(input);
  EXPECT_TRUE(std::isfinite(output.front_rate_um_y));
  EXPECT_GT(output.front_rate_um_y, 0.0);
  EXPECT_TRUE(std::isfinite(output.corrosion_rate_um_y));
  EXPECT_TRUE(std::isfinite(output.cr_diffusion_cm2_s));
  EXPECT_TRUE(std::isfinite(output.redox_acceleration_ratio));
  EXPECT_DOUBLE_EQ(output.front_depth_um, 0.0);
  EXPECT_DOUBLE_EQ(output.mass_loss_mg_cm2, 0.0);
  EXPECT_DOUBLE_EQ(output.mass_gain_mg_cm2, 0.0);
  EXPECT_DOUBLE_EQ(output.igc_depth_um, 0.0);
  EXPECT_DOUBLE_EQ(output.fe2_decrease_ppm, 0.0);
  for (unsigned int i = 0; i < 3; ++i)
  {
    EXPECT_TRUE(std::isfinite(output.affinity_log_k_over_q[i]));
    EXPECT_TRUE(std::isfinite(output.source_fraction[i]));
    EXPECT_TRUE(std::isfinite(output.saturation_activity_hot[i]));
    EXPECT_TRUE(std::isfinite(output.cold_capture_fraction[i]));
    EXPECT_TRUE(std::isfinite(output.deposit_fraction[i]));
    EXPECT_DOUBLE_EQ(output.dissolved_inventory_ppm[i], 0.0);
  }
}

TEST(MSTDBTCStandardStateCorrosionModel, RejectsInvalidRuntimeClassesAndState)
{
  ModelFixture f;
  Corrosion::MSTDBTCCorrosionFeatures input;
  input.exposure_s = 1.0;
  input.area_to_salt_mass_cm2_g = 0.25;

  input.salt_class = "unknown_salt";
  EXPECT_THROW(f.model.evaluate(input), std::exception);
  input.salt_class = "fluoride_fuel";
  input.material_class = "unknown_alloy";
  EXPECT_THROW(f.model.evaluate(input), std::exception);
  input.material_class = "hastelloy_n";
  input.redox_class = "unknown_redox";
  EXPECT_THROW(f.model.evaluate(input), std::exception);
  input.redox_class = "purified_baseline";
  input.hot_temperature_K = 0.0;
  EXPECT_THROW(f.model.evaluate(input), std::exception);
  input.hot_temperature_K = 900.0;
  input.exposure_s = -1.0;
  EXPECT_THROW(f.model.evaluate(input), std::exception);
  input.exposure_s = 1.0;

  // The core model has an explicit-geometry contract; wrappers are not the only validation layer.
  input.area_to_salt_mass_cm2_g = std::numeric_limits<Real>::quiet_NaN();
  EXPECT_THROW(f.model.evaluate(input), std::exception);
  input.area_to_salt_mass_cm2_g = 0.0;
  EXPECT_THROW(f.model.evaluate(input), std::exception);
  input.area_to_salt_mass_cm2_g = std::numeric_limits<Real>::infinity();
  EXPECT_THROW(f.model.evaluate(input), std::exception);

  input.area_to_salt_mass_cm2_g = 0.25;
  input.cold_temperature_K = 901.0;
  EXPECT_THROW(f.model.evaluate(input), std::exception);
  input.cold_temperature_K = std::numeric_limits<Real>::infinity();
  EXPECT_THROW(f.model.evaluate(input), std::exception);
}

TEST(MSTDBTCStandardStateCorrosionModel, CoreRejectsUnboundThermodynamicProvenance)
{
  Corrosion::MoltenSaltCorrosionDatabase base(
      Moose::DataFileUtils::getPath("corrosion_database.json").path);
  TestAdvancedCorrosionModelDatabase advanced(
      Moose::DataFileUtils::getPath("advanced_corrosion_models.json").path);
  const Corrosion::MSTDBTCPair pair(fixture("MSTDB-TC_V0.2_Fluorides_No_Func.dat"),
                                    fixture("MSTDB-TC_V0.2_Chlorides_No_Func.dat"),
                                    "0.2");

  // The production artifact is V3.1, so the synthetic V0.2 pair must fail at the core boundary.
  EXPECT_THROW(Corrosion::MSTDBTCStandardStateCorrosionModel(base, advanced, pair), std::exception);

  advanced.bindMSTDB("0.2", std::string(64, '0'), chloride_fixture_sha256);
  EXPECT_THROW(Corrosion::MSTDBTCStandardStateCorrosionModel(base, advanced, pair), std::exception);

  advanced.bindMSTDB("0.2", fluoride_fixture_sha256, std::string(64, '0'));
  EXPECT_THROW(Corrosion::MSTDBTCStandardStateCorrosionModel(base, advanced, pair), std::exception);
}

TEST(MSTDBTCStandardStateCorrosionModel, CoreRejectsMismatchedBaseModelSemantics)
{
  const Corrosion::MoltenSaltCorrosionDatabase mutated_base(
      advancedFixture("corrosion_database_mstdb_mutated.json"));
  TestAdvancedCorrosionModelDatabase advanced(
      Moose::DataFileUtils::getPath("advanced_corrosion_models.json").path);
  advanced.bindMSTDB("0.2", fluoride_fixture_sha256, chloride_fixture_sha256);
  const Corrosion::MSTDBTCPair pair(fixture("MSTDB-TC_V0.2_Fluorides_No_Func.dat"),
                                    fixture("MSTDB-TC_V0.2_Chlorides_No_Func.dat"),
                                    "0.2");
  EXPECT_THROW(Corrosion::MSTDBTCStandardStateCorrosionModel(mutated_base, advanced, pair),
               std::exception);
}

TEST(MSTDBTCStandardStateCorrosionModel, CoreRejectsMismatchedBaseElementSemantics)
{
  const Corrosion::MoltenSaltCorrosionDatabase mutated_base(
      advancedFixture("corrosion_database_mutated_element.json"));
  TestAdvancedCorrosionModelDatabase advanced(
      Moose::DataFileUtils::getPath("advanced_corrosion_models.json").path);
  advanced.bindMSTDB("0.2", fluoride_fixture_sha256, chloride_fixture_sha256);
  const Corrosion::MSTDBTCPair pair(fixture("MSTDB-TC_V0.2_Fluorides_No_Func.dat"),
                                    fixture("MSTDB-TC_V0.2_Chlorides_No_Func.dat"),
                                    "0.2");
  EXPECT_THROW(Corrosion::MSTDBTCStandardStateCorrosionModel(mutated_base, advanced, pair),
               std::exception);
}
