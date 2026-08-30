//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MSTDBTCData.h"
#include "MSTDBTCStandardStateCorrosionModel.h"

#include "DataFileUtils.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>

namespace
{
std::string
join(const std::string & directory, const std::string & filename)
{
  return directory.empty() || directory.back() == '/' ? directory + filename
                                                       : directory + "/" + filename;
}

constexpr const char * fluoride_sha256 =
    "a7ca53e6061c3aa16d1932e35f908fa61b328b88625e8361e63881bfc8592000";
constexpr const char * chloride_sha256 =
    "460e8d7148cc76c597064e167b4afa98ea762e59a666f9c42962f3d2d013e5c1";

void
expectPythonParity(const Real actual, const Real expected)
{
  // Cross-language tolerance: the C++ and CPython implementations evaluate the same expressions,
  // but polynomial/reaction summation and libm transcendental rounding need not be bit-identical.
  EXPECT_NEAR(actual, expected, 2.0e-10 * std::max(1.0, std::abs(expected)));
}

struct EndpointReference
{
  Real front_rate_um_y;
  Real corrosion_rate_um_y;
  Real front_depth_um;
  Real mass_loss_mg_cm2;
  Real mass_gain_mg_cm2;
  Real igc_depth_um;
  Real cr_diffusion_cm2_s;
  Real fe2_decrease_ppm;
  Real redox_log_shift;
  Real redox_acceleration_ratio;
  std::array<Real, 3> affinity;
  std::array<Real, 3> source;
  std::array<Real, 3> saturation;
  std::array<Real, 3> capture;
  std::array<Real, 3> deposit;
  std::array<Real, 3> inventory;
};

void
expectEndpointParity(const Corrosion::MSTDBTCCorrosionResult & actual,
                     const EndpointReference & expected)
{
  expectPythonParity(actual.front_rate_um_y, expected.front_rate_um_y);
  expectPythonParity(actual.corrosion_rate_um_y, expected.corrosion_rate_um_y);
  expectPythonParity(actual.front_depth_um, expected.front_depth_um);
  expectPythonParity(actual.mass_loss_mg_cm2, expected.mass_loss_mg_cm2);
  expectPythonParity(actual.mass_gain_mg_cm2, expected.mass_gain_mg_cm2);
  expectPythonParity(actual.igc_depth_um, expected.igc_depth_um);
  expectPythonParity(actual.cr_diffusion_cm2_s, expected.cr_diffusion_cm2_s);
  expectPythonParity(actual.fe2_decrease_ppm, expected.fe2_decrease_ppm);
  expectPythonParity(actual.redox_log_shift, expected.redox_log_shift);
  expectPythonParity(actual.redox_acceleration_ratio, expected.redox_acceleration_ratio);
  for (unsigned int i = 0; i < 3; ++i)
  {
    expectPythonParity(actual.affinity_log_k_over_q[i], expected.affinity[i]);
    expectPythonParity(actual.source_fraction[i], expected.source[i]);
    expectPythonParity(actual.saturation_activity_hot[i], expected.saturation[i]);
    expectPythonParity(actual.cold_capture_fraction[i], expected.capture[i]);
    expectPythonParity(actual.deposit_fraction[i], expected.deposit[i]);
    expectPythonParity(actual.dissolved_inventory_ppm[i], expected.inventory[i]);
  }
}
} // namespace

TEST(MSTDBTCV31Integration, AuthorizedDatabaseAnchors)
{
  const char * directory_environment = std::getenv("MSTDB_TC_DIR");
  if (!directory_environment || !*directory_environment)
    GTEST_SKIP() << "Set MSTDB_TC_DIR to an authorized MSTDB-TC V3.1 directory.";

  const std::string directory(directory_environment);
  const auto fluoride = join(directory, "MSTDB-TC_V3.1_Fluorides_No_Func.dat");
  const auto chloride = join(directory, "MSTDB-TC_V3.1_Chlorides_No_Func.dat");
  if (!std::ifstream(fluoride).good() || !std::ifstream(chloride).good())
    FAIL() << "MSTDB_TC_DIR is set but does not contain the expected V3.1 function-expanded pair.";

  const Corrosion::MSTDBTCPair pair(fluoride, chloride, "3.1", fluoride_sha256, chloride_sha256);
  EXPECT_EQ(pair.fluoride().recordCount(), 677);
  EXPECT_EQ(pair.chloride().recordCount(), 649);

  constexpr Real temperature_K = 923.15;
  EXPECT_NEAR(pair.fluoride().standardGibbsJMol("CrF2_L1(liq)", temperature_K),
              -881852.8502324164,
              1.0e-5);
  EXPECT_NEAR(pair.fluoride().standardGibbsJMol("Cr_S1(s)", temperature_K),
              -32479.385685162837,
              1.0e-5);
  EXPECT_NEAR(pair.chloride().standardGibbsJMol("CrCl2_L1(liq)", temperature_K),
              -528859.0468304644,
              1.0e-5);

  const std::map<std::string, Real> uranium_buffer{{"CrF2_L1(liq)", 1.0},
                                                   {"UF3_L1(liq)", 2.0},
                                                   {"Cr_S1(s)", -1.0},
                                                   {"UF4_L1(liq)", -2.0}};
  const std::map<std::string, Real> iron_buffer{{"CrF2_L1(liq)", 1.0},
                                                {"Fe_bcc(s)", 1.0},
                                                {"Cr_S1(s)", -1.0},
                                                {"FeF2_L1(liq)", -1.0}};
  EXPECT_NEAR(pair.fluoride().equilibriumLogConstant(uranium_buffer, temperature_K),
              -3.2900422614422524,
              1.0e-10);
  EXPECT_NEAR(pair.fluoride().equilibriumLogConstant(iron_buffer, temperature_K),
              9.547764169406722,
              1.0e-10);
}

TEST(MSTDBTCV31Integration, FullEndpointPythonParity)
{
  const char * directory_environment = std::getenv("MSTDB_TC_DIR");
  if (!directory_environment || !*directory_environment)
    GTEST_SKIP() << "Set MSTDB_TC_DIR to an authorized MSTDB-TC V3.1 directory.";

  const std::string directory(directory_environment);
  const auto fluoride = join(directory, "MSTDB-TC_V3.1_Fluorides_No_Func.dat");
  const auto chloride = join(directory, "MSTDB-TC_V3.1_Chlorides_No_Func.dat");
  if (!std::ifstream(fluoride).good() || !std::ifstream(chloride).good())
    FAIL() << "MSTDB_TC_DIR is set but does not contain the expected V3.1 function-expanded pair.";

  const Corrosion::MSTDBTCPair pair(fluoride, chloride, "3.1", fluoride_sha256, chloride_sha256);
  const Corrosion::MoltenSaltCorrosionDatabase base(
      Moose::DataFileUtils::getPath("corrosion_database.json").path);
  const Corrosion::AdvancedCorrosionModelDatabase advanced(
      Moose::DataFileUtils::getPath("advanced_corrosion_models.json").path);
  const Corrosion::MSTDBTCStandardStateCorrosionModel model(base, advanced, pair);

  // Frozen with Python 3.12.13 using this calibration overlay against OpenPronghorn base commit
  // bbc8b69c87bc1c90af0d1cf1484c599e9553b625.
  // Input artifact SHA-256 values:
  //   data/corrosion_database.json
  //     bfad8e09ada73900b800c4df262ae70192296e6c7b95cf5fc20943272637d0e6
  //   validation/corrosion/calibration/data/parameters.json (Python base-model reference)
  //     e7d262098ba929250e11b03fb550027a83daeaf9d5292657fe11dc67962602f8
  //   data/advanced_corrosion_models.json
  //     a7b8a3d85bc6584a42e36adfbaefe648015a7ad462357b7408de90899c002add
  // These vectors test implementation parity, not independent scientific validation or native
  // SUBQ equilibrium.
  Corrosion::MSTDBTCCorrosionFeatures fuel;
  fuel.material_class = "hastelloy_n";
  fuel.salt_class = "fluoride_fuel";
  fuel.redox_class = "msre_or_fuel_baseline";
  fuel.hot_temperature_K = 923.15;
  fuel.exposure_s = Corrosion::seconds_per_year;
  fuel.flow_factor = 1.0;
  fuel.area_to_salt_mass_cm2_g = 0.25;
  fuel.inventory_coupling_factor = 1.0;
  const EndpointReference fuel_reference{
      3.0540999703048093,
      2.0621740542864484,
      3.0540999703048093,
      1.833272734260653,
      0.0,
      10.48888490740358,
      3.9999999999999956e-14,
      99.99999999999996,
      1.116179018407099,
      6.564672081661867,
      {{19.069291071872208, 9.155522993037723, 0.17701853981425053}},
      {{0.46353330991278946, 0.11339764168974246, 0.42306904839746795}},
      {{0.6276085980934543, 0.21173939475412973, 0.028318322313217913}},
      {{0.0, 0.0, 0.0}},
      {{0.0, 0.0, 0.0}},
      {{215.7381632724302, 52.777650309905376, 196.90524388820478}}};
  const auto fuel_output = model.evaluate(fuel);
  expectEndpointParity(fuel_output, fuel_reference);

  Corrosion::MSTDBTCCorrosionFeatures flinak;
  flinak.material_class = "stainless_316h";
  flinak.salt_class = "flinak";
  flinak.redox_class = "purified_baseline";
  flinak.hot_temperature_K = 923.15;
  flinak.cold_temperature_K = 813.15;
  flinak.exposure_s = 1000.0 * 3600.0;
  flinak.flow_factor = 1.1;
  flinak.area_to_salt_mass_cm2_g = 7.767901069021746;
  flinak.inventory_coupling_factor = 1.0;
  const EndpointReference flinak_reference{
      43.26794829499689,
      15.213737670349811,
      4.9358827623770125,
      1.388431455199618,
      1.005872818676711,
      7.4470039777428525,
      3.9999999999999956e-14,
      39.711979160956204,
      0.0,
      5.311269603591153,
      {{12.827202071490573, 0.0, -8.730146882821607}},
      {{0.42125633489919373, 0.552327340154541, 0.0264163249462652}},
      {{0.6276085980934543, 0.21173939475412973, 0.028318322313217913}},
      {{0.22722708707182404, 0.533224724710303, 0.6336966207059662}},
      {{0.23520059118452535, 0.7236668571004717, 0.04113255171500289}},
      {{406.8405368923055, 414.31846951698856, 17.94530018299972}}};
  const auto flinak_output = model.evaluate(flinak);
  expectEndpointParity(flinak_output, flinak_reference);

  auto flinak_half_coupling = flinak;
  flinak_half_coupling.inventory_coupling_factor = 0.5;
  const EndpointReference flinak_half_coupling_reference{
      43.26794829499689,
      15.213737670349811,
      4.9358827623770125,
      1.388431455199618,
      1.005872818676711,
      7.4470039777428525,
      3.9999999999999956e-14,
      39.711930944041484,
      0.0,
      5.311269603591153,
      {{12.827202071490573, 0.0, -8.730146882821607}},
      {{0.42125633489919373, 0.552327340154541, 0.0264163249462652}},
      {{0.6276085980934543, 0.21173939475412973, 0.028318322313217913}},
      {{0.22722708707182404, 0.533224724710303, 0.6336966207059662}},
      {{0.23520059118452535, 0.7236668571004717, 0.04113255171500289}},
      {{406.8400429205658, 414.31796646580113, 17.945278394435533}}};
  const auto flinak_half_coupling_output = model.evaluate(flinak_half_coupling);
  expectEndpointParity(flinak_half_coupling_output, flinak_half_coupling_reference);
  // Static MSTDB v2 couples this factor only to dissolved inventory.  Mass gain is a separate
  // cold-capture closure; changing this distinction requires a new model revision and validation.
  EXPECT_DOUBLE_EQ(flinak_half_coupling_output.mass_gain_mg_cm2, flinak_output.mass_gain_mg_cm2);
  for (unsigned int i = 0; i < 3; ++i)
    EXPECT_LT(flinak_half_coupling_output.dissolved_inventory_ppm[i],
              flinak_output.dissolved_inventory_ppm[i]);

  auto flinak_small_flow = flinak;
  flinak_small_flow.flow_factor = 1.0e-4;
  const EndpointReference flinak_small_flow_reference{
      0.12868633748180805,
      0.04524827863013588,
      0.014680166265321476,
      0.004129434508796339,
      0.0029916391719220074,
      2.525801381631161,
      3.9999999999999956e-14,
      3.0907615293996726,
      0.0,
      1.002420040947854,
      {{12.827202071490573, 0.0, -8.730146882821607}},
      {{0.42125633489919373, 0.552327340154541, 0.0264163249462652}},
      {{0.6276085980934543, 0.21173939475412973, 0.028318322313217913}},
      {{0.22722708707182404, 0.533224724710303, 0.6336966207059662}},
      {{0.23520059118452535, 0.7236668571004717, 0.04113255171500289}},
      {{31.664175560993883, 32.246178950503534, 1.3966728581907577}}};
  expectEndpointParity(model.evaluate(flinak_small_flow), flinak_small_flow_reference);

  Corrosion::MSTDBTCCorrosionFeatures chloride_case;
  chloride_case.material_class = "stainless_316h";
  chloride_case.salt_class = "chloride";
  chloride_case.redox_class = "impure_moisture";
  chloride_case.hot_temperature_K = 973.15;
  chloride_case.exposure_s = 1000.0 * 3600.0;
  chloride_case.flow_factor = 0.35;
  chloride_case.area_to_salt_mass_cm2_g = 7.767901069021746;
  chloride_case.inventory_coupling_factor = 1.0;
  const EndpointReference chloride_reference{
      68.04914520009785,
      68.04914520009785,
      7.762850239573107,
      6.210280191658486,
      0.0,
      125.88738627903248,
      1.9942380658724556e-13,
      0.0,
      2.9007472534550707,
      6.17261240211636,
      {{13.248424027675343, 0.0, -3.2721956432591237}},
      {{0.4214547926230208, 0.5514672605738116, 0.02707794680316757}},
      {{0.6459663642400368, 1.0, 0.0942438634454419}},
      {{0.0, 0.0, 0.0}},
      {{0.0, 0.0, 0.0}},
      {{353.64453152446015, 462.73855329282213, 22.721221776248434}}};
  expectEndpointParity(model.evaluate(chloride_case), chloride_reference);

  auto zero_exposure = fuel;
  zero_exposure.exposure_s = 0.0;
  const EndpointReference zero_exposure_reference{
      3.0540999703048093,
      2.0621740542864484,
      0.0,
      0.0,
      0.0,
      0.0,
      3.9999999999999956e-14,
      0.0,
      1.116179018407099,
      6.564672081661867,
      {{19.069291071872208, 9.155522993037723, 0.17701853981425053}},
      {{0.46353330991278946, 0.11339764168974246, 0.42306904839746795}},
      {{0.6276085980934543, 0.21173939475412973, 0.028318322313217913}},
      {{0.0, 0.0, 0.0}},
      {{0.0, 0.0, 0.0}},
      {{0.0, 0.0, 0.0}}};
  const auto zero_exposure_output = model.evaluate(zero_exposure);
  expectEndpointParity(zero_exposure_output, zero_exposure_reference);
  EXPECT_TRUE(std::isfinite(zero_exposure_output.front_rate_um_y));
  EXPECT_TRUE(std::isfinite(zero_exposure_output.corrosion_rate_um_y));
  EXPECT_TRUE(std::isfinite(zero_exposure_output.redox_acceleration_ratio));
  EXPECT_DOUBLE_EQ(zero_exposure_output.front_depth_um, 0.0);
  EXPECT_DOUBLE_EQ(zero_exposure_output.mass_loss_mg_cm2, 0.0);
  EXPECT_DOUBLE_EQ(zero_exposure_output.mass_gain_mg_cm2, 0.0);
  EXPECT_DOUBLE_EQ(zero_exposure_output.igc_depth_um, 0.0);
  for (const auto inventory : zero_exposure_output.dissolved_inventory_ppm)
    EXPECT_DOUBLE_EQ(inventory, 0.0);

  for (const auto & salt : {"flinak", "flibe", "fluoroborate", "chloride"})
  {
    auto identity = flinak;
    identity.salt_class = salt;
    identity.redox_class = "oxidizing_fef2";
    EXPECT_DOUBLE_EQ(model.reactionLogKOverQ(
                         Corrosion::MSTDBTCStandardStateCorrosionModel::Fe,
                         identity,
                         identity.hot_temperature_K,
                         "",
                         500.0),
                     0.0);
  }
}
