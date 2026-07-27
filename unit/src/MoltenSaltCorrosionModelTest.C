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
#include "MoltenSaltCorrosionModel.h"
#include "DataFileUtils.h"

#include <cmath>
#include <limits>

// These tests verify that the C++ port of the effective Butler-Volmer correlation reproduces the
// reference predictions in msr_corrosion_plating_model/results/case_predictions_all_76_cases.csv and
// validation_predictions.csv term for term. The gold values were extracted directly from those CSVs.

namespace
{
Corrosion::MoltenSaltCorrosionDatabase
loadDatabase()
{
  return Corrosion::MoltenSaltCorrosionDatabase(
      Moose::DataFileUtils::getPath("corrosion_database.json").path);
}

// Relative tolerance for double-precision reproduction of the reference Python model.
void
expectClose(Real actual, Real expected)
{
  EXPECT_NEAR(actual, expected, std::abs(expected) * 1.0e-9 + 1.0e-15);
}

const Real NaN = std::numeric_limits<Real>::quiet_NaN();
}

TEST(MoltenSaltCorrosionModel, msreBaselineCase)
{
  // Reference case MSRE-SURV-01: Hastelloy N in fuel fluoride, MSRE baseline redox.
  const auto db = loadDatabase();
  Corrosion::MoltenSaltCorrosionModel model(db);

  Corrosion::CorrosionFeatures f;
  f.material_class = "hastelloy_n";
  f.salt_class = "fluoride_fuel";
  f.redox_class = "msre_or_fuel_baseline";
  f.position_class = "core_or_reactor";
  f.temperature_K = 918.15;
  f.flow_factor = 2.0;
  f.delta_T_C = 0.0;
  f.time_years = 0.5475701574264202;
  f.experiment_family = "MSRE reactor evidence";
  f.source_id = "ORNL-TM-1997";

  expectClose(model.corrosionRateUmY(f), 2.880225030371104);
  expectClose(model.corrosionDepthUm(f), 1.5771252733038215);
  expectClose(model.depositionRateUmY(f), 0.8928472079269775);
  expectClose(model.depositionDepthUm(f), 0.4888964862023148);
  expectClose(model.saltCrPpm(f), 54.17923675342637);
  expectClose(model.bvOverpotentialEquivalentV(f.redox_class, f.temperature_K), 0.05693796773823287);
}

TEST(MoltenSaltCorrosionModel, stainlessHotLoopCase)
{
  // Reference case ORNL-FL-01: 304L stainless, hot fluoride loop with a large thermal gradient.
  const auto db = loadDatabase();
  Corrosion::MoltenSaltCorrosionModel model(db);

  Corrosion::CorrosionFeatures f;
  f.material_class = "stainless_304l";
  f.salt_class = "fluoride_fuel";
  f.redox_class = "msre_or_fuel_baseline";
  f.position_class = "nominal";
  f.temperature_K = 961.15;
  f.flow_factor = 1.0;
  f.delta_T_C = 100.0;
  f.time_years = 9.52726443075519;
  f.experiment_family = "ORNL fluoride loops";
  f.source_id = "ORNL-TM-4286";

  expectClose(model.corrosionRateUmY(f), 16.761260038719758);
  expectClose(model.corrosionDepthUm(f), 159.68895658153312);
  expectClose(model.depositionRateUmY(f), 0.8527971023912866);
  expectClose(model.saltCrPpm(f), 71143.79021799257);
}

TEST(MoltenSaltCorrosionModel, noSaltGasControlIsGated)
{
  // Reference case MSRE-SURV-03: no-salt gas control. Corrosion is gated to a negligible floor and
  // deposition to its no-salt floor, while the effective overpotential follows the gas-control redox.
  const auto db = loadDatabase();
  Corrosion::MoltenSaltCorrosionModel model(db);

  Corrosion::CorrosionFeatures f;
  f.material_class = "hastelloy_n";
  f.salt_class = "no_salt";
  f.redox_class = "gas_control";
  f.position_class = "gas_or_offgas";
  f.temperature_K = 923.15;
  f.flow_factor = 0.02;
  f.time_years = 1.2548482774355465;

  expectClose(model.corrosionRateUmY(f), 1.0e-6);
  expectClose(model.depositionRateUmY(f), 1.0e-9);
  expectClose(model.bvOverpotentialEquivalentV(f.redox_class, f.temperature_K), -0.47730547207333357);
}

TEST(MoltenSaltCorrosionModel, telluriumRateWithUnsetTime)
{
  // Reference case MSRE-FP-03: tellurium redox with no exposure time. The rate is finite but the
  // depth is not-a-number because the time is unset.
  const auto db = loadDatabase();
  Corrosion::MoltenSaltCorrosionModel model(db);

  Corrosion::CorrosionFeatures f;
  f.material_class = "hastelloy_n";
  f.salt_class = "fluoride_fuel";
  f.redox_class = "tellurium";
  f.position_class = "core_or_reactor";
  f.temperature_K = 923.15;
  f.flow_factor = 2.0;
  f.time_years = NaN;

  expectClose(model.corrosionRateUmY(f), 0.4535068125750907);
  EXPECT_TRUE(std::isnan(model.corrosionDepthUm(f)));
  expectClose(model.depositionRateUmY(f), 1.6232890444263133);
}

TEST(MoltenSaltCorrosionModel, predictResponseDispatch)
{
  // Each branch of predict_response, checked against validation_predictions.csv.
  const auto db = loadDatabase();
  Corrosion::MoltenSaltCorrosionModel model(db);

  // M-019: igc_depth_um for oxidizing FeF2 on a hot loop.
  {
    Corrosion::CorrosionFeatures f;
    f.material_class = "hastelloy_n";
    f.salt_class = "fluoride_fuel";
    f.redox_class = "oxidizing_fef2";
    f.position_class = "nominal";
    f.temperature_K = 977.15;
    f.flow_factor = 1.0;
    f.delta_T_C = 166.0;
    f.time_years = 0.3308236367784622;
    f.experiment_family = "ORNL fluoride loops";
    f.source_id = "ORNL-TM-4188";
    f.response_kind = "igc_depth_um";
    expectClose(model.predictResponse(f), 10.213975613852675);
    f.response_kind = "redox_acceleration_ratio";
    expectClose(model.predictResponse(f), 6.583640080670622);
  }

  // M-011: mass_loss_mg_cm2 for a purified hot leg.
  {
    Corrosion::CorrosionFeatures f;
    f.material_class = "hastelloy_n";
    f.salt_class = "fluoride_fuel";
    f.redox_class = "purified_baseline";
    f.position_class = "hot_leg";
    f.temperature_K = 977.15;
    f.flow_factor = 1.0;
    f.delta_T_C = 166.0;
    f.time_years = 3.3652749258498744;
    f.experiment_family = "ORNL fluoride loops";
    f.source_id = "ORNL-TM-4188";
    f.response_kind = "mass_loss_mg_cm2";
    expectClose(model.predictResponse(f), 4.134316015558358);
  }

  // M-012: mass_gain_mg_cm2 for a purified cold leg.
  {
    Corrosion::CorrosionFeatures f;
    f.material_class = "hastelloy_n";
    f.salt_class = "fluoride_fuel";
    f.redox_class = "purified_baseline";
    f.position_class = "cold_leg";
    f.temperature_K = 811.15;
    f.flow_factor = 1.0;
    f.delta_T_C = 166.0;
    f.time_years = 3.3652749258498744;
    f.experiment_family = "ORNL fluoride loops";
    f.source_id = "ORNL-TM-4188";
    f.response_kind = "mass_gain_mg_cm2";
    expectClose(model.predictResponse(f), 1.8196552971105424);
  }

  // M-004: salt_cr_ppm for a generic metal at MSRE scale.
  {
    Corrosion::CorrosionFeatures f;
    f.material_class = "generic_metal";
    f.salt_class = "fluoride_fuel";
    f.redox_class = "msre_or_fuel_baseline";
    f.position_class = "core_or_reactor";
    f.temperature_K = 923.15;
    f.flow_factor = 2.0;
    f.time_years = 2.5704996577686514;
    f.experiment_family = "MSRE reactor evidence";
    f.source_id = "ORNL-TM-3063";
    f.response_kind = "salt_cr_ppm";
    expectClose(model.predictResponse(f), 61.39545720410084);
  }

  // M-014/M-015: NCL-16 Cr inventory needs the ORNL-TM-4188 Cr chemistry correction, while Fe
  // remains on the base inventory scale.
  {
    Corrosion::CorrosionFeatures f;
    f.material_class = "generic_metal";
    f.salt_class = "fluoride_fuel";
    f.redox_class = "purified_baseline";
    f.position_class = "nominal";
    f.temperature_K = 977.15;
    f.flow_factor = 1.0;
    f.delta_T_C = 166.0;
    f.time_years = 3.3652749258498744;
    f.experiment_family = "ORNL fluoride loops";
    f.source_id = "ORNL-TM-4188";
    f.response_kind = "salt_cr_ppm";
    expectClose(model.predictResponse(f), 493.7799732444472);
    f.response_kind = "salt_fe_decrease_ppm";
    expectClose(model.predictResponse(f), 96.0399133084755);
  }

  // M-005: cr_diffusion_cm2_s at the reference temperature.
  {
    Corrosion::CorrosionFeatures f;
    f.salt_class = "fluoride_fuel";
    f.redox_class = "msre_or_fuel_baseline";
    f.temperature_K = 923.15;
    f.response_kind = "cr_diffusion_cm2_s";
    expectClose(model.predictResponse(f), 3.9999999999999956e-14);
    // The mechanistic diffusivity is the same value converted to m^2/s.
    expectClose(model.crDiffusivityM2S(f), 3.9999999999999956e-14 * 1.0e-4);
  }

  // M-009: noble-metal deposition ranking (turbulent/graphite ratio).
  {
    Corrosion::CorrosionFeatures f;
    f.material_class = "generic_metal";
    f.salt_class = "fluoride_fuel";
    f.redox_class = "fission_product";
    f.position_class = "nominal";
    f.surface_class = "turbulent_metal";
    f.temperature_K = 923.15;
    f.flow_factor = 2.8;
    f.time_years = 0.3333607118412047;
    f.response_kind = "noble_metal_deposition_ranking";
    expectClose(model.predictResponse(f), 2.857651118063164);
  }

  // M-010 / M-039 / M-042: scalar submodels independent of the case features.
  {
    Corrosion::CorrosionFeatures f;
    f.response_kind = "offgas_fraction_percent";
    expectClose(model.predictResponse(f), 0.09999999999999996);
    f.response_kind = "te_soluble_ppm";
    expectClose(model.predictResponse(f), 2.0);
    f.response_kind = "te_redox_threshold_ratio";
    expectClose(model.predictResponse(f), 149.99999999999997);
  }
}

TEST(MoltenSaltCorrosionModel, ncl16LoopSimulation)
{
  // Reproduce the NCL-16 loop simulation in results/ncl16_simplified_loop_simulation.csv:
  // two segments (hot and cold leg), 29500 h at 250 h steps, 250 cm^3 salt.
  const auto db = loadDatabase();
  Corrosion::MoltenSaltCorrosionModel model(db);

  Corrosion::CorrosionFeatures hot;
  hot.material_class = "hastelloy_n";
  hot.salt_class = "fluoride_fuel";
  hot.redox_class = "purified_baseline";
  hot.position_class = "hot_leg";
  hot.surface_class = "metal";
  hot.temperature_K = 704.0 + 273.15;
  hot.flow_factor = 1.0;
  hot.delta_T_C = 166.0;
  hot.surface_area_cm2 = 100.0;

  Corrosion::CorrosionFeatures cold = hot;
  cold.position_class = "cold_leg";
  cold.temperature_K = 538.0 + 273.15;

  const auto sim = model.simulateLoop({hot, cold}, 29500.0, 250.0, 250.0, 0.0);

  ASSERT_EQ(sim.size(), 119u);
  EXPECT_DOUBLE_EQ(sim.front().time_h, 0.0);
  EXPECT_DOUBLE_EQ(sim.front().salt_cr_ppm, 0.0);
  EXPECT_DOUBLE_EQ(sim.back().time_h, 29500.0);
  expectClose(sim.back().salt_cr_ppm, 0.6774015746681846);
  expectClose(sim.back().total_cr_dissolved_mg, 30.548957222057208);
  expectClose(sim.back().total_deposit_mg, 30.210256434723043);
}
