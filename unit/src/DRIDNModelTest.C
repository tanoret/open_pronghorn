//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "DRIDNModel.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace Corrosion
{
namespace
{

void
expectClose(const Real actual, const Real expected, const Real tolerance = 1.0e-12)
{
  EXPECT_NEAR(actual, expected, tolerance * std::max(1.0, std::abs(expected)));
}

DRIDNModel::Parameters
historicalLegacyAlgebraParameters()
{
  // Historical pre-refit vector retained only for immutable Python/C++ RHS algebra fixtures.
  // It is neither the current AdvancedCorrosionModelDatabase vector nor an endpoint-validation
  // golden.  Production construction must load all 26 current fields from the database.
  DRIDNModel::Parameters parameters;
  parameters.log_rate_scale = 0.1543832786246943;
  parameters.affinity_feedback_scale = 2.8725723958274885e-13;
  parameters.log_surface_reservoir_um = 2.8616007511735835;
  parameters.log_surface_replenishment_y_inv = -0.36593743732588135;
  parameters.surface_availability_exponent = 1.1901514890492675;
  parameters.surface_reservoir_cr_exponent = 0.5550218064425788;
  parameters.log_dynamic_cr_exchange_bias = 0.007;
  parameters.log_dynamic_fe_capture_bias = 0.05;
  parameters.inventory_inhibition_scale = 2.50000000000004e-15;
  parameters.log_redox_relaxation_y_inv = -0.1367446316283494;
  parameters.redox_buffer_retention = 0.769179088690591;
  parameters.redox_consumption_per_um = 4.4182081345454976e-13;
  parameters.log_stress_interfacial_factor = 0.0;
  parameters.log_fluoride_impurity_interfacial_factor = 0.0;
  parameters.log_deposition_rate_y_inv_fuel = 0.8858787714715575;
  parameters.log_deposition_rate_y_inv_flinak = 1.6111042779475175;
  parameters.bulk_capture_area_multiplier_fuel = 1.0e-10;
  parameters.bulk_capture_area_multiplier_flinak = 0.0;
  parameters.log_bulk_precipitation_rate_scale = 1.761132;
  parameters.log_inventory_scale_msre = -0.19170940704877637;
  parameters.log_inventory_scale_loop = -1.4966842815300914;
  parameters.log_deposit_area_scale_fuel = 0.2672160991471665;
  parameters.log_deposit_area_scale_flinak = 1.5529376893117708;
  parameters.log_mass_loss_scale = -0.043176166193242195;
  parameters.gb_dynamic_scale = 0.9767640782979;
  parameters.damage_affinity_scale = 0.02611452873592219;
  return parameters;
}

DRIDNModel::Parameters
currentRefitLegacyParameters()
{
  // Complete v2 refit vector from AdvancedCorrosionModelDatabase.  This is mode-pinned to the
  // legacy irreversible law and used only to exercise its long-horizon boundary trajectory.
  DRIDNModel::Parameters parameters;
  parameters.log_rate_scale = 0.055365471992740514;
  parameters.affinity_feedback_scale = 2.8725723958274885e-13;
  parameters.log_surface_reservoir_um = 3.321684647830804;
  parameters.log_surface_replenishment_y_inv = -0.17253145397117015;
  parameters.surface_availability_exponent = 1.001027788010305;
  parameters.surface_reservoir_cr_exponent = 0.5550218064425788;
  parameters.log_dynamic_cr_exchange_bias = 0.007;
  parameters.log_dynamic_fe_capture_bias = 0.05;
  parameters.inventory_inhibition_scale = 2.50000000000004e-15;
  parameters.log_redox_relaxation_y_inv = 0.0003247206442954783;
  parameters.redox_buffer_retention = 0.7499538506912311;
  parameters.redox_consumption_per_um = 4.4182081345454976e-13;
  parameters.log_stress_interfacial_factor = 0.023049536841576106;
  parameters.log_fluoride_impurity_interfacial_factor = 0.17254172336369716;
  parameters.log_deposition_rate_y_inv_fuel = 0.9415843932550865;
  parameters.log_deposition_rate_y_inv_flinak = 2.9521708093322574;
  parameters.bulk_capture_area_multiplier_fuel = 1.0e-10;
  parameters.bulk_capture_area_multiplier_flinak = 0.0;
  parameters.log_bulk_precipitation_rate_scale = 1.761132;
  parameters.log_inventory_scale_msre = -0.15185382582003074;
  parameters.log_inventory_scale_loop = -0.8819669893066959;
  parameters.log_deposit_area_scale_fuel = 0.511133587761954;
  parameters.log_deposit_area_scale_flinak = 0.6632084287760844;
  parameters.log_mass_loss_scale = -1.2444239037670558e-7;
  parameters.gb_dynamic_scale = 0.9905649525756313;
  parameters.damage_affinity_scale = 0.02611452873592219;
  return parameters;
}

DRIDNModel::Context
syntheticContext()
{
  DRIDNModel::Context context;
  context.mass_fractions = {{0.17, 0.66, 0.12}};
  context.log_exchange_offsets = {{0.0, -0.5, -1.0}};
  context.affinity_baseline = {{10.0, 5.0, 0.0}};
  context.cold_capture_fraction = {{0.2, 0.4, 0.6}};
  context.initial_dissolved_ppm = {{1.0, 1.0, 1.0}};
  context.cr_fraction_ratio = 2.0;
  context.density_g_cm3 = 8.0;
  context.flow_factor = 1.21;
  context.selectivity_scale = 0.02;
  context.product_activity_floor_ppm = 1.0;
  context.redox_shift_initial = 2.0;
  context.log_charge_base_no_redox = 1.0;
  context.mass_transfer_rate_um_y = 100.0;
  context.inventory_capacity_ppm = 1000.0;
  context.area_to_salt_mass_cm2_g = 0.25;
  context.inventory_coupling_factor = 0.5;
  context.deposit_area_factor = 1.5;
  context.mass_loss_fraction = 0.4;
  context.cr_diffusion_cm2_s = 1.0e-14;
  context.front_damage_multiplier = 1.2;
  context.gb_length_multiplier = 0.8;
  context.inventory_scale = DRIDNModel::InventoryScale::Loop;
  context.deposition_closure = DRIDNModel::DepositionClosure::FLiNaK;
  context.transient_redox = true;
  return context;
}

DRIDNModel::Context
equilibriumContext()
{
  DRIDNModel::Context context;
  context.mass_fractions = {{1.0, 0.0, 0.0}};
  context.log_exchange_offsets = {{0.0, 0.0, 0.0}};
  context.affinity_baseline = {{0.0, 0.0, 0.0}};
  context.cold_capture_fraction = {{1.0, 0.0, 0.0}};
  context.initial_dissolved_ppm = {{0.0, 0.0, 0.0}};
  context.cr_fraction_ratio = 1.0;
  context.density_g_cm3 = 8.0;
  context.flow_factor = 1.0;
  context.selectivity_scale = 0.0;
  context.product_activity_floor_ppm = 1.0;
  context.redox_shift_initial = 0.0;
  context.log_charge_base_no_redox = 1.0;
  context.mass_transfer_rate_um_y = 100.0;
  context.inventory_capacity_ppm = 1.0e30;
  context.area_to_salt_mass_cm2_g = 0.25;
  context.inventory_coupling_factor = 1.0;
  context.deposit_area_factor = 1.5;
  context.mass_loss_fraction = 0.4;
  context.cr_diffusion_cm2_s = 0.0;
  context.front_damage_multiplier = 1.0;
  context.gb_length_multiplier = 0.8;
  context.inventory_scale = DRIDNModel::InventoryScale::Explicit;
  context.explicit_inventory_scale = 1.0;
  context.deposition_closure = DRIDNModel::DepositionClosure::FLiNaK;
  return context;
}

} // namespace

TEST(DRIDNModel, ExactLegacySyntheticRHS)
{
  const DRIDNModel model(historicalLegacyAlgebraParameters(),
                         DRIDNModel::ClosureConstants(),
                         DRIDNModel::ModelOptions::legacyCompatibility());
  const auto context = syntheticContext();
  const DRIDNModel::State state = {{0.8,
                                    0.6,
                                    0.4,
                                    10.0,
                                    20.0,
                                    30.0,
                                    100.0,
                                    50.0,
                                    25.0,
                                    0.1,
                                    0.2,
                                    0.3,
                                    4.0,
                                    1.0,
                                    9.0,
                                    1.0,
                                    5.0,
                                    6.0,
                                    7.0}};

  const auto rate = model.rates(context, state);
  const std::array<Real, 3> expected_affinity = {
      {8.6974149070059532, 3.0042677264460091, -2.4011973816621555}};
  const std::array<Real, 3> expected_fraction = {
      {0.38876914395060502, 0.57604104372697618, 0.035189812322418711}};
  for (std::size_t i = 0; i < 3; ++i)
  {
    expectClose(rate.affinity[i], expected_affinity[i]);
    expectClose(rate.species_fraction[i], expected_fraction[i]);
  }
  expectClose(rate.product_feedback, -2.7405262785935407);
  expectClose(rate.surface_feedback, -0.41325195178659135);
  expectClose(rate.charge_transfer_rate_um_y, 5.2727436978494691);
  expectClose(rate.mass_transfer_rate_um_y, 99.999999999999972);
  expectClose(rate.front_rate_um_y, 5.0086503995593894);
  expectClose(rate.bulk_precipitation_rate_y_inv, 0.0073319662067867003);

  const DRIDNModel::State expected_rhs = {{0.06292884032698146,
                                           0.71997126125463,
                                           0.7421705497018496,
                                           32.500735337030875,
                                           18.11135060437946,
                                           -95.43925873240984,
                                           43.592400300490056,
                                           64.59106170948823,
                                           3.9458079663156016,
                                           1.3954093308987612,
                                           5.8678139887491,
                                           12.558683978088853,
                                           5.008650399559389,
                                           1.9187992493050807,
                                           25.18342614173655,
                                           0.46955219513925317,
                                           0.073319662067867,
                                           0.146639324135734,
                                           0.21995898620360102}};
  const auto derivative = model.rhs(context, state);
  for (std::size_t i = 0; i < DRIDNModel::n_states; ++i)
    expectClose(derivative[i], expected_rhs[i]);

  // The shared source and sinks close element by element at the differential level.
  for (const auto element : {DRIDNModel::Element::Cr,
                             DRIDNModel::Element::Fe,
                             DRIDNModel::Element::Ni})
  {
    const auto i = static_cast<std::size_t>(element);
    const Real accounted =
        derivative[DRIDNModel::dissolvedIndex(element)] +
        derivative[DRIDNModel::couponDepositIndex(element)] /
            rate.deposit_conversion_mg_cm2_per_ppm +
        derivative[DRIDNModel::bulkCaptureIndex(element)];
    expectClose(accounted, derivative[DRIDNModel::cumulativeSourceIndex(element)]);
    expectClose(rate.source_rate_ppm_y[i],
                derivative[DRIDNModel::cumulativeSourceIndex(element)]);
  }
}

TEST(DRIDNModel, ExactLegacyInitialRHS)
{
  const DRIDNModel model(historicalLegacyAlgebraParameters(),
                         DRIDNModel::ClosureConstants(),
                         DRIDNModel::ModelOptions::legacyCompatibility());
  const auto context = syntheticContext();
  const auto state = model.initialState(context);
  const DRIDNModel::State expected_state = {
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0}};
  EXPECT_EQ(state, expected_state);

  const DRIDNModel::State expected_rhs = {{-0.22211578966592707,
                                           -0.4699565384010536,
                                           -0.046894112133266795,
                                           126.66945363348695,
                                           268.02413094335566,
                                           23.670146865456612,
                                           127.7713064935416,
                                           270.34080286231983,
                                           26.975668785789523,
                                           0.1395409330898761,
                                           0.293390699437455,
                                           0.4186227992696285,
                                           18.988048974965462,
                                           7.27426576271663,
                                           25.849741840850832,
                                           -0.40264073305012754,
                                           1.8329915516966752e-5,
                                           1.8329915516966752e-5,
                                           1.8329915516966752e-5}};
  const auto derivative = model.rhs(context, state);
  for (std::size_t i = 0; i < DRIDNModel::n_states; ++i)
    expectClose(derivative[i], expected_rhs[i]);
}

TEST(DRIDNModel, CorrectedDefaultIsAtEquilibriumAndDoesNotConsumeActivityFloor)
{
  // This is an analytic corrected-mode invariant, not a validation golden.  Corrected mode needs
  // a new fit before comparison with experimental data; legacy goldens above remain mode-pinned.
  const DRIDNModel corrected;
  const auto context = equilibriumContext();
  const auto state = corrected.initialState(context);
  const auto rate = corrected.rates(context, state);
  const auto derivative = corrected.rhs(context, state);

  EXPECT_DOUBLE_EQ(rate.effective_affinity, 0.0);
  EXPECT_DOUBLE_EQ(rate.charge_transfer_rate_um_y, 0.0);
  EXPECT_DOUBLE_EQ(rate.front_rate_um_y, 0.0);
  for (const auto value : derivative)
    EXPECT_DOUBLE_EQ(value, 0.0);

  const DRIDNModel legacy(DRIDNModel::Parameters(),
                          DRIDNModel::ClosureConstants(),
                          DRIDNModel::ModelOptions::legacyCompatibility());
  EXPECT_GT(legacy.rates(context, state).front_rate_um_y, 0.0);
}

TEST(DRIDNModel, AbsentAlloySpeciesAndZeroColdCaptureStayAbsent)
{
  const DRIDNModel model;
  auto context = equilibriumContext();
  context.affinity_baseline = {{2.0, 20.0, 20.0}};
  context.cold_capture_fraction = {{0.0, 0.0, 0.0}};
  const auto state = model.initialState(context);
  const auto rate = model.rates(context, state);

  EXPECT_DOUBLE_EQ(rate.species_fraction[0], 1.0);
  EXPECT_DOUBLE_EQ(rate.species_fraction[1], 0.0);
  EXPECT_DOUBLE_EQ(rate.species_fraction[2], 0.0);
  for (std::size_t i = 0; i < DRIDNModel::n_elements; ++i)
  {
    EXPECT_DOUBLE_EQ(rate.deposition_rate_y_inv[i], 0.0);
    EXPECT_DOUBLE_EQ(rate.bulk_capture_rate_y_inv[i], 0.0);
  }
}

TEST(DRIDNModel, DepositionLimitMatchesExponentialAndClosesInventory)
{
  const DRIDNModel model;
  auto context = equilibriumContext();
  context.initial_dissolved_ppm = {{10.0, 0.0, 0.0}};
  auto state = model.initialState(context);

  DRIDNModel::IntegrationOptions integration;
  integration.relative_tolerance = 1.0e-9;
  integration.absolute_tolerance = 1.0e-12;
  integration.maximum_step_y = 1.0e-3;
  const Real duration_y = 0.25;
  const auto initial_rate = model.rates(context, state);
  const Real sink_rate = initial_rate.deposition_rate_y_inv[0] +
                         initial_rate.bulk_capture_rate_y_inv[0];
  const Real expected_dissolved = 10.0 * std::exp(-sink_rate * duration_y);
  const Real expected_deposit = (10.0 - expected_dissolved) *
                                initial_rate.deposition_rate_y_inv[0] / sink_rate *
                                initial_rate.deposit_conversion_mg_cm2_per_ppm;

  const auto stats = model.advance(context, state, duration_y, integration);
  EXPECT_GT(stats.accepted_steps, 0u);
  expectClose(state[DRIDNModel::dissolvedIndex(DRIDNModel::Element::Cr)],
              expected_dissolved,
              2.0e-8);
  expectClose(state[DRIDNModel::couponDepositIndex(DRIDNModel::Element::Cr)],
              expected_deposit,
              2.0e-8);
  EXPECT_GE(state[DRIDNModel::dissolvedIndex(DRIDNModel::Element::Cr)], 0.0);
  EXPECT_LT(model.outputs(context, state, duration_y).mass_balance_relative_error, 1.0e-11);
}

TEST(DRIDNModel, LongHighSinkAdvanceRemainsPositiveAndConservative)
{
  const DRIDNModel model;
  auto context = equilibriumContext();
  context.initial_dissolved_ppm = {{10.0, 20.0, 30.0}};
  context.mass_fractions = {{0.17, 0.66, 0.12}};
  context.cold_capture_fraction = {{1.0, 1.0, 1.0}};
  auto state = model.initialState(context);

  DRIDNModel::IntegrationOptions integration;
  integration.relative_tolerance = 1.0e-8;
  integration.absolute_tolerance = 1.0e-12;
  const auto stats = model.advance(context, state, 10.0, integration);
  EXPECT_GT(stats.accepted_steps, 0u);
  for (const auto element : {DRIDNModel::Element::Cr,
                             DRIDNModel::Element::Fe,
                             DRIDNModel::Element::Ni})
  {
    EXPECT_GE(state[DRIDNModel::dissolvedIndex(element)], 0.0);
    EXPECT_GE(state[DRIDNModel::couponDepositIndex(element)], 0.0);
    EXPECT_GE(state[DRIDNModel::bulkCaptureIndex(element)], 0.0);
  }
  EXPECT_LT(model.outputs(context, state, 10.0).mass_balance_relative_error, 1.0e-10);
}

TEST(DRIDNModel, EndpointConvergesBetween120And240Intervals)
{
  const DRIDNModel model;
  const auto context = syntheticContext();
  auto coarse = model.initialState(context);
  auto fine = coarse;

  DRIDNModel::IntegrationOptions coarse_options;
  coarse_options.nominal_maximum_step_intervals = 120;
  DRIDNModel::IntegrationOptions fine_options = coarse_options;
  fine_options.nominal_maximum_step_intervals = 240;
  constexpr Real duration_y = 2.0;
  model.advance(context, coarse, duration_y, coarse_options);
  model.advance(context, fine, duration_y, fine_options);

  Real maximum_scaled_difference = 0.0;
  for (std::size_t i = 0; i < DRIDNModel::n_states; ++i)
    maximum_scaled_difference =
        std::max(maximum_scaled_difference,
                 std::abs(coarse[i] - fine[i]) / std::max(1.0, std::abs(fine[i])));
  EXPECT_LT(maximum_scaled_difference, 1.0e-4);
}

TEST(DRIDNModel, LegacyRefitLongIntervalClosesFloatingPointEndpoint)
{
  const DRIDNModel model(currentRefitLegacyParameters(),
                         DRIDNModel::ClosureConstants(),
                         DRIDNModel::ModelOptions::legacyCompatibility());
  auto context = syntheticContext();
  context.redox_shift_initial = 0.5;
  context.log_charge_base_no_redox = -3.0;
  context.inventory_scale = DRIDNModel::InventoryScale::Explicit;
  context.explicit_inventory_scale = 1.0;

  auto direct = model.initialState(context);
  auto fine_chunks = direct;
  DRIDNModel::IntegrationOptions options;
  options.maximum_step_y = 1.0;
  constexpr Real duration_y = 0.5;
  const auto direct_stats = model.advance(context, direct, duration_y, options);
  for (unsigned int chunk = 0; chunk < 100; ++chunk)
    model.advance(context, fine_chunks, duration_y / 100.0, options);

  EXPECT_GT(direct_stats.accepted_steps, 0u);
  Real maximum_scaled_difference = 0.0;
  for (const auto element : {DRIDNModel::Element::Cr,
                             DRIDNModel::Element::Fe,
                             DRIDNModel::Element::Ni})
  {
    const Real surface = direct[DRIDNModel::surfaceIndex(element)];
    EXPECT_GE(surface, model.closures().surface_availability_floor);
    EXPECT_GE(direct[DRIDNModel::dissolvedIndex(element)], 0.0);
    EXPECT_GE(direct[DRIDNModel::couponDepositIndex(element)], 0.0);
    EXPECT_GE(direct[DRIDNModel::bulkCaptureIndex(element)], 0.0);
  }
  for (std::size_t i = 0; i < DRIDNModel::n_states; ++i)
    maximum_scaled_difference =
        std::max(maximum_scaled_difference,
                 std::abs(direct[i] - fine_chunks[i]) /
                     std::max(1.0, std::abs(fine_chunks[i])));

  EXPECT_LT(maximum_scaled_difference, 1.0e-4);
  EXPECT_LT(model.outputs(context, direct, duration_y).mass_balance_relative_error, 1.0e-10);
}

TEST(DRIDNModel, ProjectedSurfaceFloorSatisfiesComplementarity)
{
  auto parameters = DRIDNModel::Parameters();
  parameters.log_surface_reservoir_um = std::log(0.5);
  parameters.log_surface_replenishment_y_inv = std::log(0.005);
  const DRIDNModel model(parameters,
                         DRIDNModel::ClosureConstants(),
                         DRIDNModel::ModelOptions::legacyCompatibility());
  auto context = equilibriumContext();
  context.affinity_baseline = {{10.0, 0.0, 0.0}};
  context.cold_capture_fraction = {{0.0, 0.0, 0.0}};
  context.log_charge_base_no_redox = 10.0;
  auto state = model.initialState(context);

  model.advance(context, state, 0.1);
  EXPECT_DOUBLE_EQ(state[DRIDNModel::surfaceIndex(DRIDNModel::Element::Cr)],
                   model.closures().surface_availability_floor);
  EXPECT_DOUBLE_EQ(model.rhs(context, state)
                       [DRIDNModel::surfaceIndex(DRIDNModel::Element::Cr)],
                   0.0);
  EXPECT_GE(state[DRIDNModel::dissolvedIndex(DRIDNModel::Element::Cr)], 0.0);
  EXPECT_LT(model.outputs(context, state, 0.1).mass_balance_relative_error, 1.0e-10);
}

TEST(DRIDNModel, RedoxRelaxationAndDamageHaveClosedFormLimits)
{
  auto parameters = DRIDNModel::Parameters();
  parameters.redox_consumption_per_um = 0.0;
  parameters.damage_affinity_scale = 0.0;
  const DRIDNModel model(parameters);

  auto context = equilibriumContext();
  context.affinity_baseline = {{-10.0, -10.0, -10.0}};
  context.redox_shift_initial = 2.0;
  context.transient_redox = true;
  context.cr_diffusion_cm2_s = 1.0e-14;
  context.gb_length_multiplier = 0.8;
  auto state = model.initialState(context);

  DRIDNModel::IntegrationOptions integration;
  integration.relative_tolerance = 1.0e-9;
  integration.absolute_tolerance = 1.0e-12;
  integration.maximum_step_y = 2.0e-3;
  const Real duration_y = 2.0;
  model.advance(context, state, duration_y, integration);

  const Real relaxation = std::exp(parameters.log_redox_relaxation_y_inv);
  const Real target = parameters.redox_buffer_retention * context.redox_shift_initial;
  const Real expected_redox =
      target + (context.redox_shift_initial - target) * std::exp(-relaxation * duration_y);
  const Real diffusion_um2_y = context.cr_diffusion_cm2_s * 1.0e8 *
                               DRIDNModel::seconds_per_year;
  constexpr Real pi = 3.141592653589793238462643383279502884;
  const Real multiplier = context.gb_length_multiplier * parameters.gb_dynamic_scale;
  const Real expected_q = 4.0 / pi * diffusion_um2_y * multiplier * multiplier * duration_y;

  expectClose(state[DRIDNModel::redox_shift], expected_redox, 2.0e-8);
  expectClose(state[DRIDNModel::grain_boundary_depth_squared], expected_q, 2.0e-8);
  EXPECT_DOUBLE_EQ(state[DRIDNModel::front_depth], 0.0);
}

TEST(DRIDNModel, RejectsInvalidParametersContextStateAndDuration)
{
  auto bad_parameters = DRIDNModel::Parameters();
  bad_parameters.log_rate_scale = 2.0;
  EXPECT_THROW(DRIDNModel(bad_parameters), std::invalid_argument);

  const DRIDNModel model;
  auto context = equilibriumContext();
  context.flow_factor = -1.0;
  EXPECT_THROW(model.initialState(context), std::invalid_argument);

  context = equilibriumContext();
  context.flow_factor = 0.0;
  EXPECT_THROW(model.initialState(context), std::invalid_argument);

  context = equilibriumContext();
  context.inventory_coupling_factor = 1.01;
  EXPECT_THROW(model.initialState(context), std::invalid_argument);

  context = equilibriumContext();
  context.mass_fractions = {{0.7, 0.3, 0.1}};
  EXPECT_THROW(model.initialState(context), std::invalid_argument);

  context = equilibriumContext();
  context.mass_loss_fraction = 1.01;
  EXPECT_THROW(model.initialState(context), std::invalid_argument);

  context = equilibriumContext();
  auto state = model.initialState(context);
  state[DRIDNModel::dissolvedIndex(DRIDNModel::Element::Cr)] = -1.0;
  EXPECT_THROW(model.rhs(context, state), std::invalid_argument);

  state = model.initialState(context);
  EXPECT_THROW(model.advance(context, state, -1.0), std::invalid_argument);

  DRIDNModel::ModelOptions bad_options;
  bad_options.charge_transfer = static_cast<DRIDNModel::ChargeTransferMode>(-1);
  EXPECT_THROW(DRIDNModel(DRIDNModel::Parameters(),
                          DRIDNModel::ClosureConstants(),
                          bad_options),
               std::invalid_argument);
}

TEST(DRIDNModel, TimeConversionIsExplicit)
{
  EXPECT_DOUBLE_EQ(DRIDNModel::seconds_per_year, 365.25 * 24.0 * 3600.0);
}

} // namespace Corrosion
