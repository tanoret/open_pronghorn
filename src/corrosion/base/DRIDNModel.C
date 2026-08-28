//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "DRIDNModel.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Corrosion
{
namespace
{

constexpr std::array<DRIDNModel::Element, DRIDNModel::n_elements> elements = {
    DRIDNModel::Element::Cr, DRIDNModel::Element::Fe, DRIDNModel::Element::Ni};

void
requireFinite(const Real value, const std::string & name)
{
  if (!std::isfinite(value))
    throw std::invalid_argument("DRIDN: '" + name + "' must be finite.");
}

void
requireRange(const Real value, const Real lower, const Real upper, const std::string & name)
{
  requireFinite(value, name);
  if (value < lower || value > upper)
  {
    std::ostringstream message;
    message << "DRIDN: '" << name << "'=" << value << " is outside [" << lower << ", "
            << upper << "].";
    throw std::invalid_argument(message.str());
  }
}

void
requirePositive(const Real value, const std::string & name)
{
  requireFinite(value, name);
  if (value <= 0.0)
    throw std::invalid_argument("DRIDN: '" + name + "' must be positive.");
}

void
requireNonnegative(const Real value, const std::string & name)
{
  requireFinite(value, name);
  if (value < 0.0)
    throw std::invalid_argument("DRIDN: '" + name + "' must be nonnegative.");
}

Real
bounded(const Real value, const Real lower, const Real upper)
{
  return std::min(std::max(value, lower), upper);
}

std::size_t
index(const DRIDNModel::Element element)
{
  return static_cast<std::size_t>(element);
}

} // namespace

void
DRIDNModel::Parameters::validate() const
{
  requireRange(log_rate_scale, -1.5, 1.5, "log_rate_scale");
  requireRange(affinity_feedback_scale, 0.0, 0.30, "affinity_feedback_scale");
  requireRange(log_surface_reservoir_um,
               std::log(0.5),
               std::log(500.0),
               "log_surface_reservoir_um");
  requireRange(log_surface_replenishment_y_inv,
               std::log(0.005),
               std::log(200.0),
               "log_surface_replenishment_y_inv");
  requireRange(surface_availability_exponent, 0.0, 3.0, "surface_availability_exponent");
  requireRange(surface_reservoir_cr_exponent,
               0.0,
               2.5,
               "surface_reservoir_cr_exponent");
  requireRange(log_dynamic_cr_exchange_bias, -0.5, 0.5, "log_dynamic_cr_exchange_bias");
  requireRange(log_dynamic_fe_capture_bias, -0.5, 0.5, "log_dynamic_fe_capture_bias");
  requireRange(inventory_inhibition_scale, 0.0, 6.0, "inventory_inhibition_scale");
  requireRange(log_redox_relaxation_y_inv,
               std::log(0.005),
               std::log(200.0),
               "log_redox_relaxation_y_inv");
  requireRange(redox_buffer_retention, 0.0, 1.0, "redox_buffer_retention");
  requireRange(redox_consumption_per_um, 0.0, 0.25, "redox_consumption_per_um");
  requireRange(log_stress_interfacial_factor,
               0.0,
               std::log(12.0),
               "log_stress_interfacial_factor");
  requireRange(log_fluoride_impurity_interfacial_factor,
               0.0,
               std::log(12.0),
               "log_fluoride_impurity_interfacial_factor");
  requireRange(log_deposition_rate_y_inv_fuel,
               std::log(1.0e-3),
               std::log(100.0),
               "log_deposition_rate_y_inv_fuel");
  requireRange(log_deposition_rate_y_inv_flinak,
               std::log(1.0e-3),
               std::log(200.0),
               "log_deposition_rate_y_inv_flinak");
  requireRange(bulk_capture_area_multiplier_fuel,
               0.0,
               12.0,
               "bulk_capture_area_multiplier_fuel");
  requireRange(bulk_capture_area_multiplier_flinak,
               0.0,
               12.0,
               "bulk_capture_area_multiplier_flinak");
  requireRange(log_bulk_precipitation_rate_scale,
               std::log(0.05),
               std::log(50.0),
               "log_bulk_precipitation_rate_scale");
  requireRange(log_inventory_scale_msre, -2.0, 2.0, "log_inventory_scale_msre");
  requireRange(log_inventory_scale_loop, -2.0, 2.0, "log_inventory_scale_loop");
  requireRange(log_deposit_area_scale_fuel,
               -std::log(10.0),
               std::log(10.0),
               "log_deposit_area_scale_fuel");
  requireRange(log_deposit_area_scale_flinak,
               -std::log(10.0),
               std::log(10.0),
               "log_deposit_area_scale_flinak");
  requireRange(log_mass_loss_scale, -1.0, 1.0, "log_mass_loss_scale");
  requireRange(gb_dynamic_scale, 0.2, 3.0, "gb_dynamic_scale");
  requireRange(damage_affinity_scale, 0.0, 0.40, "damage_affinity_scale");
}

void
DRIDNModel::ClosureConstants::validate() const
{
  requireRange(surface_availability_floor, 0.0, 1.0, "surface_availability_floor");
  if (surface_availability_floor >= 1.0)
    throw std::invalid_argument("DRIDN: surface_availability_floor must be less than one.");
  for (std::size_t i = 0; i < n_elements; ++i)
    requirePositive(replenishment_relative[i], "replenishment_relative");
  requireNonnegative(minimum_capture_flow_factor, "minimum_capture_flow_factor");
  requireNonnegative(bulk_precipitation_prefactor_y_inv,
                     "bulk_precipitation_prefactor_y_inv");
  requireNonnegative(maximum_bulk_precipitation_rate_y_inv,
                     "maximum_bulk_precipitation_rate_y_inv");
  requireRange(minimum_mass_loss_fraction, 0.0, 1.0, "minimum_mass_loss_fraction");
  requirePositive(log_rate_clip, "log_rate_clip");
  requirePositive(relative_affinity_clip, "relative_affinity_clip");
}

void
DRIDNModel::Context::validate(const ClosureConstants & closures) const
{
  closures.validate();
  Real modeled_fraction = 0.0;
  for (std::size_t i = 0; i < n_elements; ++i)
  {
    requireNonnegative(mass_fractions[i], "mass_fractions");
    modeled_fraction += mass_fractions[i];
    requireFinite(log_exchange_offsets[i], "log_exchange_offsets");
    requireFinite(affinity_baseline[i], "affinity_baseline");
    requireRange(cold_capture_fraction[i], 0.0, 1.0, "cold_capture_fraction");
    requireNonnegative(initial_dissolved_ppm[i], "initial_dissolved_ppm");
  }
  if (modeled_fraction <= 0.0)
    throw std::invalid_argument("DRIDN: at least one modeled alloy mass fraction must be positive.");
  if (modeled_fraction > 1.0 + 1.0e-12)
    throw std::invalid_argument(
        "DRIDN: modeled alloy mass fractions must not sum to more than one.");

  requirePositive(cr_fraction_ratio, "cr_fraction_ratio");
  requirePositive(density_g_cm3, "density_g_cm3");
  requirePositive(flow_factor, "flow_factor");
  requireNonnegative(selectivity_scale, "selectivity_scale");
  requirePositive(product_activity_floor_ppm, "product_activity_floor_ppm");
  requireFinite(redox_shift_initial, "redox_shift_initial");
  requireFinite(log_charge_base_no_redox, "log_charge_base_no_redox");
  requireNonnegative(mass_transfer_rate_um_y, "mass_transfer_rate_um_y");
  requirePositive(inventory_capacity_ppm, "inventory_capacity_ppm");
  requirePositive(area_to_salt_mass_cm2_g, "area_to_salt_mass_cm2_g");
  requirePositive(explicit_inventory_scale, "explicit_inventory_scale");
  requireRange(inventory_coupling_factor, 0.0, 1.0, "inventory_coupling_factor");
  requirePositive(deposit_area_factor, "deposit_area_factor");
  requireRange(mass_loss_fraction, 0.0, 1.0, "mass_loss_fraction");
  requireNonnegative(cr_diffusion_cm2_s, "cr_diffusion_cm2_s");
  requireNonnegative(front_damage_multiplier, "front_damage_multiplier");
  requireNonnegative(gb_length_multiplier, "gb_length_multiplier");

  switch (inventory_scale)
  {
    case InventoryScale::Explicit:
    case InventoryScale::MSRE:
    case InventoryScale::Loop:
      break;
    default:
      throw std::invalid_argument("DRIDN: inventory_scale is not recognized.");
  }
  switch (deposition_closure)
  {
    case DepositionClosure::FuelLike:
    case DepositionClosure::FLiNaK:
      break;
    default:
      throw std::invalid_argument("DRIDN: deposition_closure is not recognized.");
  }
}

void
DRIDNModel::ModelOptions::validate() const
{
  switch (charge_transfer)
  {
    case ChargeTransferMode::AffinityGatedDissolution:
    case ChargeTransferMode::LegacyIrreversible:
      return;
  }
  throw std::invalid_argument("DRIDN: charge_transfer mode is not recognized.");
}

void
DRIDNModel::IntegrationOptions::validate() const
{
  requirePositive(relative_tolerance, "relative_tolerance");
  requirePositive(absolute_tolerance, "absolute_tolerance");
  if (!std::isfinite(maximum_step_y) &&
      maximum_step_y != std::numeric_limits<Real>::infinity())
    throw std::invalid_argument("DRIDN: maximum_step_y must be positive or infinity.");
  if (std::isfinite(maximum_step_y) && maximum_step_y <= 0.0)
    throw std::invalid_argument("DRIDN: maximum_step_y must be positive.");
  requirePositive(minimum_step_y, "minimum_step_y");
  if (!nominal_maximum_step_intervals)
    throw std::invalid_argument("DRIDN: nominal_maximum_step_intervals must be nonzero.");
  if (!maximum_step_attempts)
    throw std::invalid_argument("DRIDN: maximum_step_attempts must be nonzero.");
}

DRIDNModel::DRIDNModel() : DRIDNModel(Parameters(), ClosureConstants(), ModelOptions()) {}

DRIDNModel::DRIDNModel(const Parameters & parameters)
  : DRIDNModel(parameters, ClosureConstants(), ModelOptions())
{
}

DRIDNModel::DRIDNModel(const Parameters & parameters, const ClosureConstants & closures)
  : DRIDNModel(parameters, closures, ModelOptions())
{
}

DRIDNModel::DRIDNModel(const Parameters & parameters,
                       const ClosureConstants & closures,
                       const ModelOptions & options)
  : _parameters(parameters), _closures(closures), _options(options)
{
  _parameters.validate();
  _closures.validate();
  _options.validate();
}

Real
DRIDNModel::inventoryScale(const Context & context) const
{
  switch (context.inventory_scale)
  {
    case InventoryScale::MSRE:
      return std::exp(_parameters.log_inventory_scale_msre);
    case InventoryScale::Loop:
      return std::exp(_parameters.log_inventory_scale_loop);
    case InventoryScale::Explicit:
      return context.explicit_inventory_scale;
  }
  throw std::logic_error("DRIDN: unrecognized inventory scale.");
}

Real
DRIDNModel::depositAreaScale(const Context & context) const
{
  return std::exp(context.deposition_closure == DepositionClosure::FLiNaK
                      ? _parameters.log_deposit_area_scale_flinak
                      : _parameters.log_deposit_area_scale_fuel);
}

Real
DRIDNModel::depositionBaseRate(const Context & context) const
{
  return std::exp(context.deposition_closure == DepositionClosure::FLiNaK
                      ? _parameters.log_deposition_rate_y_inv_flinak
                      : _parameters.log_deposition_rate_y_inv_fuel);
}

Real
DRIDNModel::bulkCaptureAreaMultiplier(const Context & context) const
{
  return context.deposition_closure == DepositionClosure::FLiNaK
             ? _parameters.bulk_capture_area_multiplier_flinak
             : _parameters.bulk_capture_area_multiplier_fuel;
}

DRIDNModel::State
DRIDNModel::initialState(const Context & context) const
{
  context.validate(_closures);
  State state{};
  for (const auto element : elements)
  {
    const auto i = index(element);
    state[surfaceIndex(element)] = 1.0;
    state[dissolvedIndex(element)] = context.initial_dissolved_ppm[i];
  }
  state[redox_shift] = context.redox_shift_initial;
  return state;
}

DRIDNModel::Rates
DRIDNModel::ratesImpl(const Context & context, const State & state) const
{
  Rates result;
  Triplet surface{};
  Triplet concentration{};
  Triplet activity_concentration{};
  Triplet log_weight{};

  const Real redox_lower = std::min(0.0, context.redox_shift_initial);
  const Real redox_upper = std::max(0.0, context.redox_shift_initial);
  const Real redox = bounded(state[redox_shift], redox_lower, redox_upper);

  for (const auto element : elements)
  {
    const auto i = index(element);
    surface[i] = bounded(state[surfaceIndex(element)],
                         _closures.surface_availability_floor,
                         1.0);
    concentration[i] = std::max(state[dissolvedIndex(element)], 0.0);
    activity_concentration[i] =
        std::max(concentration[i], context.product_activity_floor_ppm);
    const Real product_ratio = activity_concentration[i] / context.product_activity_floor_ppm;
    result.affinity[i] = context.affinity_baseline[i] + redox - std::log(product_ratio);
  }

  for (const auto element : elements)
  {
    const auto i = index(element);
    const Real relative = bounded(result.affinity[i] - result.affinity[index(Element::Cr)],
                                  -_closures.relative_affinity_clip,
                                  _closures.relative_affinity_clip);
    if (context.mass_fractions[i] == 0.0)
      log_weight[i] = -std::numeric_limits<Real>::infinity();
    else
      log_weight[i] = std::log(context.mass_fractions[i]) + context.log_exchange_offsets[i] +
                      context.selectivity_scale * relative +
                      _parameters.surface_availability_exponent *
                          std::log(std::max(surface[i], 1.0e-8));
  }
  log_weight[index(Element::Cr)] += _parameters.log_dynamic_cr_exchange_bias;

  const Real maximum_log_weight = *std::max_element(log_weight.begin(), log_weight.end());
  Real weight_sum = 0.0;
  for (std::size_t i = 0; i < n_elements; ++i)
  {
    result.species_fraction[i] = std::exp(log_weight[i] - maximum_log_weight);
    weight_sum += result.species_fraction[i];
  }
  for (std::size_t i = 0; i < n_elements; ++i)
    result.species_fraction[i] /= weight_sum;

  Real dissolved_total = 0.0;
  for (std::size_t i = 0; i < n_elements; ++i)
  {
    const Real log_product =
        std::log(activity_concentration[i] / context.product_activity_floor_ppm);
    result.product_feedback -= result.species_fraction[i] * log_product;
    result.surface_feedback += result.species_fraction[i] * std::log(surface[i]);
    result.effective_affinity += result.species_fraction[i] * result.affinity[i];
    dissolved_total += concentration[i];
  }

  Real interfacial_log_factor = 0.0;
  if (context.stress_interfacial_activation)
    interfacial_log_factor += _parameters.log_stress_interfacial_factor;
  if (context.fluoride_impurity_interfacial_activation && !context.chloride_salt)
    interfacial_log_factor += _parameters.log_fluoride_impurity_interfacial_factor;

  const Real log_charge = bounded(context.log_charge_base_no_redox + redox +
                                      _parameters.log_rate_scale + interfacial_log_factor +
                                      _parameters.affinity_feedback_scale *
                                          result.product_feedback +
                                      _parameters.surface_availability_exponent *
                                          result.surface_feedback,
                                  -_closures.log_rate_clip,
                                  _closures.log_rate_clip);
  const Real irreversible_charge_rate = std::exp(log_charge);
  if (_options.charge_transfer == ChargeTransferMode::LegacyIrreversible)
    result.charge_transfer_rate_um_y = irreversible_charge_rate;
  else if (result.effective_affinity > 0.0)
  {
    // A stable dissolution-only form of 1-exp(-A).  It approaches the calibrated irreversible
    // rate far from equilibrium while making the net branch exactly zero at A=0.
    const Real affinity = std::min(result.effective_affinity, _closures.log_rate_clip);
    result.charge_transfer_rate_um_y = irreversible_charge_rate * (-std::expm1(-affinity));
  }

  const Real inventory_ratio = dissolved_total / context.inventory_capacity_ppm;
  result.mass_transfer_rate_um_y =
      context.mass_transfer_rate_um_y * std::exp(interfacial_log_factor) /
      (1.0 + _parameters.inventory_inhibition_scale * inventory_ratio);
  if (result.charge_transfer_rate_um_y > 0.0 && result.mass_transfer_rate_um_y > 0.0)
    result.front_rate_um_y =
        1.0 / (1.0 / result.charge_transfer_rate_um_y +
               1.0 / result.mass_transfer_rate_um_y);

  const Real area_to_mass = context.area_to_salt_mass_cm2_g * inventoryScale(context);
  const Real source_conversion = context.density_g_cm3 * 100.0 * area_to_mass *
                                 context.inventory_coupling_factor;
  for (std::size_t i = 0; i < n_elements; ++i)
    result.source_rate_ppm_y[i] =
        result.front_rate_um_y * source_conversion * result.species_fraction[i];

  const Real deposition_base = depositionBaseRate(context);
  const Real capture_flow =
      std::sqrt(std::max(context.flow_factor, _closures.minimum_capture_flow_factor));
  for (std::size_t i = 0; i < n_elements; ++i)
    result.deposition_rate_y_inv[i] =
        deposition_base * context.cold_capture_fraction[i] * capture_flow;
  result.deposition_rate_y_inv[index(Element::Fe)] *=
      std::exp(_parameters.log_dynamic_fe_capture_bias);
  for (std::size_t i = 0; i < n_elements; ++i)
    result.bulk_capture_rate_y_inv[i] =
        result.deposition_rate_y_inv[i] * bulkCaptureAreaMultiplier(context);

  result.bulk_precipitation_rate_y_inv =
      std::min(_closures.maximum_bulk_precipitation_rate_y_inv,
               _closures.bulk_precipitation_prefactor_y_inv *
                   std::exp(_parameters.log_bulk_precipitation_rate_scale) *
                   inventory_ratio * inventory_ratio);
  result.deposit_conversion_mg_cm2_per_ppm =
      1.0e-3 / area_to_mass * context.deposit_area_factor * depositAreaScale(context);
  result.loss_fraction =
      bounded(context.mass_loss_fraction * std::exp(_parameters.log_mass_loss_scale),
              _closures.minimum_mass_loss_fraction,
              1.0);
  return result;
}

DRIDNModel::Rates
DRIDNModel::rates(const Context & context, const State & state) const
{
  context.validate(_closures);
  if (!admissible(context, state, 0.0))
    throw std::invalid_argument("DRIDN: rates() received a nonphysical state.");
  return ratesImpl(context, state);
}

DRIDNModel::State
DRIDNModel::rhsImpl(const Context & context, const State & state) const
{
  const Rates rate = ratesImpl(context, state);
  State derivative{};

  const Real reservoir_um =
      std::exp(_parameters.log_surface_reservoir_um) *
      std::pow(context.cr_fraction_ratio, _parameters.surface_reservoir_cr_exponent);
  const Real replenishment_y_inv = std::exp(_parameters.log_surface_replenishment_y_inv);
  for (const auto element : elements)
  {
    const auto i = index(element);
    const auto s_index = surfaceIndex(element);
    const Real surface = bounded(state[s_index], _closures.surface_availability_floor, 1.0);
    derivative[s_index] = replenishment_y_inv * _closures.replenishment_relative[i] *
                              (1.0 - surface) -
                          rate.front_rate_um_y * rate.species_fraction[i] / reservoir_um;
    if ((state[s_index] <= _closures.surface_availability_floor &&
         derivative[s_index] < 0.0) ||
        (state[s_index] >= 1.0 && derivative[s_index] > 0.0))
      derivative[s_index] = 0.0;

    const Real concentration = std::max(state[dissolvedIndex(element)], 0.0);
    derivative[dissolvedIndex(element)] =
        rate.source_rate_ppm_y[i] -
        (rate.deposition_rate_y_inv[i] + rate.bulk_capture_rate_y_inv[i] +
         rate.bulk_precipitation_rate_y_inv) *
            concentration;
    derivative[cumulativeSourceIndex(element)] = rate.source_rate_ppm_y[i];
    derivative[couponDepositIndex(element)] =
        rate.deposition_rate_y_inv[i] * concentration *
        rate.deposit_conversion_mg_cm2_per_ppm;
    derivative[bulkCaptureIndex(element)] =
        (rate.bulk_capture_rate_y_inv[i] + rate.bulk_precipitation_rate_y_inv) *
        concentration;
  }

  derivative[front_depth] = rate.front_rate_um_y;
  derivative[mass_recession] = rate.front_rate_um_y * rate.loss_fraction;
  constexpr Real pi = 3.141592653589793238462643383279502884;
  const Real diffusion_um2_y = context.cr_diffusion_cm2_s * 1.0e8 * seconds_per_year;
  const Real gb_multiplier = context.gb_length_multiplier * _parameters.gb_dynamic_scale;
  const Real redox_lower = std::min(0.0, context.redox_shift_initial);
  const Real redox_upper = std::max(0.0, context.redox_shift_initial);
  const Real redox = bounded(state[redox_shift], redox_lower, redox_upper);
  derivative[grain_boundary_depth_squared] =
      4.0 / pi * diffusion_um2_y * gb_multiplier * gb_multiplier *
      std::exp(_parameters.damage_affinity_scale * std::max(redox, 0.0));

  if (context.transient_redox)
  {
    const Real target = _parameters.redox_buffer_retention * context.redox_shift_initial;
    Real consumption = 0.0;
    if (std::abs(redox) > 1.0e-12)
      consumption = -std::copysign(_parameters.redox_consumption_per_um *
                                       rate.front_rate_um_y,
                                   redox);
    derivative[redox_shift] = std::exp(_parameters.log_redox_relaxation_y_inv) *
                                  (target - redox) +
                              consumption;
    if ((state[redox_shift] <= redox_lower && derivative[redox_shift] < 0.0) ||
        (state[redox_shift] >= redox_upper && derivative[redox_shift] > 0.0))
      derivative[redox_shift] = 0.0;
  }
  return derivative;
}

DRIDNModel::State
DRIDNModel::rhs(const Context & context, const State & state) const
{
  context.validate(_closures);
  if (!admissible(context, state, 0.0))
    throw std::invalid_argument("DRIDN: rhs() received a nonphysical state.");
  return rhsImpl(context, state);
}

bool
DRIDNModel::admissible(const Context & context, const State & state, const Real tolerance) const
{
  for (const auto value : state)
    if (!std::isfinite(value))
      return false;

  for (const auto element : elements)
  {
    if (state[surfaceIndex(element)] < _closures.surface_availability_floor - tolerance ||
        state[surfaceIndex(element)] > 1.0 + tolerance)
      return false;
    if (state[dissolvedIndex(element)] < -tolerance ||
        state[cumulativeSourceIndex(element)] < -tolerance ||
        state[couponDepositIndex(element)] < -tolerance ||
        state[bulkCaptureIndex(element)] < -tolerance)
      return false;
  }
  if (state[front_depth] < -tolerance || state[mass_recession] < -tolerance ||
      state[grain_boundary_depth_squared] < -tolerance)
    return false;
  const Real redox_lower = std::min(0.0, context.redox_shift_initial);
  const Real redox_upper = std::max(0.0, context.redox_shift_initial);
  return state[redox_shift] >= redox_lower - tolerance &&
         state[redox_shift] <= redox_upper + tolerance;
}

bool
DRIDNModel::rk4Step(const Context & context,
                    const State & state,
                    const Real step_y,
                    const Real admissibility_tolerance,
                    State & result) const
{
  // Surface availability and redox are box-constrained projected states.  Without projecting
  // RK stages, a trajectory whose unconstrained derivative points through a bound approaches the
  // bound through an infinite sequence of rejected steps (a numerical Zeno failure).  Projecting
  // only these nonconserved variables is the standard method-of-lines realization of the
  // complementarity conditions already enforced by rhsImpl().  Dissolved and accumulated
  // inventories are deliberately not projected here: a negative inventory still rejects the
  // step, and restoreConservation() closes every accepted element balance.
  const auto projectBoundaryVariables = [this, &context](State & candidate) {
    for (const auto element : elements)
      candidate[surfaceIndex(element)] =
          bounded(candidate[surfaceIndex(element)], _closures.surface_availability_floor, 1.0);
    candidate[redox_shift] = bounded(candidate[redox_shift],
                                     std::min(0.0, context.redox_shift_initial),
                                     std::max(0.0, context.redox_shift_initial));
  };

  const State k1 = rhsImpl(context, state);
  State stage{};
  for (std::size_t i = 0; i < n_states; ++i)
    stage[i] = state[i] + 0.5 * step_y * k1[i];
  projectBoundaryVariables(stage);
  if (!admissible(context, stage, admissibility_tolerance))
    return false;

  const State k2 = rhsImpl(context, stage);
  for (std::size_t i = 0; i < n_states; ++i)
    stage[i] = state[i] + 0.5 * step_y * k2[i];
  projectBoundaryVariables(stage);
  if (!admissible(context, stage, admissibility_tolerance))
    return false;

  const State k3 = rhsImpl(context, stage);
  for (std::size_t i = 0; i < n_states; ++i)
    stage[i] = state[i] + step_y * k3[i];
  projectBoundaryVariables(stage);
  if (!admissible(context, stage, admissibility_tolerance))
    return false;

  const State k4 = rhsImpl(context, stage);
  for (std::size_t i = 0; i < n_states; ++i)
    result[i] = state[i] + step_y / 6.0 *
                               (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  projectBoundaryVariables(result);
  return admissible(context, result, admissibility_tolerance);
}

void
DRIDNModel::restoreConservation(const Context & context,
                                const Triplet & inventory_invariant,
                                State & state,
                                const Real tolerance) const
{
  const Real conversion = ratesImpl(context, state).deposit_conversion_mg_cm2_per_ppm;
  for (const auto element : elements)
  {
    const auto i = index(element);
    state[surfaceIndex(element)] =
        bounded(state[surfaceIndex(element)], _closures.surface_availability_floor, 1.0);
    state[dissolvedIndex(element)] = std::max(state[dissolvedIndex(element)], 0.0);
    state[cumulativeSourceIndex(element)] =
        std::max(state[cumulativeSourceIndex(element)], 0.0);
    state[couponDepositIndex(element)] =
        std::max(state[couponDepositIndex(element)], 0.0);
    state[bulkCaptureIndex(element)] =
        inventory_invariant[i] + state[cumulativeSourceIndex(element)] -
        state[dissolvedIndex(element)] - state[couponDepositIndex(element)] / conversion;
    if (state[bulkCaptureIndex(element)] < 0.0 &&
        state[bulkCaptureIndex(element)] >= -tolerance)
      state[bulkCaptureIndex(element)] = 0.0;
  }
  state[front_depth] = std::max(state[front_depth], 0.0);
  state[mass_recession] = std::max(state[mass_recession], 0.0);
  state[grain_boundary_depth_squared] =
      std::max(state[grain_boundary_depth_squared], 0.0);
  state[redox_shift] = bounded(state[redox_shift],
                               std::min(0.0, context.redox_shift_initial),
                               std::max(0.0, context.redox_shift_initial));
}

DRIDNModel::IntegrationStats
DRIDNModel::advance(const Context & context, State & state, const Real duration_y) const
{
  return advance(context, state, duration_y, IntegrationOptions());
}

DRIDNModel::IntegrationStats
DRIDNModel::advance(const Context & context,
                    State & state,
                    const Real duration_y,
                    const IntegrationOptions & options) const
{
  context.validate(_closures);
  options.validate();
  requireNonnegative(duration_y, "duration_y");
  if (!admissible(context, state, 0.0))
    throw std::invalid_argument("DRIDN: advance() received a nonphysical initial state.");

  IntegrationStats stats;
  if (duration_y == 0.0)
    return stats;

  const Real conversion = ratesImpl(context, state).deposit_conversion_mg_cm2_per_ppm;
  Triplet inventory_invariant{};
  for (const auto element : elements)
  {
    const auto i = index(element);
    inventory_invariant[i] = state[dissolvedIndex(element)] +
                             state[couponDepositIndex(element)] / conversion +
                             state[bulkCaptureIndex(element)] -
                             state[cumulativeSourceIndex(element)];
  }

  const Real nominal_maximum_step =
      duration_y / static_cast<Real>(options.nominal_maximum_step_intervals);
  const Real maximum_step = std::min(options.maximum_step_y, nominal_maximum_step);
  Real step = std::min(duration_y, maximum_step);
  Real elapsed = 0.0;
  unsigned int attempts = 0;
  while (elapsed < duration_y)
  {
    if (++attempts > options.maximum_step_attempts)
      throw std::runtime_error("DRIDN: adaptive integrator exceeded maximum_step_attempts.");
    const Real remaining = duration_y - elapsed;
    step = std::min(step, remaining);
    const bool completes_interval = step == remaining;
    // A final interval can legitimately be shorter than minimum_step_y because it is the
    // floating-point remainder after accepted substeps (or because the caller requested a very
    // small duration).  Permit one attempt that lands exactly on the endpoint.  A rejected step
    // below the minimum is still a hard solver failure on the following iteration.
    if (step < options.minimum_step_y && !completes_interval)
    {
      std::ostringstream message;
      message << "DRIDN: adaptive integrator reached minimum_step_y=" << options.minimum_step_y
              << " at elapsed_time_y=" << elapsed << " of duration_y=" << duration_y << ".";
      throw std::runtime_error(message.str());
    }

    const Real admissibility_tolerance = 10.0 * options.absolute_tolerance;
    State full{};
    State first_half{};
    State two_half{};
    const bool valid = rk4Step(context, state, step, admissibility_tolerance, full) &&
                       rk4Step(context,
                               state,
                               0.5 * step,
                               admissibility_tolerance,
                               first_half) &&
                       rk4Step(context,
                               first_half,
                               0.5 * step,
                               admissibility_tolerance,
                               two_half);

    Real normalized_error = std::numeric_limits<Real>::infinity();
    if (valid)
    {
      normalized_error = 0.0;
      for (std::size_t i = 0; i < n_states; ++i)
      {
        const Real scale = options.absolute_tolerance +
                           options.relative_tolerance *
                               std::max(std::abs(state[i]), std::abs(two_half[i]));
        // Two half RK4 steps and one full step differ by 15 times the leading local error.
        normalized_error =
            std::max(normalized_error, std::abs(two_half[i] - full[i]) / (15.0 * scale));
      }
    }

    if (valid && normalized_error <= 1.0)
    {
      State candidate = two_half;
      restoreConservation(context,
                          inventory_invariant,
                          candidate,
                          10.0 * options.absolute_tolerance);
      if (admissible(context, candidate, 10.0 * options.absolute_tolerance))
      {
        state = candidate;
        // Assign the requested endpoint exactly when this was the final substep.  Repeated
        // floating-point addition can otherwise leave a sub-ulp remainder, which would be
        // misdiagnosed as a minimum-step failure on the next loop iteration.
        elapsed = completes_interval ? duration_y : elapsed + step;
        stats.accepted_steps++;
        stats.final_step_y = step;
        const Real factor = normalized_error > 0.0
                                ? bounded(0.9 * std::pow(1.0 / normalized_error, 0.2), 0.2, 5.0)
                                : 2.0;
        step = std::min(maximum_step, step * factor);
        continue;
      }
    }

    stats.rejected_steps++;
    const Real factor = std::isfinite(normalized_error) && normalized_error > 0.0
                            ? bounded(0.8 * std::pow(1.0 / normalized_error, 0.2), 0.1, 0.5)
                            : 0.5;
    step *= factor;
  }
  return stats;
}

DRIDNModel::Outputs
DRIDNModel::outputs(const Context & context,
                    const State & state,
                    const Real elapsed_time_y) const
{
  context.validate(_closures);
  requireNonnegative(elapsed_time_y, "elapsed_time_y");
  if (!admissible(context, state, 1.0e-10))
    throw std::invalid_argument("DRIDN: outputs() received a nonphysical state.");

  Outputs result;
  const Rates rate = ratesImpl(context, state);
  result.front_depth_um = state[front_depth];
  result.mass_recession_um = state[mass_recession];
  result.mass_loss_mg_cm2 = state[mass_recession] * context.density_g_cm3 * 0.1;
  result.igc_depth_um = context.front_damage_multiplier * state[front_depth] +
                        std::sqrt(std::max(state[grain_boundary_depth_squared], 0.0));
  result.average_corrosion_rate_um_y =
      elapsed_time_y > 0.0 ? state[mass_recession] / elapsed_time_y : 0.0;
  result.instantaneous_front_rate_um_y = rate.front_rate_um_y;
  result.redox_shift = state[redox_shift];

  for (const auto element : elements)
  {
    const auto i = index(element);
    result.surface_availability[i] = state[surfaceIndex(element)];
    result.dissolved_ppm[i] = state[dissolvedIndex(element)];
    result.cumulative_source_ppm[i] = state[cumulativeSourceIndex(element)];
    result.coupon_deposit_mg_cm2[i] = state[couponDepositIndex(element)];
    result.bulk_captured_ppm[i] = state[bulkCaptureIndex(element)];
    result.mass_gain_mg_cm2 += result.coupon_deposit_mg_cm2[i];

    const Real expected = context.initial_dissolved_ppm[i] + result.cumulative_source_ppm[i];
    const Real accounted = result.dissolved_ppm[i] +
                           result.coupon_deposit_mg_cm2[i] /
                               rate.deposit_conversion_mg_cm2_per_ppm +
                           result.bulk_captured_ppm[i];
    result.mass_balance_relative_error =
        std::max(result.mass_balance_relative_error,
                 std::abs(expected - accounted) / std::max(std::abs(expected), 1.0e-12));
  }
  return result;
}

} // namespace Corrosion
