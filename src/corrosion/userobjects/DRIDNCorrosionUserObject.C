//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "DRIDNCorrosionUserObject.h"

#include "MooseError.h"

#include <array>
#include <cmath>
#include <string>
#include <vector>

registerMooseObject("OpenPronghornApp", DRIDNCorrosionUserObject);

namespace
{
Corrosion::DRIDNModel::Triplet
triplet(const InputParameters & parameters, const std::string & name)
{
  const auto & values = parameters.get<std::vector<Real>>(name);
  if (values.size() != Corrosion::DRIDNModel::n_elements)
    mooseError(
        "DRIDNCorrosionUserObject: '", name, "' must contain exactly Cr, Fe, and Ni values.");
  return {{values[0], values[1], values[2]}};
}

Corrosion::DRIDNModel::Parameters
modelParameters(const Corrosion::AdvancedCorrosionModelDatabase & database)
{
  Corrosion::DRIDNModel::Parameters value;
  value.log_rate_scale = database.dynamicParameter("log_rate_scale");
  value.affinity_feedback_scale = database.dynamicParameter("affinity_feedback_scale");
  value.log_surface_reservoir_um = database.dynamicParameter("log_surface_reservoir_um");
  value.log_surface_replenishment_y_inv =
      database.dynamicParameter("log_surface_replenishment_y_inv");
  value.surface_availability_exponent =
      database.dynamicParameter("surface_availability_exponent");
  value.surface_reservoir_cr_exponent =
      database.dynamicParameter("surface_reservoir_cr_exponent");
  value.log_dynamic_cr_exchange_bias =
      database.dynamicParameter("log_dynamic_cr_exchange_bias");
  value.log_dynamic_fe_capture_bias =
      database.dynamicParameter("log_dynamic_fe_capture_bias");
  value.inventory_inhibition_scale = database.dynamicParameter("inventory_inhibition_scale");
  value.log_redox_relaxation_y_inv = database.dynamicParameter("log_redox_relaxation_y_inv");
  value.redox_buffer_retention = database.dynamicParameter("redox_buffer_retention");
  value.redox_consumption_per_um = database.dynamicParameter("redox_consumption_per_um");
  value.log_stress_interfacial_factor =
      database.dynamicParameter("log_stress_interfacial_factor");
  value.log_fluoride_impurity_interfacial_factor =
      database.dynamicParameter("log_fluoride_impurity_interfacial_factor");
  value.log_deposition_rate_y_inv_fuel =
      database.dynamicParameter("log_deposition_rate_y_inv_fuel");
  value.log_deposition_rate_y_inv_flinak =
      database.dynamicParameter("log_deposition_rate_y_inv_flinak");
  value.bulk_capture_area_multiplier_fuel =
      database.dynamicParameter("bulk_capture_area_multiplier_fuel");
  value.bulk_capture_area_multiplier_flinak =
      database.dynamicParameter("bulk_capture_area_multiplier_flinak");
  value.log_bulk_precipitation_rate_scale =
      database.dynamicParameter("log_bulk_precipitation_rate_scale");
  value.log_inventory_scale_msre = database.dynamicParameter("log_inventory_scale_msre");
  value.log_inventory_scale_loop = database.dynamicParameter("log_inventory_scale_loop");
  value.log_deposit_area_scale_fuel = database.dynamicParameter("log_deposit_area_scale_fuel");
  value.log_deposit_area_scale_flinak =
      database.dynamicParameter("log_deposit_area_scale_flinak");
  value.log_mass_loss_scale = database.dynamicParameter("log_mass_loss_scale");
  value.gb_dynamic_scale = database.dynamicParameter("gb_dynamic_scale");
  value.damage_affinity_scale = database.dynamicParameter("damage_affinity_scale");
  return value;
}

Corrosion::DRIDNModel::Context
modelContext(const InputParameters & parameters)
{
  Corrosion::DRIDNModel::Context value;
  value.mass_fractions = triplet(parameters, "mass_fractions");
  value.log_exchange_offsets = triplet(parameters, "log_exchange_offsets");
  value.affinity_baseline = triplet(parameters, "affinity_baseline");
  value.cold_capture_fraction = triplet(parameters, "cold_capture_fraction");
  value.initial_dissolved_ppm = triplet(parameters, "initial_dissolved_ppm");

  value.cr_fraction_ratio = parameters.get<Real>("cr_fraction_ratio");
  value.density_g_cm3 = parameters.get<Real>("density");
  value.flow_factor = parameters.get<Real>("flow_factor");
  value.selectivity_scale = parameters.get<Real>("selectivity_scale");
  value.product_activity_floor_ppm = parameters.get<Real>("product_activity_floor");
  value.redox_shift_initial = parameters.get<Real>("initial_redox_shift");
  value.log_charge_base_no_redox = parameters.get<Real>("log_charge_base_no_redox");
  value.mass_transfer_rate_um_y = parameters.get<Real>("mass_transfer_rate");
  value.inventory_capacity_ppm = parameters.get<Real>("inventory_capacity");
  value.area_to_salt_mass_cm2_g = parameters.get<Real>("area_to_salt_mass");
  value.explicit_inventory_scale = parameters.get<Real>("inventory_scale_factor");
  value.inventory_coupling_factor = parameters.get<Real>("inventory_coupling_factor");
  value.deposit_area_factor = parameters.get<Real>("deposit_area_factor");
  value.mass_loss_fraction = parameters.get<Real>("mass_loss_fraction");
  value.cr_diffusion_cm2_s = parameters.get<Real>("cr_diffusion");
  value.front_damage_multiplier = parameters.get<Real>("front_damage_multiplier");
  value.gb_length_multiplier = parameters.get<Real>("grain_boundary_length_multiplier");

  const MooseEnum & inventory = parameters.get<MooseEnum>("inventory_scale");
  value.inventory_scale = inventory == "msre"
                              ? Corrosion::DRIDNModel::InventoryScale::MSRE
                              : (inventory == "loop"
                                     ? Corrosion::DRIDNModel::InventoryScale::Loop
                                     : Corrosion::DRIDNModel::InventoryScale::Explicit);
  value.deposition_closure = parameters.get<MooseEnum>("deposition_closure") == "flinak"
                                 ? Corrosion::DRIDNModel::DepositionClosure::FLiNaK
                                 : Corrosion::DRIDNModel::DepositionClosure::FuelLike;
  value.transient_redox = parameters.get<bool>("transient_redox");
  value.stress_interfacial_activation =
      parameters.get<bool>("stress_interfacial_activation");
  value.fluoride_impurity_interfacial_activation =
      parameters.get<bool>("fluoride_impurity_interfacial_activation");
  value.chloride_salt = parameters.get<bool>("chloride_salt");
  return value;
}

Corrosion::DRIDNModel::IntegrationOptions
integrationOptions(const InputParameters & parameters)
{
  Corrosion::DRIDNModel::IntegrationOptions value;
  value.relative_tolerance = parameters.get<Real>("relative_tolerance");
  value.absolute_tolerance = parameters.get<Real>("absolute_tolerance");
  value.maximum_step_y = parameters.get<Real>("maximum_internal_step");
  value.minimum_step_y = parameters.get<Real>("minimum_internal_step");
  value.nominal_maximum_step_intervals =
      parameters.get<unsigned int>("nominal_maximum_step_intervals");
  value.maximum_step_attempts = parameters.get<unsigned int>("maximum_step_attempts");
  return value;
}

Corrosion::DRIDNModel::ModelOptions
modelOptions(const InputParameters & parameters)
{
  if (parameters.get<MooseEnum>("charge_transfer_mode") == "legacy_irreversible")
    return Corrosion::DRIDNModel::ModelOptions::legacyCompatibility();
  return Corrosion::DRIDNModel::ModelOptions();
}
} // namespace

InputParameters
DRIDNCorrosionUserObject::validParams()
{
  InputParameters params = AdvancedCorrosionModelUserObject::validParams();
  params.addClassDescription(
      "Advances the restartable 19-state Dynamic Redox-Inventory-Depletion Network at timestep "
      "end using an explicit, metadata-independent physical context.");
  params.set<ExecFlagEnum>("execute_on") = "TIMESTEP_END";

  params.addParam<DataFileName>("advanced_database",
                                "advanced_corrosion_models.json",
                                "Calibrated DRIDN parameters and model provenance.");
  params.addRequiredParam<std::vector<Real>>(
      "mass_fractions", "Cr, Fe, and Ni alloy mass fractions.");
  params.addRequiredParam<std::vector<Real>>(
      "log_exchange_offsets", "Cr, Fe, and Ni log exchange-current offsets.");
  params.addRequiredParam<std::vector<Real>>(
      "affinity_baseline", "Cr, Fe, and Ni initial dimensionless reaction affinities.");
  params.addRequiredParam<std::vector<Real>>(
      "cold_capture_fraction", "Cr, Fe, and Ni cold-side capture fractions [0,1].");
  params.addRequiredParam<std::vector<Real>>(
      "initial_dissolved_ppm", "Initial physical Cr, Fe, and Ni dissolved inventories [ppm].");

  params.addRequiredRangeCheckedParam<Real>(
      "cr_fraction_ratio", "cr_fraction_ratio > 0", "Alloy Cr fraction divided by 0.07.");
  params.addRequiredRangeCheckedParam<Real>(
      "density", "density > 0", "Alloy density [g/cm^3].");
  params.addRequiredRangeCheckedParam<Real>(
      "flow_factor", "flow_factor > 0", "Dimensionless circulation factor.");
  params.addRequiredRangeCheckedParam<Real>(
      "selectivity_scale", "selectivity_scale >= 0", "Affinity selectivity coefficient.");
  params.addRequiredRangeCheckedParam<Real>(
      "product_activity_floor",
      "product_activity_floor > 0",
      "Nernst activity floor [ppm], separate from the physical inventory.");
  params.addRequiredParam<Real>("initial_redox_shift", "Initial dimensionless redox shift.");
  params.addRequiredParam<Real>(
      "log_charge_base_no_redox", "Log charge-transfer rate before the redox term.");
  params.addRequiredRangeCheckedParam<Real>(
      "mass_transfer_rate", "mass_transfer_rate >= 0", "Mass-transfer limiting rate [um/y].");
  params.addRequiredRangeCheckedParam<Real>(
      "inventory_capacity", "inventory_capacity > 0", "Salt inventory capacity [ppm].");
  params.addRequiredRangeCheckedParam<Real>(
      "area_to_salt_mass", "area_to_salt_mass > 0", "Wetted-area/salt-mass ratio [cm^2/g].");
  params.addRequiredRangeCheckedParam<Real>(
      "inventory_scale_factor",
      "inventory_scale_factor > 0",
      "Explicit inventory multiplier, used only when inventory_scale=explicit.");
  params.addRequiredRangeCheckedParam<Real>(
      "inventory_coupling_factor",
      "inventory_coupling_factor >= 0 & inventory_coupling_factor <= 1",
      "Fractional source-to-salt inventory coupling.");
  params.addRequiredRangeCheckedParam<Real>(
      "deposit_area_factor", "deposit_area_factor > 0", "Coupon deposit-area conversion factor.");
  params.addRequiredRangeCheckedParam<Real>(
      "mass_loss_fraction",
      "mass_loss_fraction >= 0 & mass_loss_fraction <= 1",
      "Front recession contributing to mass loss.");
  params.addRequiredRangeCheckedParam<Real>(
      "cr_diffusion", "cr_diffusion >= 0", "Chromium solid diffusivity [cm^2/s].");
  params.addRequiredRangeCheckedParam<Real>(
      "front_damage_multiplier",
      "front_damage_multiplier >= 0",
      "Front-depth contribution to IGC damage.");
  params.addRequiredRangeCheckedParam<Real>(
      "grain_boundary_length_multiplier",
      "grain_boundary_length_multiplier >= 0",
      "Grain-boundary diffusion-length multiplier.");

  MooseEnum inventory_scale("explicit msre loop");
  params.addRequiredParam<MooseEnum>(
      "inventory_scale", inventory_scale, "Explicit selection of the inventory closure.");
  MooseEnum deposition_closure("fuel_like flinak");
  params.addRequiredParam<MooseEnum>(
      "deposition_closure", deposition_closure, "Explicit deposition/capture closure.");
  params.addRequiredParam<bool>("transient_redox", "Enable the DRIDN redox evolution equation.");
  params.addRequiredParam<bool>(
      "stress_interfacial_activation", "Apply the calibrated stress interfacial factor.");
  params.addRequiredParam<bool>("fluoride_impurity_interfacial_activation",
                                "Apply the calibrated fluoride-impurity interfacial factor.");
  params.addRequiredParam<bool>("chloride_salt", "Whether this context represents chloride salt.");

  MooseEnum transfer_mode("legacy_irreversible affinity_gated_dissolution", "legacy_irreversible");
  params.addParam<MooseEnum>(
      "charge_transfer_mode",
      transfer_mode,
      "Charge-transfer law. legacy_irreversible reproduces the calibrated research model; "
      "affinity_gated_dissolution is experimental and is not covered by the stored calibration "
      "metrics.");
  params.addRangeCheckedParam<Real>(
      "relative_tolerance", 2.0e-7, "relative_tolerance > 0", "Adaptive RK tolerance.");
  params.addRangeCheckedParam<Real>(
      "absolute_tolerance", 1.0e-10, "absolute_tolerance > 0", "Adaptive RK tolerance.");
  params.addRangeCheckedParam<Real>("maximum_internal_step",
                                    1.0,
                                    "maximum_internal_step > 0",
                                    "Maximum internal integrator step [y].");
  params.addRangeCheckedParam<Real>("minimum_internal_step",
                                    1.0e-14,
                                    "minimum_internal_step > 0",
                                    "Minimum internal integrator step [y].");
  params.addRangeCheckedParam<unsigned int>("nominal_maximum_step_intervals",
                                             120,
                                             "nominal_maximum_step_intervals > 0",
                                             "Minimum nominal subdivisions of every MOOSE step.");
  params.addRangeCheckedParam<unsigned int>("maximum_step_attempts",
                                            200000,
                                            "maximum_step_attempts > 0",
                                            "Maximum adaptive integration attempts per MOOSE "
                                            "step.");
  return params;
}

DRIDNCorrosionUserObject::DRIDNCorrosionUserObject(const InputParameters & parameters)
  : AdvancedCorrosionModelUserObject(parameters),
    _advanced_database(getParam<DataFileName>("advanced_database")),
    _parameters(modelParameters(_advanced_database)),
    _context(modelContext(parameters)),
    _integration_options(integrationOptions(parameters)),
    _model_options(modelOptions(parameters)),
    _model(_parameters, Corrosion::DRIDNModel::ClosureConstants(), _model_options),
    _state(declareRestartableData<Corrosion::DRIDNModel::State>("state",
                                                                _model.initialState(_context))),
    _elapsed_time_y(declareRestartableData<Real>("elapsed_time_y", 0.0)),
    _accepted_steps(declareRestartableData<unsigned int>("accepted_steps", 0)),
    _rejected_steps(declareRestartableData<unsigned int>("rejected_steps", 0))
{
  const auto & execute_on = getParam<ExecFlagEnum>("execute_on");
  if (execute_on.size() != 1 || !execute_on.contains("TIMESTEP_END"))
    paramError("execute_on",
               "DRIDNCorrosionUserObject must execute exactly once at TIMESTEP_END. Other flags "
               "could advance the full MOOSE dt zero times, more than once, or at the wrong "
               "lifecycle stage.");

  _context.validate(_model.closures());
  _integration_options.validate();
}

void
DRIDNCorrosionUserObject::execute()
{
  // The constructor fixes the public schedule. Keep this defensive lifecycle guard so framework
  // initialization can never consume a timestep.
  if (_t_step == 0)
    return;
  if (!std::isfinite(_dt) || _dt < 0.0)
    mooseError("DRIDNCorrosionUserObject '", name(), "' received invalid dt=", _dt, " s.");

  const Real duration_y = _dt / Corrosion::DRIDNModel::seconds_per_year;
  const auto stats = _model.advance(_context, _state, duration_y, _integration_options);
  _elapsed_time_y += duration_y;
  _accepted_steps += stats.accepted_steps;
  _rejected_steps += stats.rejected_steps;
}

Corrosion::DRIDNModel::Outputs
DRIDNCorrosionUserObject::outputs() const
{
  return _model.outputs(_context, _state, _elapsed_time_y);
}

Real
DRIDNCorrosionUserObject::scalarValue(const std::string & quantity) const
{
  const auto value = outputs();
  if (quantity == "front_rate_um_y" || quantity == "instantaneous_front_rate_um_y")
    return value.instantaneous_front_rate_um_y;
  if (quantity == "corrosion_rate_um_y" || quantity == "average_corrosion_rate_um_y")
    return value.average_corrosion_rate_um_y;
  if (quantity == "front_depth_um")
    return value.front_depth_um;
  if (quantity == "mass_recession_um")
    return value.mass_recession_um;
  if (quantity == "mass_loss_mg_cm2")
    return value.mass_loss_mg_cm2;
  if (quantity == "mass_gain_mg_cm2")
    return value.mass_gain_mg_cm2;
  if (quantity == "igc_depth_um")
    return value.igc_depth_um;
  if (quantity == "redox_log_shift" || quantity == "redox_shift")
    return value.redox_shift;
  if (quantity == "mass_balance_relative_error")
    return value.mass_balance_relative_error;
  if (quantity == "elapsed_time_y")
    return _elapsed_time_y;
  if (quantity == "accepted_internal_steps")
    return _accepted_steps;
  if (quantity == "rejected_internal_steps")
    return _rejected_steps;

  static constexpr std::array<const char *, 3> element_names{{"cr", "fe", "ni"}};
  for (unsigned int i = 0; i < element_names.size(); ++i)
  {
    const std::string element = element_names[i];
    if (quantity == "dissolved_inventory_" + element + "_ppm")
      return value.dissolved_ppm[i];
    if (quantity == "cumulative_source_" + element + "_ppm")
      return value.cumulative_source_ppm[i];
    if (quantity == "coupon_deposit_" + element + "_mg_cm2")
      return value.coupon_deposit_mg_cm2[i];
    if (quantity == "bulk_captured_" + element + "_ppm")
      return value.bulk_captured_ppm[i];
    if (quantity == "surface_availability_" + element)
      return value.surface_availability[i];
  }

  mooseError("DRIDNCorrosionUserObject '",
             name(),
             "' does not provide scalar quantity '",
             quantity,
             "'.");
}
