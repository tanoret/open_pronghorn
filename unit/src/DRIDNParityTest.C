//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "DRIDNModel.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace Corrosion
{
namespace
{

using CSVRow = std::map<std::string, std::string>;

struct CSVTable
{
  std::map<std::string, std::string> metadata;
  std::vector<std::string> columns;
  std::vector<CSVRow> rows;
};

std::string
fixtureFile(const std::string & name)
{
  const std::array<std::string, 5> prefixes = {"unit/data/dridn/",
                                               "../unit/data/dridn/",
                                               "data/dridn/",
                                               "../../unit/data/dridn/",
                                               ""};
  for (const auto & prefix : prefixes)
  {
    const std::string candidate = prefix + name;
    if (std::ifstream(candidate).good())
      return candidate;
  }
  throw std::runtime_error("Unable to locate DRIDN unit fixture '" + name + "'.");
}

std::vector<std::string>
splitCSVLine(const std::string & line)
{
  // These generated fixtures deliberately contain no quoted fields.  Refuse quotes instead of
  // silently pretending to implement general RFC-4180 parsing.
  if (line.find('"') != std::string::npos)
    throw std::runtime_error("DRIDN parity fixture contains an unsupported quoted CSV field.");
  std::vector<std::string> fields;
  std::size_t begin = 0;
  while (true)
  {
    const auto separator = line.find(',', begin);
    fields.push_back(line.substr(begin, separator - begin));
    if (separator == std::string::npos)
      break;
    begin = separator + 1;
  }
  return fields;
}

CSVTable
readCSV(const std::string & filename)
{
  std::ifstream stream(filename);
  if (!stream)
    throw std::runtime_error("Unable to read DRIDN parity fixture '" + filename + "'.");

  CSVTable table;
  std::string line;
  unsigned int line_number = 0;
  while (std::getline(stream, line))
  {
    ++line_number;
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    if (line.empty())
      continue;
    if (line.front() == '#')
    {
      const std::size_t begin = line.size() > 1 && line[1] == ' ' ? 2 : 1;
      const auto separator = line.find('=', begin);
      if (separator == std::string::npos)
        throw std::runtime_error("Malformed DRIDN fixture metadata at line " +
                                 std::to_string(line_number) + ".");
      const auto inserted = table.metadata.emplace(line.substr(begin, separator - begin),
                                                   line.substr(separator + 1));
      if (!inserted.second)
        throw std::runtime_error("Duplicate DRIDN fixture metadata key at line " +
                                 std::to_string(line_number) + ".");
      continue;
    }

    const auto fields = splitCSVLine(line);
    if (table.columns.empty())
    {
      table.columns = fields;
      continue;
    }
    if (fields.size() != table.columns.size())
      throw std::runtime_error("DRIDN fixture field-count mismatch at line " +
                               std::to_string(line_number) + ".");
    CSVRow row;
    for (std::size_t i = 0; i < fields.size(); ++i)
      if (!row.emplace(table.columns[i], fields[i]).second)
        throw std::runtime_error("DRIDN fixture contains a duplicate column name.");
    table.rows.push_back(std::move(row));
  }
  if (table.columns.empty() || table.rows.empty())
    throw std::runtime_error("DRIDN parity fixture is empty.");
  return table;
}

std::uint64_t
fnv1a64(const std::string & filename)
{
  std::ifstream stream(filename, std::ios::binary);
  if (!stream)
    throw std::runtime_error("Unable to hash DRIDN parity fixture '" + filename + "'.");
  std::uint64_t hash = UINT64_C(14695981039346656037);
  char character;
  while (stream.get(character))
  {
    hash ^= static_cast<unsigned char>(character);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

Real
number(const CSVRow & row, const std::string & name)
{
  const auto iterator = row.find(name);
  if (iterator == row.end())
    throw std::runtime_error("DRIDN parity fixture is missing column '" + name + "'.");
  std::size_t consumed = 0;
  const Real value = std::stod(iterator->second, &consumed);
  if (consumed != iterator->second.size() || !std::isfinite(value))
    throw std::runtime_error("DRIDN parity fixture has invalid number in column '" + name + "'.");
  return value;
}

std::string
textValue(const CSVRow & row, const std::string & name)
{
  const auto iterator = row.find(name);
  if (iterator == row.end() || iterator->second.empty())
    throw std::runtime_error("DRIDN parity fixture is missing text column '" + name + "'.");
  return iterator->second;
}

bool
booleanValue(const CSVRow & row, const std::string & name)
{
  const auto value = textValue(row, name);
  if (value == "1" || value == "true" || value == "True")
    return true;
  if (value == "0" || value == "false" || value == "False")
    return false;
  throw std::runtime_error("DRIDN parity fixture has invalid boolean in column '" + name + "'.");
}

std::map<std::string, Real>
parameterValues(const CSVTable & table)
{
  std::map<std::string, Real> values;
  for (const auto & row : table.rows)
  {
    const auto inserted = values.emplace(textValue(row, "parameter"), number(row, "value"));
    if (!inserted.second)
      throw std::runtime_error("DRIDN fitted-parameter fixture contains a duplicate name.");
  }
  if (values.size() != 26)
    throw std::runtime_error("DRIDN fitted-parameter fixture must contain exactly 26 values.");
  return values;
}

DRIDNModel::Parameters
fittedParameters(const std::map<std::string, Real> & value)
{
  DRIDNModel::Parameters result;
  result.log_rate_scale = value.at("log_rate_scale");
  result.affinity_feedback_scale = value.at("affinity_feedback_scale");
  result.log_surface_reservoir_um = value.at("log_surface_reservoir_um");
  result.log_surface_replenishment_y_inv = value.at("log_surface_replenishment_y_inv");
  result.surface_availability_exponent = value.at("surface_availability_exponent");
  result.surface_reservoir_cr_exponent = value.at("surface_reservoir_cr_exponent");
  result.log_dynamic_cr_exchange_bias = value.at("log_dynamic_cr_exchange_bias");
  result.log_dynamic_fe_capture_bias = value.at("log_dynamic_fe_capture_bias");
  result.inventory_inhibition_scale = value.at("inventory_inhibition_scale");
  result.log_redox_relaxation_y_inv = value.at("log_redox_relaxation_y_inv");
  result.redox_buffer_retention = value.at("redox_buffer_retention");
  result.redox_consumption_per_um = value.at("redox_consumption_per_um");
  result.log_stress_interfacial_factor = value.at("log_stress_interfacial_factor");
  result.log_fluoride_impurity_interfacial_factor =
      value.at("log_fluoride_impurity_interfacial_factor");
  result.log_deposition_rate_y_inv_fuel = value.at("log_deposition_rate_y_inv_fuel");
  result.log_deposition_rate_y_inv_flinak = value.at("log_deposition_rate_y_inv_flinak");
  result.bulk_capture_area_multiplier_fuel = value.at("bulk_capture_area_multiplier_fuel");
  result.bulk_capture_area_multiplier_flinak = value.at("bulk_capture_area_multiplier_flinak");
  result.log_bulk_precipitation_rate_scale = value.at("log_bulk_precipitation_rate_scale");
  result.log_inventory_scale_msre = value.at("log_inventory_scale_msre");
  result.log_inventory_scale_loop = value.at("log_inventory_scale_loop");
  result.log_deposit_area_scale_fuel = value.at("log_deposit_area_scale_fuel");
  result.log_deposit_area_scale_flinak = value.at("log_deposit_area_scale_flinak");
  result.log_mass_loss_scale = value.at("log_mass_loss_scale");
  result.gb_dynamic_scale = value.at("gb_dynamic_scale");
  result.damage_affinity_scale = value.at("damage_affinity_scale");
  result.validate();
  return result;
}

void
expectEveryParameter(const DRIDNModel::Parameters & actual,
                     const std::map<std::string, Real> & expected)
{
#define EXPECT_DRIDN_PARAMETER(name) EXPECT_DOUBLE_EQ(actual.name, expected.at(#name))
  EXPECT_DRIDN_PARAMETER(log_rate_scale);
  EXPECT_DRIDN_PARAMETER(affinity_feedback_scale);
  EXPECT_DRIDN_PARAMETER(log_surface_reservoir_um);
  EXPECT_DRIDN_PARAMETER(log_surface_replenishment_y_inv);
  EXPECT_DRIDN_PARAMETER(surface_availability_exponent);
  EXPECT_DRIDN_PARAMETER(surface_reservoir_cr_exponent);
  EXPECT_DRIDN_PARAMETER(log_dynamic_cr_exchange_bias);
  EXPECT_DRIDN_PARAMETER(log_dynamic_fe_capture_bias);
  EXPECT_DRIDN_PARAMETER(inventory_inhibition_scale);
  EXPECT_DRIDN_PARAMETER(log_redox_relaxation_y_inv);
  EXPECT_DRIDN_PARAMETER(redox_buffer_retention);
  EXPECT_DRIDN_PARAMETER(redox_consumption_per_um);
  EXPECT_DRIDN_PARAMETER(log_stress_interfacial_factor);
  EXPECT_DRIDN_PARAMETER(log_fluoride_impurity_interfacial_factor);
  EXPECT_DRIDN_PARAMETER(log_deposition_rate_y_inv_fuel);
  EXPECT_DRIDN_PARAMETER(log_deposition_rate_y_inv_flinak);
  EXPECT_DRIDN_PARAMETER(bulk_capture_area_multiplier_fuel);
  EXPECT_DRIDN_PARAMETER(bulk_capture_area_multiplier_flinak);
  EXPECT_DRIDN_PARAMETER(log_bulk_precipitation_rate_scale);
  EXPECT_DRIDN_PARAMETER(log_inventory_scale_msre);
  EXPECT_DRIDN_PARAMETER(log_inventory_scale_loop);
  EXPECT_DRIDN_PARAMETER(log_deposit_area_scale_fuel);
  EXPECT_DRIDN_PARAMETER(log_deposit_area_scale_flinak);
  EXPECT_DRIDN_PARAMETER(log_mass_loss_scale);
  EXPECT_DRIDN_PARAMETER(gb_dynamic_scale);
  EXPECT_DRIDN_PARAMETER(damage_affinity_scale);
#undef EXPECT_DRIDN_PARAMETER
}

DRIDNModel::Triplet
elementTriplet(const CSVRow & row, const std::string & prefix, const std::string & suffix)
{
  return {{number(row, prefix + "Cr" + suffix),
           number(row, prefix + "Fe" + suffix),
           number(row, prefix + "Ni" + suffix)}};
}

DRIDNModel::Context
contextFromFrozenRow(const CSVRow & row)
{
  DRIDNModel::Context context;
  context.mass_fractions = elementTriplet(row, "mass_fraction_", "");
  context.log_exchange_offsets = elementTriplet(row, "log_exchange_offset_", "");
  context.affinity_baseline = elementTriplet(row, "affinity_baseline_", "");
  context.cold_capture_fraction = elementTriplet(row, "cold_capture_fraction_", "");
  context.initial_dissolved_ppm = elementTriplet(row, "initial_dissolved_", "_ppm");
  context.product_activity_floor_ppm = number(row, "product_floor_ppm");

  context.cr_fraction_ratio = number(row, "cr_fraction_ratio");
  context.density_g_cm3 = number(row, "density_g_cm3");
  context.flow_factor = number(row, "flow_factor");
  context.selectivity_scale = number(row, "selectivity_scale");
  context.redox_shift_initial = number(row, "redox_shift_initial");
  context.log_charge_base_no_redox = number(row, "log_charge_base_no_redox");
  context.mass_transfer_rate_um_y = number(row, "mass_transfer_rate_um_y");
  context.inventory_capacity_ppm = number(row, "inventory_capacity_ppm");
  context.area_to_salt_mass_cm2_g = number(row, "area_to_salt_mass_cm2_g");
  context.explicit_inventory_scale = number(row, "explicit_inventory_scale");
  context.inventory_coupling_factor = number(row, "inventory_coupling_factor");
  context.deposit_area_factor = number(row, "deposit_area_factor");
  context.mass_loss_fraction = number(row, "mass_loss_fraction");
  context.cr_diffusion_cm2_s = number(row, "cr_diffusion_cm2_s");
  context.front_damage_multiplier = number(row, "front_damage_multiplier");
  context.gb_length_multiplier = number(row, "gb_length_multiplier");

  const auto & inventory = textValue(row, "inventory_scale");
  if (inventory == "msre")
    context.inventory_scale = DRIDNModel::InventoryScale::MSRE;
  else if (inventory == "loop")
    context.inventory_scale = DRIDNModel::InventoryScale::Loop;
  else
    throw std::runtime_error("Frozen DRIDN parity row has unknown inventory_scale.");

  const auto & deposition = textValue(row, "deposition_closure");
  if (deposition == "fuel")
    context.deposition_closure = DRIDNModel::DepositionClosure::FuelLike;
  else if (deposition == "flinak")
    context.deposition_closure = DRIDNModel::DepositionClosure::FLiNaK;
  else
    throw std::runtime_error("Frozen DRIDN parity row has unknown deposition_closure.");

  context.transient_redox = booleanValue(row, "transient_redox");
  context.stress_interfacial_activation = booleanValue(row, "stress_interfacial_activation");
  context.fluoride_impurity_interfacial_activation =
      booleanValue(row, "fluoride_impurity_interfacial_activation");
  context.chloride_salt = booleanValue(row, "chloride_salt");
  return context;
}

void
expectContextRoundTrip(const DRIDNModel::Context & context, const CSVRow & row)
{
  const auto mass = elementTriplet(row, "mass_fraction_", "");
  const auto exchange = elementTriplet(row, "log_exchange_offset_", "");
  const auto affinity = elementTriplet(row, "affinity_baseline_", "");
  const auto capture = elementTriplet(row, "cold_capture_fraction_", "");
  for (std::size_t i = 0; i < DRIDNModel::n_elements; ++i)
  {
    EXPECT_DOUBLE_EQ(context.mass_fractions[i], mass[i]);
    EXPECT_DOUBLE_EQ(context.log_exchange_offsets[i], exchange[i]);
    EXPECT_DOUBLE_EQ(context.affinity_baseline[i], affinity[i]);
    EXPECT_DOUBLE_EQ(context.cold_capture_fraction[i], capture[i]);
    static const std::array<std::string, 3> element{{"Cr", "Fe", "Ni"}};
    EXPECT_DOUBLE_EQ(context.initial_dissolved_ppm[i],
                     number(row, "initial_dissolved_" + element[i] + "_ppm"));
  }
  EXPECT_DOUBLE_EQ(context.cr_fraction_ratio, number(row, "cr_fraction_ratio"));
  EXPECT_DOUBLE_EQ(context.density_g_cm3, number(row, "density_g_cm3"));
  EXPECT_DOUBLE_EQ(context.flow_factor, number(row, "flow_factor"));
  EXPECT_DOUBLE_EQ(context.selectivity_scale, number(row, "selectivity_scale"));
  EXPECT_DOUBLE_EQ(context.product_activity_floor_ppm, number(row, "product_floor_ppm"));
  EXPECT_DOUBLE_EQ(context.redox_shift_initial, number(row, "redox_shift_initial"));
  EXPECT_DOUBLE_EQ(context.log_charge_base_no_redox, number(row, "log_charge_base_no_redox"));
  EXPECT_DOUBLE_EQ(context.mass_transfer_rate_um_y, number(row, "mass_transfer_rate_um_y"));
  EXPECT_DOUBLE_EQ(context.inventory_capacity_ppm, number(row, "inventory_capacity_ppm"));
  EXPECT_DOUBLE_EQ(context.area_to_salt_mass_cm2_g, number(row, "area_to_salt_mass_cm2_g"));
  EXPECT_DOUBLE_EQ(context.explicit_inventory_scale, number(row, "explicit_inventory_scale"));
  EXPECT_DOUBLE_EQ(context.inventory_coupling_factor, number(row, "inventory_coupling_factor"));
  EXPECT_DOUBLE_EQ(context.deposit_area_factor, number(row, "deposit_area_factor"));
  EXPECT_DOUBLE_EQ(context.mass_loss_fraction, number(row, "mass_loss_fraction"));
  EXPECT_DOUBLE_EQ(context.cr_diffusion_cm2_s, number(row, "cr_diffusion_cm2_s"));
  EXPECT_DOUBLE_EQ(context.front_damage_multiplier, number(row, "front_damage_multiplier"));
  EXPECT_DOUBLE_EQ(context.gb_length_multiplier, number(row, "gb_length_multiplier"));

  const auto & inventory = textValue(row, "inventory_scale");
  EXPECT_EQ(context.inventory_scale,
            inventory == "msre" ? DRIDNModel::InventoryScale::MSRE
                                : DRIDNModel::InventoryScale::Loop);
  const auto & deposition = textValue(row, "deposition_closure");
  EXPECT_EQ(context.deposition_closure,
            deposition == "fuel" ? DRIDNModel::DepositionClosure::FuelLike
                                  : DRIDNModel::DepositionClosure::FLiNaK);
  EXPECT_EQ(context.transient_redox, booleanValue(row, "transient_redox"));
  EXPECT_EQ(context.stress_interfacial_activation,
            booleanValue(row, "stress_interfacial_activation"));
  EXPECT_EQ(context.fluoride_impurity_interfacial_activation,
            booleanValue(row, "fluoride_impurity_interfacial_activation"));
  EXPECT_EQ(context.chloride_salt, booleanValue(row, "chloride_salt"));
}

void
expectPythonEndpoint(const Real actual, const Real expected)
{
  // Python used LSODA with rtol=2e-5 and max_step=duration/120.  The independent adaptive RK4
  // implementation is held to 2e-6 relative, comfortably above the observed 1.27e-7 maximum,
  // plus a 1e-11 absolute allowance for values analytically equal to zero.
  EXPECT_NEAR(actual, expected, 1.0e-11 + 2.0e-6 * std::abs(expected));
}

} // namespace

TEST(DRIDNParity, FixtureIntegrityAndAllFittedParameters)
{
  const auto parity_filename = fixtureFile("cpp_parity_cases.csv");
  const auto parameter_filename = fixtureFile("fitted_parameters.csv");
  EXPECT_EQ(fnv1a64(parity_filename), UINT64_C(0xad422aa9824fd6f2));
  EXPECT_EQ(fnv1a64(parameter_filename), UINT64_C(0x113b92c697d33683));

  const auto parity = readCSV(parity_filename);
  EXPECT_EQ(parity.metadata.at("source"),
            "validation/corrosion/calibration/results/advanced/cpp_parity_cases.csv");
  EXPECT_EQ(parity.metadata.at("source_sha256"),
            "7d2fc2c5e87d101ab08af2f4ff953526d6609a2bc1e5bedee3c73bd891ba2126");
  EXPECT_EQ(parity.rows.size(), 10u);

  const auto parameter_table = readCSV(parameter_filename);
  EXPECT_EQ(parameter_table.metadata.at("source"),
            "validation/corrosion/calibration/results/advanced/dynamic_network_parameters.json");
  EXPECT_EQ(parameter_table.metadata.at("source_sha256"),
            "3ed42782d7ea1a3bb2ee529abadf1b2b476e5b86d9d656ffa20b89936354c20f");
  const auto values = parameterValues(parameter_table);
  const auto parameters = fittedParameters(values);
  expectEveryParameter(parameters, values);
}

TEST(DRIDNParity, FrozenPythonLegacyEndpoints)
{
  const auto parameter_table = readCSV(fixtureFile("fitted_parameters.csv"));
  const auto model = DRIDNModel(fittedParameters(parameterValues(parameter_table)),
                                DRIDNModel::ClosureConstants(),
                                DRIDNModel::ModelOptions::legacyCompatibility());
  EXPECT_EQ(model.options().charge_transfer, DRIDNModel::ChargeTransferMode::LegacyIrreversible);

  const auto parity = readCSV(fixtureFile("cpp_parity_cases.csv"));
  const std::vector<std::string> expected_ids = {
      "M-003",
      "M-005",
      "M-014",
      "M-018",
      "M-027",
      "M-029",
      "M-030",
      "M-038",
      "M-041",
      "boundary_zero_elements"};
  std::vector<std::string> actual_ids;
  for (const auto & row : parity.rows)
  {
    const auto & measurement_id = textValue(row, "measurement_id");
    SCOPED_TRACE(measurement_id);
    actual_ids.push_back(measurement_id);

    const auto context = contextFromFrozenRow(row);
    expectContextRoundTrip(context, row);
    context.validate(model.closures());
    auto state = model.initialState(context);
    const auto initial_rate = model.rates(context, state);
    for (std::size_t i = 0; i < DRIDNModel::n_elements; ++i)
      if (context.mass_fractions[i] == 0.0)
      {
        EXPECT_DOUBLE_EQ(initial_rate.species_fraction[i], 0.0);
        EXPECT_DOUBLE_EQ(initial_rate.source_rate_ppm_y[i], 0.0);
      }
    model.advance(context, state, number(row, "time_years"));
    const auto output = model.outputs(context, state, number(row, "time_years"));

    expectPythonEndpoint(output.front_depth_um, number(row, "front_depth_um"));
    expectPythonEndpoint(output.mass_recession_um, number(row, "mass_recession_um"));
    expectPythonEndpoint(output.mass_loss_mg_cm2, number(row, "mass_loss_mg_cm2"));
    expectPythonEndpoint(output.mass_gain_mg_cm2, number(row, "mass_gain_mg_cm2"));
    expectPythonEndpoint(output.igc_depth_um, number(row, "igc_depth_um"));
    expectPythonEndpoint(output.average_corrosion_rate_um_y, number(row, "corrosion_rate_um_y"));
    expectPythonEndpoint(output.redox_shift, number(row, "redox_shift_endpoint"));
    for (std::size_t i = 0; i < DRIDNModel::n_elements; ++i)
    {
      static const std::array<std::string, 3> element{{"Cr", "Fe", "Ni"}};
      expectPythonEndpoint(output.dissolved_ppm[i],
                           number(row, "dissolved_" + element[i] + "_ppm"));
      expectPythonEndpoint(output.cumulative_source_ppm[i],
                           number(row, "cumulative_source_" + element[i] + "_ppm"));
      if (context.mass_fractions[i] == 0.0)
        EXPECT_DOUBLE_EQ(output.cumulative_source_ppm[i], 0.0);
      expectPythonEndpoint(output.coupon_deposit_mg_cm2[i],
                           number(row, "coupon_deposit_" + element[i] + "_mg_cm2"));
      expectPythonEndpoint(output.bulk_captured_ppm[i],
                           number(row, "bulk_captured_" + element[i] + "_ppm"));
      expectPythonEndpoint(output.surface_availability[i],
                           number(row, "surface_" + element[i] + "_availability"));
    }
    EXPECT_LT(number(row, "mass_balance_relative_error"), 1.0e-10);
    EXPECT_LT(output.mass_balance_relative_error, 1.0e-10);
  }
  EXPECT_EQ(actual_ids, expected_ids);
}

TEST(DRIDNParity, DescriptiveLabelsCannotChangeExplicitPhysics)
{
  const auto parameter_table = readCSV(fixtureFile("fitted_parameters.csv"));
  const auto model = DRIDNModel(fittedParameters(parameterValues(parameter_table)),
                                DRIDNModel::ClosureConstants(),
                                DRIDNModel::ModelOptions::legacyCompatibility());
  const auto parity = readCSV(fixtureFile("cpp_parity_cases.csv"));
  const auto iterator =
      std::find_if(parity.rows.begin(), parity.rows.end(), [](const CSVRow & row) {
        return textValue(row, "measurement_id") == "M-038";
      });
  if (iterator == parity.rows.end())
    throw std::runtime_error("Frozen DRIDN parity fixture is missing M-038.");

  CSVRow renamed = *iterator;
  renamed["measurement_id"] = "renamed-measurement";
  renamed["response_kind"] = "renamed-response";
  renamed["material"] = "renamed-material";
  renamed["salt_class"] = "renamed-salt-label";
  renamed["redox_class"] = "renamed-redox-label";
  renamed["position_class"] = "renamed-position";

  const auto original_context = contextFromFrozenRow(*iterator);
  const auto renamed_context = contextFromFrozenRow(renamed);
  auto original_state = model.initialState(original_context);
  auto renamed_state = model.initialState(renamed_context);
  const Real duration_y = number(*iterator, "time_years");
  model.advance(original_context, original_state, duration_y);
  model.advance(renamed_context, renamed_state, duration_y);
  EXPECT_EQ(original_state, renamed_state);
}

} // namespace Corrosion
