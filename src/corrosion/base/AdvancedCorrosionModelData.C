//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "AdvancedCorrosionModelData.h"

#include "MoltenSaltCorrosionData.h"
#include "MooseError.h"
#include "MooseUtils.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>

namespace Corrosion
{

const std::string &
AdvancedCorrosionModelDatabase::supportedModelRevision()
{
  static const std::string revision = "mstdb-nst-v2-fe-identity_dridn-v2-explicit-geometry";
  return revision;
}

AdvancedCorrosionModelDatabase::AdvancedCorrosionModelDatabase(const std::string & filename)
  : _filename(filename)
{
  MooseUtils::checkFileReadable(_filename);
  std::ifstream stream(_filename);
  try
  {
    stream >> _root;
  }
  catch (const std::exception & error)
  {
    mooseError("Advanced corrosion models: failed to parse '", _filename, "':\n", error.what());
  }

  _schema_version = _root.value("schema_version", "");
  _calibration_id = _root.value("calibration_id", "");
  if (_schema_version != "1.0")
    mooseError("Advanced corrosion models: unsupported schema_version '",
               _schema_version,
               "' in '",
               _filename,
               "' (expected 1.0).");
  if (!_root.contains("thermochemical_parameters") ||
      !_root.at("thermochemical_parameters").is_object())
    mooseError("Advanced corrosion models: '",
               _filename,
               "' has no thermochemical_parameters object.");
  if (!_root.contains("dynamic_parameters") || !_root.at("dynamic_parameters").is_object())
    mooseError("Advanced corrosion models: '", _filename, "' has no dynamic_parameters object.");

  _model_revision = _root.value("model_revision", "");
  _calibration_data_revision = _root.value("calibration_data_revision", "");
  if (_model_revision.empty() || _calibration_id.empty() || _calibration_data_revision.empty())
    mooseError("Advanced corrosion models: '",
               _filename,
               "' must define nonempty calibration_id, model_revision, and "
               "calibration_data_revision fields.");
  if (_model_revision != supportedModelRevision())
    mooseError("Advanced corrosion models: unsupported model_revision '",
               _model_revision,
               "' in '",
               _filename,
               "' (expected '",
               supportedModelRevision(),
               "'). Refit and revalidate the artifact for this model law before loading it.");
  if (!_root.contains("mstdb_provenance") || !_root.at("mstdb_provenance").is_object())
    mooseError("Advanced corrosion models: '", _filename, "' has no mstdb_provenance object.");
  const auto & provenance = _root.at("mstdb_provenance");
  _expected_mstdb_version = provenance.value("version", "");
  _expected_fluoride_sha256 = provenance.value("fluoride_sha256", "");
  _expected_chloride_sha256 = provenance.value("chloride_sha256", "");
  if (_expected_mstdb_version.empty() || _expected_fluoride_sha256.size() != 64 ||
      _expected_chloride_sha256.size() != 64)
    mooseError("Advanced corrosion models: '",
               _filename,
               "' must bind a nonempty MSTDB version and 64-character fluoride/chloride SHA-256 "
               "hashes.");

  if (!_root.contains("base_model_provenance") || !_root.at("base_model_provenance").is_object())
    mooseError("Advanced corrosion models: '",
               _filename,
               "' has no base_model_provenance object.");
  const auto & base_provenance = _root.at("base_model_provenance");
  _base_model_source_sha256 = base_provenance.value("source_sha256", "");
  const auto valid_sha256 = [](const std::string & value)
  {
    return value.size() == 64 &&
           std::all_of(value.begin(), value.end(), [](const unsigned char character)
                       {
                         return std::isxdigit(character);
                       });
  };
  if (!valid_sha256(_expected_fluoride_sha256) ||
      !valid_sha256(_expected_chloride_sha256) || !valid_sha256(_base_model_source_sha256))
    mooseError("Advanced corrosion models: '",
               _filename,
               "' must contain hexadecimal 64-character SHA-256 values for both MSTDB files and "
               "the base corrosion database source.");
  const auto normalize_sha256 = [](std::string & value)
  {
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char character)
                   {
                     return static_cast<char>(std::tolower(character));
                   });
  };
  normalize_sha256(_expected_fluoride_sha256);
  normalize_sha256(_expected_chloride_sha256);
  normalize_sha256(_base_model_source_sha256);

  const auto read_finite_map = [this, &base_provenance](const std::string & key,
                                                        std::map<std::string, Real> & destination)
  {
    if (!base_provenance.contains(key) || !base_provenance.at(key).is_object())
      mooseError("Advanced corrosion models: base_model_provenance in '",
                 _filename,
                 "' has no '",
                 key,
                 "' object.");
    for (const auto & item : base_provenance.at(key).items())
    {
      const Real value = item.value().get<Real>();
      if (!std::isfinite(value))
        mooseError("Advanced corrosion models: nonfinite base-model value '",
                   item.key(),
                   "' in '",
                   _filename,
                   "'.");
      destination[item.key()] = value;
    }
  };
  read_finite_map("required_parameters", _base_model_parameters);
  read_finite_map("required_densities_g_cm3", _base_model_densities);
  if (!base_provenance.contains("required_element_properties") ||
      !base_provenance.at("required_element_properties").is_object())
    mooseError("Advanced corrosion models: base_model_provenance in '",
               _filename,
               "' has no 'required_element_properties' object.");
  for (const auto & element : base_provenance.at("required_element_properties").items())
  {
    if (!element.value().is_object())
      mooseError("Advanced corrosion models: base-model element '",
                 element.key(),
                 "' in '",
                 _filename,
                 "' must be an object.");
    for (const auto & property : element.value().items())
    {
      const Real value = property.value().get<Real>();
      if (!std::isfinite(value))
        mooseError("Advanced corrosion models: nonfinite base-model element property '",
                   element.key(),
                   ".",
                   property.key(),
                   "' in '",
                   _filename,
                   "'.");
      _base_model_element_properties[element.key()][property.key()] = value;
    }
  }

  for (const auto & item : _root.at("thermochemical_parameters").items())
  {
    const Real value = item.value().get<Real>();
    if (!std::isfinite(value))
      mooseError("Advanced corrosion models: nonfinite thermochemical parameter '",
                 item.key(),
                 "' in '",
                 _filename,
                 "'.");
    _thermochemical_parameters[item.key()] = value;
  }
  for (const auto & item : _root.at("dynamic_parameters").items())
  {
    const Real value = item.value().get<Real>();
    if (!std::isfinite(value))
      mooseError("Advanced corrosion models: nonfinite DRIDN parameter '",
                 item.key(),
                 "' in '",
                 _filename,
                 "'.");
    _dynamic_parameters[item.key()] = value;
  }

  validateRequiredParameters();
}

Real
AdvancedCorrosionModelDatabase::lookup(const std::map<std::string, Real> & values,
                                       const std::string & name,
                                       const std::string & section) const
{
  const auto found = values.find(name);
  if (found == values.end())
    mooseError("Advanced corrosion models: parameter '",
               name,
               "' is missing from ",
               section,
               " in '",
               _filename,
               "'.");
  return found->second;
}

Real
AdvancedCorrosionModelDatabase::thermochemicalParameter(const std::string & name) const
{
  return lookup(_thermochemical_parameters, name, "thermochemical_parameters");
}

Real
AdvancedCorrosionModelDatabase::dynamicParameter(const std::string & name) const
{
  return lookup(_dynamic_parameters, name, "dynamic_parameters");
}

Real
AdvancedCorrosionModelDatabase::baseModelParameter(const std::string & name) const
{
  return lookup(_base_model_parameters, name, "base_model_provenance.required_parameters");
}

Real
AdvancedCorrosionModelDatabase::baseModelDensity(const std::string & material_class) const
{
  return lookup(_base_model_densities,
                material_class,
                "base_model_provenance.required_densities_g_cm3");
}

Real
AdvancedCorrosionModelDatabase::baseModelElementProperty(const std::string & element,
                                                          const std::string & property) const
{
  const auto found = _base_model_element_properties.find(element);
  if (found == _base_model_element_properties.end())
    mooseError("Advanced corrosion models: base-model element '",
               element,
               "' is missing from base_model_provenance.required_element_properties in '",
               _filename,
               "'.");
  return lookup(found->second,
                property,
                "base_model_provenance.required_element_properties." + element);
}

void
AdvancedCorrosionModelDatabase::validateBaseModel(
    const MoltenSaltCorrosionDatabase & database) const
{
  const auto require_equal = [this](const std::string & category,
                                    const std::string & name,
                                    const Real actual,
                                    const Real expected)
  {
    if (!std::isfinite(actual) || actual != expected)
      mooseError("Advanced corrosion models: base-model ",
                 category,
                 " '",
                 name,
                 "' does not match provenance source SHA-256 ",
                 _base_model_source_sha256,
                 ": expected ",
                 expected,
                 ", got ",
                 actual,
                 ".");
  };
  for (const auto & [name, expected] : _base_model_parameters)
    require_equal("parameter", name, database.parameter(name), expected);
  for (const auto & [name, expected] : _base_model_densities)
    require_equal("density", name, database.density(name), expected);

  struct Property
  {
    const char * name;
    Real ElementProperties::* member;
  };
  static const std::array<Property, 7> properties{{
      {"valence", &ElementProperties::valence},
      {"molar_mass_g_mol", &ElementProperties::molar_mass_g_mol},
      {"diffusivity_m2_s", &ElementProperties::diffusivity_m2_s},
      {"E0_V", &ElementProperties::E0_V},
      {"alpha_a", &ElementProperties::alpha_a},
      {"alpha_c", &ElementProperties::alpha_c},
      {"c_ref_mol_m3", &ElementProperties::c_ref_mol_m3}}};
  for (const auto & [element, expected_properties] : _base_model_element_properties)
  {
    // MoltenSaltCorrosionDatabase::element() deliberately falls back to generic_metal.  A
    // provenance check must not accept that fallback in place of an explicitly bound Cr/Fe/Ni
    // record, even if a user happens to give generic_metal identical numeric properties.
    if (!database.hasElement(element))
      mooseError("Advanced corrosion models: required explicit base-model element '",
                 element,
                 "' is missing for provenance source SHA-256 ",
                 _base_model_source_sha256,
                 ".");
    const auto actual = database.element(element);
    for (const auto & property : properties)
      require_equal("element property",
                    element + "." + property.name,
                    actual.*(property.member),
                    expected_properties.at(property.name));
  }
}

void
AdvancedCorrosionModelDatabase::validateRequiredParameters() const
{
  struct Bound
  {
    Real lower;
    Real upper;
  };
  static const std::array<const char *, 28> thermo_required = {
      "log_front_rate0_um_y",
      "cr_activity_exponent",
      "log_gamma_cr_flinak",
      "log_gamma_cr_flibe",
      "log_gamma_cr_fluoroborate",
      "log_gamma_cr_chloride",
      "selectivity_affinity_scale",
      "log_exchange_fe_relative",
      "log_exchange_ni_relative",
      "log_product_floor_ppm",
      "log_mass_transfer_cap_um_y",
      "Ea_mass_transfer_kJ_mol",
      "flow_mass_transfer_exponent",
      "gb_diffusion_multiplier",
      "multi_alloy_gb_multiplier",
      "x750_gb_multiplier",
      "chloride_reaction_multiplier",
      "chloride_gb_multiplier",
      "tellurium_gb_multiplier",
      "deposition_capture_area_factor_fuel",
      "deposition_capture_area_factor_flinak",
      "log_area_to_salt_mass_msre_cm2_g",
      "log_area_to_salt_mass_loop_cm2_g",
      "log_initial_fe2_ppm",
      "mass_loss_fraction_logit",
      "mass_loss_cr_exponent",
      "flinak_high_cr_selectivity_exponent",
      "log_inventory_capacity_ppm"};
  static const std::array<const char *, 26> dynamic_required = {
      "log_rate_scale",
      "affinity_feedback_scale",
      "log_surface_reservoir_um",
      "log_surface_replenishment_y_inv",
      "surface_availability_exponent",
      "surface_reservoir_cr_exponent",
      "log_dynamic_cr_exchange_bias",
      "log_dynamic_fe_capture_bias",
      "inventory_inhibition_scale",
      "log_redox_relaxation_y_inv",
      "redox_buffer_retention",
      "redox_consumption_per_um",
      "log_stress_interfacial_factor",
      "log_fluoride_impurity_interfacial_factor",
      "log_deposition_rate_y_inv_fuel",
      "log_deposition_rate_y_inv_flinak",
      "bulk_capture_area_multiplier_fuel",
      "bulk_capture_area_multiplier_flinak",
      "log_bulk_precipitation_rate_scale",
      "log_inventory_scale_msre",
      "log_inventory_scale_loop",
      "log_deposit_area_scale_fuel",
      "log_deposit_area_scale_flinak",
      "log_mass_loss_scale",
      "gb_dynamic_scale",
      "damage_affinity_scale"};
  static const std::array<const char *, 14> base_parameter_required = {
      "Ea_corr_kJ_mol",
      "log_Dcr_ref_cm2_s",
      "Ea_Dcr_kJ_mol",
      "redox_purified_baseline",
      "redox_msre_or_fuel_baseline",
      "redox_oxidizing_fef2",
      "redox_reducing_be",
      "redox_impure_moisture",
      "redox_chloride_unspecified",
      "redox_tellurium",
      "redox_stressed",
      "redox_multi_alloy",
      "redox_fission_product",
      "redox_gas_control"};
  static const std::array<const char *, 12> base_density_required = {
      "hastelloy_n",
      "modified_hastelloy_n",
      "gh3535",
      "stainless_304",
      "stainless_304l",
      "stainless_316",
      "stainless_316h",
      "stainless_316l",
      "alloy_x750",
      "in625",
      "ni_alloy",
      "generic_metal"};
  static const std::array<const char *, 3> base_element_required = {{"Cr", "Fe", "Ni"}};
  static const std::array<const char *, 7> base_element_property_required = {
      {"valence",
       "molar_mass_g_mol",
       "diffusivity_m2_s",
       "E0_V",
       "alpha_a",
       "alpha_c",
       "c_ref_mol_m3"}};

  for (const auto * name : thermo_required)
    if (!_thermochemical_parameters.count(name))
      mooseError("Advanced corrosion models: required thermochemical parameter '",
                 name,
                 "' is missing from '",
                 _filename,
                 "'.");
  for (const auto * name : dynamic_required)
    if (!_dynamic_parameters.count(name))
      mooseError("Advanced corrosion models: required DRIDN parameter '",
                 name,
                 "' is missing from '",
                 _filename,
                 "'.");
  for (const auto * name : base_parameter_required)
    if (!_base_model_parameters.count(name))
      mooseError("Advanced corrosion models: required base-model parameter '",
                 name,
                 "' is missing from '",
                 _filename,
                 "'.");
  for (const auto * name : base_density_required)
    if (!_base_model_densities.count(name))
      mooseError("Advanced corrosion models: required base-model density '",
                 name,
                 "' is missing from '",
                 _filename,
                 "'.");
  for (const auto * element : base_element_required)
  {
    const auto found = _base_model_element_properties.find(element);
    if (found == _base_model_element_properties.end())
      mooseError("Advanced corrosion models: required base-model element '",
                 element,
                 "' is missing from '",
                 _filename,
                 "'.");
    for (const auto * property : base_element_property_required)
      if (!found->second.count(property))
        mooseError("Advanced corrosion models: required base-model element property '",
                   element,
                   ".",
                   property,
                   "' is missing from '",
                   _filename,
                   "'.");
  }

  const std::map<std::string, Bound> thermo_bounds = {
      {"log_front_rate0_um_y", {std::log(0.01), std::log(100.0)}},
      {"cr_activity_exponent", {0.2, 4.0}},
      {"log_gamma_cr_flinak", {-12.0, 4.0}},
      {"log_gamma_cr_flibe", {-12.0, 4.0}},
      {"log_gamma_cr_fluoroborate", {-12.0, 4.0}},
      {"log_gamma_cr_chloride", {-14.0, 4.0}},
      {"selectivity_affinity_scale", {0.0, 1.0}},
      {"log_exchange_fe_relative", {-5.0, 2.0}},
      {"log_exchange_ni_relative", {-7.0, 1.0}},
      {"log_product_floor_ppm", {std::log(1.0e-4), std::log(300.0)}},
      {"log_mass_transfer_cap_um_y", {std::log(10.0), std::log(1.0e5)}},
      {"Ea_mass_transfer_kJ_mol", {0.0, 120.0}},
      {"flow_mass_transfer_exponent", {0.05, 2.0}},
      {"gb_diffusion_multiplier", {0.05, 3.0}},
      {"multi_alloy_gb_multiplier", {1.0, 40.0}},
      {"x750_gb_multiplier", {1.0, 40.0}},
      {"chloride_reaction_multiplier", {0.5, 30.0}},
      {"chloride_gb_multiplier", {0.5, 30.0}},
      {"tellurium_gb_multiplier", {1.0, 40.0}},
      {"deposition_capture_area_factor_fuel", {0.02, 8.0}},
      {"deposition_capture_area_factor_flinak", {0.02, 10.0}},
      {"log_area_to_salt_mass_msre_cm2_g", {std::log(0.002), std::log(5.0)}},
      {"log_area_to_salt_mass_loop_cm2_g", {std::log(0.03), std::log(80.0)}},
      {"log_initial_fe2_ppm", {std::log(1.0), std::log(3000.0)}},
      {"mass_loss_fraction_logit", {-5.0, 5.0}},
      {"mass_loss_cr_exponent", {-2.5, 1.5}},
      {"flinak_high_cr_selectivity_exponent", {-3.5, 0.5}},
      {"log_inventory_capacity_ppm", {std::log(100.0), std::log(2.0e4)}}};
  const std::map<std::string, Bound> dynamic_bounds = {
      {"log_rate_scale", {-1.5, 1.5}},
      {"affinity_feedback_scale", {0.0, 0.30}},
      {"log_surface_reservoir_um", {std::log(0.5), std::log(500.0)}},
      {"log_surface_replenishment_y_inv", {std::log(0.005), std::log(200.0)}},
      {"surface_availability_exponent", {0.0, 3.0}},
      {"surface_reservoir_cr_exponent", {0.0, 2.5}},
      {"log_dynamic_cr_exchange_bias", {-0.5, 0.5}},
      {"log_dynamic_fe_capture_bias", {-0.5, 0.5}},
      {"inventory_inhibition_scale", {0.0, 6.0}},
      {"log_redox_relaxation_y_inv", {std::log(0.005), std::log(200.0)}},
      {"redox_buffer_retention", {0.0, 1.0}},
      {"redox_consumption_per_um", {0.0, 0.25}},
      {"log_stress_interfacial_factor", {0.0, std::log(12.0)}},
      {"log_fluoride_impurity_interfacial_factor", {0.0, std::log(12.0)}},
      {"log_deposition_rate_y_inv_fuel", {std::log(1.0e-3), std::log(100.0)}},
      {"log_deposition_rate_y_inv_flinak", {std::log(1.0e-3), std::log(200.0)}},
      {"bulk_capture_area_multiplier_fuel", {0.0, 12.0}},
      {"bulk_capture_area_multiplier_flinak", {0.0, 12.0}},
      {"log_bulk_precipitation_rate_scale", {std::log(0.05), std::log(50.0)}},
      {"log_inventory_scale_msre", {-2.0, 2.0}},
      {"log_inventory_scale_loop", {-2.0, 2.0}},
      {"log_deposit_area_scale_fuel", {-std::log(10.0), std::log(10.0)}},
      {"log_deposit_area_scale_flinak", {-std::log(10.0), std::log(10.0)}},
      {"log_mass_loss_scale", {-1.0, 1.0}},
      {"gb_dynamic_scale", {0.2, 3.0}},
      {"damage_affinity_scale", {0.0, 0.40}}};

  auto check_bounds = [this](const std::map<std::string, Real> & values,
                             const std::map<std::string, Bound> & bounds,
                             const std::string & section)
  {
    for (const auto & [name, bound] : bounds)
    {
      const Real value = values.at(name);
      if (value < bound.lower || value > bound.upper)
        mooseError("Advanced corrosion models: parameter '",
                   name,
                   "'=",
                   value,
                   " in ",
                   section,
                   " is outside calibrated bounds [",
                   bound.lower,
                   ", ",
                   bound.upper,
                   "] in '",
                   _filename,
                   "'.");
    }
  };
  check_bounds(_thermochemical_parameters, thermo_bounds, "thermochemical_parameters");
  check_bounds(_dynamic_parameters, dynamic_bounds, "dynamic_parameters");
}

} // namespace Corrosion
