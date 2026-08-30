//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MSTDBTCCorrosionUserObject.h"

#include "MooseError.h"

#include <array>
#include <cmath>
#include <string>

registerMooseObject("OpenPronghornApp", MSTDBTCCorrosionUserObject);

namespace
{
Corrosion::MSTDBTCCorrosionFeatures
featuresFromParameters(const InputParameters & parameters)
{
  Corrosion::MSTDBTCCorrosionFeatures features;
  features.hot_temperature_K = parameters.get<Real>("hot_temperature");
  features.cold_temperature_K = parameters.get<Real>("cold_temperature");
  features.exposure_s = parameters.get<Real>("exposure_time");
  features.flow_factor = parameters.get<Real>("flow_factor");
  features.area_to_salt_mass_cm2_g = parameters.get<Real>("area_to_salt_mass");
  features.inventory_coupling_factor = parameters.get<Real>("inventory_coupling_factor");
  features.material_class = parameters.get<std::string>("material_class");
  features.salt_class = parameters.get<std::string>("salt_class");
  features.redox_class = parameters.get<std::string>("redox_class");
  return features;
}
} // namespace

InputParameters
MSTDBTCCorrosionUserObject::validParams()
{
  InputParameters params = AdvancedCorrosionModelUserObject::validParams();
  params.addClassDescription(
      "Evaluates one static Cr/Fe/Ni corrosion endpoint from MSTDB-TC standard-state Gibbs data "
      "and calibrated engineering closures. This is not a SUBQ or Gibbs-minimization solver.");

  params.addParam<DataFileName>("corrosion_database",
                                "corrosion_database.json",
                                "Base corrosion material and reduced-model database.");
  params.addParam<DataFileName>("advanced_database",
                                "advanced_corrosion_models.json",
                                "Advanced-model parameters plus mandatory MSTDB-TC and semantic "
                                "base-corrosion-database provenance.");
  params.addRequiredParam<FileName>(
      "fluoride_database",
      "Authorized external MSTDB-TC V3.1 fluoride *_No_Func.dat file.");
  params.addRequiredParam<FileName>(
      "chloride_database",
      "Authorized external MSTDB-TC V3.1 chloride *_No_Func.dat file.");
  params.addParam<bool>(
      "allow_extrapolation",
      false,
      "Permit evaluation above a species' last Gibbs interval. Keep false unless an extrapolation "
      "basis is documented for the application.");

  params.addRequiredRangeCheckedParam<Real>(
      "hot_temperature",
      "hot_temperature > 0",
      "Hot-side salt temperature [K].");
  params.addRequiredRangeCheckedParam<Real>(
      "cold_temperature",
      "cold_temperature > 0",
      "Cold-side salt temperature [K].");
  params.addRequiredRangeCheckedParam<Real>(
      "exposure_time",
      "exposure_time >= 0",
      "Physical exposure time [s].");
  params.addRequiredRangeCheckedParam<Real>(
      "flow_factor",
      "flow_factor > 0",
      "Dimensionless circulation/mass-transfer factor.");
  params.addRequiredRangeCheckedParam<Real>(
      "area_to_salt_mass",
      "area_to_salt_mass > 0",
      "Explicit wetted-area to salt-mass ratio [cm^2/g].");
  params.addRequiredRangeCheckedParam<Real>(
      "inventory_coupling_factor",
      "inventory_coupling_factor >= 0 & inventory_coupling_factor <= 1",
      "Explicit fraction applied only to dissolved inventory in the modeled salt [0,1]; it does "
      "not scale the static cold-capture mass_gain closure.");
  params.addRequiredParam<std::string>("material_class", "Explicit alloy/material class.");
  params.addRequiredParam<std::string>("salt_class", "Explicit salt chemistry class.");
  params.addRequiredParam<std::string>("redox_class", "Explicit redox condition class.");

  return params;
}

MSTDBTCCorrosionUserObject::MSTDBTCCorrosionUserObject(const InputParameters & parameters)
  : AdvancedCorrosionModelUserObject(parameters),
    _base_database(getParam<DataFileName>("corrosion_database")),
    _advanced_database(getParam<DataFileName>("advanced_database")),
    _thermodynamics(getParam<FileName>("fluoride_database"),
                    getParam<FileName>("chloride_database"),
                    _advanced_database.expectedMSTDBVersion(),
                    _advanced_database.expectedFluorideSHA256(),
                    _advanced_database.expectedChlorideSHA256(),
                    false,
                    getParam<bool>("allow_extrapolation")),
    _model(_base_database, _advanced_database, _thermodynamics),
    _features(featuresFromParameters(parameters)),
    _result(_model.evaluate(_features))
{
  if (_advanced_database.expectedMSTDBVersion() != "3.1")
    paramError("advanced_database",
               "This model release requires MSTDB-TC V3.1, but the parameter database binds "
               "version '",
               _advanced_database.expectedMSTDBVersion(),
               "'. Recalibration and revalidation are required before changing the edition.");
  if (_features.cold_temperature_K > _features.hot_temperature_K)
    paramError("cold_temperature", "cold_temperature must not exceed hot_temperature.");
}

Real
MSTDBTCCorrosionUserObject::scalarValue(const std::string & quantity) const
{
  if (quantity == "front_rate_um_y")
    return _result.front_rate_um_y;
  if (quantity == "corrosion_rate_um_y")
    return _result.corrosion_rate_um_y;
  if (quantity == "front_depth_um")
    return _result.front_depth_um;
  if (quantity == "mass_loss_mg_cm2")
    return _result.mass_loss_mg_cm2;
  if (quantity == "mass_gain_mg_cm2")
    return _result.mass_gain_mg_cm2;
  if (quantity == "igc_depth_um")
    return _result.igc_depth_um;
  if (quantity == "cr_diffusion_cm2_s")
    return _result.cr_diffusion_cm2_s;
  if (quantity == "fe2_decrease_ppm")
    return _result.fe2_decrease_ppm;
  if (quantity == "redox_log_shift")
    return _result.redox_log_shift;
  if (quantity == "redox_acceleration_ratio")
    return _result.redox_acceleration_ratio;

  static constexpr std::array<const char *, 3> element_names{{"cr", "fe", "ni"}};
  for (unsigned int i = 0; i < element_names.size(); ++i)
  {
    const std::string element = element_names[i];
    if (quantity == "affinity_" + element + "_log_k_over_q")
      return _result.affinity_log_k_over_q[i];
    if (quantity == "source_fraction_" + element)
      return _result.source_fraction[i];
    if (quantity == "saturation_activity_hot_" + element)
      return _result.saturation_activity_hot[i];
    if (quantity == "cold_capture_fraction_" + element)
      return _result.cold_capture_fraction[i];
    if (quantity == "deposit_fraction_" + element)
      return _result.deposit_fraction[i];
    if (quantity == "dissolved_inventory_" + element + "_ppm")
      return _result.dissolved_inventory_ppm[i];
  }

  mooseError("MSTDBTCCorrosionUserObject '",
             name(),
             "' does not provide scalar quantity '",
             quantity,
             "'.");
}
