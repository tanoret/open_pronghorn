//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "AdvancedCorrosionModelPostprocessor.h"

#include "AdvancedCorrosionModelUserObject.h"

registerMooseObject("OpenPronghornApp", AdvancedCorrosionModelPostprocessor);

InputParameters
AdvancedCorrosionModelPostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Reports one selected scalar result from an advanced corrosion-model UserObject.");
  params.addRequiredParam<UserObjectName>("model", "Advanced corrosion-model UserObject.");
  MooseEnum quantity(
      "front_rate_um_y corrosion_rate_um_y front_depth_um mass_loss_mg_cm2 mass_gain_mg_cm2 "
      "igc_depth_um cr_diffusion_cm2_s fe2_decrease_ppm redox_log_shift "
      "redox_acceleration_ratio "
      "affinity_cr_log_k_over_q affinity_fe_log_k_over_q affinity_ni_log_k_over_q "
      "source_fraction_cr source_fraction_fe source_fraction_ni "
      "saturation_activity_hot_cr saturation_activity_hot_fe saturation_activity_hot_ni "
      "cold_capture_fraction_cr cold_capture_fraction_fe cold_capture_fraction_ni "
      "deposit_fraction_cr deposit_fraction_fe deposit_fraction_ni "
      "dissolved_inventory_cr_ppm dissolved_inventory_fe_ppm dissolved_inventory_ni_ppm "
      "mass_recession_um average_corrosion_rate_um_y instantaneous_front_rate_um_y redox_shift "
      "mass_balance_relative_error elapsed_time_y accepted_internal_steps rejected_internal_steps "
      "cumulative_source_cr_ppm cumulative_source_fe_ppm cumulative_source_ni_ppm "
      "coupon_deposit_cr_mg_cm2 coupon_deposit_fe_mg_cm2 coupon_deposit_ni_mg_cm2 "
      "bulk_captured_cr_ppm bulk_captured_fe_ppm bulk_captured_ni_ppm "
      "surface_availability_cr surface_availability_fe surface_availability_ni");
  params.addRequiredParam<MooseEnum>("quantity", quantity, "Scalar model result to report.");
  return params;
}

AdvancedCorrosionModelPostprocessor::AdvancedCorrosionModelPostprocessor(
    const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _model(getUserObject<AdvancedCorrosionModelUserObject>("model")),
    _quantity(getParam<MooseEnum>("quantity"))
{
}

void
AdvancedCorrosionModelPostprocessor::execute()
{
  _value = _model.scalarValue(_quantity);
}
