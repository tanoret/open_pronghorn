//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "AdvancedCorrosionModelData.h"
#include "MSTDBTCData.h"
#include "MoltenSaltCorrosionModel.h"

#include <array>
#include <limits>
#include <string>

namespace Corrosion
{

/** Inputs to the MSTDB-TC standard-state/Nernst reduced corrosion model.
 *
 * The model is a lumped engineering correlation.  MSTDB supplies standard-state Gibbs energies;
 * fitted closures represent activities, charge transfer, transport, geometry, and cold-leg capture.
 * It is deliberately not described as a native Thermochimica/SUBQ equilibrium calculation.
 */
struct MSTDBTCCorrosionFeatures
{
  Real hot_temperature_K = T_ref_K;
  Real cold_temperature_K = std::numeric_limits<Real>::quiet_NaN();
  Real exposure_s = std::numeric_limits<Real>::quiet_NaN();
  Real flow_factor = 0.75;
  Real area_to_salt_mass_cm2_g = std::numeric_limits<Real>::quiet_NaN();
  // Fraction applied only to dissolved salt inventory; it does not scale static mass gain.
  Real inventory_coupling_factor = 1.0;
  // For Python calibration-context reproduction, inventory_source_material maps to this field;
  // the C++ endpoint model uses material_class as its donor/source alloy throughout.
  std::string material_class = "hastelloy_n";
  std::string salt_class = "fluoride_fuel";
  std::string redox_class = "purified_baseline";
};

/** Complete scalar endpoint and diagnostic output of one thermochemical corrosion evaluation. */
struct MSTDBTCCorrosionResult
{
  Real front_rate_um_y = 0.0;
  Real corrosion_rate_um_y = 0.0;
  Real front_depth_um = 0.0;
  Real mass_loss_mg_cm2 = 0.0;
  Real mass_gain_mg_cm2 = 0.0;
  Real igc_depth_um = 0.0;
  Real cr_diffusion_cm2_s = 0.0;
  Real fe2_decrease_ppm = 0.0;
  Real redox_log_shift = 0.0;
  Real redox_acceleration_ratio = 1.0;
  std::array<Real, 3> affinity_log_k_over_q{{0.0, 0.0, 0.0}};
  std::array<Real, 3> source_fraction{{0.0, 0.0, 0.0}};
  std::array<Real, 3> saturation_activity_hot{{0.0, 0.0, 0.0}};
  std::array<Real, 3> cold_capture_fraction{{0.0, 0.0, 0.0}};
  std::array<Real, 3> deposit_fraction{{0.0, 0.0, 0.0}};
  std::array<Real, 3> dissolved_inventory_ppm{{0.0, 0.0, 0.0}};
};

/**
 * Reduced Cr/Fe/Ni corrosion model driven by MSTDB-TC standard-state Gibbs functions.
 *
 * This class corrects the historical non-fuel Fe-buffer duplicate-key defect: the standard-state Fe
 * + FeX2 identity reaction has exactly zero net stoichiometry and affinity.  Unknown material,
 * salt, and redox labels are rejected rather than silently mapped to a generic case.
 */
class MSTDBTCStandardStateCorrosionModel
{
public:
  static constexpr unsigned int n_elements = 3;
  enum Element : unsigned int
  {
    Cr = 0,
    Fe = 1,
    Ni = 2
  };

  struct SpeciesNames
  {
    const char * metal;
    const char * dissolved;
    const char * solid;
  };

  MSTDBTCStandardStateCorrosionModel(const MoltenSaltCorrosionDatabase & base_database,
                                     const AdvancedCorrosionModelDatabase & advanced_database,
                                     const MSTDBTCPair & thermodynamics);

  MSTDBTCCorrosionResult evaluate(const MSTDBTCCorrosionFeatures & features) const;

  Real reactionLogKOverQ(Element element,
                         const MSTDBTCCorrosionFeatures & features,
                         Real temperature_K,
                         const std::string & redox_override = "",
                         Real product_ppm = -1.0) const;

  std::array<Real, 3> speciesFluxFractions(const MSTDBTCCorrosionFeatures & features,
                                           Real temperature_K,
                                           const std::string & redox_override = "") const;

  Real saturationActivity(Element element,
                          const std::string & salt_class,
                          Real temperature_K) const;
  Real coldCaptureFraction(Element element,
                           const std::string & salt_class,
                           Real hot_temperature_K,
                           Real cold_temperature_K) const;
  Real dissolutionFrontRateUmY(const MSTDBTCCorrosionFeatures & features,
                               const std::string & redox_override = "") const;
  Real massLossFraction(const MSTDBTCCorrosionFeatures & features) const;
  Real redoxLogShift(const std::string & redox_class) const;

protected:
  static const SpeciesNames & species(Element element, bool chloride);
  static const std::array<Real, 3> & alloyMassFractions(const std::string & material_class);
  static Real elementMolarMass(Element element);
  static Real halideMolarMass(Element element, bool chloride);
  static Real meanSaltMolarMass(const std::string & salt_class);
  static const char * elementName(Element element);

  const MSTDBTCData & database(const std::string & salt_class) const;
  Real parameter(const std::string & name) const;
  Real productActivity(Element element, const std::string & salt_class, Real ppm) const;
  Real metalActivity(Element element, const std::string & material_class) const;
  Real fe2BufferPpm(const std::string & salt_class, const std::string & redox_class) const;
  Real uraniumRatio(const std::string & redox_class) const;
  Real saltActivityCorrection(const std::string & salt_class) const;
  Real thermochemicalSaltLogDrive(const MSTDBTCCorrosionFeatures & features,
                                  Real temperature_K) const;
  Real areaToSaltMass(const MSTDBTCCorrosionFeatures & features) const;
  Real diffusionLengthUm(const MSTDBTCCorrosionFeatures & features) const;
  void validate(const MSTDBTCCorrosionFeatures & features) const;

  const MoltenSaltCorrosionDatabase & _base_database;
  const AdvancedCorrosionModelDatabase & _advanced_database;
  const MSTDBTCPair & _thermodynamics;
  MoltenSaltCorrosionModel _base_model;
};

} // namespace Corrosion
