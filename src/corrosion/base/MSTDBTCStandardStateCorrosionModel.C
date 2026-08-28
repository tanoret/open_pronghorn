//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MSTDBTCStandardStateCorrosionModel.h"

#include "MooseError.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace Corrosion
{
namespace
{
const std::array<MSTDBTCStandardStateCorrosionModel::SpeciesNames, 3> fluoride_species{{
    {"Cr_S1(s)", "CrF2_L1(liq)", "CrF2_P21/c_No.14(s)"},
    {"Fe_bcc(s)", "FeF2_L1(liq)", "FeF2_P42/mnm_No.136(s)"},
    {"Ni_fcc(s)", "NiF2_L1(liq)", "NiF2_P42/mnm_No.136(s)"}}};
const std::array<MSTDBTCStandardStateCorrosionModel::SpeciesNames, 3> chloride_species{{
    {"Cr_S1(s)", "CrCl2_L1(liq)", "CrCl2_Pnnm_No.58(s)"},
    {"Fe_bcc(s)", "FeCl2_L1(liq)", "FeCl2_R3M_No166(s)"},
    {"Ni_Solid_FCC(s)", "NiCl2_L1(liq)", "NiCl2_R3M_No.166(s)"}}};

const std::map<std::string, std::array<Real, 3>> alloy_mass_fractions{
    {"hastelloy_n", {0.070, 0.050, 0.710}},
    {"modified_hastelloy_n", {0.070, 0.050, 0.710}},
    {"gh3535", {0.070, 0.050, 0.700}},
    {"stainless_304", {0.190, 0.690, 0.100}},
    {"stainless_304l", {0.190, 0.690, 0.100}},
    {"stainless_316", {0.170, 0.660, 0.120}},
    {"stainless_316h", {0.170, 0.660, 0.120}},
    {"stainless_316l", {0.170, 0.660, 0.120}},
    {"alloy_x750", {0.155, 0.070, 0.700}},
    {"in625", {0.215, 0.050, 0.610}},
    {"ni_alloy", {0.150, 0.050, 0.700}},
    {"generic_metal", {0.120, 0.500, 0.300}}};

const std::set<std::string> supported_salts{
    "fluoride_fuel", "flibe", "flinak", "fluoroborate", "chloride"};
const std::set<std::string> supported_redox{
    "purified_baseline",
    "msre_or_fuel_baseline",
    "oxidizing_fef2",
    "reducing_be",
    "impure_moisture",
    "chloride_unspecified",
    "tellurium",
    "stressed",
    "multi_alloy",
    "fission_product",
    "gas_control"};

const std::map<std::string, std::string> required_base_redox_parameters{
    {"purified_baseline", "redox_purified_baseline"},
    {"msre_or_fuel_baseline", "redox_msre_or_fuel_baseline"},
    {"oxidizing_fef2", "redox_oxidizing_fef2"},
    {"reducing_be", "redox_reducing_be"},
    {"impure_moisture", "redox_impure_moisture"},
    {"chloride_unspecified", "redox_chloride_unspecified"},
    {"tellurium", "redox_tellurium"},
    {"stressed", "redox_stressed"},
    {"multi_alloy", "redox_multi_alloy"},
    {"fission_product", "redox_fission_product"},
    {"gas_control", "redox_gas_control"}};

Real
clamp(Real value, Real lower, Real upper)
{
  return std::max(lower, std::min(value, upper));
}
} // namespace

MSTDBTCStandardStateCorrosionModel::MSTDBTCStandardStateCorrosionModel(
    const MoltenSaltCorrosionDatabase & base_database,
    const AdvancedCorrosionModelDatabase & advanced_database,
    const MSTDBTCPair & thermodynamics)
  : _base_database(base_database),
    _advanced_database(advanced_database),
    _thermodynamics(thermodynamics),
    _base_model(base_database)
{
  // The core is a provenance boundary too: callers that construct the model directly must not be
  // able to bypass the UserObject/Action checks by passing an unbound MSTDBTCPair.
  if (_thermodynamics.version() != _advanced_database.expectedMSTDBVersion())
    mooseError("MSTDB-TC corrosion: thermodynamic pair version '",
               _thermodynamics.version(),
               "' does not match advanced-model calibration version '",
               _advanced_database.expectedMSTDBVersion(),
               "'.");
  if (_thermodynamics.fluoride().sha256() != _advanced_database.expectedFluorideSHA256())
    mooseError("MSTDB-TC corrosion: fluoride database SHA-256 does not match the advanced-model "
               "calibration provenance.");
  if (_thermodynamics.chloride().sha256() != _advanced_database.expectedChlorideSHA256())
    mooseError("MSTDB-TC corrosion: chloride database SHA-256 does not match the advanced-model "
               "calibration provenance.");

  // The base database class does not expose its file digest. Bind its calibrated subset
  // semantically; the advanced artifact also reports the audited source-file SHA-256.
  _advanced_database.validateBaseModel(_base_database);
  const auto require_base_value = [](const std::string & category,
                                     const std::string & name,
                                     const Real actual,
                                     const Real expected) {
    if (!std::isfinite(actual) || actual != expected)
      mooseError("MSTDB-TC corrosion: base-model ",
                 category,
                 " '",
                 name,
                 "' does not match advanced-model provenance: expected ",
                 expected,
                 ", got ",
                 actual,
                 ".");
  };
  for (const auto & [redox_class, parameter_name] : required_base_redox_parameters)
    require_base_value("redox offset",
                       redox_class,
                       _base_model.redoxOffset(redox_class),
                       _advanced_database.baseModelParameter(parameter_name));

  // Fail during construction rather than deep in the first nonlinear solve.
  for (const auto chloride : {false, true})
    for (unsigned int i = 0; i < n_elements; ++i)
    {
      const auto & names = species(static_cast<Element>(i), chloride);
      const auto & db = chloride ? _thermodynamics.chloride() : _thermodynamics.fluoride();
      db.standardGibbsJMol(names.metal, T_ref_K);
      db.standardGibbsJMol(names.dissolved, T_ref_K);
      db.standardGibbsJMol(names.solid, T_ref_K);
    }
  _thermodynamics.fluoride().standardGibbsJMol("UF3_L1(liq)", T_ref_K);
  _thermodynamics.fluoride().standardGibbsJMol("UF4_L1(liq)", T_ref_K);
  _thermodynamics.fluoride().standardGibbsJMol("Be_S1(s)", T_ref_K);
  _thermodynamics.fluoride().standardGibbsJMol("BeF2_L1(liq)", T_ref_K);
}

const MSTDBTCStandardStateCorrosionModel::SpeciesNames &
MSTDBTCStandardStateCorrosionModel::species(Element element, bool chloride)
{
  const auto index = static_cast<unsigned int>(element);
  if (index >= n_elements)
    mooseError("MSTDB-TC corrosion: unsupported element index ", index, ".");
  return chloride ? chloride_species[index] : fluoride_species[index];
}

const std::array<Real, 3> &
MSTDBTCStandardStateCorrosionModel::alloyMassFractions(const std::string & material_class)
{
  const auto found = alloy_mass_fractions.find(material_class);
  if (found == alloy_mass_fractions.end())
    mooseError("MSTDB-TC corrosion: unsupported material_class '", material_class, "'.");
  return found->second;
}

Real
MSTDBTCStandardStateCorrosionModel::elementMolarMass(Element element)
{
  static constexpr std::array<Real, 3> values{{51.9961, 55.845, 58.6934}};
  return values.at(static_cast<unsigned int>(element));
}

Real
MSTDBTCStandardStateCorrosionModel::halideMolarMass(Element element, bool chloride)
{
  static constexpr std::array<Real, 3> fluoride{{89.9929, 93.8418, 96.6902}};
  static constexpr std::array<Real, 3> chloride_values{{122.9021, 126.7510, 129.5994}};
  return (chloride ? chloride_values : fluoride).at(static_cast<unsigned int>(element));
}

Real
MSTDBTCStandardStateCorrosionModel::meanSaltMolarMass(const std::string & salt_class)
{
  static const std::map<std::string, Real> values{{"fluoride_fuel", 42.0},
                                                  {"flibe", 33.1},
                                                  {"flinak", 41.3},
                                                  {"fluoroborate", 103.0},
                                                  {"chloride", 71.1}};
  const auto found = values.find(salt_class);
  if (found == values.end())
    mooseError("MSTDB-TC corrosion: unsupported salt_class '", salt_class, "'.");
  return found->second;
}

const char *
MSTDBTCStandardStateCorrosionModel::elementName(Element element)
{
  static constexpr std::array<const char *, 3> values{{"Cr", "Fe", "Ni"}};
  return values.at(static_cast<unsigned int>(element));
}

const MSTDBTCData &
MSTDBTCStandardStateCorrosionModel::database(const std::string & salt_class) const
{
  if (!supported_salts.count(salt_class))
    mooseError("MSTDB-TC corrosion: unsupported salt_class '", salt_class, "'.");
  return salt_class == "chloride" ? _thermodynamics.chloride() : _thermodynamics.fluoride();
}

Real
MSTDBTCStandardStateCorrosionModel::parameter(const std::string & name) const
{
  return _advanced_database.thermochemicalParameter(name);
}

Real
MSTDBTCStandardStateCorrosionModel::productActivity(Element element,
                                                    const std::string & salt_class,
                                                    Real ppm) const
{
  if (!std::isfinite(ppm) || ppm < 0.0)
    mooseError("MSTDB-TC corrosion: corrosion-product concentration must be finite and "
               "nonnegative; got ",
               ppm,
               " ppm.");
  const bool chloride = salt_class == "chloride";
  return std::max(ppm * 1.0e-6 * meanSaltMolarMass(salt_class) /
                      halideMolarMass(element, chloride),
                  1.0e-30);
}

Real
MSTDBTCStandardStateCorrosionModel::metalActivity(Element element,
                                                  const std::string & material_class) const
{
  const auto & fractions = alloyMassFractions(material_class);
  std::array<Real, 3> moles{{0.0, 0.0, 0.0}};
  Real tracked = 0.0;
  Real mole_sum = 0.0;
  for (unsigned int i = 0; i < n_elements; ++i)
  {
    tracked += fractions[i];
    moles[i] = fractions[i] / elementMolarMass(static_cast<Element>(i));
    mole_sum += moles[i];
  }
  const Real remainder_moles = std::max(0.0, 1.0 - tracked) / 60.0;
  return std::max(moles.at(static_cast<unsigned int>(element)) /
                      std::max(mole_sum + remainder_moles, 1.0e-30),
                  1.0e-16);
}

Real
MSTDBTCStandardStateCorrosionModel::redoxLogShift(const std::string & redox_class) const
{
  if (!supported_redox.count(redox_class))
    mooseError("MSTDB-TC corrosion: unsupported redox_class '", redox_class, "'.");
  return _base_model.redoxOffset(redox_class) - _base_model.redoxOffset("purified_baseline");
}

Real
MSTDBTCStandardStateCorrosionModel::fe2BufferPpm(const std::string & salt_class,
                                                 const std::string & redox_class) const
{
  if (redox_class == "oxidizing_fef2")
    return 500.0;
  Real baseline = 100.0;
  if (salt_class == "flibe")
    baseline = 30.0;
  else if (salt_class == "chloride")
    baseline = 250.0;
  return clamp(baseline * std::exp(redoxLogShift(redox_class)), 1.0e-3, 1.0e6);
}

Real
MSTDBTCStandardStateCorrosionModel::uraniumRatio(const std::string & redox_class) const
{
  return 100.0 * std::exp(0.5 * redoxLogShift(redox_class));
}

Real
MSTDBTCStandardStateCorrosionModel::reactionLogKOverQ(
    Element element,
    const MSTDBTCCorrosionFeatures & features,
    Real temperature_K,
    const std::string & redox_override,
    Real product_ppm) const
{
  if (!std::isfinite(temperature_K) || temperature_K <= 0.0)
    mooseError("MSTDB-TC corrosion: temperature must be finite and positive; got ",
               temperature_K,
               " K.");
  const auto & salt = features.salt_class;
  const auto redox = redox_override.empty() ? features.redox_class : redox_override;
  if (!supported_redox.count(redox))
    mooseError("MSTDB-TC corrosion: unsupported redox_class '", redox, "'.");
  const bool chloride = salt == "chloride";
  const auto & names = species(element, chloride);
  const Real ppm = product_ppm < 0.0 ? std::exp(parameter("log_product_floor_ppm")) : product_ppm;
  const Real product = productActivity(element, salt, ppm);
  const Real metal = metalActivity(element, features.material_class);
  const auto & db = database(salt);

  if (salt == "fluoride_fuel")
  {
    const std::map<std::string, Real> reaction{{names.dissolved, 1.0},
                                               {"UF3_L1(liq)", 2.0},
                                               {names.metal, -1.0},
                                               {"UF4_L1(liq)", -2.0}};
    const Real ratio = std::max(uraniumRatio(redox), 1.0e-30);
    return db.equilibriumLogConstant(reaction, temperature_K) - std::log(product) +
           std::log(metal) + 2.0 * std::log(ratio);
  }

  if (salt == "flibe" && redox == "reducing_be")
  {
    const std::map<std::string, Real> reaction{{names.dissolved, 1.0},
                                               {"Be_S1(s)", 1.0},
                                               {names.metal, -1.0},
                                               {"BeF2_L1(liq)", -1.0}};
    return db.equilibriumLogConstant(reaction, temperature_K) - std::log(product) +
           std::log(metal);
  }

  const auto & fe = species(Fe, chloride);
  if (element == Fe)
    // Fe + FeX2 <=> FeX2 + Fe is the identity reaction.  Constructing it as a map in the original
    // Python implementation overwrote duplicate keys and introduced a spurious Gibbs term.
    return 0.0;

  const Real fe2 = productActivity(Fe, salt, fe2BufferPpm(salt, redox));
  const std::map<std::string, Real> reaction{{names.dissolved, 1.0},
                                             {fe.metal, 1.0},
                                             {names.metal, -1.0},
                                             {fe.dissolved, -1.0}};
  return db.equilibriumLogConstant(reaction, temperature_K) - std::log(product) -
         std::log(metalActivity(Fe, features.material_class)) + std::log(metal) + std::log(fe2);
}

Real
MSTDBTCStandardStateCorrosionModel::saltActivityCorrection(const std::string & salt_class) const
{
  if (salt_class == "flinak")
    return parameter("log_gamma_cr_flinak");
  if (salt_class == "flibe")
    return parameter("log_gamma_cr_flibe");
  if (salt_class == "fluoroborate")
    return parameter("log_gamma_cr_fluoroborate");
  if (salt_class == "chloride")
    return parameter("log_gamma_cr_chloride");
  return 0.0;
}

Real
MSTDBTCStandardStateCorrosionModel::thermochemicalSaltLogDrive(
    const MSTDBTCCorrosionFeatures & features, Real temperature_K) const
{
  if (features.salt_class == "fluoride_fuel")
    return 0.0;
  auto target = features;
  target.redox_class = "purified_baseline";
  auto reference = target;
  reference.salt_class = "fluoride_fuel";
  return reactionLogKOverQ(Cr, target, temperature_K) -
         reactionLogKOverQ(Cr, reference, temperature_K) -
         saltActivityCorrection(features.salt_class);
}

std::array<Real, 3>
MSTDBTCStandardStateCorrosionModel::speciesFluxFractions(
    const MSTDBTCCorrosionFeatures & features,
    Real temperature_K,
    const std::string & redox_override) const
{
  std::array<Real, 3> affinity;
  std::array<Real, 3> logs;
  const auto & fractions = alloyMassFractions(features.material_class);
  for (unsigned int i = 0; i < n_elements; ++i)
    affinity[i] = reactionLogKOverQ(
        static_cast<Element>(i), features, temperature_K, redox_override);
  const std::array<Real, 3> offsets{{0.0,
                                     parameter("log_exchange_fe_relative"),
                                     parameter("log_exchange_ni_relative")}};
  for (unsigned int i = 0; i < n_elements; ++i)
    logs[i] = std::log(std::max(fractions[i], 1.0e-16)) + offsets[i] +
              parameter("selectivity_affinity_scale") *
                  clamp(affinity[i] - affinity[Cr], -80.0, 80.0);
  const Real maximum = *std::max_element(logs.begin(), logs.end());
  Real total = 0.0;
  std::array<Real, 3> result;
  for (unsigned int i = 0; i < n_elements; ++i)
  {
    result[i] = std::exp(logs[i] - maximum);
    total += result[i];
  }
  for (auto & value : result)
    value /= std::max(total, 1.0e-30);
  return result;
}

Real
MSTDBTCStandardStateCorrosionModel::dissolutionFrontRateUmY(
    const MSTDBTCCorrosionFeatures & features, const std::string & redox_override) const
{
  validate(features);
  const auto redox = redox_override.empty() ? features.redox_class : redox_override;
  const auto & fractions = alloyMassFractions(features.material_class);
  const Real cr_ratio = std::max(fractions[Cr] / 0.07, 1.0e-6);
  const Real thermal = _base_database.parameter("Ea_corr_kJ_mol") * 1000.0 / R_gas *
                       (1.0 / T_ref_K - 1.0 / features.hot_temperature_K);
  const Real log_charge = parameter("log_front_rate0_um_y") + thermal +
                          parameter("cr_activity_exponent") * std::log(cr_ratio) +
                          thermochemicalSaltLogDrive(features, features.hot_temperature_K) +
                          redoxLogShift(redox);
  const Real charge = std::exp(clamp(log_charge, -60.0, 60.0));
  const Real thermal_mt = parameter("Ea_mass_transfer_kJ_mol") * 1000.0 / R_gas *
                          (1.0 / T_ref_K - 1.0 / features.hot_temperature_K);
  const Real log_mt = parameter("log_mass_transfer_cap_um_y") + thermal_mt +
                      parameter("flow_mass_transfer_exponent") * std::log(features.flow_factor);
  const Real mass_transfer = std::exp(clamp(log_mt, -60.0, 60.0));
  return 1.0 / (1.0 / std::max(charge, 1.0e-30) +
                1.0 / std::max(mass_transfer, 1.0e-30));
}

Real
MSTDBTCStandardStateCorrosionModel::massLossFraction(
    const MSTDBTCCorrosionFeatures & features) const
{
  const Real cr_ratio = std::max(alloyMassFractions(features.material_class)[Cr] / 0.07, 1.0e-6);
  const Real logit = parameter("mass_loss_fraction_logit");
  Real value = 1.0 / (1.0 + std::exp(-logit));
  value *= std::pow(cr_ratio, parameter("mass_loss_cr_exponent"));
  if (features.salt_class == "flinak" && cr_ratio > 1.0)
    value *= std::pow(cr_ratio, parameter("flinak_high_cr_selectivity_exponent"));
  return clamp(value, 0.03, 1.0);
}

Real
MSTDBTCStandardStateCorrosionModel::saturationActivity(Element element,
                                                       const std::string & salt_class,
                                                       Real temperature_K) const
{
  if (!std::isfinite(temperature_K) || temperature_K <= 0.0)
    mooseError("MSTDB-TC corrosion: saturation temperature must be finite and positive.");
  const auto & names = species(element, salt_class == "chloride");
  const auto & db = database(salt_class);
  const Real delta_g = db.standardGibbsJMol(names.solid, temperature_K) -
                       db.standardGibbsJMol(names.dissolved, temperature_K);
  return clamp(std::exp(clamp(delta_g / (R_gas * temperature_K), -80.0, 0.0)), 1.0e-30, 1.0);
}

Real
MSTDBTCStandardStateCorrosionModel::coldCaptureFraction(Element element,
                                                        const std::string & salt_class,
                                                        Real hot_temperature_K,
                                                        Real cold_temperature_K) const
{
  if (hot_temperature_K <= cold_temperature_K)
    return 0.0;
  const Real hot = saturationActivity(element, salt_class, hot_temperature_K);
  const Real cold = saturationActivity(element, salt_class, cold_temperature_K);
  return clamp(1.0 - cold / std::max(hot, 1.0e-30), 0.0, 1.0);
}

Real
MSTDBTCStandardStateCorrosionModel::areaToSaltMass(
    const MSTDBTCCorrosionFeatures & features) const
{
  if (!std::isfinite(features.area_to_salt_mass_cm2_g) ||
      features.area_to_salt_mass_cm2_g <= 0.0)
    mooseError(
        "MSTDB-TC corrosion: area_to_salt_mass_cm2_g must be supplied as a finite positive value.");
  return features.area_to_salt_mass_cm2_g;
}

Real
MSTDBTCStandardStateCorrosionModel::diffusionLengthUm(
    const MSTDBTCCorrosionFeatures & features) const
{
  if (features.exposure_s <= 0.0)
    return 0.0;
  CorrosionFeatures base;
  base.temperature_K = features.hot_temperature_K;
  base.material_class = features.material_class;
  base.salt_class = features.salt_class;
  base.redox_class = features.redox_class;
  const Real diffusivity = _base_model.crDiffusionCm2S(base);
  const Real pi = std::acos(-1.0);
  return 2.0 * std::sqrt(std::max(diffusivity, 0.0) * features.exposure_s / pi) * 1.0e4;
}

void
MSTDBTCStandardStateCorrosionModel::validate(const MSTDBTCCorrosionFeatures & features) const
{
  if (!std::isfinite(features.hot_temperature_K) || features.hot_temperature_K <= 0.0)
    mooseError("MSTDB-TC corrosion: hot_temperature_K must be finite and positive.");
  if (!std::isnan(features.cold_temperature_K))
  {
    if (!std::isfinite(features.cold_temperature_K) || features.cold_temperature_K <= 0.0)
      mooseError(
          "MSTDB-TC corrosion: cold_temperature_K must be finite and positive when supplied.");
    if (features.cold_temperature_K > features.hot_temperature_K)
      mooseError("MSTDB-TC corrosion: cold_temperature_K must not exceed hot_temperature_K.");
  }
  if (!std::isfinite(features.exposure_s) || features.exposure_s < 0.0)
    mooseError("MSTDB-TC corrosion: exposure_s must be finite and nonnegative.");
  if (!std::isfinite(features.flow_factor) || features.flow_factor <= 0.0)
    mooseError("MSTDB-TC corrosion: flow_factor must be finite and positive.");
  if (!std::isfinite(features.inventory_coupling_factor) ||
      features.inventory_coupling_factor < 0.0 || features.inventory_coupling_factor > 1.0)
    mooseError("MSTDB-TC corrosion: inventory_coupling_factor must be in [0,1].");
  if (!std::isfinite(features.area_to_salt_mass_cm2_g) ||
      features.area_to_salt_mass_cm2_g <= 0.0)
    mooseError(
        "MSTDB-TC corrosion: area_to_salt_mass_cm2_g must be supplied as a finite positive value.");
  if (!supported_salts.count(features.salt_class))
    mooseError("MSTDB-TC corrosion: unsupported salt_class '", features.salt_class, "'.");
  if (!supported_redox.count(features.redox_class))
    mooseError("MSTDB-TC corrosion: unsupported redox_class '", features.redox_class, "'.");
  alloyMassFractions(features.material_class);
}

MSTDBTCCorrosionResult
MSTDBTCStandardStateCorrosionModel::evaluate(const MSTDBTCCorrosionFeatures & features) const
{
  validate(features);
  MSTDBTCCorrosionResult result;
  const Real time_y = features.exposure_s / seconds_per_year;
  const Real hot = features.hot_temperature_K;
  const Real cold = std::isfinite(features.cold_temperature_K) ? features.cold_temperature_K : hot;
  const bool has_cold_leg = cold < hot;
  const Real density = _base_database.density(features.material_class);

  result.redox_log_shift = redoxLogShift(features.redox_class);
  result.front_rate_um_y = dissolutionFrontRateUmY(features);
  result.redox_acceleration_ratio =
      dissolutionFrontRateUmY(features, "oxidizing_fef2") /
      std::max(dissolutionFrontRateUmY(features, "purified_baseline"), 1.0e-30);
  result.corrosion_rate_um_y = result.front_rate_um_y * massLossFraction(features);
  result.front_depth_um = result.front_rate_um_y * time_y;
  result.mass_loss_mg_cm2 = result.corrosion_rate_um_y * time_y * density * 0.1;
  result.source_fraction = speciesFluxFractions(features, hot);

  Real deposit_weight_sum = 0.0;
  Real weighted_capture = 0.0;
  for (unsigned int i = 0; i < n_elements; ++i)
  {
    const auto element = static_cast<Element>(i);
    result.affinity_log_k_over_q[i] = reactionLogKOverQ(element, features, hot);
    result.saturation_activity_hot[i] = saturationActivity(element, features.salt_class, hot);
    result.cold_capture_fraction[i] =
        has_cold_leg ? coldCaptureFraction(element, features.salt_class, hot, cold) : 0.0;
    result.deposit_fraction[i] = result.source_fraction[i] * result.cold_capture_fraction[i];
    deposit_weight_sum += result.deposit_fraction[i];
    weighted_capture += result.deposit_fraction[i];
  }
  if (deposit_weight_sum > 0.0)
    for (auto & value : result.deposit_fraction)
      value /= deposit_weight_sum;

  const Real capture_area_factor = features.salt_class == "flinak"
                                       ? parameter("deposition_capture_area_factor_flinak")
                                       : parameter("deposition_capture_area_factor_fuel");
  const Real donor_mass_mg_cm2 = result.front_depth_um * density * 0.1;
  result.mass_gain_mg_cm2 = donor_mass_mg_cm2 * capture_area_factor * weighted_capture;

  const Real total_mass_g_cm2 = result.front_depth_um * 1.0e-4 * density;
  const Real area_to_mass = areaToSaltMass(features);
  Real raw_total = 0.0;
  // Model-law invariant: the static closure couples this factor only to dissolved inventory.
  // Coupon mass gain is the independent cold-capture closure above. Changing that distinction
  // requires a new model revision, refit, and validation.
  for (unsigned int i = 0; i < n_elements; ++i)
  {
    const Real retained = has_cold_leg
                              ? std::max(0.0,
                                         1.0 - std::min(1.0,
                                                        capture_area_factor *
                                                            result.cold_capture_fraction[i]))
                              : 1.0;
    result.dissolved_inventory_ppm[i] =
        total_mass_g_cm2 * result.source_fraction[i] * area_to_mass *
        features.inventory_coupling_factor * retained * 1.0e6;
    raw_total += result.dissolved_inventory_ppm[i];
  }
  if (raw_total > 0.0)
  {
    const Real capacity = std::exp(parameter("log_inventory_capacity_ppm"));
    const Real dissolved_total = capacity * (1.0 - std::exp(-raw_total / capacity));
    const Real scale = dissolved_total / raw_total;
    for (auto & value : result.dissolved_inventory_ppm)
      value *= scale;
  }

  const Real cementation_moles =
      std::max(0.0,
               result.dissolved_inventory_ppm[Cr] / elementMolarMass(Cr) +
                   result.dissolved_inventory_ppm[Ni] / elementMolarMass(Ni) -
                   result.dissolved_inventory_ppm[Fe] / elementMolarMass(Fe));
  result.fe2_decrease_ppm = std::min(std::exp(parameter("log_initial_fe2_ppm")),
                                     cementation_moles * elementMolarMass(Fe));

  CorrosionFeatures base;
  base.temperature_K = hot;
  base.material_class = features.material_class;
  base.salt_class = features.salt_class;
  base.redox_class = features.redox_class;
  result.cr_diffusion_cm2_s = _base_model.crDiffusionCm2S(base);
  const Real diffusion = diffusionLengthUm(features);
  Real front_multiplier = 1.0;
  Real gb_multiplier = parameter("gb_diffusion_multiplier");
  if (features.redox_class == "multi_alloy")
  {
    gb_multiplier = parameter("multi_alloy_gb_multiplier");
    if (features.material_class == "alloy_x750")
      gb_multiplier *= parameter("x750_gb_multiplier");
  }
  else if (features.salt_class == "chloride")
  {
    front_multiplier = parameter("chloride_reaction_multiplier");
    gb_multiplier = parameter("chloride_gb_multiplier");
  }
  else if (features.redox_class == "tellurium")
    gb_multiplier = parameter("tellurium_gb_multiplier");
  result.igc_depth_um = front_multiplier * result.front_depth_um + gb_multiplier * diffusion;
  return result;
}

} // namespace Corrosion
