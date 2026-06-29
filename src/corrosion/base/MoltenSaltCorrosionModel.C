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

#include "MoltenSaltCorrosionModel.h"
#include "MooseError.h"

#include <algorithm>
#include <cmath>

namespace Corrosion
{

namespace
{
// Reproduce Python's "value or fallback" idiom: a numeric value of exactly zero (falsy) selects the
// fallback. The reference model relies on this for the temperature and flow defaults.
Real
pyOr(Real value, Real fallback)
{
  return (value == 0.0) ? fallback : value;
}
}

MoltenSaltCorrosionModel::MoltenSaltCorrosionModel(const MoltenSaltCorrosionDatabase & db) : _db(db)
{
}

Real
MoltenSaltCorrosionModel::materialOffset(const std::string & material_class) const
{
  // Hastelloy N is the reference alloy (zero offset). Stainless and graphite map onto shared keys.
  if (material_class == "hastelloy_n")
    return 0.0;
  if (material_class == "modified_hastelloy_n")
    return param("mat_modified_hastelloy_n");
  if (material_class == "gh3535")
    return param("mat_gh3535");
  if (material_class == "stainless_304" || material_class == "stainless_304l")
    return param("mat_stainless_304");
  if (material_class == "stainless_316" || material_class == "stainless_316h" ||
      material_class == "stainless_316l")
    return param("mat_stainless_316");
  if (material_class == "alloy_x750")
    return param("mat_alloy_x750");
  if (material_class == "in625")
    return param("mat_in625");
  if (material_class == "ni_alloy")
    return param("mat_ni_alloy");
  // generic_metal and graphite both fall onto the generic-metal offset.
  return param("mat_generic_metal");
}

Real
MoltenSaltCorrosionModel::saltOffset(const std::string & salt_class) const
{
  // Fuel fluoride is the reference salt (zero offset); unknown salts also default to zero.
  if (salt_class == "flinak")
    return param("salt_flinak");
  if (salt_class == "flibe")
    return param("salt_flibe");
  if (salt_class == "fluoroborate")
    return param("salt_fluoroborate");
  if (salt_class == "chloride")
    return param("salt_chloride");
  if (salt_class == "generic_salt")
    return param("salt_generic_salt");
  if (salt_class == "no_salt")
    return param("salt_no_salt");
  return 0.0;
}

Real
MoltenSaltCorrosionModel::redoxOffset(const std::string & redox_class) const
{
  if (redox_class == "msre_or_fuel_baseline")
    return param("redox_msre_or_fuel_baseline");
  if (redox_class == "oxidizing_fef2")
    return param("redox_oxidizing_fef2");
  if (redox_class == "reducing_be")
    return param("redox_reducing_be");
  if (redox_class == "impure_moisture")
    return param("redox_impure_moisture");
  if (redox_class == "chloride_unspecified")
    return param("redox_chloride_unspecified");
  if (redox_class == "tellurium")
    return param("redox_tellurium");
  if (redox_class == "stressed")
    return param("redox_stressed");
  if (redox_class == "multi_alloy")
    return param("redox_multi_alloy");
  if (redox_class == "fission_product")
    return param("redox_fission_product");
  if (redox_class == "gas_control")
    return param("redox_gas_control");
  // purified_baseline and any unrecognized class default to the purified baseline offset.
  return param("redox_purified_baseline");
}

Real
MoltenSaltCorrosionModel::positionOffset(const std::string & position_class) const
{
  if (position_class == "hot_leg")
    return param("hot_leg_bonus");
  if (position_class == "cold_leg")
    return param("cold_leg_corrosion_penalty");
  return 0.0;
}

Real
MoltenSaltCorrosionModel::depositionSaltOffset(const std::string & salt_class) const
{
  if (salt_class == "flinak")
    return param("dep_salt_flinak");
  if (salt_class == "flibe")
    return param("dep_salt_flibe");
  if (salt_class == "fluoroborate")
    return param("dep_salt_fluoroborate");
  if (salt_class == "chloride")
    return param("dep_salt_chloride");
  return 0.0;
}

Real
MoltenSaltCorrosionModel::thermalTerm(Real temperature_K, Real Ea_kJ_mol) const
{
  if (!std::isfinite(temperature_K) || temperature_K <= 0.0)
    temperature_K = T_ref_K;
  return (Ea_kJ_mol * 1000.0 / R_gas) * (1.0 / T_ref_K - 1.0 / temperature_K);
}

Real
MoltenSaltCorrosionModel::bvOverpotentialEquivalentV(const std::string & redox_class,
                                                     Real temperature_K,
                                                     Real alpha_n) const
{
  return redoxOffset(redox_class) * R_gas * temperature_K / (alpha_n * faraday);
}

Real
MoltenSaltCorrosionModel::corrosionRateUmY(const CorrosionFeatures & feat,
                                           const std::string & redox_override) const
{
  const Real T = pyOr(feat.temperature_K, T_ref_K);
  const std::string & material = feat.material_class;
  const std::string & salt = feat.salt_class;
  const std::string redox = redox_override.empty() ? feat.redox_class : redox_override;
  const std::string & position = feat.position_class;
  const Real flow = std::max(pyOr(feat.flow_factor, 0.75), 1.0e-3);
  const Real dT = std::max(feat.delta_T_C, 0.0);

  if (salt == "no_salt" || flow < 0.05)
    // Keep finite for log plots while indicating negligible molten-salt corrosion.
    return 1.0e-6;

  const Real log_kin = param("log_rate0_um_y") + thermalTerm(T, param("Ea_corr_kJ_mol")) +
                       materialOffset(material) + saltOffset(salt) + redoxOffset(redox) +
                       positionOffset(position) + param("gamma_flow_corr") * std::log(flow) +
                       param("theta_dT_corr") * std::log1p(dT / 100.0);
  const Real rate_kin = expClip(log_kin);

  const Real log_mt = param("log_mt_cap_um_y") + thermalTerm(T, param("Ea_mt_kJ_mol")) +
                      0.35 * saltOffset(salt) + param("gamma_flow_mt") * std::log(flow) +
                      param("theta_dT_mt") * std::log1p(dT / 100.0);
  const Real rate_mt = expClip(log_mt);

  // Harmonic mean: kinetic-limited at low current, transport-limited at high driving force.
  return 1.0 / (1.0 / std::max(rate_kin, 1.0e-12) + 1.0 / std::max(rate_mt, 1.0e-12));
}

Real
MoltenSaltCorrosionModel::damageMultiplier(const CorrosionFeatures & feat) const
{
  const std::string & redox = feat.redox_class;
  const std::string & salt = feat.salt_class;
  Real value = param("damage_base_log");
  if (salt == "chloride" || redox == "impure_moisture" || redox == "chloride_unspecified")
    value += param("damage_chloride_log");
  if (redox == "tellurium")
    value += param("damage_tellurium_log");
  if (redox == "stressed")
    value += param("damage_stress_log");
  if (redox == "multi_alloy")
    value += param("damage_multi_alloy_log");
  return expClip(value, -10.0, 10.0);
}

Real
MoltenSaltCorrosionModel::corrosionDepthUm(const CorrosionFeatures & feat) const
{
  Real time_y = (feat.time_years == 0.0) ? std::numeric_limits<Real>::quiet_NaN() : feat.time_years;
  if (!std::isfinite(time_y))
    return std::numeric_limits<Real>::quiet_NaN();
  return corrosionRateUmY(feat) * std::max(time_y, 0.0);
}

Real
MoltenSaltCorrosionModel::igcDepthUm(const CorrosionFeatures & feat) const
{
  Real time_y = (feat.time_years == 0.0) ? std::numeric_limits<Real>::quiet_NaN() : feat.time_years;
  if (!std::isfinite(time_y))
    return std::numeric_limits<Real>::quiet_NaN();
  const Real uniform = corrosionRateUmY(feat) * std::max(time_y, 0.0);
  // Mixed linear/parabolic morphology: the square-root term allows deep IGC/void penetration with
  // limited mass loss while remaining tied to the Butler-Volmer corrosion drive.
  const Real dmg = damageMultiplier(feat);
  return uniform * dmg + std::sqrt(std::max(uniform, 0.0)) * dmg;
}

Real
MoltenSaltCorrosionModel::massLossMgCm2(const CorrosionFeatures & feat) const
{
  return umToMgCm2(corrosionDepthUm(feat), _db.density(feat.material_class));
}

Real
MoltenSaltCorrosionModel::depositionRateUmY(const CorrosionFeatures & feat,
                                            const std::string & surface_override) const
{
  const Real T = pyOr(feat.temperature_K, T_ref_K);
  const std::string & salt = feat.salt_class;
  const std::string & redox = feat.redox_class;
  const std::string & position = feat.position_class;
  const std::string surface = surface_override.empty() ? feat.surface_class : surface_override;
  const Real flow = std::max(pyOr(feat.flow_factor, 0.75), 1.0e-3);
  const Real dT = std::max(feat.delta_T_C, 0.0);
  if (salt == "no_salt")
    return 1.0e-9;

  Real surface_offset = 0.0;
  if (surface == "graphite")
    surface_offset += param("dep_graphite_offset");
  else if (surface == "turbulent_metal")
    surface_offset += param("dep_turbulent_bonus");
  else if (surface == "laminar_metal")
    surface_offset += param("dep_laminar_bonus");
  if (redox == "multi_alloy")
    surface_offset += param("dep_multi_alloy_bonus");

  // Cathodic/plating branch: the redox drive enters with a reduced sign because noble metals and
  // corrosion products deposit as cathodic/precipitation processes.
  const Real redox_cathodic = -0.30 * redoxOffset(redox);
  const Real cold_bonus = (position == "cold_leg") ? param("dep_cold_bonus") : 0.0;
  const Real log_dep = param("log_dep0_um_y") + thermalTerm(T, param("Ea_dep_kJ_mol")) +
                       0.20 * saltOffset(salt) + depositionSaltOffset(salt) + redox_cathodic +
                       param("gamma_flow_dep") * std::log(flow) +
                       param("theta_dT_dep") * std::log1p(dT / 100.0) + cold_bonus + surface_offset;
  return expClip(log_dep);
}

Real
MoltenSaltCorrosionModel::depositionDepthUm(const CorrosionFeatures & feat) const
{
  Real time_y = (feat.time_years == 0.0) ? std::numeric_limits<Real>::quiet_NaN() : feat.time_years;
  if (!std::isfinite(time_y))
    return std::numeric_limits<Real>::quiet_NaN();
  return depositionRateUmY(feat) * std::max(time_y, 0.0);
}

Real
MoltenSaltCorrosionModel::massGainMgCm2(const CorrosionFeatures & feat) const
{
  return umToMgCm2(depositionDepthUm(feat), _db.density(feat.material_class));
}

DepositionRanking
MoltenSaltCorrosionModel::depositionRanking(const CorrosionFeatures & feat) const
{
  DepositionRanking r;
  r.graphite = depositionRateUmY(feat, "graphite");
  r.laminar_metal = depositionRateUmY(feat, "laminar_metal");
  r.metal = depositionRateUmY(feat, "metal");
  r.turbulent_metal = depositionRateUmY(feat, "turbulent_metal");
  return r;
}

Real
MoltenSaltCorrosionModel::redoxAccelerationRatio(const CorrosionFeatures & feat) const
{
  const Real oxidized = corrosionRateUmY(feat, "oxidizing_fef2");
  const Real baseline = corrosionRateUmY(feat, "purified_baseline");
  return oxidized / std::max(baseline, 1.0e-12);
}

Real
MoltenSaltCorrosionModel::crDiffusionCm2S(const CorrosionFeatures & feat) const
{
  const Real T = pyOr(feat.temperature_K, T_ref_K);
  const Real logD = param("log_Dcr_ref_cm2_s") + thermalTerm(T, param("Ea_Dcr_kJ_mol"));
  return expClip(logD, -80.0, -20.0);
}

Real
MoltenSaltCorrosionModel::saltCrPpmBase(const CorrosionFeatures & feat) const
{
  const Real depth = std::max(corrosionDepthUm(feat), 1.0e-9);
  const Real cr_factor = _db.crWeightFraction(feat.material_class) / 0.07;
  const bool msre = feat.experiment_family.find("MSRE") != std::string::npos ||
                    feat.source_id.rfind("ORNL-TM-3", 0) == 0;
  const Real scale =
      msre ? expClip(param("log_ppm_scale_msre")) : expClip(param("log_ppm_scale_loop"));
  return scale * depth * cr_factor;
}

Real
MoltenSaltCorrosionModel::saltCrPpm(const CorrosionFeatures & feat) const
{
  Real ppm = saltCrPpmBase(feat);

  // NCL-16 reports +500 ppm Cr after the 29,500 h baseline exposure. Treat this as a
  // source-specific Cr inventory correction; the generic scale already matches the Fe decrease.
  if (feat.source_id == "ORNL-TM-4188" && feat.salt_class == "fluoride_fuel" &&
      feat.redox_class == "purified_baseline")
    ppm *= expClip(param("log_ncl16_cr_inventory_bonus"));

  return ppm;
}

Real
MoltenSaltCorrosionModel::saltFeDecreasePpm(const CorrosionFeatures & feat) const
{
  // Do not propagate the NCL-16 Cr inventory correction to Fe.
  const Real cr_ppm = (feat.source_id == "ORNL-TM-4188" && feat.salt_class == "fluoride_fuel" &&
                       feat.redox_class == "purified_baseline")
                          ? saltCrPpmBase(feat)
                          : saltCrPpm(feat);
  return cr_ppm * expClip(param("log_fe_to_cr_ppm_ratio"));
}

Real
MoltenSaltCorrosionModel::offgasFractionPercent() const
{
  const Real logit = param("logit_offgas_fraction");
  const Real frac = 1.0 / (1.0 + std::exp(-logit));
  return 100.0 * frac;
}

Real
MoltenSaltCorrosionModel::teSolublePpm() const
{
  return expClip(param("log_te_soluble_ppm"));
}

Real
MoltenSaltCorrosionModel::teRedoxThresholdRatio() const
{
  return expClip(param("log_te_threshold_ratio"));
}

Real
MoltenSaltCorrosionModel::predictResponse(const CorrosionFeatures & feat) const
{
  const std::string & kind = feat.response_kind;
  if (kind == "corrosion_rate_um_y")
    return corrosionRateUmY(feat);
  if (kind == "corrosion_depth_um")
    return corrosionDepthUm(feat);
  if (kind == "igc_depth_um")
    return igcDepthUm(feat);
  if (kind == "mass_loss_mg_cm2")
    return massLossMgCm2(feat);
  if (kind == "mass_gain_mg_cm2")
    return massGainMgCm2(feat);
  if (kind == "redox_acceleration_ratio" || kind == "redox_acceleration_qualitative")
    return redoxAccelerationRatio(feat);
  if (kind == "salt_cr_ppm")
    return saltCrPpm(feat);
  if (kind == "salt_fe_decrease_ppm")
    return saltFeDecreasePpm(feat);
  if (kind == "cr_diffusion_cm2_s")
    return crDiffusionCm2S(feat);
  if (kind == "offgas_fraction_percent")
    return offgasFractionPercent();
  if (kind == "te_soluble_ppm")
    return teSolublePpm();
  if (kind == "te_redox_threshold_ratio")
    return teRedoxThresholdRatio();
  if (kind == "noble_metal_deposition_ranking")
  {
    // Return the turbulent-metal/graphite ratio as a scalar diagnostic.
    const Real turbulent = depositionRateUmY(feat, "turbulent_metal");
    const Real graphite = depositionRateUmY(feat, "graphite");
    return turbulent / std::max(graphite, 1.0e-12);
  }
  return std::numeric_limits<Real>::quiet_NaN();
}

Real
MoltenSaltCorrosionModel::crDiffusivityM2S(const CorrosionFeatures & feat) const
{
  // cr_diffusion_cm2_s is in cm^2/s; convert to m^2/s (1 cm^2 = 1e-4 m^2).
  return crDiffusionCm2S(feat) * 1.0e-4;
}

Real
MoltenSaltCorrosionModel::exchangeCurrentDensity(const CorrosionFeatures & feat,
                                                 const std::string & element) const
{
  // Seed the mechanistic exchange current from the calibrated baseline dissolution rate by the
  // Faradaic conversion: i0 [A/m^2] = 1e4 * i0 [A/cm^2]. The element supplies valence and molar
  // mass; the material class supplies density.
  const ElementProperties props = _db.element(element);
  const Real rate_um_y = corrosionRateUmY(feat);
  const Real i0_a_cm2 =
      umYToCorrosionCurrent(rate_um_y, props.valence, props.molar_mass_g_mol, _db.density(feat.material_class));
  return i0_a_cm2 * 1.0e4;
}

std::vector<LoopState>
MoltenSaltCorrosionModel::simulateLoop(const std::vector<CorrosionFeatures> & segments,
                                       Real duration_h,
                                       Real dt_h,
                                       Real salt_volume_cm3,
                                       Real initial_cr_ppm) const
{
  if (duration_h <= 0.0)
    mooseError("Corrosion: simulateLoop requires duration_h > 0.");
  if (dt_h <= 0.0)
    mooseError("Corrosion: simulateLoop requires dt_h > 0.");

  const std::size_t n_seg = segments.size();
  const int n_steps = static_cast<int>(std::ceil(duration_h / dt_h));
  const Real salt_mass_mg = salt_volume_cm3 * _db.saltDensity() * 1000.0;

  Real cr_ppm = initial_cr_ppm;
  std::vector<Real> cum_diss_mg(n_seg, 0.0);
  std::vector<Real> cum_dep_mg(n_seg, 0.0);
  std::vector<LoopState> rows;

  for (int step = 0; step <= n_steps; ++step)
  {
    const Real time_h = std::min(static_cast<Real>(step) * dt_h, duration_h);
    Real total_diss = 0.0;
    Real total_dep = 0.0;
    for (std::size_t i = 0; i < n_seg; ++i)
    {
      total_diss += cum_diss_mg[i];
      total_dep += cum_dep_mg[i];
    }
    rows.push_back({time_h, cr_ppm, total_diss, total_dep});
    if (step == n_steps)
      break;

    const Real actual_dt_h = std::min(dt_h, duration_h - time_h);
    const Real dt_y = actual_dt_h / hours_per_year;
    for (std::size_t i = 0; i < n_seg; ++i)
    {
      const CorrosionFeatures & seg = segments[i];
      const Real area = (seg.surface_area_cm2 == 0.0) ? 1.0 : seg.surface_area_cm2;
      const Real rho = _db.density(seg.material_class);

      const Real corr_depth = corrosionRateUmY(seg) * dt_y;
      const Real total_mg = umToMgCm2(corr_depth, rho) * area;
      const Real cr_mg = total_mg * _db.crWeightFraction(seg.material_class);
      cum_diss_mg[i] += cr_mg;
      cr_ppm += cr_mg / std::max(salt_mass_mg, 1.0) * 1.0e6;

      Real dep_depth = depositionRateUmY(seg) * dt_y;
      Real dep_mg = umToMgCm2(dep_depth, rho) * area;
      // Deposition cannot exceed a relaxed available inventory in this simplified closure.
      const Real available_mg = std::max(cr_ppm, 0.0) / 1.0e6 * salt_mass_mg;
      dep_mg = std::min(dep_mg, 0.25 * available_mg);
      cum_dep_mg[i] += dep_mg;
      cr_ppm -= dep_mg / std::max(salt_mass_mg, 1.0) * 1.0e6;
      cr_ppm = std::max(cr_ppm, 0.0);
    }
  }
  return rows;
}

} // namespace Corrosion
