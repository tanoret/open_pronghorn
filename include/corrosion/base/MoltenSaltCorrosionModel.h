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

#pragma once

#include "MoltenSaltCorrosionData.h"
#include "CorrosionChemistry.h"

#include <limits>
#include <string>
#include <vector>

namespace Corrosion
{

/**
 * The feature set of a single corrosion/plating case, mirroring one row of the reference model's
 * case_features.csv / targets.csv. Defaults match the per-feature defaults used by
 * msr_corrosion_bv/model.py so that an unset field reproduces the reference behaviour.
 */
struct CorrosionFeatures
{
  /// Operating temperature [K].
  Real temperature_K = T_ref_K;
  /// Alloy class (e.g. "hastelloy_n", "stainless_316").
  std::string material_class = "generic_metal";
  /// Salt class (e.g. "fluoride_fuel", "flinak", "chloride").
  std::string salt_class = "generic_salt";
  /// Redox/effective-overpotential class (e.g. "purified_baseline", "oxidizing_fef2").
  std::string redox_class = "purified_baseline";
  /// Loop position class ("hot_leg", "cold_leg" or nominal).
  std::string position_class = "nominal";
  /// Deposition surface class ("metal", "graphite", "turbulent_metal", "laminar_metal").
  std::string surface_class = "metal";
  /// Circulation/mass-transfer factor (dimensionless).
  Real flow_factor = 0.75;
  /// Loop thermal gradient [C].
  Real delta_T_C = 0.0;
  /// Exposure time [years]; not-a-number marks an unset time (rate-only response).
  Real time_years = std::numeric_limits<Real>::quiet_NaN();
  /// Experiment family string (used to select the Cr inventory scale).
  std::string experiment_family = "";
  /// Source identifier string (used to select the Cr inventory scale).
  std::string source_id = "";
  /// The response the case reports (selects the predict_response branch).
  std::string response_kind = "corrosion_rate_um_y";
  /// Wetted area of a loop segment [cm^2] (loop simulation only).
  Real surface_area_cm2 = 1.0;
};

/// Plating rates by deposition surface, returned by depositionRanking.
struct DepositionRanking
{
  Real graphite = 0.0;
  Real laminar_metal = 0.0;
  Real metal = 0.0;
  Real turbulent_metal = 0.0;
};

/// One recorded state of the simplified multi-segment loop simulation.
struct LoopState
{
  Real time_h = 0.0;
  Real salt_cr_ppm = 0.0;
  Real total_cr_dissolved_mg = 0.0;
  Real total_deposit_mg = 0.0;
};

/**
 * Effective Butler-Volmer corrosion, mass-transfer and plating correlation, a faithful C++ port of
 * msr_corrosion_bv/model.py::MoltenSaltBVModel. The calibrated parameters and the engineering
 * material tables are read from MoltenSaltCorrosionDatabase, so this class reproduces the reference
 * predictions term for term while also supplying the per-element exchange current, effective
 * overpotential and solid-state chromium diffusivity that parameterize the mechanistic spatial
 * Butler-Volmer objects.
 */
class MoltenSaltCorrosionModel
{
public:
  /// Construct from a loaded database. The database must outlive the model.
  MoltenSaltCorrosionModel(const MoltenSaltCorrosionDatabase & db);

  // --- Log-rate offsets (port of the model.py *_offset helpers) -------------------------------
  Real materialOffset(const std::string & material_class) const;
  Real saltOffset(const std::string & salt_class) const;
  Real redoxOffset(const std::string & redox_class) const;
  Real positionOffset(const std::string & position_class) const;
  Real depositionSaltOffset(const std::string & salt_class) const;

  /// Arrhenius temperature term (Ea/R)(1/Tref - 1/T) with Ea in kJ/mol.
  Real thermalTerm(Real temperature_K, Real Ea_kJ_mol) const;

  /// Effective overpotential [V] implied by the fitted redox offset (offset = alpha*n*F*eta/RT).
  Real bvOverpotentialEquivalentV(const std::string & redox_class,
                                  Real temperature_K,
                                  Real alpha_n = 1.0) const;

  // --- Corrosion branch -----------------------------------------------------------------------
  /// Anodic dissolution rate [um/y] (harmonic mean of the kinetic and transport branches).
  /// An empty redox_override uses the feature's redox class.
  Real corrosionRateUmY(const CorrosionFeatures & feat,
                        const std::string & redox_override = "") const;
  Real damageMultiplier(const CorrosionFeatures & feat) const;
  Real corrosionDepthUm(const CorrosionFeatures & feat) const;
  Real igcDepthUm(const CorrosionFeatures & feat) const;
  Real massLossMgCm2(const CorrosionFeatures & feat) const;

  // --- Deposition / plating branch ------------------------------------------------------------
  Real depositionRateUmY(const CorrosionFeatures & feat,
                         const std::string & surface_override = "") const;
  Real depositionDepthUm(const CorrosionFeatures & feat) const;
  Real massGainMgCm2(const CorrosionFeatures & feat) const;
  DepositionRanking depositionRanking(const CorrosionFeatures & feat) const;

  // --- Auxiliary submodels --------------------------------------------------------------------
  Real redoxAccelerationRatio(const CorrosionFeatures & feat) const;
  Real crDiffusionCm2S(const CorrosionFeatures & feat) const;
  Real saltCrPpmBase(const CorrosionFeatures & feat) const;
  Real saltCrPpm(const CorrosionFeatures & feat) const;
  Real saltFeDecreasePpm(const CorrosionFeatures & feat) const;
  Real offgasFractionPercent() const;
  Real teSolublePpm() const;
  Real teRedoxThresholdRatio() const;

  /// Dispatch to the predictor for feat.response_kind (port of model.py::predict_response).
  Real predictResponse(const CorrosionFeatures & feat) const;

  // --- Faradaic bridge to the mechanistic Butler-Volmer ---------------------------------------
  /// Solid-state chromium diffusivity [m^2/s] (cr_diffusion_cm2_s converted from cm^2/s).
  Real crDiffusivityM2S(const CorrosionFeatures & feat) const;
  /// Exchange current density [A/m^2] for a tracked element, seeded from the baseline corrosion
  /// rate via the Faradaic conversion. element gives valence/molar mass; material gives density.
  Real exchangeCurrentDensity(const CorrosionFeatures & feat, const std::string & element) const;

  /// Simulate the simplified multi-segment loop (port of model.py::simulate_loop).
  std::vector<LoopState> simulateLoop(const std::vector<CorrosionFeatures> & segments,
                                      Real duration_h,
                                      Real dt_h = 24.0,
                                      Real salt_volume_cm3 = 1000.0,
                                      Real initial_cr_ppm = 0.0) const;

protected:
  /// One calibrated parameter by name.
  Real param(const std::string & name) const { return _db.parameter(name); }

  const MoltenSaltCorrosionDatabase & _db;
};

} // namespace Corrosion
