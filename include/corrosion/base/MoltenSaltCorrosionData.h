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

#include "MooseTypes.h"
#include "json.h"

#include <map>
#include <string>

namespace Corrosion
{

/// Mechanistic transport and electrode-kinetics properties of one metal cation, read from the
/// 'elements' section of the corrosion database.
struct ElementProperties
{
  /// Charge number z of the dissolved cation M^{z+}.
  Real valence = 2.0;
  /// Molar mass of the metal [g/mol].
  Real molar_mass_g_mol = 55.0;
  /// Molecular diffusivity of the cation in the salt [m^2/s].
  Real diffusivity_m2_s = 1.0e-9;
  /// Standard electrode potential E0 [V].
  Real E0_V = 0.0;
  /// Anodic charge-transfer coefficient.
  Real alpha_a = 0.5;
  /// Cathodic charge-transfer coefficient.
  Real alpha_c = 0.5;
  /// Reference concentration for the Nernst activity term [mol/m^3].
  Real c_ref_mol_m3 = 1.0;
};

/// Solid-state Arrhenius diffusivity properties for one element in one alloy.

struct SolidDiffusivityProperties
{
  /// Arrhenius pre-exponential factor [cm^2/s].
  Real D0_cm2_s = 0.0;

  /// Arrhenius activation energy [kJ/mol].
  Real Q_kJ_mol = 0.0;

  /// Minimum temperature represented by the experimental measurements [K].
  Real measurement_temperature_min_K = 0.0;

  /// Maximum temperature represented by the experimental measurements [K].
  Real measurement_temperature_max_K = 0.0;

  /// Human-readable source description.
  std::string source;

  /// Provenance/status classification, e.g. "direct_tracer_measurement".
  std::string status;
};

/**
 * Reads and queries the molten salt corrosion and plating database from a JSON file.
 *
 * The database bundles three kinds of information:
 *   - 'elements': the per-cation mechanistic properties (valence, molar mass, diffusivity, E0,
 *     transfer coefficients, reference concentration) consumed by the Nernst-Planck and
 *     Butler-Volmer objects.
 *   - 'densities_g_cm3' and 'alloy_cr_wt_frac': the engineering material tables from the reference
 *     model (msr_corrosion_bv/chemistry.py) used by the Faradaic conversions and the ported
 *     correlation.
 *   - 'calibrated_parameters': the verbatim fitted parameter set of the reference effective
 *     Butler-Volmer correlation (results/parameters.json), keyed exactly as in that model so the
 *     C++ port reproduces it term by term.
 *
 * The default file is data/corrosion_database.json; users can supply their own to change the
 * elements, materials or parameters without recompiling.
 */
class MoltenSaltCorrosionDatabase
{
public:
  /// Load and parse the database from a JSON file.
  MoltenSaltCorrosionDatabase(const std::string & filename);

  /// Mechanistic properties of a tracked element (case-insensitive); falls back to 'generic_metal'.
  ElementProperties element(const std::string & name) const;

  /// True when the database lists the element explicitly (no fall-back applied).
  bool hasElement(const std::string & name) const;

  /// Metal density [g/cm^3] for a material class; falls back to 'generic_metal'
  /// (port of chemistry.py::density).
  Real density(const std::string & material_class) const;

  /// Chromium mass fraction of a material class; falls back to 'generic_metal'
  /// (port of chemistry.py::cr_weight_fraction).
  Real crWeightFraction(const std::string & material_class) const;

/// True when a material/element pair has an exact or explicitly configured
  /// solid-state diffusivity. The generic_metal legacy fallback is not applied here.
  bool hasSolidDiffusivity(const std::string & material_class,
                          const std::string & element) const;

  /// Solid-state Arrhenius diffusivity properties for an exact or explicitly configured
  /// fallback material/element pair. For configured fallbacks, the returned provenance
  /// describes the source-material correlation. Errors if neither is present.
  SolidDiffusivityProperties
  solidDiffusivity(const std::string & material_class,
                  const std::string & element) const;

  /// Effective salt mass density [g/cm^3] used by the loop inventory closure.
  Real saltDensity() const;

  /// One calibrated correlation parameter by name (e.g. "log_rate0_um_y"); errors if absent.
  Real parameter(const std::string & name) const;

  /// The complete calibrated parameter set, keyed as in results/parameters.json.
  const std::map<std::string, Real> & calibratedParameters() const { return _parameters; }

protected:
  /// Look up a numeric member of an object with a 'generic_metal' fall-back key.
  Real lookupWithFallback(const nlohmann::json & object,
                          const std::string & name,
                          const std::string & fallback_key,
                          const std::string & what) const;

  /// Find an object member whose key matches name case-insensitively; returns nullptr if absent.
  const nlohmann::json * findCaseInsensitive(const nlohmann::json & object,
                                             const std::string & name) const;

  /// File the database was read from (for error messages).
  const std::string _filename;

  /// Parsed JSON document.
  nlohmann::json _root;

  /// Flattened calibrated correlation parameters.
  std::map<std::string, Real> _parameters;
};

} // namespace Corrosion
