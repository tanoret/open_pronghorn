//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "MooseTypes.h"
#include "nlohmann/json.h"

#include <map>
#include <string>

namespace Corrosion
{

class MoltenSaltCorrosionDatabase;

/**
 * Calibrated closure parameters and provenance for the MSTDB-TC and DRIDN model layers.
 *
 * Thermodynamic standard-state data are deliberately not stored here. They are loaded from an
 * authorized, version-matched pair of MSTDB-TC *_No_Func.dat files at runtime. This JSON database
 * contains fitted corrosion/transport closures, solver defaults, and semantic provenance bindings.
 */
class AdvancedCorrosionModelDatabase
{
public:
  explicit AdvancedCorrosionModelDatabase(const std::string & filename);

  /// Exact model-law revision understood by this implementation.
  static const std::string & supportedModelRevision();

  /// A calibrated MSTDB-TC thermochemical-model parameter.
  Real thermochemicalParameter(const std::string & name) const;

  /// A calibrated DRIDN parameter.
  Real dynamicParameter(const std::string & name) const;

  const std::map<std::string, Real> & thermochemicalParameters() const
  {
    return _thermochemical_parameters;
  }
  const std::map<std::string, Real> & dynamicParameters() const { return _dynamic_parameters; }

  /// Schema identifier recorded in the JSON input.
  const std::string & schemaVersion() const { return _schema_version; }

  /// Human-readable calibration identifier recorded in the JSON input.
  const std::string & calibrationId() const { return _calibration_id; }

  /// Model-law revision, verified exactly against supportedModelRevision() while loading.
  const std::string & modelRevision() const { return _model_revision; }

  /// Required MSTDB-TC edition and raw-file hashes bound to this calibration.
  const std::string & expectedMSTDBVersion() const { return _expected_mstdb_version; }
  const std::string & expectedFluorideSHA256() const { return _expected_fluoride_sha256; }
  const std::string & expectedChlorideSHA256() const { return _expected_chloride_sha256; }

  /// Provenance and exact semantic subset of the base corrosion database used by the static model.
  const std::string & baseModelSourceSHA256() const { return _base_model_source_sha256; }
  Real baseModelParameter(const std::string & name) const;
  Real baseModelDensity(const std::string & material_class) const;
  Real baseModelElementProperty(const std::string & element, const std::string & property) const;
  /// Reject any base database whose bound parameters, densities, or Cr/Fe/Ni properties differ.
  void validateBaseModel(const MoltenSaltCorrosionDatabase & database) const;
  const std::map<std::string, Real> & baseModelParameters() const
  {
    return _base_model_parameters;
  }
  const std::map<std::string, Real> & baseModelDensities() const
  {
    return _base_model_densities;
  }
  const std::map<std::string, std::map<std::string, Real>> & baseModelElementProperties() const
  {
    return _base_model_element_properties;
  }

  /// Commit/revision of the authoritative normalized validation data.
  const std::string & calibrationDataRevision() const { return _calibration_data_revision; }

  /// Full parsed JSON document, used by validation/reporting objects.
  const nlohmann::json & json() const { return _root; }

protected:
  Real lookup(const std::map<std::string, Real> & values,
              const std::string & name,
              const std::string & section) const;
  void validateRequiredParameters() const;

  const std::string _filename;
  nlohmann::json _root;
  std::string _schema_version;
  std::string _calibration_id;
  std::string _model_revision;
  std::string _expected_mstdb_version;
  std::string _expected_fluoride_sha256;
  std::string _expected_chloride_sha256;
  std::string _base_model_source_sha256;
  std::string _calibration_data_revision;
  std::map<std::string, Real> _base_model_parameters;
  std::map<std::string, Real> _base_model_densities;
  std::map<std::string, std::map<std::string, Real>> _base_model_element_properties;
  std::map<std::string, Real> _thermochemical_parameters;
  std::map<std::string, Real> _dynamic_parameters;
};

} // namespace Corrosion
