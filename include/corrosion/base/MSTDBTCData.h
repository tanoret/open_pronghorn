//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "MooseTypes.h"
#include "nlohmann/json.h"

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Corrosion
{

/** One additional term in a ChemSage standard-state Gibbs function. */
struct MSTDBTCExtraTerm
{
  Real coefficient;
  Real exponent;
};

/** One temperature interval in a function-expanded ChemSage species record. */
struct MSTDBTCGibbsInterval
{
  Real upper_temperature_K;
  std::array<Real, 6> coefficients;
  std::vector<MSTDBTCExtraTerm> additional_terms;

  Real evaluate(Real temperature_K) const;
};

/** A parsed ChemSage standard-state species record. */
struct MSTDBTCSpecies
{
  std::string name;
  unsigned int equation_type;
  std::vector<Real> stoichiometry;
  std::vector<MSTDBTCGibbsInterval> intervals;
  std::vector<Real> magnetic_parameters;
  unsigned int source_line;
};

/**
 * Strict reader for function-expanded MSTDB-TC/ChemSage standard-state data.
 *
 * This class evaluates standard-state species Gibbs functions. It deliberately does not parse
 * ChemSage solution-model declarations, perform Gibbs-energy minimization, or evaluate SUBQ
 * activities. Database files remain external runtime inputs and are never distributed by this
 * class.
 */
class MSTDBTCData
{
public:
  static constexpr Real gas_constant_J_mol_K = 8.31446261815324;

  MSTDBTCData(const std::string & filename,
              bool allow_extrapolation = false,
              const std::string & expected_sha256 = "",
              bool allow_hash_mismatch = false);

  /// Evaluate G^0 in J/mol using the database-wide extrapolation policy.
  Real standardGibbsJMol(const std::string & species_name, Real temperature_K) const;

  /// Evaluate G^0 in J/mol with an explicit per-call extrapolation policy.
  Real standardGibbsJMol(const std::string & species_name,
                         Real temperature_K,
                         bool allow_extrapolation) const;

  /// Evaluate a particular zero-based occurrence; -1 selects the last occurrence.
  Real standardGibbsJMol(const std::string & species_name,
                         Real temperature_K,
                         std::ptrdiff_t occurrence,
                         bool allow_extrapolation) const;

  /// Return sum(nu_i G_i^0), with products carrying positive coefficients.
  Real reactionGibbsJMol(const std::map<std::string, Real> & reaction,
                         Real temperature_K) const;
  Real reactionGibbsJMol(const std::map<std::string, Real> & reaction,
                         Real temperature_K,
                         bool allow_extrapolation) const;

  /// Return natural-log K = -Delta G^0/(R T).
  Real equilibriumLogConstant(const std::map<std::string, Real> & reaction,
                              Real temperature_K) const;
  Real equilibriumLogConstant(const std::map<std::string, Real> & reaction,
                              Real temperature_K,
                              bool allow_extrapolation) const;

  /// Exact, case-sensitive species lookup. -1 selects the last occurrence.
  const MSTDBTCSpecies & species(const std::string & name,
                                 std::ptrdiff_t occurrence = -1) const;
  std::size_t occurrenceCount(const std::string & name) const;

  const std::string & filename() const { return _filename; }
  const std::string & version() const { return _version; }
  const std::string & sha256() const { return _sha256; }
  const std::string & expectedSha256() const { return _expected_sha256; }
  const std::string & systemLine() const { return _system_line; }
  unsigned int elementCount() const { return _n_elements; }
  std::size_t recordCount() const { return _records.size(); }
  bool allowExtrapolation() const { return _allow_extrapolation; }
  bool hashMatchesExpected() const { return _hash_matches_expected; }

  /// Machine-readable database provenance for result files and reports.
  nlohmann::json provenance() const;

  /// Public for independent unit testing of the SGTE expression.
  static Real magneticGibbsJMol(Real temperature_K,
                                Real critical_temperature_K,
                                Real average_magnetic_moment,
                                Real structure_factor,
                                Real p);

protected:
  void parse();
  Real evaluateSpecies(const MSTDBTCSpecies & record,
                       Real temperature_K,
                       bool allow_extrapolation) const;

  const std::string _filename;
  const bool _allow_extrapolation;
  const std::string _expected_sha256;
  const bool _allow_hash_mismatch;

  std::string _sha256;
  bool _hash_matches_expected = true;
  std::string _version;
  std::string _system_line;
  unsigned int _n_elements = 0;
  std::vector<MSTDBTCSpecies> _records;
  std::map<std::string, std::vector<std::size_t>> _record_indices;
};

/** A version- and provenance-checked fluoride/chloride MSTDB-TC pair. */
class MSTDBTCPair
{
public:
  MSTDBTCPair(const std::string & fluoride_file,
              const std::string & chloride_file,
              const std::string & expected_version = "",
              const std::string & expected_fluoride_sha256 = "",
              const std::string & expected_chloride_sha256 = "",
              bool allow_hash_mismatch = false,
              bool allow_extrapolation = false);

  const MSTDBTCData & fluoride() const { return _fluoride; }
  const MSTDBTCData & chloride() const { return _chloride; }
  const std::string & version() const { return _version; }

  nlohmann::json provenance() const;

protected:
  MSTDBTCData _fluoride;
  MSTDBTCData _chloride;
  std::string _version;
  std::string _expected_version;
};

} // namespace Corrosion
