//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MSTDBTCData.h"

#include "MooseError.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace Corrosion
{
namespace
{

std::string
trim(const std::string & value)
{
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
    return "";
  const auto last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::vector<std::string>
tokens(const std::string & line)
{
  std::istringstream stream(line);
  return {std::istream_iterator<std::string>(stream), std::istream_iterator<std::string>()};
}

bool
tryInteger(const std::string & token, long long & value)
{
  try
  {
    std::size_t consumed = 0;
    value = std::stoll(token, &consumed, 10);
    return consumed == token.size();
  }
  catch (const std::exception &)
  {
    return false;
  }
}

long long
strictInteger(const std::string & token, const std::string & context)
{
  long long value = 0;
  if (!tryInteger(token, value))
    mooseError("MSTDB-TC: expected an integer for ", context, ", got '", token, "'.");
  return value;
}

Real
strictReal(std::string token, const std::string & context)
{
  std::replace(token.begin(), token.end(), 'D', 'E');
  std::replace(token.begin(), token.end(), 'd', 'e');
  try
  {
    std::size_t consumed = 0;
    const Real value = std::stod(token, &consumed);
    if (consumed != token.size() || !std::isfinite(value))
      mooseError("MSTDB-TC: expected a finite number for ", context, ", got '", token, "'.");
    return value;
  }
  catch (const std::exception & error)
  {
    mooseError("MSTDB-TC: could not parse '", token, "' for ", context, ": ", error.what());
  }
}

bool
supportedEquationType(const long long equation_type)
{
  return equation_type == 1 || equation_type == 4 || equation_type == 10 ||
         equation_type == 13 || equation_type == 16;
}

bool
looksLikeSpeciesName(const std::string & name)
{
  return !name.empty() && std::isalpha(static_cast<unsigned char>(name.front())) &&
         std::none_of(name.begin(), name.end(), [](const unsigned char character)
         {
           return std::isspace(character);
         });
}

std::string
basename(const std::string & path)
{
  const auto separator = path.find_last_of("/\\");
  return separator == std::string::npos ? path : path.substr(separator + 1);
}

bool
endsWith(const std::string & value, const std::string & suffix)
{
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string
databaseVersion(const std::string & filename)
{
  static const std::regex pattern("V([0-9]+(?:\\.[0-9]+)*)", std::regex::icase);
  std::smatch match;
  return std::regex_search(filename, match, pattern) ? match[1].str() : "";
}

std::string
normalizeHash(const std::string & hash)
{
  if (hash.empty())
    return "";
  if (hash.size() != 64 ||
      !std::all_of(hash.begin(), hash.end(), [](const unsigned char character)
      {
        return std::isxdigit(character);
      }))
    mooseError("MSTDB-TC: expected SHA-256 must contain exactly 64 hexadecimal characters; got '",
               hash,
               "'.");
  std::string result = hash;
  std::transform(result.begin(), result.end(), result.begin(), [](const unsigned char character)
  {
    return static_cast<char>(std::tolower(character));
  });
  return result;
}

std::uint32_t
rotateRight(const std::uint32_t value, const unsigned int shift)
{
  return (value >> shift) | (value << (32U - shift));
}

std::string
sha256File(const std::string & filename)
{
  std::ifstream stream(filename, std::ios::binary);
  if (!stream)
    mooseError("MSTDB-TC: unable to read database file '", filename, "'.");

  std::vector<std::uint8_t> bytes;
  for (std::istreambuf_iterator<char> iterator(stream), end; iterator != end; ++iterator)
    bytes.push_back(static_cast<std::uint8_t>(static_cast<unsigned char>(*iterator)));

  if (bytes.size() > std::numeric_limits<std::uint64_t>::max() / 8U)
    mooseError("MSTDB-TC: database file '", filename, "' is too large to hash.");
  const std::uint64_t bit_length = static_cast<std::uint64_t>(bytes.size()) * 8U;
  bytes.push_back(0x80U);
  while (bytes.size() % 64U != 56U)
    bytes.push_back(0U);
  for (int shift = 56; shift >= 0; shift -= 8)
    bytes.push_back(static_cast<std::uint8_t>((bit_length >> shift) & 0xffU));

  static constexpr std::array<std::uint32_t, 64> constants = {
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U,
      0x923f82a4U, 0xab1c5ed5U, 0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
      0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U, 0xe49b69c1U, 0xefbe4786U,
      0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
      0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U,
      0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
      0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U, 0xa2bfe8a1U, 0xa81a664bU,
      0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
      0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU,
      0x5b9cca4fU, 0x682e6ff3U, 0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
      0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};
  std::array<std::uint32_t, 8> state = {0x6a09e667U,
                                        0xbb67ae85U,
                                        0x3c6ef372U,
                                        0xa54ff53aU,
                                        0x510e527fU,
                                        0x9b05688cU,
                                        0x1f83d9abU,
                                        0x5be0cd19U};

  for (std::size_t offset = 0; offset < bytes.size(); offset += 64U)
  {
    std::array<std::uint32_t, 64> schedule{};
    for (unsigned int index = 0; index < 16; ++index)
    {
      const std::size_t byte = offset + 4U * index;
      schedule[index] = (static_cast<std::uint32_t>(bytes[byte]) << 24U) |
                        (static_cast<std::uint32_t>(bytes[byte + 1]) << 16U) |
                        (static_cast<std::uint32_t>(bytes[byte + 2]) << 8U) |
                        static_cast<std::uint32_t>(bytes[byte + 3]);
    }
    for (unsigned int index = 16; index < 64; ++index)
    {
      const std::uint32_t s0 = rotateRight(schedule[index - 15], 7U) ^
                               rotateRight(schedule[index - 15], 18U) ^
                               (schedule[index - 15] >> 3U);
      const std::uint32_t s1 = rotateRight(schedule[index - 2], 17U) ^
                               rotateRight(schedule[index - 2], 19U) ^
                               (schedule[index - 2] >> 10U);
      schedule[index] = schedule[index - 16] + s0 + schedule[index - 7] + s1;
    }

    std::uint32_t a = state[0];
    std::uint32_t b = state[1];
    std::uint32_t c = state[2];
    std::uint32_t d = state[3];
    std::uint32_t e = state[4];
    std::uint32_t f = state[5];
    std::uint32_t g = state[6];
    std::uint32_t h = state[7];
    for (unsigned int index = 0; index < 64; ++index)
    {
      const std::uint32_t sum1 = rotateRight(e, 6U) ^ rotateRight(e, 11U) ^
                                 rotateRight(e, 25U);
      const std::uint32_t choice = (e & f) ^ ((~e) & g);
      const std::uint32_t temp1 = h + sum1 + choice + constants[index] + schedule[index];
      const std::uint32_t sum0 = rotateRight(a, 2U) ^ rotateRight(a, 13U) ^
                                 rotateRight(a, 22U);
      const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temp2 = sum0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + temp1;
      d = c;
      c = b;
      b = a;
      a = temp1 + temp2;
    }
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    state[5] += f;
    state[6] += g;
    state[7] += h;
  }

  std::ostringstream digest;
  digest << std::hex << std::setfill('0');
  for (const auto value : state)
    digest << std::setw(8) << value;
  return digest.str();
}

std::vector<Real>
numericValues(const std::vector<std::string> & line_tokens,
              const std::size_t first_token,
              const std::string & context)
{
  std::vector<Real> values;
  values.reserve(line_tokens.size() - std::min(first_token, line_tokens.size()));
  for (std::size_t index = first_token; index < line_tokens.size(); ++index)
    values.push_back(strictReal(line_tokens[index], context));
  return values;
}

void
appendNumericLine(const std::vector<std::string> & lines,
                  std::size_t & cursor,
                  std::vector<Real> & values,
                  const std::string & context)
{
  while (cursor < lines.size() && trim(lines[cursor]).empty())
    ++cursor;
  if (cursor >= lines.size())
    mooseError("MSTDB-TC: unexpected end of file while reading ", context, ".");
  const auto line_tokens = tokens(lines[cursor]);
  if (line_tokens.empty())
    mooseError("MSTDB-TC: empty numeric record while reading ", context, ".");
  const auto parsed = numericValues(line_tokens, 0, context);
  values.insert(values.end(), parsed.begin(), parsed.end());
  ++cursor;
}

} // namespace

Real
MSTDBTCGibbsInterval::evaluate(const Real temperature_K) const
{
  if (!(temperature_K > 0.0) || !std::isfinite(temperature_K))
    mooseError("MSTDB-TC: temperature must be finite and positive; got ", temperature_K, " K.");
  const Real T = temperature_K;
  const Real value = coefficients[0] + coefficients[1] * T +
                     coefficients[2] * T * std::log(T) + coefficients[3] * T * T +
                     coefficients[4] * T * T * T + coefficients[5] / T;
  Real result = value;
  for (const auto & term : additional_terms)
    result += term.coefficient *
              (term.exponent == 99.0 ? std::log(T) : std::pow(T, term.exponent));
  if (!std::isfinite(result))
    mooseError("MSTDB-TC: Gibbs function evaluated to a non-finite value at ", T, " K.");
  return result;
}

MSTDBTCData::MSTDBTCData(const std::string & filename,
                         const bool allow_extrapolation,
                         const std::string & expected_sha256,
                         const bool allow_hash_mismatch)
  : _filename(filename),
    _allow_extrapolation(allow_extrapolation),
    _expected_sha256(normalizeHash(expected_sha256)),
    _allow_hash_mismatch(allow_hash_mismatch)
{
  _sha256 = sha256File(_filename);
  _hash_matches_expected = _expected_sha256.empty() || _sha256 == _expected_sha256;
  if (!_hash_matches_expected && !_allow_hash_mismatch)
    mooseError("MSTDB-TC: SHA-256 mismatch for '",
               _filename,
               "': expected ",
               _expected_sha256,
               ", calculated ",
               _sha256,
               ". Set allow_hash_mismatch only for an intentional, documented database change.");
  parse();
}

void
MSTDBTCData::parse()
{
  std::ifstream stream(_filename);
  if (!stream)
    mooseError("MSTDB-TC: unable to open database file '", _filename, "'.");
  std::vector<std::string> lines;
  for (std::string line; std::getline(stream, line);)
  {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    lines.push_back(line);
  }
  if (lines.size() < 3)
    mooseError("MSTDB-TC: '", _filename, "' is too short to be a ChemSage database.");

  _system_line = trim(lines.front());
  if (_system_line.rfind("System", 0) != 0)
    mooseError("MSTDB-TC: first line of '", _filename, "' must begin with 'System'.");
  const auto header_tokens = tokens(lines[1]);
  if (header_tokens.empty())
    mooseError("MSTDB-TC: missing element count on line 2 of '", _filename, "'.");
  const long long n_elements = strictInteger(header_tokens[0], "the ChemSage element count");
  if (n_elements <= 0 || n_elements > 10000)
    mooseError("MSTDB-TC: invalid element count ", n_elements, " in '", _filename, "'.");
  _n_elements = static_cast<unsigned int>(n_elements);
  _version = databaseVersion(basename(_filename));

  std::size_t index = 2;
  while (index + 1 < lines.size())
  {
    const std::string raw_name =
        lines[index].substr(0, std::min<std::size_t>(26, lines[index].size()));
    const std::string name = trim(raw_name);
    if (!looksLikeSpeciesName(name))
    {
      ++index;
      continue;
    }

    const auto equation_tokens = tokens(lines[index + 1]);
    if (equation_tokens.empty())
    {
      ++index;
      continue;
    }
    long long equation_type = 0;
    if (!tryInteger(equation_tokens[0], equation_type))
    {
      ++index;
      continue;
    }
    // SUBL/SUBQ solution-model declarations can contain an endmember label followed by a bare
    // two-integer sublattice record (for example, "F" followed by "1 1"). A function-expanded
    // standard-state species record always begins its stoichiometry on the equation line.
    if (supportedEquationType(equation_type) && equation_tokens.size() < 3)
    {
      ++index;
      continue;
    }
    if (!supportedEquationType(equation_type))
    {
      // SUBQ is a solution-model declaration rather than a species record.
      if (name == "SUBQ")
      {
        ++index;
        continue;
      }
      mooseError("MSTDB-TC: unsupported equation type ",
                 equation_type,
                 " following candidate species '",
                 name,
                 "' on line ",
                 index + 1,
                 " of '",
                 _filename,
                 "'. Supported types are 1, 4, 10, 13, and 16.");
    }
    const long long n_intervals =
        strictInteger(equation_tokens[1], "the interval count for species '" + name + "'");
    if (n_intervals < 1 || n_intervals > 50)
      mooseError("MSTDB-TC: invalid interval count ",
                 n_intervals,
                 " for species '",
                 name,
                 "' on line ",
                 index + 1,
                 ".");

    MSTDBTCSpecies record;
    record.name = name;
    record.equation_type = static_cast<unsigned int>(equation_type);
    record.source_line = static_cast<unsigned int>(index + 1);

    std::vector<Real> stoichiometry =
        numericValues(equation_tokens, 2, "stoichiometry for species '" + name + "'");
    std::size_t cursor = index + 2;
    while (stoichiometry.size() < _n_elements)
      appendNumericLine(lines,
                        cursor,
                        stoichiometry,
                        "stoichiometry for species '" + name + "'");
    if (stoichiometry.size() != _n_elements)
      mooseError("MSTDB-TC: species '",
                 name,
                 "' has ",
                 stoichiometry.size(),
                 " stoichiometric coefficients; expected exactly ",
                 _n_elements,
                 ".");
    record.stoichiometry = std::move(stoichiometry);

    Real previous_upper_temperature = 0.0;
    for (long long interval_index = 0; interval_index < n_intervals; ++interval_index)
    {
      std::vector<Real> interval_values;
      while (interval_values.size() < 7)
        appendNumericLine(lines,
                          cursor,
                          interval_values,
                          "Gibbs interval for species '" + name + "'");
      if (interval_values.size() != 7)
        mooseError("MSTDB-TC: Gibbs interval ",
                   interval_index,
                   " for species '",
                   name,
                   "' contains ",
                   interval_values.size(),
                   " values; expected exactly 7.");
      if (!(interval_values[0] > previous_upper_temperature))
        mooseError("MSTDB-TC: interval upper temperatures for species '",
                   name,
                   "' must be finite, positive, and strictly increasing; got ",
                   interval_values[0],
                   " K after ",
                   previous_upper_temperature,
                   " K.");

      MSTDBTCGibbsInterval interval;
      interval.upper_temperature_K = interval_values[0];
      std::copy_n(interval_values.begin() + 1, 6, interval.coefficients.begin());
      previous_upper_temperature = interval.upper_temperature_K;

      if (equation_type == 4 || equation_type == 16)
      {
        while (cursor < lines.size() && trim(lines[cursor]).empty())
          ++cursor;
        if (cursor >= lines.size())
          mooseError("MSTDB-TC: missing additional-term record for species '", name, "'.");
        const auto extra_tokens = tokens(lines[cursor]);
        if (extra_tokens.empty())
          mooseError("MSTDB-TC: empty additional-term record for species '", name, "'.");
        const long long n_extra =
            strictInteger(extra_tokens[0], "the additional-term count for species '" + name + "'");
        if (n_extra < 0 || n_extra > 1000)
          mooseError("MSTDB-TC: invalid additional-term count ",
                     n_extra,
                     " for species '",
                     name,
                     "'.");
        std::vector<Real> extra_values =
            numericValues(extra_tokens, 1, "additional terms for species '" + name + "'");
        ++cursor;
        const std::size_t expected_values = 2U * static_cast<std::size_t>(n_extra);
        while (extra_values.size() < expected_values)
          appendNumericLine(lines,
                            cursor,
                            extra_values,
                            "additional terms for species '" + name + "'");
        if (extra_values.size() != expected_values)
          mooseError("MSTDB-TC: additional-term record for species '",
                     name,
                     "' contains ",
                     extra_values.size(),
                     " values; expected exactly ",
                     expected_values,
                     ".");
        for (std::size_t term = 0; term < static_cast<std::size_t>(n_extra); ++term)
          interval.additional_terms.push_back(
              {extra_values[2U * term], extra_values[2U * term + 1U]});
      }
      record.intervals.push_back(std::move(interval));
    }

    if (equation_type == 13 || equation_type == 16)
    {
      while (cursor < lines.size() && trim(lines[cursor]).empty())
        ++cursor;
      if (cursor >= lines.size())
        mooseError("MSTDB-TC: missing SGTE magnetic parameters for species '", name, "'.");
      const auto magnetic_tokens = tokens(lines[cursor]);
      record.magnetic_parameters =
          numericValues(magnetic_tokens, 0, "SGTE magnetic parameters for species '" + name + "'");
      ++cursor;
      if (record.magnetic_parameters.size() < 4)
        mooseError("MSTDB-TC: species '",
                   name,
                   "' has ",
                   record.magnetic_parameters.size(),
                   " SGTE magnetic parameters; expected at least 4.");
    }

    const std::size_t record_index = _records.size();
    _records.push_back(std::move(record));
    _record_indices[name].push_back(record_index);
    index = cursor;
  }

  if (_records.empty())
    mooseError("MSTDB-TC: no standard-state species records were parsed from '", _filename, "'.");
}

const MSTDBTCSpecies &
MSTDBTCData::species(const std::string & name, const std::ptrdiff_t occurrence) const
{
  const auto found = _record_indices.find(name);
  if (found == _record_indices.end())
    mooseError("MSTDB-TC: exact species '", name, "' was not found in '", _filename, "'.");
  const auto & indices = found->second;
  std::ptrdiff_t selected = occurrence;
  if (selected == -1)
    selected = static_cast<std::ptrdiff_t>(indices.size()) - 1;
  if (selected < 0 || static_cast<std::size_t>(selected) >= indices.size())
    mooseError("MSTDB-TC: occurrence ",
               occurrence,
               " is invalid for species '",
               name,
               "', which has ",
               indices.size(),
               " occurrence(s).");
  return _records[indices[static_cast<std::size_t>(selected)]];
}

std::size_t
MSTDBTCData::occurrenceCount(const std::string & name) const
{
  const auto found = _record_indices.find(name);
  return found == _record_indices.end() ? 0U : found->second.size();
}

Real
MSTDBTCData::evaluateSpecies(const MSTDBTCSpecies & record,
                             const Real temperature_K,
                             const bool allow_extrapolation) const
{
  if (!(temperature_K > 0.0) || !std::isfinite(temperature_K))
    mooseError("MSTDB-TC: temperature must be finite and positive; got ", temperature_K, " K.");
  const auto interval = std::find_if(record.intervals.begin(),
                                     record.intervals.end(),
                                     [temperature_K](const MSTDBTCGibbsInterval & candidate)
                                     {
                                       return temperature_K <= candidate.upper_temperature_K;
                                     });
  if (interval == record.intervals.end() && !allow_extrapolation)
    mooseError("MSTDB-TC: temperature ",
               temperature_K,
               " K exceeds the last standard-state interval (",
               record.intervals.back().upper_temperature_K,
               " K) for species '",
               record.name,
               "'. Set allow_extrapolation=true only with a documented basis.");
  const auto & selected = interval == record.intervals.end() ? record.intervals.back() : *interval;
  Real value = selected.evaluate(temperature_K);
  if (record.magnetic_parameters.size() >= 4)
    value += magneticGibbsJMol(temperature_K,
                               record.magnetic_parameters[0],
                               record.magnetic_parameters[1],
                               record.magnetic_parameters[2],
                               record.magnetic_parameters[3]);
  if (!std::isfinite(value))
    mooseError("MSTDB-TC: standard Gibbs energy for species '",
               record.name,
               "' is non-finite at ",
               temperature_K,
               " K.");
  return value;
}

Real
MSTDBTCData::standardGibbsJMol(const std::string & species_name,
                               const Real temperature_K) const
{
  return standardGibbsJMol(species_name, temperature_K, -1, _allow_extrapolation);
}

Real
MSTDBTCData::standardGibbsJMol(const std::string & species_name,
                               const Real temperature_K,
                               const bool allow_extrapolation) const
{
  return standardGibbsJMol(species_name, temperature_K, -1, allow_extrapolation);
}

Real
MSTDBTCData::standardGibbsJMol(const std::string & species_name,
                               const Real temperature_K,
                               const std::ptrdiff_t occurrence,
                               const bool allow_extrapolation) const
{
  return evaluateSpecies(species(species_name, occurrence), temperature_K, allow_extrapolation);
}

Real
MSTDBTCData::reactionGibbsJMol(const std::map<std::string, Real> & reaction,
                               const Real temperature_K) const
{
  return reactionGibbsJMol(reaction, temperature_K, _allow_extrapolation);
}

Real
MSTDBTCData::reactionGibbsJMol(const std::map<std::string, Real> & reaction,
                               const Real temperature_K,
                               const bool allow_extrapolation) const
{
  if (reaction.empty())
    mooseError("MSTDB-TC: a reaction must contain at least one species.");
  Real result = 0.0;
  for (const auto & [name, coefficient] : reaction)
  {
    if (!std::isfinite(coefficient))
      mooseError("MSTDB-TC: non-finite reaction coefficient for species '", name, "'.");
    result += coefficient * standardGibbsJMol(name, temperature_K, allow_extrapolation);
  }
  if (!std::isfinite(result))
    mooseError("MSTDB-TC: reaction Gibbs energy is non-finite at ", temperature_K, " K.");
  return result;
}

Real
MSTDBTCData::equilibriumLogConstant(const std::map<std::string, Real> & reaction,
                                    const Real temperature_K) const
{
  return equilibriumLogConstant(reaction, temperature_K, _allow_extrapolation);
}

Real
MSTDBTCData::equilibriumLogConstant(const std::map<std::string, Real> & reaction,
                                    const Real temperature_K,
                                    const bool allow_extrapolation) const
{
  if (!(temperature_K > 0.0) || !std::isfinite(temperature_K))
    mooseError("MSTDB-TC: temperature must be finite and positive; got ", temperature_K, " K.");
  return -reactionGibbsJMol(reaction, temperature_K, allow_extrapolation) /
         (gas_constant_J_mol_K * temperature_K);
}

Real
MSTDBTCData::magneticGibbsJMol(const Real temperature_K,
                               const Real critical_temperature_K,
                               const Real average_magnetic_moment,
                               const Real structure_factor,
                               const Real p)
{
  // Adapted from Thermochimica src/gem/CompGibbsMagnetic.f90 at
  // 0c35c8d7d1cf2084b4e2ca5d6608f7dcdf60adad; see THIRD_PARTY_NOTICES.md.
  if (!(temperature_K > 0.0) || !std::isfinite(temperature_K))
    mooseError("MSTDB-TC: magnetic Gibbs temperature must be finite and positive; got ",
               temperature_K,
               " K.");
  if (!std::isfinite(critical_temperature_K) || !std::isfinite(average_magnetic_moment) ||
      !std::isfinite(structure_factor) || !std::isfinite(p))
    mooseError("MSTDB-TC: SGTE magnetic parameters must be finite.");
  if (critical_temperature_K == 0.0 || p == 0.0)
    return 0.0;

  Real critical = critical_temperature_K;
  Real moment = average_magnetic_moment;
  if (critical < 0.0)
  {
    critical = -critical * structure_factor;
    moment = -moment * structure_factor;
  }
  if (!(critical > 0.0) || !(moment + 1.0 > 0.0))
    return 0.0;

  const Real tau = temperature_K / critical;
  const Real inv_p_minus_one = 1.0 / p - 1.0;
  const Real denominator = 518.0 / 1125.0 + 11692.0 / 15975.0 * inv_p_minus_one;
  if (denominator == 0.0 || !std::isfinite(denominator))
    mooseError("MSTDB-TC: invalid SGTE magnetic denominator for p=", p, ".");

  Real g_tau = 0.0;
  if (tau > 1.0)
    g_tau = -(std::pow(tau, -5.0) / 10.0 + std::pow(tau, -15.0) / 315.0 +
              std::pow(tau, -25.0) / 1500.0) /
            denominator;
  else
    g_tau =
        1.0 -
        (79.0 / (140.0 * p * tau) +
         474.0 / 497.0 * inv_p_minus_one *
             (std::pow(tau, 3.0) / 6.0 + std::pow(tau, 9.0) / 135.0 +
              std::pow(tau, 15.0) / 600.0)) /
            denominator;
  const Real result =
      gas_constant_J_mol_K * temperature_K * std::log(moment + 1.0) * g_tau;
  if (!std::isfinite(result))
    mooseError("MSTDB-TC: SGTE magnetic Gibbs contribution is non-finite.");
  return result;
}

nlohmann::json
MSTDBTCData::provenance() const
{
  return {{"filename", basename(_filename)},
          {"version", _version},
          {"sha256", _sha256},
          {"expected_sha256", _expected_sha256.empty() ? nlohmann::json(nullptr)
                                                        : nlohmann::json(_expected_sha256)},
          {"sha256_matches_expected", _hash_matches_expected},
          {"system", _system_line},
          {"n_elements", _n_elements},
          {"n_species_records", _records.size()},
          {"allow_hash_mismatch", _allow_hash_mismatch},
          {"allow_extrapolation", _allow_extrapolation}};
}

MSTDBTCPair::MSTDBTCPair(const std::string & fluoride_file,
                         const std::string & chloride_file,
                         const std::string & expected_version,
                         const std::string & expected_fluoride_sha256,
                         const std::string & expected_chloride_sha256,
                         const bool allow_hash_mismatch,
                         const bool allow_extrapolation)
  : _fluoride(fluoride_file,
              allow_extrapolation,
              expected_fluoride_sha256,
              allow_hash_mismatch),
    _chloride(chloride_file,
              allow_extrapolation,
              expected_chloride_sha256,
              allow_hash_mismatch),
    _expected_version(expected_version)
{
  if (!endsWith(basename(fluoride_file), "Fluorides_No_Func.dat"))
    mooseError("MSTDB-TC: fluoride database filename must end in 'Fluorides_No_Func.dat'; got '",
               basename(fluoride_file),
               "'.");
  if (!endsWith(basename(chloride_file), "Chlorides_No_Func.dat"))
    mooseError("MSTDB-TC: chloride database filename must end in 'Chlorides_No_Func.dat'; got '",
               basename(chloride_file),
               "'.");
  if (_fluoride.version().empty() || _chloride.version().empty())
    mooseError("MSTDB-TC: both database filenames must contain an explicit V<version> token.");
  if (_fluoride.version() != _chloride.version())
    mooseError("MSTDB-TC: fluoride version ",
               _fluoride.version(),
               " does not match chloride version ",
               _chloride.version(),
               ".");
  _version = _fluoride.version();
  if (!_expected_version.empty() && _version != _expected_version)
    mooseError("MSTDB-TC: paired database version ",
               _version,
               " does not match calibrated expected version ",
               _expected_version,
               ".");
}

nlohmann::json
MSTDBTCPair::provenance() const
{
  return {{"version", _version},
          {"expected_version", _expected_version.empty() ? nlohmann::json(nullptr)
                                                          : nlohmann::json(_expected_version)},
          {"fluoride", _fluoride.provenance()},
          {"chloride", _chloride.provenance()}};
}

} // namespace Corrosion
