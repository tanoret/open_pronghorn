//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MoltenSaltRadiolysisData.h"
#include "MooseError.h"
#include "MooseUtils.h"

#include <cmath>
#include <cctype>
#include <fstream>

namespace MSR
{

namespace
{
// Case-insensitive string comparison so callers can pass salt/metal/gas tokens in any case
// (MultiMooseEnum, for example, upper-cases its values).
bool
iequals(const std::string & a, const std::string & b)
{
  if (a.size() != b.size())
    return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (std::tolower(static_cast<unsigned char>(a[i])) !=
        std::tolower(static_cast<unsigned char>(b[i])))
      return false;
  return true;
}

// Convert a JSON object of {species: stoichiometric_coefficient} into SpeciesTerm entries, applying
// any kinetic-order overrides from the optional 'orders' object (orders default to the coefficient).
std::vector<SpeciesTerm>
toSpeciesTerms(const nlohmann::json & terms, const nlohmann::json & orders)
{
  std::vector<SpeciesTerm> result;
  for (const auto & item : terms.items())
  {
    const Real coeff = item.value().get<Real>();
    Real order = coeff;
    if (!orders.is_null() && orders.contains(item.key()))
      order = orders.at(item.key()).get<Real>();
    result.push_back({item.key(), coeff, order});
  }
  return result;
}
}

Real
gToSource(Real g_value, Real dose_rate)
{
  // S = G * dose_rate * 1 / (100 eV) * 1 / N_A
  const Real factor = 1.0 / (100.0 * eV_to_J) / avogadro;
  return g_value * dose_rate * factor;
}

Real
ReactionData::rate(Real temperature) const
{
  // Arrhenius forms mirror msr_radiolysis/core.py::Reaction.k_forward
  if (A > 0.0 && Ea > 0.0)
    return A * std::exp(-Ea / (R_gas * temperature));
  if (k_ref > 0.0 && T_ref > 0.0 && Ea > 0.0)
    return k_ref * std::exp(Ea / R_gas * (1.0 / T_ref - 1.0 / temperature));
  if (k_ref > 0.0)
    return k_ref;
  return 0.0;
}

MoltenSaltRadiolysisDatabase::MoltenSaltRadiolysisDatabase(const std::string & filename)
  : _filename(filename)
{
  MooseUtils::checkFileReadable(_filename);
  std::ifstream stream(_filename);
  try
  {
    stream >> _root;
  }
  catch (std::exception & e)
  {
    mooseError("MSR: failed to parse the chemistry database '", _filename, "':\n", e.what());
  }
}

const nlohmann::json &
MoltenSaltRadiolysisDatabase::findCaseInsensitive(const nlohmann::json & object,
                                                  const std::string & name,
                                                  const std::string & what) const
{
  for (const auto & item : object.items())
    if (iequals(item.key(), name))
      return item.value();
  mooseError("MSR: ", what, " '", name, "' is not present in the database '", _filename, "'.");
}

const nlohmann::json &
MoltenSaltRadiolysisDatabase::kernel(const std::string & salt) const
{
  if (!_root.contains("kernels"))
    mooseError("MSR: the database '", _filename, "' has no 'kernels' section.");
  return findCaseInsensitive(_root.at("kernels"), salt, "salt kernel");
}

std::vector<ReactionData>
MoltenSaltRadiolysisDatabase::parseReactions(const nlohmann::json & reactions) const
{
  std::vector<ReactionData> result;
  for (const auto & rxn : reactions)
  {
    ReactionData data;
    data.name = rxn.value("name", std::string("reaction"));
    const nlohmann::json no_orders = nlohmann::json();
    data.reactants =
        toSpeciesTerms(rxn.at("reactants"), rxn.contains("orders") ? rxn.at("orders") : no_orders);
    data.products = toSpeciesTerms(rxn.at("products"), no_orders);
    data.A = rxn.value("A", 0.0);
    data.k_ref = rxn.value("k_ref", 0.0);
    data.T_ref = rxn.value("T_ref", 0.0);
    data.Ea = rxn.value("Ea_J_mol", 0.0);
    result.push_back(data);
  }
  return result;
}

std::vector<std::string>
MoltenSaltRadiolysisDatabase::coreSpecies(const std::string & salt) const
{
  return kernel(salt).at("species").get<std::vector<std::string>>();
}

std::vector<ReactionData>
MoltenSaltRadiolysisDatabase::coreReactions(const std::string & salt) const
{
  return parseReactions(kernel(salt).at("reactions"));
}

std::vector<std::string>
MoltenSaltRadiolysisDatabase::metalSpecies(const std::string & metal) const
{
  if (!_root.contains("metals"))
    mooseError("MSR: the database '", _filename, "' has no 'metals' section.");
  return findCaseInsensitive(_root.at("metals"), metal, "metal").at("species").get<std::vector<std::string>>();
}

std::vector<ReactionData>
MoltenSaltRadiolysisDatabase::metalReactions(const std::string & metal,
                                             const std::string & salt) const
{
  if (!_root.contains("metals"))
    mooseError("MSR: the database '", _filename, "' has no 'metals' section.");
  const auto & metal_db = findCaseInsensitive(_root.at("metals"), metal, "metal");
  const std::string key = salt + "_reactions";
  // A metal with no reactions for this salt kernel contributes its species but no chemistry.
  if (!metal_db.contains(key))
    return {};
  return parseReactions(metal_db.at(key));
}

std::map<std::string, Real>
MoltenSaltRadiolysisDatabase::defaultGValues(const std::string & radiation,
                                             const std::string & salt) const
{
  std::map<std::string, Real> g_values;
  if (!_root.contains("G_values") || !_root.at("G_values").contains(radiation))
    mooseError("MSR: the database '", _filename, "' has no G values for radiation '", radiation, "'.");
  const auto & by_salt = findCaseInsensitive(_root.at("G_values").at(radiation), salt, "salt kernel");
  for (const auto & item : by_salt.items())
    g_values[item.key()] = item.value().get<Real>();
  return g_values;
}

Real
MoltenSaltRadiolysisDatabase::henryCoefficient(const std::string & gas) const
{
  if (!_root.contains("henry"))
    mooseError("MSR: the database '", _filename, "' has no 'henry' section.");
  return findCaseInsensitive(_root.at("henry"), gas, "gas").get<Real>();
}

std::string
MoltenSaltRadiolysisDatabase::gasLiquidSpecies(const std::string & gas) const
{
  if (!_root.contains("gas_pairs"))
    mooseError("MSR: the database '", _filename, "' has no 'gas_pairs' section.");
  return findCaseInsensitive(_root.at("gas_pairs"), gas, "gas").at("liquid").get<std::string>();
}

std::string
MoltenSaltRadiolysisDatabase::gasPhaseSpecies(const std::string & gas) const
{
  if (!_root.contains("gas_pairs"))
    mooseError("MSR: the database '", _filename, "' has no 'gas_pairs' section.");
  return findCaseInsensitive(_root.at("gas_pairs"), gas, "gas").at("gas").get<std::string>();
}

} // namespace MSR
