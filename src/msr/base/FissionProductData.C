//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FissionProductData.h"
#include "MooseError.h"
#include "MooseUtils.h"

#include <cmath>
#include <fstream>

namespace FissionProducts
{

Real
NuclideData::decayConstant() const
{
  return half_life > 0.0 ? std::log(2.0) / half_life : 0.0;
}

FissionProductDatabase::FissionProductDatabase(const std::string & filename) : _filename(filename)
{
  MooseUtils::checkFileReadable(_filename);
  std::ifstream stream(_filename);
  nlohmann::json root;
  try
  {
    stream >> root;
  }
  catch (std::exception & e)
  {
    mooseError("FissionProducts: failed to parse the database '", _filename, "':\n", e.what());
  }

  if (!root.contains("nuclides"))
    mooseError("FissionProducts: the database '", _filename, "' has no 'nuclides' section.");

  for (const auto & item : root.at("nuclides").items())
  {
    const auto & entry = item.value();
    NuclideData data;
    data.name = item.key();
    data.group = entry.at("group").get<std::string>();
    data.half_life = entry.at("half_life").is_null() ? 0.0 : entry.at("half_life").get<Real>();
    data.molar_mass = entry.value("molar_mass", 0.0);
    for (const auto & branch : entry.at("decay"))
      data.decay.push_back(
          {branch.at("daughter").get<std::string>(), branch.at("branching").get<Real>()});
    if (entry.contains("capture"))
    {
      data.capture_cross_section = entry.at("capture").value("cross_section", 0.0);
      data.capture_product = entry.at("capture").value("product", std::string());
    }
    for (const auto & ys : entry.at("yields").items())
    {
      data.yields[ys.key()] = {ys.value().value("independent", 0.0),
                               ys.value().value("cumulative", 0.0)};
      if (std::find(_yield_sets.begin(), _yield_sets.end(), ys.key()) == _yield_sets.end())
        _yield_sets.push_back(ys.key());
    }
    _nuclides[data.name] = data;
  }
}

const NuclideData &
FissionProductDatabase::nuclide(const std::string & name) const
{
  const auto it = _nuclides.find(name);
  if (it == _nuclides.end())
    mooseError(
        "FissionProducts: nuclide '", name, "' is not present in the database '", _filename, "'.");
  return it->second;
}

std::vector<std::string>
FissionProductDatabase::nuclideNames() const
{
  std::vector<std::string> names;
  for (const auto & item : _nuclides)
    names.push_back(item.first);
  return names;
}

} // namespace FissionProducts
