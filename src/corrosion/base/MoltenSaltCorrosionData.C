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

#include "MoltenSaltCorrosionData.h"
#include "MooseError.h"
#include "MooseUtils.h"

#include <cctype>
#include <fstream>

namespace Corrosion
{

namespace
{
// Case-insensitive string comparison so callers can pass element/material tokens in any case
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
}

MoltenSaltCorrosionDatabase::MoltenSaltCorrosionDatabase(const std::string & filename)
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
    mooseError("Corrosion: failed to parse the database '", _filename, "':\n", e.what());
  }

  if (!_root.contains("calibrated_parameters"))
    mooseError("Corrosion: the database '", _filename, "' has no 'calibrated_parameters' section.");
  for (const auto & item : _root.at("calibrated_parameters").items())
    _parameters[item.key()] = item.value().get<Real>();
}

const nlohmann::json *
MoltenSaltCorrosionDatabase::findCaseInsensitive(const nlohmann::json & object,
                                                 const std::string & name) const
{
  for (const auto & item : object.items())
    if (iequals(item.key(), name))
      return &item.value();
  return nullptr;
}

bool
MoltenSaltCorrosionDatabase::hasElement(const std::string & name) const
{
  if (!_root.contains("elements"))
    return false;
  return findCaseInsensitive(_root.at("elements"), name) != nullptr;
}

ElementProperties
MoltenSaltCorrosionDatabase::element(const std::string & name) const
{
  if (!_root.contains("elements"))
    mooseError("Corrosion: the database '", _filename, "' has no 'elements' section.");

  const auto & elements = _root.at("elements");
  const nlohmann::json * entry = findCaseInsensitive(elements, name);
  if (!entry)
    entry = findCaseInsensitive(elements, "generic_metal");
  if (!entry)
    mooseError("Corrosion: element '",
               name,
               "' is not in the database '",
               _filename,
               "' and no 'generic_metal' fall-back is provided.");

  ElementProperties props;
  props.valence = entry->value("valence", props.valence);
  props.molar_mass_g_mol = entry->value("molar_mass_g_mol", props.molar_mass_g_mol);
  props.diffusivity_m2_s = entry->value("diffusivity_m2_s", props.diffusivity_m2_s);
  props.E0_V = entry->value("E0_V", props.E0_V);
  props.alpha_a = entry->value("alpha_a", props.alpha_a);
  props.alpha_c = entry->value("alpha_c", props.alpha_c);
  props.c_ref_mol_m3 = entry->value("c_ref_mol_m3", props.c_ref_mol_m3);
  return props;
}

Real
MoltenSaltCorrosionDatabase::lookupWithFallback(const nlohmann::json & object,
                                                const std::string & name,
                                                const std::string & fallback_key,
                                                const std::string & what) const
{
  if (const nlohmann::json * entry = findCaseInsensitive(object, name))
    return entry->get<Real>();
  if (const nlohmann::json * fallback = findCaseInsensitive(object, fallback_key))
    return fallback->get<Real>();
  mooseError("Corrosion: ",
             what,
             " '",
             name,
             "' is not in the database '",
             _filename,
             "' and no '",
             fallback_key,
             "' fall-back is provided.");
}

Real
MoltenSaltCorrosionDatabase::density(const std::string & material_class) const
{
  if (!_root.contains("densities_g_cm3"))
    mooseError("Corrosion: the database '", _filename, "' has no 'densities_g_cm3' section.");
  return lookupWithFallback(
      _root.at("densities_g_cm3"), material_class, "generic_metal", "material density");
}

Real
MoltenSaltCorrosionDatabase::crWeightFraction(const std::string & material_class) const
{
  if (!_root.contains("alloy_cr_wt_frac"))
    mooseError("Corrosion: the database '", _filename, "' has no 'alloy_cr_wt_frac' section.");
  return lookupWithFallback(
      _root.at("alloy_cr_wt_frac"), material_class, "generic_metal", "chromium weight fraction");
}

Real
MoltenSaltCorrosionDatabase::saltDensity() const
{
  return _root.value("salt_density_g_cm3", 2.0);
}

Real
MoltenSaltCorrosionDatabase::parameter(const std::string & name) const
{
  const auto it = _parameters.find(name);
  if (it == _parameters.end())
    mooseError("Corrosion: calibrated parameter '",
               name,
               "' is not present in the database '",
               _filename,
               "'.");
  return it->second;
}

} // namespace Corrosion
