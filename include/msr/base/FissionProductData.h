//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"
#include "json.h"

#include <map>
#include <string>
#include <vector>

/**
 * Fission-product nuclide database.
 *
 * The data (half-lives, decay branchings, 2200 m/s capture cross sections, U-235 independent and
 * cumulative fission yields, chemical groups) are read from a JSON file at run time. The bundled
 * default, data/fission_product_database.json, is generated from the ENDF/B-VIII.0 decay, fission
 * yield and neutron sub-libraries by scripts/fission_products/build_fission_product_database.py.
 */
namespace FissionProducts
{

/// One decay branch to a daughter nuclide
struct DecayBranch
{
  std::string daughter;
  Real branching;
};

/// Data for one nuclide
struct NuclideData
{
  std::string name;
  /// Chemical group: noble_gas, noble_metal, halogen or salt_seeking
  std::string group;
  /// Half-life [s]; zero for stable nuclides
  Real half_life = 0.0;
  /// Molar mass [g/mol]
  Real molar_mass = 0.0;
  /// Decay branches (daughters may be untracked, in which case the decay is a pure sink)
  std::vector<DecayBranch> decay;
  /// 2200 m/s (n,gamma) cross section [barn]
  Real capture_cross_section = 0.0;
  /// Capture product name
  std::string capture_product;
  /// Independent and cumulative yields per yield set (e.g. "U235_thermal")
  std::map<std::string, std::pair<Real, Real>> yields;

  /// Decay constant [1/s] (zero for stable nuclides)
  Real decayConstant() const;
};

class FissionProductDatabase
{
public:
  FissionProductDatabase(const std::string & filename);

  /// Whether the database contains a nuclide
  bool hasNuclide(const std::string & name) const { return _nuclides.count(name); }

  /// Data for a nuclide (errors if absent)
  const NuclideData & nuclide(const std::string & name) const;

  /// All nuclide names in the database
  std::vector<std::string> nuclideNames() const;

  /// The yield sets available (e.g. U235_thermal, U235_fast)
  const std::vector<std::string> & yieldSets() const { return _yield_sets; }

  /// Name of the file the database was read from
  const std::string & filename() const { return _filename; }

private:
  const std::string _filename;
  std::map<std::string, NuclideData> _nuclides;
  std::vector<std::string> _yield_sets;
};

/// Barn to m^2
const Real barn = 1.0e-28;
/// Avogadro number [1/mol]
const Real avogadro = 6.02214076e23;
/// MeV to Joule
const Real MeV_to_J = 1.602176634e-13;

} // namespace FissionProducts
