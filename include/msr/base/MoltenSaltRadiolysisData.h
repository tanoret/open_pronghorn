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

#include <string>
#include <vector>
#include <map>

/**
 * Molten salt radiolysis chemistry database.
 *
 * The chemistry data (species lists, mass-action reaction networks, default radiolytic yields and
 * Henry coefficients for the chloride and fluoride salt kernels, plus the optional metal redox
 * templates) is read from a JSON file at run time. The bundled default is data/msr_database.json,
 * a translation of the validated 0D Python model in MoltenSaltRadiolysis
 * (msr_radiolysis/data/database.yaml); users can supply their own file to change the network
 * without recompiling.
 *
 * Species are named with plain identifiers (no symbols) so that they are valid MOOSE variable
 * names, e.g. e_sol (e_s-), Cl_ion (Cl-), Cl_rad (Cl.), Cl2m_rad (Cl2.-), Cl3_ion (Cl3-),
 * Cl2_diss, Cl2_gas, F_ion, F_rad, F2_diss, F2_gas, Zn_II/Zn_I, Cr_II/Cr_III/Cr_I, U_IV/U_III.
 */
namespace MSR
{

/// Universal gas constant [J/(mol K)] (matches msr_radiolysis core.py R_GAS)
const Real R_gas = 8.314462618;
/// Avogadro number [1/mol]
const Real avogadro = 6.02214076e23;
/// Electron volt to Joule conversion [J/eV]
const Real eV_to_J = 1.602176634e-19;

/**
 * Convert a radiolytic yield (G value, molecules per 100 eV) and a volumetric dose rate [J/m^3/s]
 * into a molar source term [mol/m^3/s]. Mirrors msr_radiolysis/utils.py::g_to_source:
 *
 *   S = G * dose_rate * 1 / (100 eV) * 1 / N_A
 */
Real gToSource(Real g_value, Real dose_rate);

/// A single reactant or product entry: a species, its stoichiometric coefficient and its kinetic
/// order (the order defaults to the stoichiometric coefficient).
struct SpeciesTerm
{
  std::string species;
  Real coeff;
  Real order;
};

/**
 * A mass-action reaction with Arrhenius kinetics. Mirrors msr_radiolysis/core.py::Reaction.
 *
 * The forward rate constant is selected, in order of precedence:
 *   k(T) = A * exp(-Ea / (R T))                    when A and Ea are supplied
 *   k(T) = k_ref * exp(Ea / R * (1/T_ref - 1/T))   when k_ref, T_ref and Ea are supplied
 *   k(T) = k_ref                                    when only k_ref is supplied
 *
 * A value of zero marks an unsupplied parameter; all supplied physical values are strictly positive.
 */
struct ReactionData
{
  std::string name;
  std::vector<SpeciesTerm> reactants;
  std::vector<SpeciesTerm> products;
  Real A = 0.0;
  Real k_ref = 0.0;
  Real T_ref = 0.0;
  Real Ea = 0.0;

  /// Forward rate constant k(T) in SI units (m^3/mol/s for bimolecular, 1/s for unimolecular).
  Real rate(Real temperature) const;
};

/**
 * Reads and queries the molten salt radiolysis chemistry database from a JSON file.
 */
class MoltenSaltRadiolysisDatabase
{
public:
  /// Load and parse the database from a JSON file.
  MoltenSaltRadiolysisDatabase(const std::string & filename);

  /// Names of the core species for a salt kernel ("chloride" or "fluoride").
  std::vector<std::string> coreSpecies(const std::string & salt) const;

  /// Base (non metal) reaction network for a salt kernel.
  std::vector<ReactionData> coreReactions(const std::string & salt) const;

  /// Species introduced by a metal ("Zn", "Cr" or "U"); case-insensitive.
  std::vector<std::string> metalSpecies(const std::string & metal) const;

  /// Templated metal reactions for a given salt kernel; case-insensitive metal name.
  std::vector<ReactionData> metalReactions(const std::string & metal,
                                           const std::string & salt) const;

  /// Default G values [molecules/100eV] keyed by species, for a radiation type and salt kernel.
  std::map<std::string, Real> defaultGValues(const std::string & radiation,
                                             const std::string & salt) const;

  /// Henry coefficient [mol/(m^3 Pa)] for a dissolved diatomic gas ("Cl2" or "F2"); case-insensitive.
  Real henryCoefficient(const std::string & gas) const;

  /// Dissolved-liquid species name partnered with a gas ("Cl2" -> "Cl2_diss"); case-insensitive.
  std::string gasLiquidSpecies(const std::string & gas) const;

  /// Gas-phase species name for a gas ("Cl2" -> "Cl2_gas"); case-insensitive.
  std::string gasPhaseSpecies(const std::string & gas) const;

protected:
  /// Parse a JSON array of reactions into ReactionData objects.
  std::vector<ReactionData> parseReactions(const nlohmann::json & reactions) const;

  /// Look up a salt kernel object, erroring with context if absent.
  const nlohmann::json & kernel(const std::string & salt) const;

  /// Find an object member whose key matches name case-insensitively; throws if absent.
  const nlohmann::json & findCaseInsensitive(const nlohmann::json & object,
                                             const std::string & name,
                                             const std::string & what) const;

  /// File the database was read from (for error messages).
  const std::string _filename;

  /// Parsed JSON document.
  nlohmann::json _root;
};

} // namespace MSR
