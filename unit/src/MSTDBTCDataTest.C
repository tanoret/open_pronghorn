//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MSTDBTCData.h"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

using namespace Corrosion;

namespace
{

std::string
dataFile(const std::string & name)
{
  const std::string source_path(__FILE__);
  const auto separator = source_path.find_last_of("/\\");
  const std::string source_relative =
      separator == std::string::npos
          ? ""
          : source_path.substr(0, separator) + "/../data/mstdb/" + name;
  const std::array<std::string, 4> candidates = {
      source_relative, "unit/data/mstdb/" + name, "../unit/data/mstdb/" + name, name};
  for (const auto & candidate : candidates)
  {
    if (!candidate.empty() && std::ifstream(candidate).good())
      return candidate;
  }
  throw std::runtime_error("Unable to locate MSTDB-TC unit fixture '" + name + "'.");
}

const std::string fixture_sha256 =
    "3c5b66f569552b67d775307fcc558bd888e0b1de34c288afb719eb566be5543a";

} // namespace

TEST(MSTDBTCData, ParsesAllSupportedEquationTypesAndMetadata)
{
  const MSTDBTCData database(dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat"),
                             false,
                             fixture_sha256);
  EXPECT_EQ(database.version(), "0.0");
  EXPECT_EQ(database.sha256(), fixture_sha256);
  EXPECT_TRUE(database.hashMatchesExpected());
  EXPECT_EQ(database.elementCount(), 3);
  EXPECT_EQ(database.recordCount(), 7);
  EXPECT_EQ(database.species("X_type1").equation_type, 1);
  EXPECT_EQ(database.species("X_poly").equation_type, 4);
  EXPECT_EQ(database.species("X_type10").equation_type, 10);
  EXPECT_EQ(database.species("X_mag13").equation_type, 13);
  EXPECT_EQ(database.species("X_mag16").equation_type, 16);

  const auto provenance = database.provenance();
  EXPECT_EQ(provenance.at("version").get<std::string>(), "0.0");
  EXPECT_EQ(provenance.at("sha256").get<std::string>(), fixture_sha256);
  EXPECT_EQ(provenance.at("n_species_records").get<std::size_t>(), 7);
  EXPECT_FALSE(provenance.at("allow_extrapolation").get<bool>());
}

TEST(MSTDBTCData, GibbsPolynomialIntervalsAndAdditionalTerms)
{
  const MSTDBTCData database(dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat"));
  EXPECT_NEAR(database.standardGibbsJMol("X_poly", 400.0), 320648032.71270835, 1e-6);
  EXPECT_NEAR(database.standardGibbsJMol("X_poly", 800.0), 25625776440.756466, 1e-5);

  const Real at_boundary = 1.0 + 2.0 * 500.0 + 3.0 * 500.0 * std::log(500.0) +
                           4.0 * 500.0 * 500.0 + 5.0 * std::pow(500.0, 3.0) +
                           6.0 / 500.0 + 7.0 * std::log(500.0);
  EXPECT_NEAR(database.standardGibbsJMol("X_poly", 500.0), at_boundary, 1e-6);
  EXPECT_DOUBLE_EQ(database.standardGibbsJMol("X_type1", 923.15), 10.0);
  EXPECT_DOUBLE_EQ(database.standardGibbsJMol("X_type10", 923.15), 20.0);
}

TEST(MSTDBTCData, SGTEBranchesAndMagneticRecords)
{
  const MSTDBTCData database(dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat"));
  const Real expected = -1179.7062552844277;
  EXPECT_NEAR(MSTDBTCData::magneticGibbsJMol(923.15, 1043.0, 2.22, 1.0, 0.4),
              expected,
              1e-10);
  EXPECT_NEAR(database.standardGibbsJMol("X_mag13", 923.15), expected, 1e-10);
  EXPECT_NEAR(database.standardGibbsJMol("X_mag16", 923.15), 243.04989984118538, 1e-10);
  EXPECT_DOUBLE_EQ(MSTDBTCData::magneticGibbsJMol(923.15, 0.0, 2.22, 1.0, 0.4),
                   0.0);
}

TEST(MSTDBTCData, ExactCaseSensitiveOccurrenceLookup)
{
  const MSTDBTCData database(dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat"));
  EXPECT_EQ(database.occurrenceCount("X_duplicate"), 2);
  EXPECT_DOUBLE_EQ(database.standardGibbsJMol("X_duplicate", 900.0), 40.0);
  EXPECT_DOUBLE_EQ(database.standardGibbsJMol("X_duplicate", 900.0, 0, false), 30.0);
  EXPECT_DOUBLE_EQ(database.standardGibbsJMol("X_duplicate", 900.0, 1, false), 40.0);
  EXPECT_THROW(database.species("x_duplicate"), std::exception);
  EXPECT_THROW(database.species("X_duplicate", 2), std::exception);
  EXPECT_THROW(database.species("X_duplicate", -2), std::exception);
}

TEST(MSTDBTCData, ReactionGibbsAndEquilibriumConstant)
{
  const MSTDBTCData database(dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat"));
  const std::map<std::string, Real> reaction = {{"X_type10", 1.0}, {"X_type1", -1.0}};
  EXPECT_DOUBLE_EQ(database.reactionGibbsJMol(reaction, 900.0), 10.0);
  EXPECT_NEAR(database.equilibriumLogConstant(reaction, 900.0),
              -10.0 / (MSTDBTCData::gas_constant_J_mol_K * 900.0),
              1e-15);
  EXPECT_THROW(database.reactionGibbsJMol({}, 900.0), std::exception);
}

TEST(MSTDBTCData, ExtrapolationAndInvalidStateAreExplicit)
{
  const MSTDBTCData database(dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat"));
  EXPECT_THROW(database.standardGibbsJMol("X_poly", 2100.0), std::exception);
  EXPECT_TRUE(std::isfinite(database.standardGibbsJMol("X_poly", 2100.0, true)));
  EXPECT_THROW(database.standardGibbsJMol("X_poly", 0.0), std::exception);
  EXPECT_THROW(database.standardGibbsJMol(
                   "X_poly", std::numeric_limits<Real>::quiet_NaN()),
               std::exception);

  const MSTDBTCData compatibility_database(
      dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat"), true);
  EXPECT_TRUE(std::isfinite(compatibility_database.standardGibbsJMol("X_poly", 2100.0)));
}

TEST(MSTDBTCData, HashBindingAndMismatchOverride)
{
  const std::string filename = dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat");
  const std::string wrong_hash(64, '0');
  EXPECT_THROW(MSTDBTCData(filename, false, wrong_hash, false), std::exception);

  const MSTDBTCData intentionally_unbound(filename, false, wrong_hash, true);
  EXPECT_FALSE(intentionally_unbound.hashMatchesExpected());
  EXPECT_FALSE(
      intentionally_unbound.provenance().at("sha256_matches_expected").get<bool>());
}

TEST(MSTDBTCData, PairBindsVersionAndBothHashes)
{
  const MSTDBTCPair pair(dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat"),
                         dataFile("MSTDB-TC_V0.0_Chlorides_No_Func.dat"),
                         "0.0",
                         fixture_sha256,
                         fixture_sha256);
  EXPECT_EQ(pair.version(), "0.0");
  EXPECT_EQ(pair.fluoride().recordCount(), 7);
  EXPECT_EQ(pair.chloride().recordCount(), 7);
  EXPECT_EQ(pair.provenance().at("expected_version").get<std::string>(), "0.0");

  EXPECT_THROW(MSTDBTCPair(dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat"),
                           dataFile("MSTDB-TC_V0.1_Chlorides_No_Func.dat")),
               std::exception);
  EXPECT_THROW(MSTDBTCPair(dataFile("MSTDB-TC_V0.0_Fluorides_No_Func.dat"),
                           dataFile("MSTDB-TC_V0.0_Chlorides_No_Func.dat"),
                           "3.1"),
               std::exception);
}

TEST(MSTDBTCData, MalformedCandidateFailsFast)
{
  EXPECT_THROW(MSTDBTCData(dataFile("MSTDB-TC_V0.0_Malformed_No_Func.dat")),
               std::exception);
  EXPECT_THROW(MSTDBTCData(dataFile("MSTDB-TC_V0.0_Truncated_No_Func.dat")),
               std::exception);
}

TEST(MSTDBTCData, SkipsSolutionHeaderWithoutStoichiometry)
{
  // Real V3.1 solution declarations contain an element label such as F followed by a bare
  // equation-type/interval-count pair ("1 1"). It is not a standard-state species record.
  const MSTDBTCData database(dataFile("MSTDB-TC_V0.0_SolutionHeader_No_Func.dat"));
  EXPECT_EQ(database.recordCount(), 1);
  EXPECT_EQ(database.occurrenceCount("F"), 0);
  EXPECT_DOUBLE_EQ(database.standardGibbsJMol("X_valid", 923.15), 10.0);
}
