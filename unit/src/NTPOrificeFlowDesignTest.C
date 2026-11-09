//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Copyright 2025, Battelle Energy Alliance, LLC
//* ALL RIGHTS RESERVED
//*

#include "gtest/gtest.h"
#include "NTPOrificeFlowDesignAux.h"

#include <cmath>

TEST(NTPOrificeFlowDesign, overfedSliceIsThrottledAtFixedMean)
{
  // flows 3 and 1 for equal shares: target 2 each; slice 0 throttles by ln(3/2), slice 1 opens by
  // ln(1/2), then the profile is shifted to a mean log of 2
  std::vector<Real> log_chi = {2, 2};
  const Real mismatch =
      NTPOrificeFlowDesignAux::matchFlows(log_chi, {3, 1}, {0.5, 0.5}, 1, 10, 2);
  EXPECT_NEAR(mismatch, std::log(2.0), 1e-12);
  EXPECT_NEAR(log_chi[0] - log_chi[1], std::log(3.0), 1e-12);
  EXPECT_NEAR(0.5 * (log_chi[0] + log_chi[1]), 2, 1e-12);
}

TEST(NTPOrificeFlowDesign, matchedFlowsAreAFixedPoint)
{
  std::vector<Real> log_chi = {0.5, 1.5, 4.0};
  const Real mismatch =
      NTPOrificeFlowDesignAux::matchFlows(log_chi, {0.2, 0.3, 0.5}, {0.2, 0.3, 0.5}, 0.7, 2, 2);
  EXPECT_NEAR(mismatch, 0, 1e-12);
  EXPECT_NEAR(log_chi[0], 0.5, 1e-12);
  EXPECT_NEAR(log_chi[1], 1.5, 1e-12);
  EXPECT_NEAR(log_chi[2], 4.0, 1e-12);
}

TEST(NTPOrificeFlowDesign, reversedFlowIsOpenedByTheLimit)
{
  // slice 1 has a reversed flow: it opens by the largest step ln 2, slice 0 throttles by ln 2
  std::vector<Real> log_chi = {0, 0};
  NTPOrificeFlowDesignAux::matchFlows(log_chi, {1.1, -0.1}, {0.5, 0.5}, 1, 2, 0);
  EXPECT_NEAR(log_chi[0] - log_chi[1], 2 * std::log(2.0), 1e-12);
}

TEST(NTPOrificeFlowDesign, hotSliceGetsALargerShare)
{
  // rises 200 and 100 K, mean 150: the shares move by (rise / mean)^1 and are renormalized
  std::vector<Real> shares = {0.5, 0.5};
  EXPECT_TRUE(NTPOrificeFlowDesignAux::correctShares(shares, {230, 130}, 30, 1, 10, 1));
  EXPECT_NEAR(shares[0] / shares[1], 2, 1e-12);
  EXPECT_NEAR(shares[0] + shares[1], 1, 1e-12);
  // below the minimum mean rise, or with gamma = 0, the shares stay
  std::vector<Real> noise = {0.4, 0.6};
  EXPECT_FALSE(NTPOrificeFlowDesignAux::correctShares(noise, {30.2, 29.9}, 30, 1, 2, 1));
  EXPECT_FALSE(NTPOrificeFlowDesignAux::correctShares(noise, {230, 130}, 30, 0, 2, 1));
  EXPECT_DOUBLE_EQ(noise[0], 0.4);
  EXPECT_DOUBLE_EQ(noise[1], 0.6);
}
