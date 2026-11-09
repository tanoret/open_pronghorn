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
#include "NTPOrificeDesignAux.h"

#include <cmath>

TEST(NTPOrificeDesign, hotSliceIsOpenedAtFixedMean)
{
  // rises 200 and 100 K, mean 150: slice 0 opens by (200/150)^-1, slice 1 throttles by
  // (100/150)^-1, then the profile is shifted to a mean log of 2
  std::vector<Real> log_chi = {2, 2};
  EXPECT_TRUE(NTPOrificeDesignAux::updateProfile(log_chi, {230, 130}, 30, 1, 10, 2, 1));
  EXPECT_NEAR(log_chi[1] - log_chi[0], std::log(200.0 / 100.0), 1e-12);
  EXPECT_NEAR(0.5 * (log_chi[0] + log_chi[1]), 2, 1e-12);
}

TEST(NTPOrificeDesign, uniformRiseIsAFixedPoint)
{
  std::vector<Real> log_chi = {0.5, 1.5, 4.0};
  EXPECT_TRUE(NTPOrificeDesignAux::updateProfile(log_chi, {500, 500, 500}, 30, 0.7, 2, 2, 1));
  EXPECT_NEAR(log_chi[0], 0.5, 1e-12);
  EXPECT_NEAR(log_chi[1], 1.5, 1e-12);
  EXPECT_NEAR(log_chi[2], 4.0, 1e-12);
}

TEST(NTPOrificeDesign, changeIsLimitedAndNoiseIgnored)
{
  // rises 3000 and 30 K, mean 1515: slice 0 opens by ln(3000/1515) (inside the limit ln 2), slice 1
  // would throttle by ln(1515/30) but is limited to ln 2
  std::vector<Real> log_chi = {0, 0};
  NTPOrificeDesignAux::updateProfile(log_chi, {3030, 60}, 30, 1, 2, 0, 1);
  EXPECT_NEAR(log_chi[1] - log_chi[0], std::log(2.0) + std::log(3000.0 / 1515.0), 1e-12);
  // a mean rise below the threshold leaves the profile
  std::vector<Real> noise = {1, 2};
  EXPECT_FALSE(NTPOrificeDesignAux::updateProfile(noise, {30.2, 29.9}, 30, 1, 2, 0, 1));
  EXPECT_DOUBLE_EQ(noise[0], 1.0);
  EXPECT_DOUBLE_EQ(noise[1], 2.0);
}
