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

#include "NTPElementColumnsDivision.h"

// 12 slices over [0, 1.2] times 4 sectors over the half element [0, pi]
static unsigned int
column(Real axial, Real theta)
{
  return NTPElementColumnsDivision::columnIndex(axial, theta, 0, 1.2, 12, 0, libMesh::pi, 4);
}

TEST(NTPElementColumnsDivision, sliceMajorNumbering)
{
  EXPECT_EQ(column(0.05, 0.1), 0u);
  EXPECT_EQ(column(0.05, 3.0), 3u);
  EXPECT_EQ(column(0.15, 0.1), 4u);
  EXPECT_EQ(column(1.15, 1.0), 45u);
}

TEST(NTPElementColumnsDivision, outOfRangeGoesToNearestColumn)
{
  // round-off below the symmetry plane theta = 0, and just past theta = pi (atan2 gives -pi)
  EXPECT_EQ(column(0.05, -1e-14), 0u);
  EXPECT_EQ(column(0.05, -libMesh::pi + 1e-14), 3u);
  // faceted cells slightly outside the axial range
  EXPECT_EQ(column(-1e-9, 0.1), 0u);
  EXPECT_EQ(column(1.2 + 1e-9, 0.1), 44u);
}

TEST(NTPElementColumnsDivision, singleSectorIsAxialSlicing)
{
  for (const auto j : make_range(12))
    EXPECT_EQ(NTPElementColumnsDivision::columnIndex(0.1 * j + 0.05, 2.0, 0, 1.2, 12, 0, 3.2, 1),
              static_cast<unsigned int>(j));
}
