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

#pragma once

#include "FaceInfo.h"
#include "MooseTypes.h"
#include "libmesh/vector_value.h"

/**
 * Face interpolation of a transport coefficient between the two cells of a face, for coefficients
 * acting in series along the line between the cell centers (a mobility, a conductance).
 *
 * The element-side half segment has length d_Cf = (1 - gC) d and the neighbor-side one
 * d_fN = gC d, d = |x_N - x_C|, gC the linear interpolation weight of the element. A flux through
 * both halves in series sees the resistance d_Cf / a_C + d_fN / a_N, so the face coefficient is
 *
 *   a_f = d / (d_Cf / a_C + d_fN / a_N) = 1 / ((1 - gC) / a_C + gC / a_N).
 *
 * Moose::FV::harmonicInterpolation uses the linear interpolation weights instead,
 * 1 / (gC / a_C + (1 - gC) / a_N), which attributes each half segment to the wrong cell. The two
 * agree on uniform meshes only. Across a thin porous layer next to a thick one, the difference is
 * the resistance of the thick cell's half segment counted at the thin layer's coefficient.
 *
 * A source term of the momentum balance acting along the segment, of force density F_C and F_N in
 * the two cells, enters with the matching segment average (1 - gC) F_C + gC F_N.
 */
namespace NTPSeriesInterpolation
{
/// Series (resistance-weighted harmonic) face value, elem value first
template <typename T>
T
coefficient(const T & elem_value, const T & neighbor_value, const FaceInfo & fi)
{
  const Real gc = fi.gC();
  return 1.0 / ((1.0 - gc) / elem_value + gc / neighbor_value);
}

/// Component-wise series face value of a diagonal coefficient
template <typename T>
libMesh::VectorValue<T>
coefficient(const libMesh::VectorValue<T> & elem_value,
            const libMesh::VectorValue<T> & neighbor_value,
            const FaceInfo & fi)
{
  libMesh::VectorValue<T> result;
  for (const auto i : make_range(Moose::dim))
    result(i) = coefficient(elem_value(i), neighbor_value(i), fi);
  return result;
}

/// Segment average of a field acting along the line between the cell centers, elem value first
template <typename T>
T
segmentAverage(const T & elem_value, const T & neighbor_value, const FaceInfo & fi)
{
  const Real gc = fi.gC();
  return (1.0 - gc) * elem_value + gc * neighbor_value;
}
}
