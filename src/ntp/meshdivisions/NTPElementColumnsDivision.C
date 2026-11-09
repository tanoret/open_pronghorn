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

#include "NTPElementColumnsDivision.h"
#include "libmesh/elem.h"

#include <cmath>

registerMooseObject("OpenPronghornApp", NTPElementColumnsDivision);

InputParameters
NTPElementColumnsDivision::validParams()
{
  InputParameters params = MeshDivision::validParams();
  params.addClassDescription(
      "Columns of a cylindrical fuel element: equal axial slices times equal azimuthal sectors "
      "of an angle range, numbered slice-major (j n_azimuthal + k).");
  params.addParam<Point>("axis_direction", Point(0, 1, 0), "Direction of the element axis");
  params.addParam<Point>("center", Point(0, 0, 0), "Point on the element axis");
  params.addParam<Point>(
      "azimuthal_start",
      Point(1, 0, 0),
      "Direction e1 of zero angle, normal to the axis; the angle grows towards e1 x axis");
  params.addRequiredParam<Real>("axial_min", "Axial coordinate of the start of the first slice");
  params.addRequiredParam<Real>("axial_max", "Axial coordinate of the end of the last slice");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "n_axial", "n_axial > 0", "Number of axial slices");
  params.addParam<Real>("theta_min", 0, "Start of the angle range [rad]");
  params.addParam<Real>("theta_max", libMesh::pi, "End of the angle range [rad]");
  params.addRangeCheckedParam<unsigned int>(
      "n_azimuthal", 1, "n_azimuthal > 0", "Number of azimuthal sectors over the angle range");
  return params;
}

NTPElementColumnsDivision::NTPElementColumnsDivision(const InputParameters & parameters)
  : MeshDivision(parameters),
    _axis(getParam<Point>("axis_direction").unit()),
    _center(getParam<Point>("center")),
    _axial_min(getParam<Real>("axial_min")),
    _axial_max(getParam<Real>("axial_max")),
    _n_axial(getParam<unsigned int>("n_axial")),
    _theta_min(getParam<Real>("theta_min")),
    _theta_max(getParam<Real>("theta_max")),
    _n_azimuthal(getParam<unsigned int>("n_azimuthal"))
{
  if (!(_axial_max > _axial_min))
    paramError("axial_max", "Must be larger than 'axial_min'");
  if (!(_theta_max > _theta_min) || _theta_max - _theta_min > 2 * libMesh::pi + 1e-12)
    paramError("theta_max", "The angle range must be positive and at most 2 pi");
  Point start = getParam<Point>("azimuthal_start");
  start -= (start * _axis) * _axis;
  if (start.norm() < 1e-12)
    paramError("azimuthal_start", "Must not be parallel to 'axis_direction'");
  _e1 = start.unit();
  _e2 = _e1.cross(_axis);
  NTPElementColumnsDivision::initialize();
}

void
NTPElementColumnsDivision::initialize()
{
  setNumDivisions(_n_axial * _n_azimuthal);
  // every point is assigned to the nearest column
  _mesh_fully_indexed = true;
}

unsigned int
NTPElementColumnsDivision::columnIndex(const Real axial,
                                       const Real theta,
                                       const Real axial_min,
                                       const Real axial_max,
                                       const unsigned int n_axial,
                                       const Real theta_min,
                                       const Real theta_max,
                                       const unsigned int n_azimuthal)
{
  const auto bin = [](const Real value, const Real lo, const Real hi, const unsigned int n)
  {
    const Real f = std::floor((value - lo) / (hi - lo) * n);
    return static_cast<unsigned int>(std::min(std::max(f, Real(0)), Real(n - 1)));
  };
  // the angle in [theta_min, theta_min + 2 pi)
  Real t = theta;
  while (t < theta_min)
    t += 2 * libMesh::pi;
  while (t >= theta_min + 2 * libMesh::pi)
    t -= 2 * libMesh::pi;
  // outside the range: the nearer end
  if (t > theta_max)
    t = (t - theta_max < theta_min + 2 * libMesh::pi - t) ? theta_max : theta_min;
  return bin(axial, axial_min, axial_max, n_axial) * n_azimuthal +
         bin(t, theta_min, theta_max, n_azimuthal);
}

unsigned int
NTPElementColumnsDivision::divisionIndex(const Point & pt) const
{
  const Point d = pt - _center;
  const Real theta = std::atan2(d * _e2, d * _e1);
  return columnIndex(
      d * _axis, theta, _axial_min, _axial_max, _n_axial, _theta_min, _theta_max, _n_azimuthal);
}

unsigned int
NTPElementColumnsDivision::divisionIndex(const Elem & elem) const
{
  return divisionIndex(elem.vertex_average());
}
