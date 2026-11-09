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

#include "MeshDivision.h"

/**
 * Columns of a cylindrical fuel element: n_axial equal slices along the element axis times
 * n_azimuthal equal sectors of the angle range [theta_min, theta_max] about the axis. The index of
 * slice j and sector k is j n_azimuthal + k. With one sector this is the axial slicing of the
 * orifice design; with several, the design can shape the flow around the element as well.
 *
 * The angle is measured about 'axis_direction' through 'center', from 'azimuthal_start' (e1)
 * towards e1 x axis: for the axis y and the start x, theta = atan2(z, x). Points
 * slightly outside the axial or angular range (the faceted mesh, round-off on the symmetry
 * planes) are assigned to the nearest slice or sector.
 */
class NTPElementColumnsDivision : public MeshDivision
{
public:
  static InputParameters validParams();

  NTPElementColumnsDivision(const InputParameters & parameters);

  virtual void initialize() override;
  virtual unsigned int divisionIndex(const Point & pt) const override;
  virtual unsigned int divisionIndex(const Elem & elem) const override;

  /// Column index of a point given its axial coordinate and angle (unit-testable)
  static unsigned int columnIndex(Real axial,
                                  Real theta,
                                  Real axial_min,
                                  Real axial_max,
                                  unsigned int n_axial,
                                  Real theta_min,
                                  Real theta_max,
                                  unsigned int n_azimuthal);

protected:
  /// Unit vector of the element axis
  const Point _axis;
  /// Point on the axis
  const Point _center;
  /// Unit vector of zero angle, normal to the axis
  Point _e1;
  /// Unit vector of angle pi/2: e1 x axis (z for the axis y and the start x)
  Point _e2;
  const Real _axial_min;
  const Real _axial_max;
  const unsigned int _n_axial;
  const Real _theta_min;
  const Real _theta_max;
  const unsigned int _n_azimuthal;
};
