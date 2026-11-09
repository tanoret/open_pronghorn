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

#pragma once

#include "MooseTypes.h"
#include "MooseFunctor.h"
#include "MathFVUtils.h"
#include "MooseEnum.h"
#include "MooseVariableFV.h"

#include <set>

class MeshDivision;
class FaceInfo;
class MooseVariableScalar;

/**
 * Face mass flux of the Darcy-type flow of a particle-bed element, as a plain calculator shared by
 * the kernels (mass, energy) and by the cell-velocity reconstruction, so that all of them see the
 * same face fluxes:
 *
 *   m_f = -(rho M)_f [ grad(p) + grad(J) ] . n_f            (positive along the face normal)
 *
 * The mobility rho M of each side is that side's region law evaluated at the face, with the
 * compact two-point face pressure gradient. By default the two sides are combined as series
 * resistances of the two half segments between the cell centers (NTPSeriesInterpolation); the
 * MOOSE harmonic and linear averages are available for comparison. The optional orifice divides the
 * mobility of the cells of the orifice blocks in slice j by exp(lambda_j); the optional momentum
 * flux J acts on the axial faces between two channel cells.
 */
class NTPDarcyFaceFlux
{
public:
  /// Face interpolation of the mobility
  enum class Interpolation
  {
    SERIES,
    HARMONIC,
    AVERAGE
  };

  /// Parameter enum of the face interpolation, default series
  static MooseEnum interpolationEnum();
  /// The interpolation selected by a parameter value of interpolationEnum()
  static Interpolation interpolation(const MooseEnum & value);

  NTPDarcyFaceFlux(const MooseVariableFV<Real> & pressure,
                   const Moose::Functor<ADRealVectorValue> & rho_mobility,
                   Interpolation interp_method,
                   bool face_mobility);

  /// Add a momentum flux J acting between two cells of the given subdomains
  void setMomentumFlux(const Moose::Functor<ADReal> & J, const std::set<SubdomainID> & blocks);

  /// Add an orifice exp(lambda_j) on the drag of the cells of the given subdomains
  void setOrifice(const MeshDivision & division,
                  const MooseVariableScalar & log_multiplier,
                  const std::set<SubdomainID> & blocks);

  /**
   * Component j of a scalar variable as an AD number, from the current solution vector with the
   * derivative seeded at the component's dof. Scalar AD values cached by the variable are not
   * always seeded in finite-volume face loops; this evaluation is.
   */
  static ADReal scalarComponent(const MooseVariableScalar & var, unsigned int j);

  /// Mass flux through the face, positive along the face normal (element to neighbor)
  ADReal massFlux(const FaceInfo & fi, const Moose::StateArg & state) const;

private:
  /// Mobility of one side of the face, divided by the orifice multiplier if applicable
  ADRealVectorValue
  sideMobility(const FaceInfo & fi, const Elem & elem, const Moose::StateArg & state) const;

  const MooseVariableFV<Real> & _pressure;
  const Moose::Functor<ADRealVectorValue> & _rho_mobility;
  const Interpolation _interp_method;
  const bool _face_mobility;

  const Moose::Functor<ADReal> * _momentum_flux = nullptr;
  std::set<SubdomainID> _momentum_blocks;

  const MeshDivision * _division = nullptr;
  const MooseVariableScalar * _log_orifice = nullptr;
  std::set<SubdomainID> _orifice_blocks;
};
