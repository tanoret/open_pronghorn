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

#include "LinearFVFluxKernel.h"

#include "MooseLinearVariableFV.h"

/**
 * The part of the Newtonian viscous stress of a compressible gas that the Laplacian momentum
 * flux omits, on the internal faces between cells of the kernel's blocks:
 *
 *   tau = mu (grad u + grad u^T - 2/3 (div u) I),   tau_extra . n = mu (grad u^T - 2/3 div u I) . n
 *
 * for the momentum component i: mu (d_i u_j - 2/3 delta_ij div u) n_j. The velocity gradients are
 * the cell gradients of the last iteration, interpolated linearly to the face (explicit). In RZ
 * the divergence includes u_r / r. The kernel is meant for free channels (the plenum and the
 * chamber of a particle-bed element): across a porosity jump the superficial velocity gradient
 * is not a strain rate, and in the beds the drag dominates the stress by many orders of magnitude.
 * Faces on the edge of the blocks and boundary faces get no contribution.
 */
class NTPLinearFVViscousStressCorrection : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();

  NTPLinearFVViscousStressCorrection(const InputParameters & params);

  virtual bool hasFaceSide(const FaceInfo & fi, bool fi_elem_side) const override;

  virtual Real computeElemMatrixContribution() override { return 0; }
  virtual Real computeNeighborMatrixContribution() override { return 0; }
  virtual Real computeElemRightHandSideContribution() override;
  virtual Real computeNeighborRightHandSideContribution() override;
  virtual Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition &) override
  {
    return 0;
  }
  virtual Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition &) override
  {
    return 0;
  }

protected:
  /// The extra stress flux through the current face, mu (grad u^T - 2/3 div u I)_i . n
  Real faceStress();

  /// Dynamic viscosity
  const Moose::Functor<Real> & _mu;
  /// Momentum component of this kernel
  const unsigned int _index;
  /// Mesh dimension
  const unsigned int _dim;
  /// Velocity components
  std::vector<MooseLinearVariableFVReal *> _vel;
  /// Functors whose cell gradients are used (the velocity by default, e.g. the channel velocity)
  std::vector<const Moose::Functor<Real> *> _grad_functors;
  /// Coordinate system of the blocks and the radial coordinate in RZ
  const Moose::CoordinateSystemType _coord_type;
  const unsigned int _rz_radial_coord;
  /// Cached face stress
  Real _face_stress;
};
