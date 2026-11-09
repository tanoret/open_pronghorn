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

#include "FVFluxKernel.h"
#include "FVDiffusionInterpolationInterface.h"
#include "NTPDarcyFaceFlux.h"

/**
 * Face mass flux of a Darcy-type flow with an optional momentum-flux correction, shared by the
 * mass and energy kernels so that both see exactly the same fluxes:
 *
 *   m_f = -(rho M)_f [ grad(p) + grad(J) ] . n,
 *
 * rho M being the diagonal mass mobility (series face combination by default, see
 * NTPDarcyFaceFlux) and J = rho u^2 a scalar
 * momentum flux. With J = 0 this is Darcy-Forchheimer flow. In a slender free channel, J = rho
 * u_z^2 turns the relation into the steady one-dimensional momentum balance d(rho u_z^2)/dz + dp/dz
 * = -(friction) u_z, which carries the acceleration of heated gas and the pressure recovery of a
 * dead-ended manifold. J only acts on faces whose two cells are both in 'momentum_blocks' (the free
 * channels).
 *
 * Orifice solved with the flow: with a scalar 'orifice_log_multiplier' lambda holding one value per
 * slice of 'mesh_division', the mobility of a cell in 'orifice_blocks' (the cold frit) and slice j
 * is divided by exp(lambda_j) before the face interpolation, i.e. its drag is multiplied by
 * chi_j = exp(lambda_j). The scalar is read here, in the kernel, so that its derivatives enter the
 * Jacobian.
 */
class NTPDarcyFluxInterface : public FVFluxKernel, public FVDiffusionInterpolationInterface
{
public:
  static InputParameters validParams();

  NTPDarcyFluxInterface(const InputParameters & params);

protected:
  /// Mass flux through the current face, positive along the face normal (element to neighbor)
  ADReal faceMassFlux(const Moose::StateArg & state) const;

  /// The face-flux calculator shared with the other Darcy objects
  NTPDarcyFaceFlux _flux;
};
