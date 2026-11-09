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

#include "NTPDarcyFluxInterface.h"

/**
 * Upwind advection of a functor q with the face mass flux m_f of NTPDarcyFluxInterface, the same
 * flux NTPDarcyMassFlux uses for the mass equation, so mass and energy see identical fluxes.
 * Used with 'advected_quantity = h' for the gas energy of a Darcy-Forchheimer flow.
 *
 * The term is assembled in convective form, rho u . grad(q) = div(rho u q) - q div(rho u): each
 * cell P receives  sum_f m_f (q_up - q_P).  At a converged flow (sum_f m_f = 0) this equals the
 * conservative divergence exactly, but during the iterations it keeps a positive diagonal,
 * sum over the inflow faces of |m_f|, even where the current pressure iterate does not conserve
 * mass. The conservative form has a zero diagonal in a cell with inflow and no outflow, which
 * makes the Newton Jacobian singular.
 *
 * On a boundary the flux is computed from the boundary pressure gradient. 'boundary_value' sets q
 * on the inflow boundaries listed in 'inflow_boundaries' (for example the inlet enthalpy).
 *
 * Limitations. On a boundary with an imposed mass flux (a Neumann pressure condition) the inflow
 * term uses the flux of the extrapolated boundary pressure gradient, not the imposed one; the
 * difference only multiplies (q_boundary - q_P), and the balances of the tutorial inputs close to
 * 1e-7. The boundary term is assembled only where the advected variable has a boundary condition,
 * and assumes the domain cell is the face's element.
 */
class NTPDarcyAdvection : public NTPDarcyFluxInterface
{
public:
  static InputParameters validParams();

  NTPDarcyAdvection(const InputParameters & params);

  void computeResidual(const FaceInfo & fi) override;
  void computeJacobian(const FaceInfo & fi) override;

protected:
  virtual ADReal computeQpResidual() override;

  /// Residuals of the current face on the element and on the neighbor side
  std::pair<ADReal, ADReal> sideResiduals();

  /// Advected quantity, e.g. specific enthalpy
  const Moose::Functor<ADReal> & _q;
  /// Advected quantity on the inflow boundaries
  const Moose::Functor<ADReal> * const _q_boundary;
  /// Boundaries on which the inflow value is prescribed
  std::set<BoundaryID> _inflow_boundaries;
};
