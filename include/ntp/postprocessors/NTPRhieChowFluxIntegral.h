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

#include "SideIntegralPostprocessor.h"

class RhieChowMassFlux;

/**
 * Boundary integral of the face mass fluxes of the segregated solvers' Rhie-Chow object, positive
 * out of the domain, optionally times a quantity carried by the flux, evaluated at the boundary
 * face from the inside cell: for a variable, its boundary value (the Dirichlet value of an inflow,
 * the extrapolated value of an outflow), which is what the advection kernels transport. These are
 * the fluxes the discrete balances see, so the mass flow rate and, with the enthalpy as carried
 * quantity, the enthalpy flow rate of a boundary are exact.
 *
 * With 'normal_kinetic_energy' the carried quantity is the kinetic energy of the face velocity
 * normal to the boundary, u_n^2 / 2 with u_n = mdot_f / rho_f, from the exact face flux and the
 * density 'rho' at the face: the kinetic energy flow through an outlet. The cell velocities of the
 * segregated solvers next to a porosity jump are not a reliable alternative.
 */
class NTPRhieChowFluxIntegral : public SideIntegralPostprocessor
{
public:
  static InputParameters validParams();

  NTPRhieChowFluxIntegral(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual Real computeQpIntegral() override;
  virtual Real computeFaceInfoIntegral(const FaceInfo * fi) override;

  /// The Rhie-Chow object holding the face fluxes
  const RhieChowMassFlux * _rc;
  /// Quantity carried by the mass flux
  const Moose::Functor<Real> & _q;
  /// Carry the kinetic energy of the normal face velocity instead, with this density
  const bool _normal_ke;
  const Moose::Functor<Real> * const _rho;
};
