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

#include "FunctorMaterial.h"
#include "NTPDarcyFaceFlux.h"

/**
 * Cell superficial mass flux G = rho u [kg/(m^2 s)] of the Darcy-type flow of the Newton element
 * model, reconstructed from its face mass fluxes,
 *
 *   G_cell = (1 / V) sum_f m_f A_f (x_f - x_C),
 *
 * as an AD functor of the solution. The interphase heat transfer depends on G; through this functor
 * its dependence on the pressure and the temperatures enters the Newton Jacobian, where an
 * auxiliary variable would be a constant, lagged by one residual evaluation. A face of the domain
 * boundary carries the flux of the kernels where the pressure is imposed (Dirichlet), the imposed
 * flux on 'inflow_boundaries', and nothing elsewhere.
 */
class NTPDarcyMassFluxFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  NTPDarcyMassFluxFunctorMaterial(const InputParameters & parameters);

protected:
  /// The cell mass flux vector
  ADRealVectorValue cellMassFlux(const Elem & elem, const Moose::StateArg & state) const;

  /// The pressure variable
  const MooseVariableFV<Real> & _pressure;
  /// The face fluxes
  NTPDarcyFaceFlux _flux;
  /// Imposed inflow mass flux of the inflow boundaries
  std::map<BoundaryID, const Moose::Functor<ADReal> *> _inflow;
};
