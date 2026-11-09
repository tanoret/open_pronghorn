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

#include "SideIntegralVariablePostprocessor.h"

/**
 * Integral over a boundary of the Darcy mass flux -(rho M) grad p . n times an optional functor q:
 * with q = 1 the mass flow through the boundary [kg/s], with q = h the enthalpy flow [W].
 * The mass flux uses the same face gradient and mobility as NTPDarcyAdvection, and q is taken from
 * the cell inside the domain, the upwind value of an outflow, so the result is exactly the flow the
 * discrete conservation equations see leave the domain. 'variable' is the pressure.
 *
 * Use it on boundaries where the pressure is prescribed (an outlet). On a boundary with an
 * imposed mass flux the reconstructed boundary gradient is not that flux, and the imposed value
 * itself should be reported instead.
 */
class NTPDarcyFluxIntegral : public SideIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  NTPDarcyFluxIntegral(const InputParameters & parameters);

protected:
  Real computeQpIntegral() override;
  Real computeFaceInfoIntegral(const FaceInfo * fi) override;

  /// Diagonal mass mobility rho M [s]
  const Moose::Functor<ADRealVectorValue> & _rho_mobility;
  /// Quantity carried by the flux
  const Moose::Functor<ADReal> & _q;
};
