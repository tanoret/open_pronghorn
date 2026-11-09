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

#include "AuxKernel.h"

class RhieChowMassFlux;

/**
 * One component of the cell superficial mass flux G = rho u [kg/(m^2 s)] of the segregated
 * solvers, reconstructed from the Rhie-Chow face mass fluxes the balances use,
 *
 *   G_cell = (1 / V) sum_f m_f A_f (x_f - x_C),
 *
 * as NTPDarcyMassFluxFunctorMaterial does for the Newton Darcy model. In porous layers one cell
 * thick, the cell velocity of the segregated solvers comes from a pressure gradient reconstructed
 * across porosity jumps and can be far from the flux through the layer; this is the flow for the
 * interphase heat transfer and the inertial drag.
 */
class NTPRhieChowMassFluxAux : public AuxKernel
{
public:
  static InputParameters validParams();

  NTPRhieChowMassFluxAux(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual Real computeValue() override;

  /// Component of the mass flux
  const unsigned int _component;
  /// The Rhie-Chow object holding the face fluxes
  const RhieChowMassFlux * _rc;
};
