//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVAdvectionDiffusionBC.h"

/**
 * Surface degassing (outgassing) boundary condition for a dissolved radiolysis species at a free
 * surface / cover-gas interface. The dissolved species leaves the melt at the convective
 * mass-transfer rate
 *
 *   J = k_surf * (C - kH * p_cover)
 *
 * where k_surf is the surface mass-transfer coefficient [m/s], kH the Henry coefficient
 * [mol/(m^3 Pa)] and p_cover the cover-gas partial pressure [Pa]. This is the mass-transfer analogue
 * of a convective heat-transfer boundary condition and represents direct escape of the volatile
 * species to a swept cover gas (the off-gas pathway), distinct from the volumetric bubble exchange
 * handled by MSRGasExchange. The flux already carries the full transfer coefficient, so the
 * diffusion kernel does not multiply it by the molecular diffusivity
 * (includesMaterialPropertyMultiplier returns true).
 *
 * The species must have a diffusion (dispersion) kernel for this boundary flux to be assembled.
 */
class MSRDegassingBC : public LinearFVAdvectionDiffusionBC
{
public:
  static InputParameters validParams();

  MSRDegassingBC(const InputParameters & params);

  virtual Real computeBoundaryValue() const override;
  virtual Real computeBoundaryNormalGradient() const override;
  virtual Real computeBoundaryValueMatrixContribution() const override;
  virtual Real computeBoundaryValueRHSContribution() const override;
  virtual Real computeBoundaryGradientMatrixContribution() const override;
  virtual Real computeBoundaryGradientRHSContribution() const override;

  /// The flux already includes the transfer coefficient, so do not multiply by the diffusivity.
  virtual bool includesMaterialPropertyMultiplier() const override { return true; }

protected:
  /// Equilibrium dissolved concentration C_eq = kH * p_cover at the current face
  Real equilibriumConcentration() const;

  /// Surface mass-transfer coefficient k_surf [m/s]
  const Moose::Functor<Real> & _k_surf;
  /// Henry coefficient kH [mol/(m^3 Pa)]
  const Real _kH;
  /// Cover-gas partial pressure p_cover [Pa]
  const Moose::Functor<Real> & _p_cover;
};
