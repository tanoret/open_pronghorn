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

#include "LinearFVAdvectionDiffusionBC.h"

/**
 * Butler-Volmer electrode reaction as a linear finite-volume boundary condition, for corrosion or
 * plating of a wall in a flowing molten salt solved on the segregated (SIMPLE/PIMPLE) basis.
 *
 * With the concentration-explicit Butler-Volmer kinetics and a prescribed (metal minus salt)
 * overpotential, the molar flux of the cation into the salt is linear in the cell concentration,
 *   J_in = A - B c,
 *   A = (i0 / zF) exp( alpha_a z F eta / R T),     (anodic dissolution, metal activity = 1)
 *   B = (i0 / (zF c_ref)) exp(-alpha_c z F eta / R T),  (cathodic, dissolved-ion activity c/c_ref)
 * so it assembles directly into the linear system without lagging (the outward boundary flux is
 * B c - A). The exchange current is Arrhenius in the local temperature, so the wall corrodes faster
 * where the salt is hotter; the temperature, metal potential and salt potential are read as functors,
 * which lets the boundary couple one-way to a flow / energy / radiolysis solution.
 */
class CorrosionLinearFVButlerVolmerBC : public LinearFVAdvectionDiffusionBC
{
public:
  static InputParameters validParams();

  CorrosionLinearFVButlerVolmerBC(const InputParameters & parameters);

  virtual Real computeBoundaryValue() const override;
  virtual Real computeBoundaryNormalGradient() const override;
  virtual Real computeBoundaryValueMatrixContribution() const override;
  virtual Real computeBoundaryValueRHSContribution() const override;
  virtual Real computeBoundaryGradientMatrixContribution() const override;
  virtual Real computeBoundaryGradientRHSContribution() const override;
  /// The Butler-Volmer flux is a complete molar flux, not a diffusive D*grad(c) term.
  virtual bool includesMaterialPropertyMultiplier() const override { return true; }

protected:
  /// Compute the linear flux coefficients A (anodic source) and B (cathodic, multiplies c) at the
  /// current boundary face, including the Arrhenius temperature dependence of the exchange current.
  void computeFluxCoefficients(Real & A, Real & B) const;

  /// Concentration of the reacting cation in the boundary cell.
  Real cellConcentration() const;

  /// Charge number z of the cation.
  const Real _valence;
  /// Exchange current density i0 at the reference temperature [A/m^2].
  const Real _i0;
  /// Standard electrode potential E0 [V].
  const Real _E0;
  /// Anodic charge-transfer coefficient.
  const Real _alpha_a;
  /// Cathodic charge-transfer coefficient.
  const Real _alpha_c;
  /// Reference concentration at which i0 is defined [mol/m^3].
  const Real _c_ref;
  /// Reference temperature for the Arrhenius exchange current [K].
  const Real _reference_temperature;
  /// Activation energy of the exchange current [kJ/mol] (0 disables the temperature dependence).
  const Real _activation_energy;

  /// Local temperature functor [K].
  const Moose::Functor<Real> & _temperature;
  /// Metal-phase potential functor [V].
  const Moose::Functor<Real> & _metal_potential;
  /// Salt-phase potential functor [V].
  const Moose::Functor<Real> & _salt_potential;
};
