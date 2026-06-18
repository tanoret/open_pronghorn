//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVElementalKernel.h"
#include "MoltenSaltRadiolysisData.h"

/**
 * Linear finite-volume kernel for a single mass-action reaction's contribution to one radiolysis
 * species. The reaction rate is r = k(T) * product over reactants of [S]^order, with the forward
 * rate constant k(T) following the Arrhenius forms in the molten salt radiolysis database.
 *
 * A species in a given reaction is either consumed (a reactant) or produced (a product); the two
 * modes are split following the linear finite-volume convention used by LinearFVReaction and
 * LinearFVSource:
 *
 *  - Consumption: an implicit diagonal sink is added to the matrix,
 *      coeff = stoich * k(T) * [self]^(self_order - 1) * product(partner reactants, lagged)
 *    so that coeff * [self] reproduces the consumption rate of the species. Treating one power of
 *    the variable implicitly keeps the diagonal positive (segregated fixed-point convergence).
 *
 *  - Production: an explicit source is added to the right hand side,
 *      source = stoich * k(T) * product(all reactants, lagged)
 *
 * Coupled (partner) species concentrations are evaluated as functors at the current element. Within
 * the segregated fixed-point loop these return the latest available iterate, which is the standard
 * convergent linearization for mass-action coupling. Lagged concentrations are clamped at zero
 * inside the rate products to keep the coefficients well behaved.
 */
class MSRReaction : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  MSRReaction(const InputParameters & params);

  /// Volumetric contribution to the system matrix (implicit consumption sink)
  virtual Real computeMatrixContribution() override;

  /// Volumetric contribution to the right-hand side (explicit production source)
  virtual Real computeRightHandSideContribution() override;

protected:
  /// Whether the species is consumed (reactant) or produced (product) in this reaction
  enum class Mode
  {
    Consumption,
    Production
  };

  /// Evaluate the partner-reactant mass-action product at the current element and state
  Real partnerProduct(const Moose::ElemArg & elem_arg, const Moose::StateArg & state) const;

  /// Temperature functor [K]
  const Moose::Functor<Real> & _temperature;

  /// Selected mode (consumption or production)
  const Mode _mode;

  /// Net stoichiometric coefficient of the species in this reaction (always positive)
  const Real _stoichiometric_coefficient;

  /// Kinetic order of the species itself (used only in consumption mode)
  const Real _self_order;

  /// First partner reactant functor (may be null)
  const Moose::Functor<Real> * const _partner1;
  /// Kinetic order of the first partner reactant
  const Real _partner1_order;

  /// Second partner reactant functor (may be null)
  const Moose::Functor<Real> * const _partner2;
  /// Kinetic order of the second partner reactant
  const Real _partner2_order;

  /// Reaction kinetics holder, used to evaluate the Arrhenius rate constant k(T)
  MSR::ReactionData _kinetics;
};
