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

/**
 * Linear finite-volume kernel for gas-liquid exchange of a dissolved diatomic (Cl2, F2) between the
 * melt and a well-mixed headspace, following the lumped two-volume model of MoltenSaltRadiolysis
 * (core.py::GasExchange):
 *
 *   R = kLa * (C_liq - kH * p_gas),     p_gas = C_gas * R_gas * T
 *
 * The gas inventory is carried as a gas-phase concentration C_gas [mol/m^3 of headspace], so that
 * n_gas = C_gas * V_gas and the ideal-gas pressure reduces to the purely local form
 * p_gas = C_gas * R_gas * T (no postprocessor needed). The dissolved and gas equations are
 *
 *   dC_liq/dt = -R
 *   dC_gas/dt = (V_liq / V_gas) * R
 *
 * which conserve total moles. Each instance contributes one of the two equations (selected by
 * "mode"). The variable's own term is implicit (positive diagonal); the partner concentration is
 * lagged on the right-hand side, matching the segregated fixed-point linearization.
 */
class MSRGasExchange : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  MSRGasExchange(const InputParameters & params);

  virtual Real computeMatrixContribution() override;

  virtual Real computeRightHandSideContribution() override;

protected:
  /// Which side of the exchange this instance represents
  enum class Mode
  {
    Liquid,
    Gas
  };

  /// Selected mode (liquid or gas equation)
  const Mode _mode;

  /// Temperature functor [K]
  const Moose::Functor<Real> & _temperature;

  /// Overall mass-transfer coefficient kLa [1/s]
  const Real _kLa;

  /// Henry coefficient kH [mol/(m^3 Pa)]
  const Real _kH;

  /// Ratio of liquid to gas (headspace) volume V_liq / V_gas (used in gas mode)
  const Real _volume_ratio;

  /// Partner concentration functor (gas concentration in liquid mode, liquid concentration in gas mode)
  const Moose::Functor<Real> & _partner;
};
