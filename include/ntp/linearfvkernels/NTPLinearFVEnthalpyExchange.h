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

#include "LinearFVElementalKernel.h"

/**
 * Volumetric heat exchange alpha (T_s - T) between a porous solid and a gas whose energy equation
 * is solved for the specific enthalpy h (linear finite volume, segregated solvers).
 *
 * The exchange depends on the gas temperature T(h); it is linearized about the current state
 * (h0, T0, cp0) with dT/dh = 1/cp:
 *
 *   alpha (T_s - T) ~ alpha (T_s - T0) - (alpha / cp0) (h - h0),
 *
 * which adds alpha / cp0 to the matrix diagonal of h and alpha (T_s - T0 + h0 / cp0) to the
 * right-hand side, both per unit volume. At convergence (h = h0) the source is exactly
 * alpha (T_s - T). LinearFVVolumetricHeatTransfer on the enthalpy variable would instead give
 * alpha (T_s - h), with inconsistent units.
 *
 * Partial elimination (optional). The segregated solvers solve the gas and the particle energy one
 * after the other. Where the exchange alpha V dominates the particle balance (a packed bed: the
 * particle conduction conductance G of a cell is about 1 % of alpha V), this Gauss-Seidel coupling
 * contracts by only ~1 - G / (alpha V) per iteration. The particles follow the gas:
 * dT_s = beta dT with beta = alpha V / (alpha V + G) for neighbours held fixed. Accounting for that
 * response in the gas equation,
 *
 *   alpha (T_s - T) ~ alpha (T_s - T0) - alpha (1 - beta) (T - T0),
 *
 * makes the coupled pair converge at the rate G / (alpha V) instead ('elimination_floor' bounds
 * 1 - beta from below, so that a cell with little flow keeps a gas diagonal). The source is still
 * exactly alpha (T_s - T) at convergence, whatever beta: G (from 'solid_conductivity', the
 * conductance of the particle conduction of the cell to its neighbours in the same blocks) only
 * sets the convergence rate.
 */
class NTPLinearFVEnthalpyExchange : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  NTPLinearFVEnthalpyExchange(const InputParameters & params);

  Real computeMatrixContribution() override;
  Real computeRightHandSideContribution() override;

protected:
  /// Volumetric exchange coefficient [W/(m^3 K)]
  const Moose::Functor<Real> & _alpha;
  /// Solid temperature
  const Moose::Functor<Real> & _T_solid;
  /// Gas temperature T(h)
  const Moose::Functor<Real> & _T_fluid;
  /// Gas isobaric specific heat
  const Moose::Functor<Real> & _cp;
  /// Particle conductivity (diagonal tensor) for the partial elimination, or null
  const Moose::Functor<RealVectorValue> * const _kappa_s;
  /// Smallest fraction 1 - beta kept on the gas diagonal
  const Real _floor;

  /// The fraction 1 - beta of the exchange coefficient on the gas diagonal (1 without elimination)
  Real diagonalFraction(const Moose::ElemArg & elem, const Moose::StateArg & state) const;
};
