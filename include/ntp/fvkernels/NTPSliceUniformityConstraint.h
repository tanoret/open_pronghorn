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

#include "FVElementalKernel.h"

class MeshDivision;

/**
 * Equations that make a field uniform across the slices of a mesh division, for a scalar unknown
 * with one component per slice (for example the logarithm of the orifice multiplier of each axial
 * slice of a fuel element, solved together with the flow).
 *
 * With <q>_j the volume average of the field q over slice j (within the kernel's blocks) and r the
 * reference slice, the scalar receives the residuals
 *
 *   R_j = w (<q>_j - <q>_r) / s + (1 - w + delta) lambda_j     for j != r,     R_r = lambda_r,
 *
 * so lambda_r = 0 fixes the scale and every other component adjusts until its slice average
 * matches the reference one. The scale s (e.g. the slice temperature rise) makes the equations
 * dimensionless, so their conditioning does not depend on the power level, and the small
 * regularization delta keeps lambda defined where the field does not depend on it (a uniform field
 * at zero power), where it pins lambda to 0. The activation w(t) switches the design on smoothly
 * once the field differs between slices: before the gas is heated the orifice cannot be seen at
 * all.
 *
 * With 'reference = mean' the slice averages are equalized with their volume average q_mean
 * instead, and the mean of the components is fixed:
 *
 *   R_j = w (<q>_j - q_mean) / s + (1 - w + delta) (lambda_j - lambda_mean)   for j < n - 1,
 *   R_(n-1) = (1 / n) sum_k lambda_k - lambda_mean.
 *
 * No slice is then pinned: a slice that would need its drag removed altogether (lambda -> -inf,
 * where the slice average no longer depends on lambda) with a fixed reference keeps a finite drag,
 * the level exp(lambda_mean) setting the margin.
 *
 * Nothing is added to the field's own equation. The field-to-scalar Jacobian comes from automatic
 * differentiation; the scalar must enter the field equations (here through the face mobility of
 * NTPDarcyFaceFlux, 'orifice_log_multiplier') for the system to be well posed.
 */
class NTPSliceUniformityConstraint : public FVElementalKernel
{
public:
  static InputParameters validParams();

  NTPSliceUniformityConstraint(const InputParameters & parameters);

  void initialSetup() override;

protected:
  void computeResidual() override;
  void computeJacobian() override;
  void computeResidualAndJacobian() override;
  ADReal computeQpResidual() override;

  /// Residual contributions of the current cell: (scalar component, value)
  std::vector<std::pair<unsigned int, ADReal>> cellContributions();

  /// Slices
  const MeshDivision & _division;
  /// The scalar unknown, one component per slice
  const MooseVariableScalar & _lambda_var;
  /// Reference slice
  /// Whether the slices are equalized with their mean (true) or with the reference slice
  const bool _mean_reference;
  const unsigned int _ref;
  /// Mean of the components with the mean reference
  const Real _lambda_mean;
  /// Field made uniform
  const Moose::Functor<ADReal> & _q;
  /// Regularization coefficient
  const Real _delta;
  /// Weight of the uniformity equations, w(t) in [0, 1]
  const Function & _activation;
  /// Scale of the field differences, fixed and optionally from a postprocessor
  const Real _scale;
  const PostprocessorValue * const _scale_pp;
  /// Volume of each slice within the blocks of this kernel
  std::vector<Real> _slice_volume;
  /// Volume of all slices
  Real _total_volume = 0;
};
