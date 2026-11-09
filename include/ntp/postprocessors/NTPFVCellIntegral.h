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

#include "ElementPostprocessor.h"

/**
 * Volume integral of a functor with the finite-volume cell volumes: the sum over the cells of the
 * functor at the cell centroid times the cell volume of the finite-volume system (ElemInfo, times
 * the coordinate factor). This is the integral the linear finite-volume kernels assemble, e.g. the
 * heat a LinearFVSource deposits. The element quadrature of ElementIntegralFunctorPostprocessor
 * integrates the cell volume differently on cells with non-planar faces (a 3-D mesh mapped onto
 * cones), by a few 1e-4 here, and a balance against it does not close to round-off.
 */
class NTPFVCellIntegral : public ElementPostprocessor
{
public:
  static InputParameters validParams();

  NTPFVCellIntegral(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;
  virtual Real getValue() const override { return _integral; }

protected:
  /// Integrand
  const Moose::Functor<Real> & _functor;
  /// Integral over the cells of this rank, then over all ranks
  Real _integral;
};
