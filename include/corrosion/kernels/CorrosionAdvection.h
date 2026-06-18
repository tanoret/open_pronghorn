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

#include "ADKernel.h"

/**
 * Advection of a dissolved species by a prescribed velocity field, in non-conservative form
 * (test, u . grad c). The velocity components are read as functors, so the field can be supplied by
 * a separate flow solve (e.g. the segregated Navier-Stokes velocity vel_x/vel_y) as a one-way
 * coupling. For strictly mass-conserving transport use the linear finite-volume advection objects of
 * the flow framework; this kernel is intended for the nonlinear Newton corrosion system where the
 * advective term is a secondary transport contribution.
 */
class CorrosionAdvection : public ADKernel
{
public:
  static InputParameters validParams();

  CorrosionAdvection(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// Velocity component functors [m/s].
  const Moose::Functor<ADReal> & _vel_x;
  const Moose::Functor<ADReal> & _vel_y;
  const Moose::Functor<ADReal> & _vel_z;
};
