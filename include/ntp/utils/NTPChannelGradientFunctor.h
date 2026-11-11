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

#include "MooseFunctor.h"
#include "MooseMesh.h"

#include <set>

/**
 * A velocity component of a free channel next to porous regions, for closures that need its strain
 * rate (turbulence production, viscous stress): the values are those of the wrapped functor, the
 * cell gradient is a Green-Gauss reconstruction over the cell's faces in which a face towards a
 * cell outside the channel blocks takes the cell's own value (no normal gradient through a porous
 * interface). The sum is (1/V) sum_f (phi_f - phi_C) A_f n_f, exact for a constant field in
 * Cartesian and RZ coordinates. The superficial velocity jumps by orders of magnitude between a free channel and a
 * packed bed; that jump is not a strain rate of the channel flow. Boundary faces use the wrapped
 * functor's face value, internal channel faces the linear interpolation of the two cells.
 */
class NTPChannelGradientFunctor : public Moose::FunctorBase<Real>
{
public:
  NTPChannelGradientFunctor(const MooseFunctorName & name,
                            const Moose::Functor<Real> & velocity,
                            const MooseMesh & mesh,
                            const std::set<SubdomainID> & channel_blocks);

  virtual bool hasBlocks(const SubdomainID id) const override { return _velocity.hasBlocks(id); }
  virtual bool supportsFaceArg() const override { return true; }
  virtual bool supportsElemSideQpArg() const override { return true; }

  /// Green-Gauss gradient, unit-testable: outward area vectors with the face value minus the
  /// cell value, and the cell volume
  static RealVectorValue
  greenGauss(const std::vector<std::pair<Point, Real>> & area_vectors_and_values, Real volume);

protected:
  using typename Moose::FunctorBase<Real>::ValueType;
  using typename Moose::FunctorBase<Real>::GradientType;

  ValueType evaluate(const Moose::ElemArg & a, const Moose::StateArg & s) const override
  {
    return _velocity(a, s);
  }
  ValueType evaluate(const Moose::FaceArg & a, const Moose::StateArg & s) const override
  {
    return _velocity(a, s);
  }
  ValueType evaluate(const Moose::ElemQpArg & a, const Moose::StateArg & s) const override
  {
    return _velocity(a, s);
  }
  ValueType evaluate(const Moose::ElemSideQpArg & a, const Moose::StateArg & s) const override
  {
    return _velocity(a, s);
  }
  ValueType evaluate(const Moose::ElemPointArg & a, const Moose::StateArg & s) const override
  {
    return _velocity(a, s);
  }
  ValueType evaluate(const Moose::NodeArg & a, const Moose::StateArg & s) const override
  {
    return _velocity(a, s);
  }

  GradientType evaluateGradient(const Moose::ElemArg & elem_arg,
                                const Moose::StateArg & state) const override;

  const Moose::Functor<Real> & _velocity;
  const MooseMesh & _mesh;
  const std::set<SubdomainID> _channel_blocks;
};
