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

#include "CorrosionAdvection.h"

registerMooseObject("OpenPronghornApp", CorrosionAdvection);

InputParameters
CorrosionAdvection::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Non-conservative advection (u . grad c) of a dissolved species by a "
                             "velocity field read as functors.");
  params.addParam<MooseFunctorName>("vel_x", "0", "x-velocity functor [m/s].");
  params.addParam<MooseFunctorName>("vel_y", "0", "y-velocity functor [m/s].");
  params.addParam<MooseFunctorName>("vel_z", "0", "z-velocity functor [m/s].");
  return params;
}

CorrosionAdvection::CorrosionAdvection(const InputParameters & parameters)
  : ADKernel(parameters),
    _vel_x(getFunctor<ADReal>("vel_x")),
    _vel_y(getFunctor<ADReal>("vel_y")),
    _vel_z(getFunctor<ADReal>("vel_z"))
{
}

ADReal
CorrosionAdvection::computeQpResidual()
{
  const Moose::ElemQpArg qp_arg = {_current_elem, _qp, _qrule, _q_point[_qp]};
  const auto state = Moose::currentState();
  const VectorValue<ADReal> velocity(
      _vel_x(qp_arg, state), _vel_y(qp_arg, state), _vel_z(qp_arg, state));
  return _test[_i][_qp] * (velocity * _grad_u[_qp]);
}
