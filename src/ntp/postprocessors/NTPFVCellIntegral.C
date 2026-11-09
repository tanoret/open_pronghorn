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

#include "NTPFVCellIntegral.h"
#include "MooseMesh.h"

registerMooseObject("OpenPronghornApp", NTPFVCellIntegral);

InputParameters
NTPFVCellIntegral::validParams()
{
  InputParameters params = ElementPostprocessor::validParams();
  params.addClassDescription(
      "Volume integral of a functor with the finite-volume cell volumes: the sum of the functor "
      "at the cell centroids times the cell volumes, as the linear finite-volume kernels assemble "
      "it.");
  params.addParam<MooseFunctorName>("functor", 1, "Integrand (1: the volume of the blocks)");
  return params;
}

NTPFVCellIntegral::NTPFVCellIntegral(const InputParameters & parameters)
  : ElementPostprocessor(parameters),
    _functor(getFunctor<Real>("functor")),
    _integral(0)
{
}

void
NTPFVCellIntegral::initialize()
{
  _integral = 0;
}

void
NTPFVCellIntegral::execute()
{
  const auto & info = _fe_problem.mesh().elemInfo(_current_elem->id());
  _integral += _functor(makeElemArg(_current_elem), determineState()) * info.volume() *
               info.coordFactor();
}

void
NTPFVCellIntegral::threadJoin(const UserObject & y)
{
  _integral += static_cast<const NTPFVCellIntegral &>(y)._integral;
}

void
NTPFVCellIntegral::finalize()
{
  gatherSum(_integral);
}
