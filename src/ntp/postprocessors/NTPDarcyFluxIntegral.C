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

#include "NTPDarcyFluxIntegral.h"

registerMooseObject("OpenPronghornApp", NTPDarcyFluxIntegral);

InputParameters
NTPDarcyFluxIntegral::validParams()
{
  InputParameters params = SideIntegralVariablePostprocessor::validParams();
  params.addClassDescription("Boundary integral of the Darcy mass flux -(rho M) grad p . n, "
                             "optionally times a carried quantity such as the enthalpy.");
  params.addRequiredParam<MooseFunctorName>("rho_mobility", "Diagonal mass mobility rho M [s]");
  params.addParam<MooseFunctorName>("carried_quantity", 1, "Quantity carried by the mass flux");
  return params;
}

NTPDarcyFluxIntegral::NTPDarcyFluxIntegral(const InputParameters & parameters)
  : SideIntegralVariablePostprocessor(parameters),
    _rho_mobility(getFunctor<ADRealVectorValue>("rho_mobility")),
    _q(getFunctor<ADReal>("carried_quantity"))
{
  if (!_fv)
    paramError("variable", "Must be a finite volume pressure variable");
}

Real
NTPDarcyFluxIntegral::computeQpIntegral()
{
  mooseError("Only finite volume variables are supported");
}

Real
NTPDarcyFluxIntegral::computeFaceInfoIntegral(const FaceInfo * const fi)
{
  // Same face gradient (including the boundary-condition treatment) and the same one-sided
  // mobility as NTPDarcyAdvection on a boundary face, so this is the flux the equations see.
  const auto state = determineState();
  const auto grad_p = MetaPhysicL::raw_value(_fv_variable->adGradSln(*fi, state));
  const auto face = makeFace(*fi, Moose::FV::LimiterType::CentralDifference, true);
  const auto rhoM = MetaPhysicL::raw_value(_rho_mobility(face, state));
  Real mdot = 0;
  for (const auto i : make_range(Moose::dim))
    mdot -= rhoM(i) * grad_p(i) * _normals[_qp](i);
  // the carried quantity is taken from the cell inside the domain, which is the upwind value of
  // an outflow and the value NTPDarcyAdvection transports out of the domain
  const bool elem_inside = _fv_variable->hasBlocks(fi->elem().subdomain_id());
  const Elem * const inside = elem_inside ? fi->elemPtr() : fi->neighborPtr();
  const Moose::ElemArg cell{inside, false};
  return mdot * MetaPhysicL::raw_value(_q(cell, state));
}
