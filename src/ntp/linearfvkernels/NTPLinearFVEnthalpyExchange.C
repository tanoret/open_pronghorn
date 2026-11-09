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

#include "NTPLinearFVEnthalpyExchange.h"
#include "NS.h"
#include "MooseMesh.h"
#include "FaceInfo.h"

registerMooseObject("OpenPronghornApp", NTPLinearFVEnthalpyExchange);

InputParameters
NTPLinearFVEnthalpyExchange::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription("Volumetric solid-to-gas heat exchange alpha (T_s - T(h)) in a gas "
                             "energy equation solved for the specific enthalpy.");
  params.addRequiredParam<MooseFunctorName>("alpha", "Volumetric exchange coefficient [W/(m^3 K)]");
  params.addRequiredParam<MooseFunctorName>(NS::T_solid, "Solid temperature");
  params.addRequiredParam<MooseFunctorName>(NS::T_fluid, "Gas temperature T(h)");
  params.addRequiredParam<MooseFunctorName>(NS::cp, "Gas isobaric specific heat");
  params.addParam<MooseFunctorName>(
      "solid_conductivity",
      "Particle conductivity (diagonal tensor). If given, the gas equation accounts for the "
      "response of the particles of the cell (partial elimination), which makes the segregated "
      "gas-particle iterations converge much faster; the converged solution does not change.");
  params.addRangeCheckedParam<Real>(
      "elimination_floor",
      0.05,
      "elimination_floor > 0 & elimination_floor <= 1",
      "Smallest fraction 1 - beta of the exchange kept on the gas diagonal with the partial "
      "elimination. A cell with little flow keeps a diagonal this way; the convergence rate of "
      "the gas-particle iterations is about this fraction where the exchange dominates.");
  return params;
}

NTPLinearFVEnthalpyExchange::NTPLinearFVEnthalpyExchange(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _alpha(getFunctor<Real>("alpha")),
    _T_solid(getFunctor<Real>(NS::T_solid)),
    _T_fluid(getFunctor<Real>(NS::T_fluid)),
    _cp(getFunctor<Real>(NS::cp)),
    _kappa_s(isParamValid("solid_conductivity") ? &getFunctor<RealVectorValue>("solid_conductivity")
                                                : nullptr),
    _floor(getParam<Real>("elimination_floor"))
{
}

Real
NTPLinearFVEnthalpyExchange::diagonalFraction(const Moose::ElemArg & elem,
                                              const Moose::StateArg & state) const
{
  if (!_kappa_s)
    return 1;

  // conduction conductance G of the particles of this cell to its neighbours in the same blocks,
  // two-point, as the particle diffusion kernel assembles it
  const Elem * const cell = _current_elem_info->elem();
  const auto & mesh = _fe_problem.mesh();
  Real G = 0;
  for (const auto side : make_range(cell->n_sides()))
  {
    const Elem * const neighbor = cell->neighbor_ptr(side);
    if (!neighbor || !hasBlocks(neighbor->subdomain_id()))
      continue;
    const FaceInfo * fi = mesh.faceInfo(cell, side);
    if (!fi)
      fi = mesh.faceInfo(neighbor, neighbor->which_neighbor_am_i(cell));
    if (!fi)
      continue;
    const auto & n = fi->normal();
    const auto k_c = (*_kappa_s)(elem, state);
    const auto k_n = (*_kappa_s)(Moose::ElemArg{neighbor, false}, state);
    Real k_f = 0;
    for (const auto i : make_range(Moose::dim))
      k_f += n(i) * n(i) * 0.5 * (k_c(i) + k_n(i));
    G += k_f * fi->faceArea() * fi->faceCoord() / fi->dCNMag();
  }
  const Real alphaV = _alpha(elem, state) * _current_elem_volume;
  if (alphaV <= 0)
    return 1;
  // 1 - beta = G / (alpha V + G), kept above the floor
  return std::max(G / (alphaV + G), _floor);
}

Real
NTPLinearFVEnthalpyExchange::computeMatrixContribution()
{
  const auto elem = makeElemArg(_current_elem_info->elem());
  const auto state = determineState();
  return diagonalFraction(elem, state) * _alpha(elem, state) / _cp(elem, state) *
         _current_elem_volume;
}

Real
NTPLinearFVEnthalpyExchange::computeRightHandSideContribution()
{
  const auto elem = makeElemArg(_current_elem_info->elem());
  const auto state = determineState();
  const Real alpha = _alpha(elem, state);
  const Real cp = _cp(elem, state);
  const Real h = _var.getElemValue(*_current_elem_info, state);
  const Real T0 = _T_fluid(elem, state);
  // alpha (T_s - T0) - alpha (1 - beta) (T - T0), T - T0 = (h - h0) / cp: the h0 part goes here
  return alpha * (_T_solid(elem, state) - T0 + diagonalFraction(elem, state) * h / cp) *
         _current_elem_volume;
}
