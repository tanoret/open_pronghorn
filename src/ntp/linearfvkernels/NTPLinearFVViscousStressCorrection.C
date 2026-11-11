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

#include "NTPLinearFVViscousStressCorrection.h"
#include "MooseLinearVariableFV.h"
#include "SubProblem.h"
#include "NS.h"

registerMooseObject("OpenPronghornApp", NTPLinearFVViscousStressCorrection);

InputParameters
NTPLinearFVViscousStressCorrection::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription(
      "Transpose and dilatation part of the compressible Newtonian viscous stress, "
      "mu (grad u^T - 2/3 div(u) I) . n, on the internal faces of the kernel's blocks (explicit).");
  params.addRequiredParam<MooseFunctorName>(NS::mu, "Dynamic viscosity (effective if turbulent)");
  MooseEnum component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>("momentum_component", component, "Momentum component");
  params.addRequiredParam<SolverVariableName>("u", "Velocity in the x direction");
  params.addParam<SolverVariableName>("v", "Velocity in the y direction");
  params.addParam<SolverVariableName>("w", "Velocity in the z direction");
  params.addParam<std::vector<MooseFunctorName>>(
      "gradient_functors",
      "Functors, one per velocity component, whose cell gradients form the stress (e.g. the "
      "'<component>_channel' functors of NTPChannelVelocity); by default the velocity's own cell "
      "gradients");
  return params;
}

NTPLinearFVViscousStressCorrection::NTPLinearFVViscousStressCorrection(
    const InputParameters & params)
  : LinearFVFluxKernel(params),
    _mu(getFunctor<Real>(NS::mu)),
    _index(getParam<MooseEnum>("momentum_component")),
    _dim(_subproblem.mesh().dimension()),
    _coord_type(getBlockCoordSystem()),
    _rz_radial_coord(_coord_type == Moose::COORD_RZ ? _subproblem.getAxisymmetricRadialCoord() : 0),
    _face_stress(0)
{
  const std::vector<std::string> names = {"u", "v", "w"};
  for (const auto d : make_range(_dim))
  {
    if (!isParamValid(names[d]))
      paramError(names[d], "Needed in ", _dim, " dimensions");
    auto & var = _fe_problem.getVariable(_tid, getParam<SolverVariableName>(names[d]));
    auto * ptr = dynamic_cast<MooseLinearVariableFVReal *>(&var);
    if (!ptr)
      paramError(names[d], "Must be a linear finite-volume variable");
    ptr->computeCellGradients();
    _vel.push_back(ptr);
  }
  if (_index >= _dim)
    paramError("momentum_component", "Not a component of the mesh dimension");
  if (isParamValid("gradient_functors"))
  {
    const auto & names = getParam<std::vector<MooseFunctorName>>("gradient_functors");
    if (names.size() != _dim)
      paramError("gradient_functors", "Give one functor per velocity component");
    for (const auto & n : names)
      _grad_functors.push_back(&getFunctor<Real>(n));
  }
}

bool
NTPLinearFVViscousStressCorrection::hasFaceSide(const FaceInfo & fi, const bool fi_elem_side) const
{
  // only faces with both cells in the blocks: no stress across a porosity jump or a boundary
  if (!fi.neighborPtr() || !hasBlocks(fi.elemPtr()->subdomain_id()) ||
      !hasBlocks(fi.neighborPtr()->subdomain_id()))
    return false;
  return LinearFVFluxKernel::hasFaceSide(fi, fi_elem_side);
}

Real
NTPLinearFVViscousStressCorrection::faceStress()
{
  if (_cached_rhs_contribution)
    return _face_stress;

  const auto state = determineState();
  const auto & fi = *_current_face_info;
  const Real gc = fi.gC();

  // d_i u_j of both cells and the divergence
  RealVectorValue col_elem, col_neighbor;
  Real div_elem = 0, div_neighbor = 0;
  for (const auto j : make_range(_dim))
  {
    const auto grad_elem = _grad_functors.empty()
                               ? _vel[j]->gradSln(*fi.elemInfo(), state)
                               : _grad_functors[j]->gradient(makeElemArg(fi.elemPtr()), state);
    const auto grad_neighbor =
        _grad_functors.empty() ? _vel[j]->gradSln(*fi.neighborInfo(), state)
                               : _grad_functors[j]->gradient(makeElemArg(fi.neighborPtr()), state);
    col_elem(j) = grad_elem(_index);
    col_neighbor(j) = grad_neighbor(_index);
    div_elem += grad_elem(j);
    div_neighbor += grad_neighbor(j);
  }
  if (_coord_type == Moose::COORD_RZ)
  {
    const auto & ur = *_vel[_rz_radial_coord];
    div_elem += ur.getElemValue(*fi.elemInfo(), state) / fi.elemInfo()->centroid()(_rz_radial_coord);
    div_neighbor +=
        ur.getElemValue(*fi.neighborInfo(), state) / fi.neighborInfo()->centroid()(_rz_radial_coord);
  }
  col_elem(_index) -= 2. / 3 * div_elem;
  col_neighbor(_index) -= 2. / 3 * div_neighbor;

  const auto face_vector = gc * col_elem + (1 - gc) * col_neighbor;
  // both cells are in the blocks: the face viscosity is the linear average of the two
  const Real mu_face = gc * _mu(makeElemArg(fi.elemPtr()), state) +
                       (1 - gc) * _mu(makeElemArg(fi.neighborPtr()), state);
  _face_stress = mu_face * (face_vector * fi.normal());
  _cached_rhs_contribution = true;
  return _face_stress;
}

Real
NTPLinearFVViscousStressCorrection::computeElemRightHandSideContribution()
{
  return faceStress() * _current_face_area;
}

Real
NTPLinearFVViscousStressCorrection::computeNeighborRightHandSideContribution()
{
  return -faceStress() * _current_face_area;
}
