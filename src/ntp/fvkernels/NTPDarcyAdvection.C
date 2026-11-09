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

#include "NTPDarcyAdvection.h"

registerMooseObject("OpenPronghornApp", NTPDarcyAdvection);

InputParameters
NTPDarcyAdvection::validParams()
{
  InputParameters params = NTPDarcyFluxInterface::validParams();
  params.addClassDescription("Upwind advection of a functor with the Darcy-Forchheimer face mass "
                             "flux of a pressure field, in convective form.");
  params.addRequiredParam<MooseFunctorName>("advected_quantity",
                                            "The advected functor, e.g. the specific enthalpy");
  params.addParam<MooseFunctorName>("boundary_value",
                                    "Advected quantity entering through 'inflow_boundaries'");
  params.addParam<std::vector<BoundaryName>>(
      "inflow_boundaries", {}, "Boundaries on which 'boundary_value' is advected in");
  return params;
}

NTPDarcyAdvection::NTPDarcyAdvection(const InputParameters & params)
  : NTPDarcyFluxInterface(params),
    _q(getFunctor<ADReal>("advected_quantity")),
    _q_boundary(isParamValid("boundary_value") ? &getFunctor<ADReal>("boundary_value") : nullptr)
{
  const auto & inflow = getParam<std::vector<BoundaryName>>("inflow_boundaries");
  if (inflow.size() && !_q_boundary)
    paramError("boundary_value", "Required when 'inflow_boundaries' is set");
  for (const auto & id : _mesh.getBoundaryIDs(inflow))
    _inflow_boundaries.insert(id);
}

ADReal
NTPDarcyAdvection::computeQpResidual()
{
  mooseError("NTPDarcyAdvection assembles the element and neighbor sides separately");
}

std::pair<ADReal, ADReal>
NTPDarcyAdvection::sideResiduals()
{
  const auto state = determineState();
  const auto mdot = faceMassFlux(state);
  const auto area = _face_info->faceArea() * _face_info->faceCoord();
  const auto q_elem = _q(elemArg(), state);

  if (_var.isInternalFace(*_face_info))
  {
    const auto q_neighbor = _q(neighborArg(), state);
    const auto q_up = mdot > 0 ? q_elem : q_neighbor;
    // element: +m (q_up - q_E); neighbor (normal reversed): -m (q_up - q_N)
    return {area * mdot * (q_up - q_elem), -area * mdot * (q_up - q_neighbor)};
  }

  // Boundary face: outflow carries the cell value (no contribution in convective form); inflow on
  // a listed boundary brings the inlet value
  bool inflow_boundary = false;
  for (const auto id : _face_info->boundaryIDs())
    if (_inflow_boundaries.count(id))
      inflow_boundary = true;
  ADReal r = 0;
  if (mdot < 0 && inflow_boundary)
    r = area * mdot * ((*_q_boundary)(singleSidedFaceArg(), state) - q_elem);
  return {r, 0};
}

void
NTPDarcyAdvection::computeResidual(const FaceInfo & fi)
{
  if (skipForBoundary(fi))
    return;
  _face_info = &fi;
  _normal = fi.normal();
  _face_type = fi.faceType(std::make_pair(_var.number(), _var.sys().number()));
  const auto [r_elem, r_neighbor] = sideResiduals();

  if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    prepareVectorTag(_assembly, _var.number());
    _local_re(0) = MetaPhysicL::raw_value(r_elem);
    accumulateTaggedLocalResidual();
  }
  if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    prepareVectorTagNeighbor(_assembly, _var.number());
    _local_re(0) = MetaPhysicL::raw_value(r_neighbor);
    accumulateTaggedLocalResidual();
  }
}

void
NTPDarcyAdvection::computeJacobian(const FaceInfo & fi)
{
  if (skipForBoundary(fi))
    return;
  _face_info = &fi;
  _normal = fi.normal();
  _face_type = fi.faceType(std::make_pair(_var.number(), _var.sys().number()));
  const auto [r_elem, r_neighbor] = sideResiduals();

  if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    addResidualsAndJacobian(
        _assembly, std::array<ADReal, 1>{{r_elem}}, _var.dofIndices(), _var.scalingFactor());
  if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    addResidualsAndJacobian(_assembly,
                            std::array<ADReal, 1>{{r_neighbor}},
                            _var.dofIndicesNeighbor(),
                            _var.scalingFactor());
}
