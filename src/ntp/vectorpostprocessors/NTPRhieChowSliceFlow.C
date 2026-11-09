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

#include "NTPRhieChowSliceFlow.h"
#include "RhieChowMassFlux.h"
#include "MeshDivision.h"
#include "FEProblemBase.h"
#include "MooseMesh.h"

registerMooseObject("OpenPronghornApp", NTPRhieChowSliceFlow);

InputParameters
NTPRhieChowSliceFlow::validParams()
{
  InputParameters params = GeneralVectorPostprocessor::validParams();
  params.addClassDescription(
      "Mass flow rate through the interface between two blocks per division of a mesh division, "
      "from the Rhie-Chow face mass fluxes of the segregated solvers.");
  params.addRequiredParam<UserObjectName>("rhie_chow_user_object", "The Rhie-Chow object");
  params.addRequiredParam<MeshDivisionName>("mesh_division",
                                            "Divisions (e.g. axial slices) of the interface");
  params.addRequiredParam<std::vector<SubdomainName>>("from_block",
                                                      "Blocks the flow leaves (e.g. the plenum)");
  params.addRequiredParam<std::vector<SubdomainName>>(
      "to_block", "Blocks the flow enters (e.g. the cold frit); a face is attributed to the "
                  "division of its cell in these blocks");
  return params;
}

NTPRhieChowSliceFlow::NTPRhieChowSliceFlow(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _division(_fe_problem.getMeshDivision(getParam<MeshDivisionName>("mesh_division"))),
    _rc(nullptr),
    _flow(declareVector("flow")),
    _division_index(declareVector("division"))
{
  for (const auto id : _fe_problem.mesh().getSubdomainIDs(
           getParam<std::vector<SubdomainName>>("from_block")))
    _from_ids.insert(id);
  for (const auto id :
       _fe_problem.mesh().getSubdomainIDs(getParam<std::vector<SubdomainName>>("to_block")))
    _to_ids.insert(id);
  for (const auto id : _from_ids)
    if (_to_ids.count(id))
      paramError("to_block", "A block cannot be both in 'from_block' and 'to_block'");
}

void
NTPRhieChowSliceFlow::initialSetup()
{
  _rc = &getUserObject<RhieChowMassFlux>("rhie_chow_user_object");
}

void
NTPRhieChowSliceFlow::initialize()
{
  const auto n = _division.getNumDivisions();
  _flow.assign(n, 0);
  _division_index.resize(n);
  for (const auto j : make_range(n))
    _division_index[j] = j;
}

void
NTPRhieChowSliceFlow::execute()
{
  for (const auto * const fi : _fe_problem.mesh().faceInfo())
  {
    if (!fi->neighborPtr())
      continue;
    const auto elem_sub = fi->elem().subdomain_id();
    const auto neighbor_sub = fi->neighbor().subdomain_id();
    // the stored flux is along the face normal, from the face's element to its neighbor
    Real sign;
    const Elem * to_elem;
    if (_from_ids.count(elem_sub) && _to_ids.count(neighbor_sub))
    {
      sign = 1;
      to_elem = fi->neighborPtr();
    }
    else if (_to_ids.count(elem_sub) && _from_ids.count(neighbor_sub))
    {
      sign = -1;
      to_elem = fi->elemPtr();
    }
    else
      continue;
    // faces on processor boundaries are on both processors: count each on the owner of its
    // 'to_block' cell
    if (to_elem->processor_id() != processor_id())
      continue;
    const auto j = _division.divisionIndex(*to_elem);
    if (j == MooseMeshDivision::INVALID_DIVISION_INDEX || j >= _flow.size())
      mooseError(name(),
                 ": element ",
                 to_elem->id(),
                 " (centroid ",
                 to_elem->vertex_average(),
                 ") is outside the divisions of 'mesh_division'");
    _flow[j] += sign * _rc->getMassFlux(*fi) * fi->faceArea() * fi->faceCoord();
  }
}

void
NTPRhieChowSliceFlow::finalize()
{
  gatherSum(_flow);
}
