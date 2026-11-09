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

#include "NTPRhieChowMassFluxAux.h"
#include "RhieChowMassFlux.h"
#include "FaceInfo.h"
#include "MooseMesh.h"

registerMooseObject("OpenPronghornApp", NTPRhieChowMassFluxAux);

InputParameters
NTPRhieChowMassFluxAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Cell superficial mass flux of the segregated solvers, reconstructed "
                             "from the Rhie-Chow face mass fluxes.");
  params.addRequiredParam<UserObjectName>("rhie_chow_user_object", "The Rhie-Chow object");
  params.addRequiredParam<MooseEnum>(
      "component", MooseEnum("x=0 y=1 z=2"), "Component of the mass flux");
  return params;
}

NTPRhieChowMassFluxAux::NTPRhieChowMassFluxAux(const InputParameters & parameters)
  : AuxKernel(parameters), _component(getParam<MooseEnum>("component")), _rc(nullptr)
{
  if (isNodal())
    paramError("variable", "Must be an elemental variable");
}

void
NTPRhieChowMassFluxAux::initialSetup()
{
  AuxKernel::initialSetup();
  _rc = &getUserObject<RhieChowMassFlux>("rhie_chow_user_object");
}

Real
NTPRhieChowMassFluxAux::computeValue()
{
  const auto & elem_info = _mesh.elemInfo(_current_elem->id());
  const Point & xc = elem_info.centroid();

  Real sum = 0;
  for (const auto side : make_range(_current_elem->n_sides()))
  {
    // a face is stored on one of its two cells; its normal points away from that cell
    const FaceInfo * face = _mesh.faceInfo(_current_elem, side);
    bool outward = true;
    if (!face)
    {
      const Elem * const neighbor = _current_elem->neighbor_ptr(side);
      if (!neighbor)
        continue;
      face = _mesh.faceInfo(neighbor, neighbor->which_neighbor_am_i(_current_elem));
      outward = false;
    }
    if (!face)
      continue;
    const Real mdot = _rc->getMassFlux(*face);
    const Real area = face->faceArea() * face->faceCoord();
    sum += (outward ? mdot : -mdot) * area * (face->faceCentroid()(_component) - xc(_component));
  }
  return sum / (elem_info.volume() * elem_info.coordFactor());
}
