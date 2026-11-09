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

#include "NTPDarcyMassFluxFunctorMaterial.h"
#include "FaceInfo.h"
#include "MooseMesh.h"
#include "MeshDivision.h"
#include "MooseVariableScalar.h"
#include "FEProblemBase.h"

registerMooseObject("OpenPronghornApp", NTPDarcyMassFluxFunctorMaterial);

InputParameters
NTPDarcyMassFluxFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription("Cell superficial mass flux of a Darcy-type flow, reconstructed from "
                             "its face mass fluxes, as an AD functor.");
  params.addRequiredParam<VariableName>("pressure", "The pressure variable");
  params.addRequiredParam<MooseFunctorName>("rho_mobility",
                                            "Diagonal mass mobility rho M [s], as in the kernels");
  params.addParam<bool>("face_mobility", true, "As in NTPDarcyMassFlux");
  params.addParam<MooseEnum>(
      "coeff_interp_method", NTPDarcyFaceFlux::interpolationEnum(), "As in NTPDarcyMassFlux");
  params.addCoupledVar("orifice_log_multiplier", "As in NTPDarcyMassFlux");
  params.addParam<MeshDivisionName>("mesh_division", "As in NTPDarcyMassFlux");
  params.addParam<std::vector<SubdomainName>>("orifice_blocks", {}, "As in NTPDarcyMassFlux");
  params.addParam<std::vector<BoundaryName>>(
      "inflow_boundaries", {}, "Boundaries with an imposed inflow mass flux (Neumann condition)");
  params.addParam<std::vector<MooseFunctorName>>(
      "inflow_mass_fluxes",
      {},
      "Imposed inflow mass flux of each inflow boundary [kg/(m^2 s)] (functors: numbers, "
      "postprocessors, ...)");
  params.addParam<MooseFunctorName>("mass_flux_x_name", "G_x", "Name of the radial component");
  params.addParam<MooseFunctorName>("mass_flux_y_name", "G_y", "Name of the axial component");
  params.addParam<MooseFunctorName>("mass_flux_z_name", "G_z", "Name of the third component");
  return params;
}

NTPDarcyMassFluxFunctorMaterial::NTPDarcyMassFluxFunctorMaterial(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _pressure(dynamic_cast<const MooseVariableFV<Real> &>(
        _subproblem.getVariable(_tid, getParam<VariableName>("pressure")))),
    _flux(_pressure,
          getFunctor<ADRealVectorValue>("rho_mobility"),
          NTPDarcyFaceFlux::interpolation(getParam<MooseEnum>("coeff_interp_method")),
          getParam<bool>("face_mobility"))
{
  const auto & inflow = getParam<std::vector<BoundaryName>>("inflow_boundaries");
  const auto & fluxes = getParam<std::vector<MooseFunctorName>>("inflow_mass_fluxes");
  if (inflow.size() != fluxes.size())
    paramError("inflow_mass_fluxes", "One mass flux per inflow boundary");
  for (const auto i : index_range(inflow))
    _inflow[_mesh.getBoundaryID(inflow[i])] = &getFunctor<ADReal>(fluxes[i]);

  if (isCoupledScalar("orifice_log_multiplier"))
  {
    if (!isParamValid("mesh_division"))
      paramError("mesh_division", "Required with 'orifice_log_multiplier'");
    const auto ids = _mesh.getSubdomainIDs(getParam<std::vector<SubdomainName>>("orifice_blocks"));
    _flux.setOrifice(_fe_problem.getMeshDivision(getParam<MeshDivisionName>("mesh_division")),
                     *getScalarVar("orifice_log_multiplier", 0),
                     std::set<SubdomainID>(ids.begin(), ids.end()));
  }

  const std::vector<std::string> names = {
      "mass_flux_x_name", "mass_flux_y_name", "mass_flux_z_name"};
  for (const auto component : make_range(3u))
    addFunctorProperty<ADReal>(
        getParam<MooseFunctorName>(names[component]),
        [this, component](const auto & r, const auto & t) -> ADReal
        {
          const Elem * elem = nullptr;
          if constexpr (std::is_same_v<std::decay_t<decltype(r)>, Moose::ElemArg>)
            elem = r.elem;
          else if constexpr (std::is_same_v<std::decay_t<decltype(r)>, Moose::FaceArg>)
            elem = r.face_side ? r.face_side : r.fi->elemPtr();
          else
            mooseError("NTPDarcyMassFluxFunctorMaterial is evaluated on cells (or faces) only");
          return cellMassFlux(*elem, t)(component);
        });
}

ADRealVectorValue
NTPDarcyMassFluxFunctorMaterial::cellMassFlux(const Elem & elem,
                                              const Moose::StateArg & state) const
{
  const auto & elem_info = _mesh.elemInfo(elem.id());
  const Point & xc = elem_info.centroid();

  ADRealVectorValue sum;
  for (const auto side : make_range(elem.n_sides()))
  {
    // a face is stored on one of its two cells; its normal points away from that cell
    const FaceInfo * face = _mesh.faceInfo(&elem, side);
    bool outward = true;
    const Elem * const neighbor = elem.neighbor_ptr(side);
    if (!face)
    {
      if (!neighbor)
        continue;
      face = _mesh.faceInfo(neighbor, neighbor->which_neighbor_am_i(&elem));
      outward = false;
    }
    if (!face)
      continue;
    // a domain boundary carries the kernel flux where the pressure is imposed, the imposed flux of
    // an inflow boundary, and nothing elsewhere
    const Real area = face->faceArea() * face->faceCoord();
    ADReal out_flux;
    if (_pressure.isInternalFace(*face))
    {
      const ADReal mdot = _flux.massFlux(*face, state);
      out_flux = outward ? mdot : -mdot;
    }
    else
    {
      const auto inflow = std::find_if(face->boundaryIDs().begin(),
                                       face->boundaryIDs().end(),
                                       [this](const BoundaryID id) { return _inflow.count(id); });
      if (inflow != face->boundaryIDs().end())
        out_flux = -(*_inflow.at(*inflow))(Moose::ElemArg{&elem, false}, state);
      else if (_pressure.isDirichletBoundaryFace(*face, &elem, state))
      {
        // the flux is along the stored normal, which points away from the face's own element
        const ADReal mdot = _flux.massFlux(*face, state);
        out_flux = outward ? mdot : -mdot;
      }
      else
        continue;
    }
    for (const auto i : make_range(Moose::dim))
      sum(i) += out_flux * area * (face->faceCentroid()(i) - xc(i));
  }
  return sum / (elem_info.volume() * elem_info.coordFactor());
}
