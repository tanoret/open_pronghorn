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

#include "NTPDarcyFluxInterface.h"
#include "RelationshipManager.h"
#include "MeshDivision.h"
#include "MooseVariableScalar.h"
#include "FEProblemBase.h"

InputParameters
NTPDarcyFluxInterface::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params += FVDiffusionInterpolationInterface::validParams();
  params.addRequiredParam<VariableName>("pressure", "The pressure variable");
  params.addRequiredParam<MooseFunctorName>(
      "rho_mobility", "Diagonal mass mobility rho M [s], the coefficient of the pressure equation");
  params.addParam<bool>(
      "face_mobility",
      true,
      "Evaluate the mobility of each side at the face, with the face pressure gradient (true), "
      "or in the cells with the reconstructed cell gradient (false). The face evaluation is exact "
      "for the flow across a one-cell porous layer bounded by porosity jumps.");
  params.addParam<MooseFunctorName>(
      "momentum_flux",
      "Momentum flux J = rho u^2 [Pa] added to the driving pressure in 'momentum_blocks'; omit "
      "for pure Darcy-Forchheimer flow");
  params.addParam<std::vector<SubdomainName>>(
      "momentum_blocks", {}, "Subdomains (free channels) in which the momentum flux acts");
  params.addCoupledVar("orifice_log_multiplier",
                       "Scalar with one component per slice of 'mesh_division': the drag of the "
                       "cells of 'orifice_blocks' in slice j is multiplied by exp(component j)");
  params.addParam<MeshDivisionName>("mesh_division", "Slices of 'orifice_log_multiplier'");
  params.addParam<std::vector<SubdomainName>>(
      "orifice_blocks", {}, "Subdomains whose drag 'orifice_log_multiplier' scales");
  params.addParam<MooseEnum>(
      "coeff_interp_method",
      NTPDarcyFaceFlux::interpolationEnum(),
      "Face interpolation of the mobility: series resistances of the two half segments between "
      "the cell centers, or the MOOSE harmonic or linear average");
  params.set<unsigned short>("ghost_layers") = 2;
  params.addRelationshipManager(
      "ElementSideNeighborLayers",
      Moose::RelationshipManagerType::GEOMETRIC | Moose::RelationshipManagerType::ALGEBRAIC |
          Moose::RelationshipManagerType::COUPLING,
      [](const InputParameters & obj_params, InputParameters & rm_params)
      { FVRelationshipManagerInterface::setRMParamsDiffusion(obj_params, rm_params, 3); });
  return params;
}

NTPDarcyFluxInterface::NTPDarcyFluxInterface(const InputParameters & params)
  : FVFluxKernel(params),
    FVDiffusionInterpolationInterface(params),
    _flux(dynamic_cast<const MooseVariableFV<Real> &>(
              _subproblem.getVariable(_tid, getParam<VariableName>("pressure"))),
          getFunctor<ADRealVectorValue>("rho_mobility"),
          NTPDarcyFaceFlux::interpolation(getParam<MooseEnum>("coeff_interp_method")),
          getParam<bool>("face_mobility"))
{
  const bool has_division = isParamValid("mesh_division");
  const bool has_orifice = isCoupledScalar("orifice_log_multiplier");
  if (has_division != has_orifice)
    paramError("mesh_division", "'orifice_log_multiplier' and 'mesh_division' go together");
  if (has_orifice)
  {
    const auto & division =
        _fe_problem.getMeshDivision(getParam<MeshDivisionName>("mesh_division"));
    if (getScalarVar("orifice_log_multiplier", 0)->order() != division.getNumDivisions())
      paramError("orifice_log_multiplier", "Needs one component per slice of 'mesh_division'");
    const auto & blocks = getParam<std::vector<SubdomainName>>("orifice_blocks");
    if (blocks.empty())
      paramError("orifice_blocks", "Required with 'orifice_log_multiplier'");
    const auto ids = _mesh.getSubdomainIDs(blocks);
    _flux.setOrifice(division,
                     *getScalarVar("orifice_log_multiplier", 0),
                     std::set<SubdomainID>(ids.begin(), ids.end()));
  }
  if (isParamValid("momentum_flux"))
  {
    const auto & blocks = getParam<std::vector<SubdomainName>>("momentum_blocks");
    if (blocks.empty())
      paramError("momentum_blocks", "Required with 'momentum_flux'");
    const auto ids = _mesh.getSubdomainIDs(blocks);
    _flux.setMomentumFlux(getFunctor<ADReal>("momentum_flux"),
                          std::set<SubdomainID>(ids.begin(), ids.end()));
  }
}

ADReal
NTPDarcyFluxInterface::faceMassFlux(const Moose::StateArg & state) const
{
  return _flux.massFlux(*_face_info, state);
}
