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

#include "NTPChannelVelocity.h"
#include "FEProblemBase.h"
#include "MooseMesh.h"

registerMooseObject("OpenPronghornApp", NTPChannelVelocity);

InputParameters
NTPChannelVelocity::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Adds, for each velocity component, the functor '<component>_channel' with a cell gradient "
      "that ignores the porous neighbours of the free-channel blocks.");
  params.addRequiredParam<std::vector<MooseFunctorName>>("velocity", "Velocity components");
  params.addRequiredParam<std::vector<SubdomainName>>("channel_blocks",
                                                      "Free-channel blocks (plenum, chamber)");
  params.addParam<std::string>("suffix", "_channel", "Suffix of the added functor names");
  return params;
}

NTPChannelVelocity::NTPChannelVelocity(const InputParameters & parameters)
  : GeneralUserObject(parameters)
{
  const auto & mesh = _fe_problem.mesh();
  const auto ids_vector = mesh.getSubdomainIDs(getParam<std::vector<SubdomainName>>("channel_blocks"));
  const std::set<SubdomainID> ids(ids_vector.begin(), ids_vector.end());
  const auto & suffix = getParam<std::string>("suffix");
  for (const auto & name : getParam<std::vector<MooseFunctorName>>("velocity"))
    for (const auto tid : make_range(libMesh::n_threads()))
    {
      const auto & velocity = _fe_problem.getFunctor<Real>(name, tid, this->name(), true);
      _functors.push_back(
          std::make_unique<NTPChannelGradientFunctor>(name + suffix, velocity, mesh, ids));
      _fe_problem.addFunctor(name + suffix, *_functors.back(), tid);
    }
}
