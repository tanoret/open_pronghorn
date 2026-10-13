//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* OpenPronghorn is powered by the MOOSE Framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVAdvectionUtils.h"
#include "FEProblemBase.h"
#include "Factory.h"
#include "MooseEnum.h"

namespace OpenPronghorn
{
void
addAdvectedInterpolationParams(InputParameters & params)
{
  MooseEnum interp("average upwind", "upwind");
  params.addParam<MooseEnum>(
      "advected_interp_method",
      interp,
      "Shorthand for the advected interpolation method of the species. An FVInterpolationMethod "
      "of the matching type is created automatically unless 'advected_interp_method_name' is "
      "supplied.");
  params.addParam<InterpolationMethodName>(
      "advected_interp_method_name",
      "Name of an existing [FVInterpolationMethods] object used for the advected species "
      "(overrides 'advected_interp_method').");
}

InterpolationMethodName
advectedInterpolationMethodName(const InputParameters & params)
{
  if (params.isParamValid("advected_interp_method_name"))
    return params.get<InterpolationMethodName>("advected_interp_method_name");
  return InterpolationMethodName(std::string(params.get<MooseEnum>("advected_interp_method")));
}

InterpolationMethodName
addAdvectedInterpolationMethod(const InputParameters & params,
                               FEProblemBase & problem,
                               Factory & factory)
{
  const auto name = advectedInterpolationMethodName(params);
  if (params.isParamValid("advected_interp_method_name") || problem.hasFVInterpolationMethod(name))
    return name;

  const std::string type = (name == "average") ? "FVGeometricAverage" : "FVAdvectedUpwind";
  auto method_params = factory.getValidParams(type);
  problem.addFVInterpolationMethod(type, name, method_params);
  return name;
}
}
