//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* OpenPronghorn is powered by the MOOSE Framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InputParameters.h"
#include "MooseTypes.h"

class FEProblemBase;
class Factory;

namespace OpenPronghorn
{
/**
 * Add the advected-interpolation parameters shared by the chemistry actions that create linear
 * finite-volume advection kernels: a shorthand 'advected_interp_method' enum and an optional
 * 'advected_interp_method_name' that selects an existing [FVInterpolationMethods] object.
 */
void addAdvectedInterpolationParams(InputParameters & params);

/**
 * Return the FVInterpolationMethod name the action's advection kernels should use. When the user
 * did not name one, an object of the enum's type is created (once, shared by all actions) under the
 * enum value's name. Call on the 'add_interpolation_method_physics' task, which precedes the linear
 * FV kernel task.
 */
InterpolationMethodName addAdvectedInterpolationMethod(const InputParameters & params,
                                                       FEProblemBase & problem,
                                                       Factory & factory);

/// The FVInterpolationMethod name the advection kernels should reference (no object creation).
InterpolationMethodName advectedInterpolationMethodName(const InputParameters & params);
}
