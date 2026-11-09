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

#include "NTPSliceUniformityConstraint.h"
#include "MeshDivision.h"
#include "MooseVariableScalar.h"
#include "MooseVariableFV.h"
#include "Assembly.h"
#include "FEProblemBase.h"
#include "NTPDarcyFaceFlux.h"

#include <numeric>
#include "Function.h"

registerMooseObject("OpenPronghornApp", NTPSliceUniformityConstraint);

InputParameters
NTPSliceUniformityConstraint::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription("Equations for a scalar with one component per slice of a mesh "
                             "division that make a field's slice averages equal.");
  params.addRequiredParam<MeshDivisionName>("mesh_division", "The slices");
  params.addRequiredCoupledVar(
      "lambda", "Scalar unknown with one component per slice (order = number of slices)");
  params.addRequiredParam<MooseFunctorName>("field", "Field whose slice averages are equalized");
  params.addParam<MooseEnum>(
      "reference",
      MooseEnum("slice mean", "slice"),
      "Equalize the slice averages with the one of 'reference_slice', whose component is fixed "
      "to 0, or with the volume average of all slices, the mean of the components being fixed to "
      "'mean_log_multiplier'. With 'mean' no slice needs a drag below the level exp(lambda_mean) "
      "while it is far enough above 1, which keeps every slice controllable");
  params.addParam<unsigned int>("reference_slice", 0, "Slice whose scalar component is fixed to 0");
  params.addParam<Real>("mean_log_multiplier",
                        0,
                        "Mean of the components with 'reference = mean' (the logarithm of the "
                        "geometric mean of the orifice multiplier)");
  params.addRangeCheckedParam<Real>(
      "field_scale",
      1,
      "field_scale > 0",
      "Scale of the field differences, e.g. the expected slice temperature rise [field units]; "
      "the equations are written for the dimensionless mismatch");
  params.addParam<PostprocessorName>(
      "field_scale_pp",
      "Postprocessor giving the scale instead of 'field_scale' (e.g. the current mean temperature "
      "rise, so the equations stay well scaled along a power ramp); floored by 'field_scale'");
  params.addParam<FunctionName>(
      "activation",
      1,
      "Weight w(t) in [0, 1] of the uniformity equations: each row is w (scaled mismatch) + (1 - "
      "w) "
      "lambda_j, so w = 0 pins lambda to 0. Ramp it up once the field differs between slices (the "
      "orifice cannot be seen before the gas is heated)");
  params.addRangeCheckedParam<Real>(
      "regularization",
      1e-3,
      "regularization >= 0",
      "Dimensionless coefficient delta of a term delta * lambda_j added to every equation. It "
      "keeps "
      "the scalar defined where the field does not depend on it (e.g. before any heating) and is "
      "negligible once the scaled field differences respond to lambda.");
  return params;
}

NTPSliceUniformityConstraint::NTPSliceUniformityConstraint(const InputParameters & parameters)
  : FVElementalKernel(parameters),
    _division(_fe_problem.getMeshDivision(getParam<MeshDivisionName>("mesh_division"))),
    _lambda_var(*getScalarVar("lambda", 0)),
    _mean_reference(getParam<MooseEnum>("reference") == "mean"),
    _ref(_mean_reference ? _division.getNumDivisions() - 1
                         : getParam<unsigned int>("reference_slice")),
    _lambda_mean(getParam<Real>("mean_log_multiplier")),
    _q(getFunctor<ADReal>("field")),
    _delta(getParam<Real>("regularization")),
    _activation(getFunction("activation")),
    _scale(getParam<Real>("field_scale")),
    _scale_pp(isParamValid("field_scale_pp") ? &getPostprocessorValue("field_scale_pp") : nullptr),
    _slice_volume(_division.getNumDivisions(), 0)
{
  if (_lambda_var.order() != _division.getNumDivisions())
    paramError("lambda",
               "Needs one component per slice: order ",
               _division.getNumDivisions(),
               ", not ",
               _lambda_var.order());
  if (_ref >= _division.getNumDivisions())
    paramError("reference_slice", "Exceeds the number of slices");
}

void
NTPSliceUniformityConstraint::initialSetup()
{
  FVElementalKernel::initialSetup();
  std::fill(_slice_volume.begin(), _slice_volume.end(), 0);
  for (const auto * const elem : *_mesh.getActiveLocalElementRange())
    if (hasBlocks(elem->subdomain_id()))
    {
      const auto j = _division.divisionIndex(*elem);
      if (j < _slice_volume.size())
        _slice_volume[j] += _assembly.elementVolume(elem);
    }
  _communicator.sum(_slice_volume);
  for (const auto j : index_range(_slice_volume))
    if (_slice_volume[j] <= 0)
      mooseError(name(), ": slice ", j, " contains no cell of the kernel's blocks");
  _total_volume = std::accumulate(_slice_volume.begin(), _slice_volume.end(), 0.0);
}

std::vector<std::pair<unsigned int, ADReal>>
NTPSliceUniformityConstraint::cellContributions()
{
  std::vector<std::pair<unsigned int, ADReal>> contributions;
  const auto j = _division.divisionIndex(*_current_elem);
  if (j >= _slice_volume.size())
    return contributions;

  const auto dV = _assembly.elemVolume();
  const Real scale = _scale_pp ? std::max(*_scale_pp, _scale) : _scale;
  const Real w = std::min(std::max(_activation.value(_t, Point()), 0.0), 1.0);
  const auto q = w * _q(makeElemArg(_current_elem), determineState()) / scale;
  // lambda from the solution vector with seeded derivatives (cached scalar AD values are not
  // always seeded), see NTPDarcyFaceFlux::scalarComponent
  const auto lambda = [this](unsigned int k)
  { return NTPDarcyFaceFlux::scalarComponent(_lambda_var, k); };

  if (_mean_reference)
  {
    // rows k < n - 1: w (<q>_k - q_mean) / s + (1 - w + delta) (lambda_k - lambda_mean);
    // row n - 1: mean_k lambda_k - lambda_mean
    const auto n = _slice_volume.size();
    const Real v_slice = dV / _slice_volume[j];
    const Real v_total = dV / _total_volume;
    for (const auto k : make_range(n - 1))
    {
      ADReal r = -q * v_total;
      if (k == j)
        r += (q + ((1 - w) + _delta) * (lambda(j) - _lambda_mean)) * v_slice;
      contributions.emplace_back(k, r);
    }
    contributions.emplace_back(n - 1, lambda(j) / n * v_slice - _lambda_mean * v_total);
    return contributions;
  }

  if (j != _ref)
    contributions.emplace_back(j, (q + ((1 - w) + _delta) * lambda(j)) * dV / _slice_volume[j]);
  else
  {
    const auto v = dV / _slice_volume[_ref];
    for (const auto k : index_range(_slice_volume))
      if (k != _ref)
        contributions.emplace_back(k, -q * v);
    // spread over the reference slice, this sums to R_r = lambda_r
    contributions.emplace_back(_ref, lambda(_ref) * v);
  }
  return contributions;
}

void
NTPSliceUniformityConstraint::computeResidual()
{
  const auto & dofs = _lambda_var.dofIndices();
  for (const auto & [k, r] : cellContributions())
    addResiduals(_assembly,
                 std::array<Real, 1>{{MetaPhysicL::raw_value(r)}},
                 std::array<dof_id_type, 1>{{dofs[k]}},
                 _lambda_var.scalingFactor());
}

void
NTPSliceUniformityConstraint::computeResidualAndJacobian()
{
  const auto & dofs = _lambda_var.dofIndices();
  for (const auto & [k, r] : cellContributions())
    addResidualsAndJacobian(_assembly,
                            std::array<ADReal, 1>{{r}},
                            std::array<dof_id_type, 1>{{dofs[k]}},
                            _lambda_var.scalingFactor());
}

void
NTPSliceUniformityConstraint::computeJacobian()
{
  // the only Jacobian path: the base class forwards computeOffDiagJacobian() here, so diagonal and
  // full coupling (and the scaling Jacobian) all assemble the scalar rows
  const auto & dofs = _lambda_var.dofIndices();
  for (const auto & [k, r] : cellContributions())
    addJacobian(_assembly,
                std::array<ADReal, 1>{{r}},
                std::array<dof_id_type, 1>{{dofs[k]}},
                _lambda_var.scalingFactor());
}

ADReal
NTPSliceUniformityConstraint::computeQpResidual()
{
  mooseError("NTPSliceUniformityConstraint assembles its residuals directly");
}
