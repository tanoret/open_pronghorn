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

#include "NTPDarcyFaceFlux.h"
#include "MeshDivision.h"
#include "FaceInfo.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"
#include "ADUtils.h"
#include "NTPSeriesInterpolation.h"

MooseEnum
NTPDarcyFaceFlux::interpolationEnum()
{
  return MooseEnum("series harmonic average", "series");
}

NTPDarcyFaceFlux::Interpolation
NTPDarcyFaceFlux::interpolation(const MooseEnum & value)
{
  if (value == "series")
    return Interpolation::SERIES;
  if (value == "harmonic")
    return Interpolation::HARMONIC;
  return Interpolation::AVERAGE;
}

NTPDarcyFaceFlux::NTPDarcyFaceFlux(const MooseVariableFV<Real> & pressure,
                                   const Moose::Functor<ADRealVectorValue> & rho_mobility,
                                   const Interpolation interp_method,
                                   const bool face_mobility)
  : _pressure(pressure),
    _rho_mobility(rho_mobility),
    _interp_method(interp_method),
    _face_mobility(face_mobility)
{
}

void
NTPDarcyFaceFlux::setMomentumFlux(const Moose::Functor<ADReal> & J,
                                  const std::set<SubdomainID> & blocks)
{
  _momentum_flux = &J;
  _momentum_blocks = blocks;
}

ADReal
NTPDarcyFaceFlux::scalarComponent(const MooseVariableScalar & var, const unsigned int j)
{
  const auto & dofs = var.dofIndices();
  if (j >= dofs.size())
    mooseError("Scalar variable '", var.name(), "' has no component ", j);
  ADReal value = (*var.sys().currentSolution())(dofs[j]);
  // seed only for a variable of the system whose Jacobian is assembled, as MOOSE does for its own
  // scalar variables: the dof index is a column of that system
  if (ADReal::do_derivatives &&
      var.sys().number() == var.sys().feProblem().currentNonlinearSystem().number())
    Moose::derivInsert(value.derivatives(), dofs[j], 1.);
  return value;
}

void
NTPDarcyFaceFlux::setOrifice(const MeshDivision & division,
                             const MooseVariableScalar & log_multiplier,
                             const std::set<SubdomainID> & blocks)
{
  _division = &division;
  _log_orifice = &log_multiplier;
  _orifice_blocks = blocks;
}

ADRealVectorValue
NTPDarcyFaceFlux::sideMobility(const FaceInfo & fi,
                               const Elem & elem,
                               const Moose::StateArg & state) const
{
  // The mobility law of the side's region evaluated at this face: with a face argument the
  // pressure gradient is the compact two-point face gradient, exact for the flow across the face,
  // instead of a cell gradient reconstructed across the porosity jumps of the neighbouring faces.
  auto coeff =
      _face_mobility
          ? _rho_mobility(
                Moose::FaceArg{
                    &fi, Moose::FV::LimiterType::CentralDifference, true, false, &elem, nullptr},
                state)
          : _rho_mobility(Moose::ElemArg{&elem, false}, state);
  if (_log_orifice && _orifice_blocks.count(elem.subdomain_id()))
  {
    const auto j = _division->divisionIndex(elem);
    if (j == MooseMeshDivision::INVALID_DIVISION_INDEX)
      mooseError("An orifice cell (element ",
                 elem.id(),
                 ") lies outside the slices of the orifice's mesh division");
    using std::exp;
    coeff /= exp(scalarComponent(*_log_orifice, j));
  }
  return coeff;
}

ADReal
NTPDarcyFaceFlux::massFlux(const FaceInfo & fi, const Moose::StateArg & state) const
{
  const auto & grad_p = _pressure.adGradSln(fi, state, false);
  const bool internal = _pressure.isInternalFace(fi);

  ADRealVectorValue coeff;
  if (internal)
  {
    const auto elem_coeff = sideMobility(fi, fi.elem(), state);
    const auto neighbor_coeff = sideMobility(fi, *fi.neighborPtr(), state);
    switch (_interp_method)
    {
      case Interpolation::SERIES:
        coeff = NTPSeriesInterpolation::coefficient(elem_coeff, neighbor_coeff, fi);
        break;
      case Interpolation::HARMONIC:
        Moose::FV::interpolate(
            Moose::FV::InterpMethod::HarmonicAverage, coeff, elem_coeff, neighbor_coeff, fi, true);
        break;
      case Interpolation::AVERAGE:
        Moose::FV::interpolate(
            Moose::FV::InterpMethod::Average, coeff, elem_coeff, neighbor_coeff, fi, true);
        break;
    }
  }
  else
  {
    // a boundary face takes the mobility of the cell inside the domain
    const Elem & cell = _pressure.hasBlocks(fi.elem().subdomain_id()) ? fi.elem() : fi.neighbor();
    coeff = sideMobility(fi, cell, state);
  }

  const auto & normal = fi.normal();
  // Darcy part: the diagonal mobility tensor times the face pressure gradient
  ADReal mdot = 0;
  for (const auto i : make_range(Moose::dim))
    mdot -= normal(i) * coeff(i) * grad_p(i);

  // Momentum part, between two channel cells: the axial gradient of J from its jump over the
  // axial distance between the cell centers (y is the channel axis)
  if (_momentum_flux && internal && _momentum_blocks.count(fi.elem().subdomain_id()) &&
      _momentum_blocks.count(fi.neighbor().subdomain_id()))
  {
    const Real dz = fi.dCN()(1);
    if (std::abs(dz) > 0)
      mdot -= normal(1) * coeff(1) *
              ((*_momentum_flux)(Moose::ElemArg{fi.neighborPtr(), false}, state) -
               (*_momentum_flux)(Moose::ElemArg{fi.elemPtr(), false}, state)) /
              dz;
  }
  return mdot;
}
