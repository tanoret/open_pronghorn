//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MSRDegassingBC.h"

registerMooseObject("OpenPronghornApp", MSRDegassingBC);

InputParameters
MSRDegassingBC::validParams()
{
  InputParameters params = LinearFVAdvectionDiffusionBC::validParams();
  params.addClassDescription("Surface degassing of a dissolved species to a cover gas: the species "
                             "leaves the melt at the convective mass-transfer rate "
                             "J = k_surf * (C - kH * p_cover).");
  params.addRequiredParam<MooseFunctorName>(
      "mass_transfer_coefficient", "Surface mass-transfer coefficient k_surf [m/s].");
  params.addParam<Real>("henry_coefficient", 0.0, "Henry coefficient kH [mol/(m^3 Pa)].");
  params.addParam<MooseFunctorName>(
      "cover_gas_pressure",
      "0",
      "Partial pressure of the species in the cover gas p_cover [Pa] (functor; default 0 for a "
      "continuously swept clean cover gas).");
  return params;
}

MSRDegassingBC::MSRDegassingBC(const InputParameters & params)
  : LinearFVAdvectionDiffusionBC(params),
    _k_surf(getFunctor<Real>("mass_transfer_coefficient")),
    _kH(getParam<Real>("henry_coefficient")),
    _p_cover(getFunctor<Real>("cover_gas_pressure"))
{
  _var.computeCellGradients();
}

Real
MSRDegassingBC::equilibriumConcentration() const
{
  const auto face = singleSidedFaceArg(_current_face_info);
  return _kH * _p_cover(face, determineState());
}

Real
MSRDegassingBC::computeBoundaryValue() const
{
  // The free surface is a no-penetration boundary for the melt, so the advected face value is
  // approximated by the adjacent cell value.
  const auto elem_info = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM)
                             ? _current_face_info->elemInfo()
                             : _current_face_info->neighborInfo();
  return _var.getElemValue(*elem_info, determineState());
}

Real
MSRDegassingBC::computeBoundaryNormalGradient() const
{
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto elem_info = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM)
                             ? _current_face_info->elemInfo()
                             : _current_face_info->neighborInfo();
  // Outward degassing flux k_surf * (C - C_eq); used for cell-gradient reconstruction.
  return _k_surf(face, determineState()) *
         (_var.getElemValue(*elem_info, determineState()) - equilibriumConcentration());
}

Real
MSRDegassingBC::computeBoundaryValueMatrixContribution() const
{
  // Face value approximated by the cell value.
  return 1.0;
}

Real
MSRDegassingBC::computeBoundaryValueRHSContribution() const
{
  return 0.0;
}

Real
MSRDegassingBC::computeBoundaryGradientMatrixContribution() const
{
  // Coefficient of the boundary value in the degassing flux k_surf * (C - C_eq).
  const auto face = singleSidedFaceArg(_current_face_info);
  return _k_surf(face, determineState());
}

Real
MSRDegassingBC::computeBoundaryGradientRHSContribution() const
{
  // Constant part of the degassing flux: k_surf * C_eq.
  const auto face = singleSidedFaceArg(_current_face_info);
  return _k_surf(face, determineState()) * equilibriumConcentration();
}
