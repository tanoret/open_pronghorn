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

#include "NTPChannelGradientFunctor.h"
#include "FaceInfo.h"
#include "libmesh/elem.h"

NTPChannelGradientFunctor::NTPChannelGradientFunctor(const MooseFunctorName & name,
                                                     const Moose::Functor<Real> & velocity,
                                                     const MooseMesh & mesh,
                                                     const std::set<SubdomainID> & channel_blocks)
  : Moose::FunctorBase<Real>(name), _velocity(velocity), _mesh(mesh), _channel_blocks(channel_blocks)
{
}

RealVectorValue
NTPChannelGradientFunctor::greenGauss(
    const std::vector<std::pair<Point, Real>> & area_vectors_and_values, const Real volume)
{
  RealVectorValue grad;
  for (const auto & [area_vector, value] : area_vectors_and_values)
    grad += value * area_vector;
  return grad / volume;
}

NTPChannelGradientFunctor::GradientType
NTPChannelGradientFunctor::evaluateGradient(const Moose::ElemArg & elem_arg,
                                            const Moose::StateArg & state) const
{
  const Elem * const elem = elem_arg.elem;
  // outside the channels: the wrapped functor's own gradient
  if (!_channel_blocks.count(elem->subdomain_id()))
    return _velocity.gradient(elem_arg, state);

  const Real value = _velocity(elem_arg, state);
  std::vector<std::pair<Point, Real>> faces;
  for (const auto side : make_range(elem->n_sides()))
  {
    const Elem * const neighbor = elem->neighbor_ptr(side);
    const FaceInfo * fi = _mesh.faceInfo(elem, side);
    bool elem_is_fi_elem = true;
    if (!fi)
    {
      fi = _mesh.faceInfo(neighbor, neighbor->which_neighbor_am_i(elem));
      elem_is_fi_elem = false;
    }
    // outward area vector of this cell (the coordinate factor makes it consistent in RZ)
    const Point area = (elem_is_fi_elem ? 1. : -1.) * fi->normal() * fi->faceArea() * fi->faceCoord();
    Real face_value;
    if (!neighbor)
      face_value = _velocity(Moose::FaceArg{fi, Moose::FV::LimiterType::CentralDifference, true,
                                            false, elem, nullptr},
                             state);
    else if (!_channel_blocks.count(neighbor->subdomain_id()))
      face_value = value; // porous interface: no normal gradient
    else
    {
      const Real gc = elem_is_fi_elem ? fi->gC() : 1 - fi->gC();
      face_value = gc * value + (1 - gc) * _velocity(Moose::ElemArg{neighbor, false}, state);
    }
    // (phi_f - phi_C): exact for a constant field, also on RZ cells whose area vectors do not sum
    // to zero
    faces.emplace_back(area, face_value - value);
  }
  const auto & info = _mesh.elemInfo(elem->id());
  return greenGauss(faces, info.volume() * info.coordFactor());
}
