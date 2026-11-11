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

#include "NTPRhieChowFluxIntegral.h"
#include "RhieChowMassFlux.h"

registerMooseObject("OpenPronghornApp", NTPRhieChowFluxIntegral);

InputParameters
NTPRhieChowFluxIntegral::validParams()
{
  InputParameters params = SideIntegralPostprocessor::validParams();
  params.addClassDescription("Boundary integral of the Rhie-Chow face mass fluxes of the "
                             "segregated solvers, positive out of the domain, optionally times a "
                             "carried quantity such as the enthalpy.");
  params.addRequiredParam<UserObjectName>("rhie_chow_user_object", "The Rhie-Chow object");
  params.addParam<MooseFunctorName>("carried_quantity", 1, "Quantity carried by the mass flux");
  params.addParam<bool>("normal_kinetic_energy",
                        false,
                        "Carry the kinetic energy of the normal face velocity (mdot_f / rho)^2 / 2 "
                        "instead of 'carried_quantity'");
  params.addParam<MooseFunctorName>("rho", "rho", "Density for 'normal_kinetic_energy'");
  return params;
}

NTPRhieChowFluxIntegral::NTPRhieChowFluxIntegral(const InputParameters & parameters)
  : SideIntegralPostprocessor(parameters),
    _rc(nullptr),
    _q(getFunctor<Real>("carried_quantity")),
    _normal_ke(getParam<bool>("normal_kinetic_energy")),
    _rho(_normal_ke ? &getFunctor<Real>("rho") : nullptr)
{
  _qp_integration = false;
}

void
NTPRhieChowFluxIntegral::initialSetup()
{
  SideIntegralPostprocessor::initialSetup();
  _rc = &getUserObject<RhieChowMassFlux>("rhie_chow_user_object");
}

Real
NTPRhieChowFluxIntegral::computeQpIntegral()
{
  mooseError("Integrates over faces only");
}

Real
NTPRhieChowFluxIntegral::computeFaceInfoIntegral(const FaceInfo * const fi)
{
  // the stored flux is along the face normal, from the face's element to its neighbor
  const bool elem_inside = _rc->hasBlocks(fi->elem().subdomain_id());
  const Real outward = elem_inside ? _rc->getMassFlux(*fi) : -_rc->getMassFlux(*fi);
  const Elem * const inside = elem_inside ? fi->elemPtr() : fi->neighborPtr();
  // the carried quantity at the face, from the inside: the boundary value where the variable has
  // one (an inflow), its boundary extrapolation otherwise, as the advection kernels use it
  const Moose::FaceArg face{
      fi, Moose::FV::LimiterType::CentralDifference, true, false, inside, nullptr};
  if (_normal_ke)
  {
    const Real un = outward / (*_rho)(face, determineState());
    return outward * 0.5 * un * un;
  }
  return outward * _q(face, determineState());
}
