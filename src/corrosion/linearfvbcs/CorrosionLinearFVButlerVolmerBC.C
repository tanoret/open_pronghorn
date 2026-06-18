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

#include "CorrosionLinearFVButlerVolmerBC.h"
#include "CorrosionChemistry.h"

#include <algorithm>
#include <cmath>

registerMooseObject("OpenPronghornApp", CorrosionLinearFVButlerVolmerBC);

InputParameters
CorrosionLinearFVButlerVolmerBC::validParams()
{
  InputParameters params = LinearFVAdvectionDiffusionBC::validParams();
  params.addClassDescription("Butler-Volmer corrosion/plating electrode reaction as a linear "
                             "finite-volume wall boundary condition for a flowing molten salt.");

  params.addRequiredParam<Real>("valence", "Charge number z of the cation.");
  params.addRequiredParam<Real>("exchange_current_density",
                                "Exchange current density i0 at the reference temperature [A/m^2].");
  params.addParam<Real>("E0", 0.0, "Standard electrode potential [V].");
  params.addParam<Real>("alpha_a", 0.5, "Anodic charge-transfer coefficient.");
  params.addParam<Real>("alpha_c", 0.5, "Cathodic charge-transfer coefficient.");
  params.addParam<Real>("c_ref", 1.0, "Reference concentration at which i0 is defined [mol/m^3].");
  params.addParam<Real>(
      "reference_temperature", 923.15, "Reference temperature for the Arrhenius i0 [K].");
  params.addParam<Real>("activation_energy",
                        0.0,
                        "Activation energy of the exchange current [kJ/mol]; 0 disables the "
                        "temperature dependence.");

  params.addRequiredParam<MooseFunctorName>("temperature", "Local temperature functor [K].");
  params.addParam<MooseFunctorName>("metal_potential", "0", "Metal-phase potential functor [V].");
  params.addParam<MooseFunctorName>("salt_potential", "0", "Salt-phase potential functor [V].");
  return params;
}

CorrosionLinearFVButlerVolmerBC::CorrosionLinearFVButlerVolmerBC(const InputParameters & parameters)
  : LinearFVAdvectionDiffusionBC(parameters),
    _valence(getParam<Real>("valence")),
    _i0(getParam<Real>("exchange_current_density")),
    _E0(getParam<Real>("E0")),
    _alpha_a(getParam<Real>("alpha_a")),
    _alpha_c(getParam<Real>("alpha_c")),
    _c_ref(getParam<Real>("c_ref")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _activation_energy(getParam<Real>("activation_energy")),
    _temperature(getFunctor<Real>("temperature")),
    _metal_potential(getFunctor<Real>("metal_potential")),
    _salt_potential(getFunctor<Real>("salt_potential"))
{
  _var.computeCellGradients();
}

void
CorrosionLinearFVButlerVolmerBC::computeFluxCoefficients(Real & A, Real & B) const
{
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto state = determineState();

  const Real temperature = _temperature(face, state);
  const Real eta = _metal_potential(face, state) - _salt_potential(face, state) - _E0;
  const Real f = Corrosion::faraday / (Corrosion::R_gas * temperature);

  // Arrhenius exchange current, anchored so that i0(reference_temperature) reproduces the calibrated
  // value.
  Real i0 = _i0;
  if (_activation_energy != 0.0)
    i0 *= std::exp(_activation_energy * 1000.0 / Corrosion::R_gas *
                   (1.0 / _reference_temperature - 1.0 / temperature));

  // Clamp the Butler-Volmer exponents (matching the reference exp_clip) for numerical robustness.
  const Real anodic_exp = std::min(std::max(_alpha_a * _valence * f * eta, -60.0), 60.0);
  const Real cathodic_exp = std::min(std::max(-_alpha_c * _valence * f * eta, -60.0), 60.0);

  const Real zF = _valence * Corrosion::faraday;
  A = i0 / zF * std::exp(anodic_exp);
  B = i0 / (zF * _c_ref) * std::exp(cathodic_exp);
}

Real
CorrosionLinearFVButlerVolmerBC::cellConcentration() const
{
  const auto elem_info = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM)
                             ? _current_face_info->elemInfo()
                             : _current_face_info->neighborInfo();
  return _var.getElemValue(*elem_info, determineState());
}

Real
CorrosionLinearFVButlerVolmerBC::computeBoundaryValue() const
{
  // No advective penetration at the wall, so the face value is the adjacent cell value.
  return cellConcentration();
}

Real
CorrosionLinearFVButlerVolmerBC::computeBoundaryNormalGradient() const
{
  // Outward molar flux B c - A (the inward Butler-Volmer source is A - B c); used for cell-gradient
  // reconstruction.
  Real A, B;
  computeFluxCoefficients(A, B);
  return B * cellConcentration() - A;
}

Real
CorrosionLinearFVButlerVolmerBC::computeBoundaryValueMatrixContribution() const
{
  return 1.0;
}

Real
CorrosionLinearFVButlerVolmerBC::computeBoundaryValueRHSContribution() const
{
  return 0.0;
}

Real
CorrosionLinearFVButlerVolmerBC::computeBoundaryGradientMatrixContribution() const
{
  // Coefficient of the cell concentration in the outward flux (the cathodic, plating term).
  Real A, B;
  computeFluxCoefficients(A, B);
  return B;
}

Real
CorrosionLinearFVButlerVolmerBC::computeBoundaryGradientRHSContribution() const
{
  // Constant part of the outward flux (the anodic dissolution source).
  Real A, B;
  computeFluxCoefficients(A, B);
  return A;
}
