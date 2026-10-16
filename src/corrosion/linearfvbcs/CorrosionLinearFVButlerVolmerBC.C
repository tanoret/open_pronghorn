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
  params.addRequiredParam<Real>(
      "exchange_current_density",
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

  params.addRangeCheckedParam<Real>(
      "equilibrium_concentration",
      0.0,
      "equilibrium_concentration>=0",
      "Equilibrium salt concentration of the cation at the reference temperature [mol/m^3]. When "
      "positive, the overpotential follows from c_eq(T) instead of the potential functors, so the "
      "net flux vanishes at c = c_eq(T). 0 (default) uses metal_potential - salt_potential - E0.");
  params.addParam<Real>(
      "dissolution_enthalpy",
      0.0,
      "Dissolution enthalpy [kJ/mol] of the equilibrium concentration: c_eq(T) = c_eq,ref "
      "exp[-dH / R (1/T - 1/T_ref)]. Positive values make the equilibrium rise with temperature "
      "(hot walls dissolve, cold walls plate). Used only with equilibrium_concentration.");

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
    _equilibrium_concentration(getParam<Real>("equilibrium_concentration")),
    _dissolution_enthalpy(getParam<Real>("dissolution_enthalpy")),
    _temperature(getFunctor<Real>("temperature")),
    _metal_potential(getFunctor<Real>("metal_potential")),
    _salt_potential(getFunctor<Real>("salt_potential"))
{
  if (_equilibrium_concentration > 0.0 && !(_alpha_a + _alpha_c > 0.0))
    paramError("alpha_a", "alpha_a + alpha_c must be positive with an equilibrium_concentration.");
  if (_equilibrium_concentration > 0.0 &&
      (isParamSetByUser("metal_potential") || isParamSetByUser("salt_potential")))
    paramError("equilibrium_concentration",
               "Give either the equilibrium concentration or the metal/salt potentials, not both.");
  if (_dissolution_enthalpy != 0.0 && !(_equilibrium_concentration > 0.0))
    paramError("dissolution_enthalpy", "Requires a positive equilibrium_concentration.");
  _var.computeCellGradients();
}

void
CorrosionLinearFVButlerVolmerBC::computeFluxCoefficients(Real & A, Real & B) const
{
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto state = determineState();

  const Real temperature = _temperature(face, state);
  const Real f = Corrosion::faraday / (Corrosion::R_gas * temperature);
  Real eta;
  if (_equilibrium_concentration > 0.0)
  {
    // Overpotential at which the net Butler-Volmer flux vanishes at c = c_eq(T).
    const Real c_eq =
        _equilibrium_concentration * std::exp(-_dissolution_enthalpy * 1000.0 / Corrosion::R_gas *
                                              (1.0 / temperature - 1.0 / _reference_temperature));
    eta = std::log(c_eq / _c_ref) / ((_alpha_a + _alpha_c) * _valence * f);
  }
  else
    eta = _metal_potential(face, state) - _salt_potential(face, state) - _E0;

  // Arrhenius exchange current, anchored so that i0(reference_temperature) reproduces the
  // calibrated value.
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
