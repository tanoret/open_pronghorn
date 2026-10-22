//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVTurbulentWallCHTBC.h"
#include "NavierStokesMethods.h"

#include <cmath>

registerMooseObject("OpenPronghornApp", LinearFVTurbulentWallCHTBC);

InputParameters
LinearFVTurbulentWallCHTBC::validParams()
{
  InputParameters params = LinearFVAdvectionDiffusionBC::validParams();
  params.addClassDescription(
      "Conjugate heat transfer condition for a turbulent wall with the Jayatilleke thermal law "
      "of the wall; applied on both the fluid and the solid side of the interface.");
  MooseEnum side("fluid solid");
  params.addRequiredParam<MooseEnum>("side", side, "The side of the interface this BC acts on.");
  params.addRequiredParam<MooseFunctorName>("u", "The fluid velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The fluid velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The fluid velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Fluid density.");
  params.addRequiredParam<MooseFunctorName>(NS::mu, "Fluid dynamic (molecular) viscosity.");
  params.addRequiredParam<MooseFunctorName>(NS::cp, "Fluid specific heat capacity.");
  params.addRequiredParam<MooseFunctorName>(
      "k_fluid",
      "Fluid molecular thermal conductivity (sets the Prandtl number and the sublayer "
      "coefficient).");
  params.addRequiredParam<MooseFunctorName>(
      "thermal_conductivity",
      "The conductivity used by the conduction kernel of this side: k + mu_t cp / Pr_t on the "
      "fluid side, the solid conductivity on the solid side.");
  params.addParam<MooseFunctorName>(NS::TKE, "Turbulent kinetic energy (wall_treatment = neq).");
  params.addParam<Real>("Pr_t", 0.85, "Turbulent Prandtl number.");
  params.addParam<Real>("C_mu", 0.09, "k-epsilon closure coefficient.");
  MooseEnum wall_treatment("eq_newton eq_incremental eq_linearized neq", "eq_newton");
  params.addParam<MooseEnum>(
      "wall_treatment", wall_treatment, "The method used to obtain the friction velocity.");
  params.addRequiredParam<MooseFunctorName>(
      "other_temperature",
      "The temperature handed over by the other side: interface_temperature_solid_<boundary> "
      "on the fluid side, interface_temperature_fluid_<boundary> on the solid side.");
  return params;
}

LinearFVTurbulentWallCHTBC::LinearFVTurbulentWallCHTBC(const InputParameters & parameters)
  : LinearFVAdvectionDiffusionBC(parameters),
    LinearFVCHTBCInterface(),
    _is_fluid_side(getParam<MooseEnum>("side") == "fluid"),
    _u(getFunctor<Real>("u")),
    _v(isParamValid("v") ? &getFunctor<Real>("v") : nullptr),
    _w(isParamValid("w") ? &getFunctor<Real>("w") : nullptr),
    _rho(getFunctor<Real>(NS::density)),
    _mu(getFunctor<Real>(NS::mu)),
    _cp(getFunctor<Real>(NS::cp)),
    _k_fluid(getFunctor<Real>("k_fluid")),
    _k_kernel(getFunctor<Real>("thermal_conductivity")),
    _tke(isParamValid(NS::TKE) ? &getFunctor<Real>(NS::TKE) : nullptr),
    _Pr_t(getParam<Real>("Pr_t")),
    _C_mu(getParam<Real>("C_mu")),
    _wall_treatment(getParam<MooseEnum>("wall_treatment").getEnum<NS::WallTreatmentEnum>()),
    _other_temperature(getFunctor<Real>("other_temperature"))
{
  if (_wall_treatment == NS::WallTreatmentEnum::NEQ && !_tke)
    paramError(NS::TKE, "The non-equilibrium wall treatment needs the turbulent kinetic energy.");
}

Real
LinearFVTurbulentWallCHTBC::jayatillekeP(const Real Pr, const Real Pr_t)
{
  const Real r = Pr / Pr_t;
  return 9.24 * (std::pow(r, 0.75) - 1.0) * (1.0 + 0.28 * std::exp(-0.007 * r));
}

Real
LinearFVTurbulentWallCHTBC::wallFunctionHTC() const
{
  // The fluid cell is on this side for the fluid condition and on the other side for the solid
  const bool this_is_elem = _current_face_type == FaceInfo::VarFaceNeighbors::ELEM;
  const bool fluid_is_elem = _is_fluid_side ? this_is_elem : !this_is_elem;
  const auto & fi = *_current_face_info;
  const Elem * const fluid_elem = fluid_is_elem ? fi.elemPtr() : fi.neighborPtr();
  const Point & fluid_centroid = fluid_is_elem ? fi.elemCentroid() : fi.neighborCentroid();
  const Real y = std::abs((fluid_centroid - fi.faceCentroid()) * fi.normal());

  const auto re = makeElemArg(fluid_elem);
  // Lag the flow quantities by one iteration, as the momentum wall function does
  const auto old_state = Moose::StateArg(1, Moose::SolutionIterationType::Nonlinear);
  const Real rho = _rho(re, old_state);
  const Real mu = _mu(re, old_state);
  const Real cp = _cp(re, old_state);
  const Real k = _k_fluid(re, old_state);

  RealVectorValue velocity(_u(re, old_state));
  if (_v)
    velocity(1) = (*_v)(re, old_state);
  if (_w)
    velocity(2) = (*_w)(re, old_state);
  const auto & n = fi.normal();
  const Real u_par = NS::computeSpeed<Real>(velocity - (velocity * n) * n);

  Real u_tau = 0.0;
  if (_wall_treatment == NS::WallTreatmentEnum::NEQ)
    u_tau = std::pow(_C_mu, 0.25) * std::sqrt(std::max((*_tke)(re, old_state), 0.0));
  else if (_wall_treatment == NS::WallTreatmentEnum::EQ_INCREMENTAL)
  {
    const Real y_plus = NS::findyPlus<Real>(mu, rho, std::max(u_par, 1e-10), y);
    u_tau = u_par * NS::von_karman_constant /
            std::log(std::max(NS::E_turb_constant * y_plus, 1.0 + 1e-4));
  }
  else
    // eq_newton and eq_linearized share the Newton solve of the log law here
    u_tau = u_par > 0 ? NS::findUStar<Real>(mu, rho, u_par, y) : 0.0;

  const Real y_plus = rho * u_tau * y / mu;
  const Real h_lam = k / y;
  if (y_plus <= 5.0)
    return h_lam;

  const Real Pr = cp * mu / k;
  const Real T_plus =
      _Pr_t * (std::log(NS::E_turb_constant * std::max(y_plus, 30.0)) / NS::von_karman_constant +
               jayatillekeP(Pr, _Pr_t));
  const Real h_log = rho * cp * u_tau / T_plus;
  if (y_plus >= 30.0)
    return h_log;

  // buffer layer: blend between the sublayer conduction and the y+ = 30 log-law coefficient
  const Real blend = (y_plus - 5.0) / 25.0;
  return blend * h_log + (1.0 - blend) * h_lam;
}

Real
LinearFVTurbulentWallCHTBC::couplingCoefficient() const
{
  const Real h = wallFunctionHTC();
  if (_is_fluid_side)
    return h;
  // solid: the wall function in series with the half solid cell
  const Real k_s = _k_kernel(singleSidedFaceArg(_current_face_info), determineState());
  return 1.0 / (1.0 / h + computeCellToFaceDistance() / k_s);
}

Real
LinearFVTurbulentWallCHTBC::cellTemperature() const
{
  const auto * const elem_info = _current_face_type == FaceInfo::VarFaceNeighbors::ELEM
                                     ? _current_face_info->elemInfo()
                                     : _current_face_info->neighborInfo();
  return _var.getElemValue(*elem_info, determineState());
}

Real
LinearFVTurbulentWallCHTBC::otherTemperature() const
{
  return _other_temperature(functorFaceArg(_other_temperature, _current_face_info),
                            determineState());
}

Real
LinearFVTurbulentWallCHTBC::computeBoundaryValue() const
{
  // The value handed to the other side by the CHT handler
  if (_is_fluid_side)
    // the near-wall fluid cell temperature, the reference of the solid-side coupling
    return cellTemperature();

  // the solid surface temperature: T_w = T_S + U d_S / k_S (T_P - T_S)
  const Real k_s = _k_kernel(singleSidedFaceArg(_current_face_info), determineState());
  const Real T_S = cellTemperature();
  return T_S +
         couplingCoefficient() * computeCellToFaceDistance() / k_s * (otherTemperature() - T_S);
}

Real
LinearFVTurbulentWallCHTBC::computeBoundaryNormalGradient() const
{
  // outward-normal gradient such that -k dT/dn equals the outward flux C (T_cell - T_other)
  return couplingCoefficient() * (otherTemperature() - cellTemperature()) /
         _k_kernel(singleSidedFaceArg(_current_face_info), determineState());
}

Real
LinearFVTurbulentWallCHTBC::computeBoundaryValueMatrixContribution() const
{
  // only used by advection kernels, which do not act on a wall
  return 0.0;
}

Real
LinearFVTurbulentWallCHTBC::computeBoundaryValueRHSContribution() const
{
  return computeBoundaryValue();
}

Real
LinearFVTurbulentWallCHTBC::computeBoundaryGradientMatrixContribution() const
{
  // the diffusion kernel multiplies by the face area only (includesMaterialPropertyMultiplier)
  return couplingCoefficient();
}

Real
LinearFVTurbulentWallCHTBC::computeBoundaryGradientRHSContribution() const
{
  return couplingCoefficient() * otherTemperature();
}
