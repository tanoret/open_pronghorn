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

#include "CorrosionTotalRateAux.h"
#include "ButlerVolmerKinetics.h"
#include "CorrosionChemistry.h"

#include <cmath>
#include <string>

registerMooseObject("OpenPronghornApp", CorrosionTotalRateAux);

InputParameters
CorrosionTotalRateAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Sums species-resolved Butler-Volmer penetration rates or integrates their total recession.");

  MooseEnum mode("penetration_rate recession", "penetration_rate");
  params.addParam<MooseEnum>("mode", mode, "The total-front diagnostic quantity to output.");

  params.addCoupledVar(
      "concentrations", "One salt-side cation concentration variable per species [mol/m^3].");
  params.addParam<MooseFunctorName>(
      "concentration_functor",
      "One external salt concentration functor [mol/m^3] shared by all species. Supply exactly "
      "one of this parameter and the coupled concentration-variable vector.");

  params.addRequiredParam<std::vector<Real>>("valences", "Charge number z of each cation.");
  params.addRequiredParam<std::vector<Real>>(
      "molar_masses", "Molar mass of each metal [g/mol].");
  params.addRequiredParam<std::vector<Real>>(
      "exchange_current_densities", "Exchange current density i0 of each species [A/m^2].");
  params.addRequiredParam<std::vector<Real>>("E0_values", "Standard electrode potentials [V].");
  params.addRequiredParam<std::vector<Real>>("alpha_a_values", "Anodic transfer coefficients.");
  params.addRequiredParam<std::vector<Real>>("alpha_c_values", "Cathodic transfer coefficients.");
  params.addRequiredParam<std::vector<Real>>(
      "c_ref_values", "Reference concentrations of the species [mol/m^3].");
  params.addRequiredParam<Real>("density", "Metal density [g/cm^3].");
  params.addRequiredParam<Real>("temperature", "Butler-Volmer temperature [K].");

  params.addCoupledVar("salt_potential", "Salt-phase electric potential variable [V].");
  params.addParam<MooseFunctorName>("salt_potential_functor",
                                    "Salt-phase electric potential functor [V].");
  params.addCoupledVar("metal_potential", "Metal-phase electric potential variable [V].");
  params.addParam<MooseFunctorName>("metal_potential_functor",
                                    "Metal-phase electric potential functor [V].");
  params.addParam<Real>(
      "salt_potential_value", 0.0, "Fixed salt-phase potential [V] when not coupled or a functor.");
  params.addParam<Real>("metal_potential_value",
                        0.0,
                        "Fixed metal-phase potential [V] when not coupled or a functor.");
  return params;
}

CorrosionTotalRateAux::CorrosionTotalRateAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _mode(getParam<MooseEnum>("mode") == "recession" ? Mode::Recession : Mode::PenetrationRate),
    _valences(getParam<std::vector<Real>>("valences")),
    _molar_masses(getParam<std::vector<Real>>("molar_masses")),
    _exchange_current_densities(getParam<std::vector<Real>>("exchange_current_densities")),
    _E0_values(getParam<std::vector<Real>>("E0_values")),
    _alpha_a_values(getParam<std::vector<Real>>("alpha_a_values")),
    _alpha_c_values(getParam<std::vector<Real>>("alpha_c_values")),
    _c_ref_values(getParam<std::vector<Real>>("c_ref_values")),
    _density(getParam<Real>("density")),
    _temperature(getParam<Real>("temperature")),
    _concentration_functor(isParamValid("concentration_functor")
                               ? &getFunctor<Real>("concentration_functor")
                               : nullptr),
    _has_salt_potential(isCoupled("salt_potential")),
    _salt_potential_var(_has_salt_potential ? coupledValue("salt_potential") : _zero),
    _salt_potential_functor(isParamValid("salt_potential_functor")
                                ? &getFunctor<Real>("salt_potential_functor")
                                : nullptr),
    _salt_potential(getParam<Real>("salt_potential_value")),
    _has_metal_potential(isCoupled("metal_potential")),
    _metal_potential_var(_has_metal_potential ? coupledValue("metal_potential") : _zero),
    _metal_potential_functor(isParamValid("metal_potential_functor")
                                 ? &getFunctor<Real>("metal_potential_functor")
                                 : nullptr),
    _metal_potential(getParam<Real>("metal_potential_value")),
    _u_old(uOld())
{
  const auto n_species = _valences.size();
  if (n_species == 0)
    paramError("valences", "At least one species is required for a total-front diagnostic.");

  const auto require_size = [this, n_species](const std::string & name, std::size_t size)
  {
    if (size != n_species)
      paramError(name,
                 "Must contain one entry per species (",
                 n_species,
                 " entries, matching 'valences').");
  };
  require_size("molar_masses", _molar_masses.size());
  require_size("exchange_current_densities", _exchange_current_densities.size());
  require_size("E0_values", _E0_values.size());
  require_size("alpha_a_values", _alpha_a_values.size());
  require_size("alpha_c_values", _alpha_c_values.size());
  require_size("c_ref_values", _c_ref_values.size());

  const auto n_concentrations = coupledComponents("concentrations");
  const bool has_concentrations = n_concentrations != 0;
  if (has_concentrations == static_cast<bool>(_concentration_functor))
    paramError(has_concentrations ? "concentration_functor" : "concentrations",
               "Supply exactly one salt concentration source: one coupled 'concentrations' "
               "variable per species or the shared 'concentration_functor'.");
  if (has_concentrations && n_concentrations != n_species)
    paramError("concentrations",
               "Must contain one coupled concentration per species (",
               n_species,
               " entries, matching 'valences').");
  _concentrations.reserve(n_concentrations);
  for (unsigned int i = 0; i < n_concentrations; ++i)
    _concentrations.push_back(&coupledValue("concentrations", i));

  const auto validate_potential = [this](const std::string & coupled_name,
                                         bool has_coupled,
                                         const std::string & functor_name,
                                         const Moose::Functor<Real> * functor,
                                         const std::string & value_name)
  {
    const unsigned int source_count = static_cast<unsigned int>(has_coupled) +
                                      static_cast<unsigned int>(functor != nullptr) +
                                      static_cast<unsigned int>(isParamSetByUser(value_name));
    if (source_count > 1)
      paramError(functor ? functor_name : coupled_name,
                 "Supply at most one ",
                 coupled_name,
                 " source: the coupled variable, '",
                 functor_name,
                 "', or the fixed '",
                 value_name,
                 "'.");
  };
  validate_potential("salt_potential",
                     _has_salt_potential,
                     "salt_potential_functor",
                     _salt_potential_functor,
                     "salt_potential_value");
  validate_potential("metal_potential",
                     _has_metal_potential,
                     "metal_potential_functor",
                     _metal_potential_functor,
                     "metal_potential_value");

  if (!std::isfinite(_density) || _density <= 0.0)
    paramError("density", "Must be finite and positive.");
  if (!std::isfinite(_temperature) || _temperature <= 0.0)
    paramError("temperature", "Must be finite and positive.");
  for (std::size_t i = 0; i < n_species; ++i)
  {
    if (!std::isfinite(_valences[i]) || _valences[i] <= 0.0)
      paramError("valences", "Every entry must be finite and positive.");
    if (!std::isfinite(_molar_masses[i]) || _molar_masses[i] <= 0.0)
      paramError("molar_masses", "Every entry must be finite and positive.");
    if (!std::isfinite(_exchange_current_densities[i]) ||
        _exchange_current_densities[i] < 0.0)
      paramError("exchange_current_densities", "Every entry must be finite and nonnegative.");
    if (!std::isfinite(_E0_values[i]))
      paramError("E0_values", "Every entry must be finite.");
    if (!std::isfinite(_alpha_a_values[i]) || !std::isfinite(_alpha_c_values[i]))
      paramError("alpha_a_values", "All transfer coefficients must be finite.");
    if (!std::isfinite(_c_ref_values[i]) || _c_ref_values[i] <= 0.0)
      paramError("c_ref_values", "Every entry must be finite and positive.");
  }

  if (_mode == Mode::Recession && !isNodal())
    paramError("mode", "Recession accumulation requires a nodal auxiliary variable.");
}

Real
CorrosionTotalRateAux::functorValue(const Moose::Functor<Real> & functor) const
{
  const auto state = determineState();
  if (isNodal())
  {
    const Moose::NodeArg node_arg = {_current_node,
                                     &Moose::NodeArg::undefined_subdomain_connection};
    return functor(node_arg, state);
  }
  if (_bnd)
  {
    const Moose::ElemSideQpArg side_qp_arg = {
        _current_elem, _current_side, _qp, _qrule, _q_point[_qp]};
    return functor(side_qp_arg, state);
  }
  const Moose::ElemQpArg qp_arg = {_current_elem, _qp, _qrule, _q_point[_qp]};
  return functor(qp_arg, state);
}

Real
CorrosionTotalRateAux::computeValue()
{
  const Real phi_salt = _has_salt_potential
                            ? _salt_potential_var[_qp]
                            : (_salt_potential_functor ? functorValue(*_salt_potential_functor)
                                                       : _salt_potential);
  const Real phi_metal = _has_metal_potential
                             ? _metal_potential_var[_qp]
                             : (_metal_potential_functor ? functorValue(*_metal_potential_functor)
                                                        : _metal_potential);
  const Real shared_concentration =
      _concentration_functor ? functorValue(*_concentration_functor) : 0.0;

  Real total_rate_um_y = 0.0;
  for (std::size_t i = 0; i < _valences.size(); ++i)
  {
    const Real concentration =
        _concentration_functor ? shared_concentration : (*_concentrations[i])[_qp];
    const Real current = Corrosion::bvCurrent(phi_metal,
                                               phi_salt,
                                               concentration,
                                               _E0_values[i],
                                               _valences[i],
                                               _exchange_current_densities[i],
                                               _alpha_a_values[i],
                                               _alpha_c_values[i],
                                               _temperature,
                                               _c_ref_values[i]);
    total_rate_um_y += Corrosion::corrosionCurrentToUmY(
        current * 1.0e-4, _valences[i], _molar_masses[i], _density);
  }

  if (_mode == Mode::PenetrationRate)
    return total_rate_um_y;

  return _u_old[_qp] + total_rate_um_y * _dt / Corrosion::seconds_per_year;
}
