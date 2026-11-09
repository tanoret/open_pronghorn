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

#include "NTPHydrogenFluidProperties.h"

#include <cmath>

registerMooseObject("OpenPronghornApp", NTPHydrogenFluidProperties);

InputParameters
NTPHydrogenFluidProperties::validParams()
{
  InputParameters params = SinglePhaseFluidProperties::validParams();
  params.addClassDescription("Hot hydrogen gas for nuclear thermal propulsion: ideal gas with NASA "
                             "polynomial heat capacity, Sutherland viscosity and Eucken "
                             "conductivity.");
  params.addParam<Real>("molar_mass", 2.01588e-3, "Molar mass [kg/mol]");
  params.addParam<Real>(
      "universal_gas_constant",
      8.314462618,
      "Universal gas constant [J/(mol K)]; the default is the exact SI value, "
      "also used by the reference channel model (MOOSE's FluidProperties uses 8.3144598)");
  // GRI-Mech 3.0 coefficients for H2 (a1..a7 of the NASA 7-coefficient form)
  params.addParam<std::vector<Real>>("coefficients_low",
                                     {2.34433112,
                                      7.98052075e-03,
                                      -1.94781510e-05,
                                      2.01572094e-08,
                                      -7.37611761e-12,
                                      -9.17935173e02,
                                      6.83010238e-01},
                                     "NASA coefficients a1..a7 below the switch temperature");
  params.addParam<std::vector<Real>>("coefficients_high",
                                     {3.33727920,
                                      -4.94024731e-05,
                                      4.99456778e-07,
                                      -1.79566394e-10,
                                      2.00255376e-14,
                                      -9.50158922e02,
                                      -3.20502331},
                                     "NASA coefficients a1..a7 above the switch temperature");
  params.addParam<Real>("switch_temperature", 1000, "Temperature between the two branches [K]");
  params.addParam<Real>("mu_0", 8.76e-6, "Sutherland reference viscosity [Pa s]");
  params.addParam<Real>("T_mu_0", 300, "Sutherland reference temperature [K]");
  params.addParam<Real>("sutherland_constant", 72, "Sutherland constant [K]");
  params.addParam<Real>(
      "eucken_factor", 1.25, "Eucken factor f in k = mu (cp + f R) [-]; 5/4 for a diatomic gas");
  params.addRangeCheckedParam<Real>(
      "T_min", 20, "T_min > 0", "Lowest temperature at which properties are evaluated [K]");
  params.addParam<Real>(
      "T_max",
      3500,
      "Highest temperature at which properties are evaluated [K]: the upper limit of the "
      "GRI-Mech 3.0 fit; above it cp, mu and k are held and h continues linearly");
  return params;
}

NTPHydrogenFluidProperties::NTPHydrogenFluidProperties(const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters),
    _molar_mass(getParam<Real>("molar_mass")),
    _R_specific(getParam<Real>("universal_gas_constant") / _molar_mass),
    _a_low(getParam<std::vector<Real>>("coefficients_low")),
    _a_high(getParam<std::vector<Real>>("coefficients_high")),
    _T_switch(getParam<Real>("switch_temperature")),
    _mu_0(getParam<Real>("mu_0")),
    _T_mu_0(getParam<Real>("T_mu_0")),
    _S_mu(getParam<Real>("sutherland_constant")),
    _eucken(getParam<Real>("eucken_factor")),
    _T_min(getParam<Real>("T_min")),
    _T_max(getParam<Real>("T_max"))
{
  if (_a_low.size() != 7 || _a_high.size() != 7)
    paramError("coefficients_low", "Both coefficient sets need the seven NASA coefficients");
  if (_T_max <= _T_min)
    paramError("T_max", "Must exceed T_min");
}

Real
NTPHydrogenFluidProperties::clamp(Real T) const
{
  return std::min(std::max(T, _T_min), _T_max);
}

Real
NTPHydrogenFluidProperties::cp_from_T(Real T) const
{
  T = clamp(T);
  const auto & a = coeffs(T);
  return _R_specific * (a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4]))));
}

Real
NTPHydrogenFluidProperties::dcp_dT_from_T(Real T) const
{
  if (T < _T_min || T > _T_max)
    return 0;
  const auto & a = coeffs(T);
  return _R_specific * (a[1] + T * (2 * a[2] + T * (3 * a[3] + T * 4 * a[4])));
}

Real
NTPHydrogenFluidProperties::h_from_T(Real T) const
{
  // Linear continuation with the end-point cp outside the range keeps h monotone and smooth.
  const Real Tc = clamp(T);
  const auto & a = coeffs(Tc);
  const Real h =
      _R_specific *
      (Tc * (a[0] + Tc * (a[1] / 2 + Tc * (a[2] / 3 + Tc * (a[3] / 4 + Tc * a[4] / 5)))) + a[5]);
  return h + cp_from_T(Tc) * (T - Tc);
}

Real
NTPHydrogenFluidProperties::T_from_h(Real h) const
{
  // Newton iteration started from a constant-cp estimate; h(T) is monotone with cp > 0.
  Real T = std::max(_T_min, (h - h_from_T(300)) / (3.5 * _R_specific) + 300);
  for (unsigned int it = 0; it < 50; ++it)
  {
    const Real dT = (h_from_T(T) - h) / cp_from_T(T);
    T -= dT;
    if (std::abs(dT) < 1e-10 * std::max(T, 1.0))
      return T;
  }
  mooseException(name(), ": T_from_h did not converge for h = ", h);
}

void
NTPHydrogenFluidProperties::muT(Real T, Real & mu, Real & dmu_dT) const
{
  const Real Tc = clamp(T);
  mu = _mu_0 * std::pow(Tc / _T_mu_0, 1.5) * (_T_mu_0 + _S_mu) / (Tc + _S_mu);
  dmu_dT = (T < _T_min || T > _T_max) ? 0 : mu * (1.5 / Tc - 1 / (Tc + _S_mu));
}

// ---------------------------------------------------------------------------- (p, T)

Real
NTPHydrogenFluidProperties::rho_from_p_T(Real p, Real T) const
{
  // exact ideal gas, not clamped: intermediate Newton or segregated iterates may pass through
  // non-physical temperatures and recover
  return p / (_R_specific * T);
}

void
NTPHydrogenFluidProperties::rho_from_p_T(
    Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = rho_from_p_T(p, T);
  drho_dp = 1 / (_R_specific * T);
  drho_dT = -rho / T;
}

Real
NTPHydrogenFluidProperties::v_from_p_T(Real p, Real T) const
{
  return _R_specific * T / p;
}

void
NTPHydrogenFluidProperties::v_from_p_T(Real p, Real T, Real & v, Real & dv_dp, Real & dv_dT) const
{
  v = v_from_p_T(p, T);
  dv_dp = -v / p;
  dv_dT = _R_specific / p;
}

Real
NTPHydrogenFluidProperties::h_from_p_T(Real /*p*/, Real T) const
{
  return h_from_T(T);
}

void
NTPHydrogenFluidProperties::h_from_p_T(Real p, Real T, Real & h, Real & dh_dp, Real & dh_dT) const
{
  h = h_from_p_T(p, T);
  dh_dp = 0;
  dh_dT = cp_from_T(T);
}

Real
NTPHydrogenFluidProperties::e_from_p_T(Real /*p*/, Real T) const
{
  return h_from_T(T) - _R_specific * T;
}

void
NTPHydrogenFluidProperties::e_from_p_T(Real p, Real T, Real & e, Real & de_dp, Real & de_dT) const
{
  e = e_from_p_T(p, T);
  de_dp = 0;
  de_dT = cp_from_T(T) - _R_specific;
}

Real
NTPHydrogenFluidProperties::s_from_p_T(Real p, Real T) const
{
  const Real Tc = clamp(T);
  const auto & a = coeffs(Tc);
  const Real s0 =
      a[0] * std::log(Tc) + Tc * (a[1] + Tc * (a[2] / 2 + Tc * (a[3] / 3 + Tc * a[4] / 4))) + a[6];
  return _R_specific * (s0 - std::log(p / 1e5)) + cp_from_T(Tc) * std::log(T / Tc);
}

void
NTPHydrogenFluidProperties::s_from_p_T(Real p, Real T, Real & s, Real & ds_dp, Real & ds_dT) const
{
  s = s_from_p_T(p, T);
  ds_dp = -_R_specific / p;
  ds_dT = cp_from_T(T) / T;
}

Real
NTPHydrogenFluidProperties::cp_from_p_T(Real /*p*/, Real T) const
{
  return cp_from_T(T);
}

void
NTPHydrogenFluidProperties::cp_from_p_T(
    Real p, Real T, Real & cp, Real & dcp_dp, Real & dcp_dT) const
{
  cp = cp_from_p_T(p, T);
  dcp_dp = 0;
  dcp_dT = dcp_dT_from_T(T);
}

Real
NTPHydrogenFluidProperties::cv_from_p_T(Real /*p*/, Real T) const
{
  return cp_from_T(T) - _R_specific;
}

void
NTPHydrogenFluidProperties::cv_from_p_T(
    Real p, Real T, Real & cv, Real & dcv_dp, Real & dcv_dT) const
{
  cv = cv_from_p_T(p, T);
  dcv_dp = 0;
  dcv_dT = dcp_dT_from_T(T);
}

Real
NTPHydrogenFluidProperties::gamma_from_p_T(Real /*p*/, Real T) const
{
  const Real cp = cp_from_T(T);
  return cp / (cp - _R_specific);
}

void
NTPHydrogenFluidProperties::gamma_from_p_T(
    Real p, Real T, Real & gamma, Real & dgamma_dp, Real & dgamma_dT) const
{
  const Real cp = cp_from_T(T);
  const Real cv = cp - _R_specific;
  gamma = gamma_from_p_T(p, T);
  dgamma_dp = 0;
  dgamma_dT = -_R_specific * dcp_dT_from_T(T) / (cv * cv);
}

Real
NTPHydrogenFluidProperties::c_from_p_T(Real p, Real T) const
{
  return std::sqrt(gamma_from_p_T(p, T) * _R_specific * T);
}

void
NTPHydrogenFluidProperties::c_from_p_T(Real p, Real T, Real & c, Real & dc_dp, Real & dc_dT) const
{
  Real gamma, dgamma_dp, dgamma_dT;
  gamma_from_p_T(p, T, gamma, dgamma_dp, dgamma_dT);
  c = std::sqrt(gamma * _R_specific * T);
  dc_dp = 0;
  dc_dT = _R_specific * (gamma + T * dgamma_dT) / (2 * c);
}

Real
NTPHydrogenFluidProperties::mu_from_p_T(Real /*p*/, Real T) const
{
  Real mu, dmu_dT;
  muT(T, mu, dmu_dT);
  return mu;
}

void
NTPHydrogenFluidProperties::mu_from_p_T(
    Real /*p*/, Real T, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  muT(T, mu, dmu_dT);
  dmu_dp = 0;
}

Real
NTPHydrogenFluidProperties::k_from_p_T(Real p, Real T) const
{
  return mu_from_p_T(p, T) * (cp_from_T(T) + _eucken * _R_specific);
}

void
NTPHydrogenFluidProperties::k_from_p_T(
    Real /*p*/, Real T, Real & k, Real & dk_dp, Real & dk_dT) const
{
  Real mu, dmu_dT;
  muT(T, mu, dmu_dT);
  const Real cp = cp_from_T(T);
  k = mu * (cp + _eucken * _R_specific);
  dk_dp = 0;
  dk_dT = dmu_dT * (cp + _eucken * _R_specific) + mu * dcp_dT_from_T(T);
}

Real
NTPHydrogenFluidProperties::beta_from_p_T(Real /*p*/, Real T) const
{
  return 1 / T;
}

void
NTPHydrogenFluidProperties::beta_from_p_T(
    Real p, Real T, Real & beta, Real & dbeta_dp, Real & dbeta_dT) const
{
  beta = beta_from_p_T(p, T);
  dbeta_dp = 0;
  dbeta_dT = -1 / (T * T);
}

// ---------------------------------------------------------------------------- (p, h)

Real
NTPHydrogenFluidProperties::T_from_p_h(Real /*p*/, Real h) const
{
  return T_from_h(h);
}

void
NTPHydrogenFluidProperties::T_from_p_h(Real p, Real h, Real & T, Real & dT_dp, Real & dT_dh) const
{
  T = T_from_p_h(p, h);
  dT_dp = 0;
  dT_dh = 1 / cp_from_T(T);
}

// ---------------------------------------------------------------------------- (v, e)

Real
NTPHydrogenFluidProperties::T_from_v_e(Real /*v*/, Real e) const
{
  // e = h - R T, so solve h(T) - R T = e with the same Newton scheme
  Real T = std::max(_T_min, (e - e_from_p_T(0, 300)) / (2.5 * _R_specific) + 300);
  for (unsigned int it = 0; it < 50; ++it)
  {
    const Real dT = (e_from_p_T(0, T) - e) / (cp_from_T(T) - _R_specific);
    T -= dT;
    if (std::abs(dT) < 1e-10 * std::max(T, 1.0))
      return T;
  }
  mooseException(name(), ": T_from_v_e did not converge for e = ", e);
}

void
NTPHydrogenFluidProperties::T_from_v_e(Real v, Real e, Real & T, Real & dT_dv, Real & dT_de) const
{
  T = T_from_v_e(v, e);
  dT_dv = 0;
  dT_de = 1 / (cp_from_T(T) - _R_specific);
}

Real
NTPHydrogenFluidProperties::p_from_v_e(Real v, Real e) const
{
  return _R_specific * T_from_v_e(v, e) / v;
}

void
NTPHydrogenFluidProperties::p_from_v_e(Real v, Real e, Real & p, Real & dp_dv, Real & dp_de) const
{
  const Real T = T_from_v_e(v, e);
  p = _R_specific * T / v;
  dp_dv = -p / v;
  dp_de = _R_specific / (v * (cp_from_T(T) - _R_specific));
}

Real
NTPHydrogenFluidProperties::c_from_v_e(Real v, Real e) const
{
  return c_from_p_T(0, T_from_v_e(v, e));
}

void
NTPHydrogenFluidProperties::c_from_v_e(Real v, Real e, Real & c, Real & dc_dv, Real & dc_de) const
{
  const Real T = T_from_v_e(v, e);
  Real dc_dp, dc_dT;
  c_from_p_T(0, T, c, dc_dp, dc_dT);
  dc_dv = 0;
  dc_de = dc_dT / (cp_from_T(T) - _R_specific);
}

Real
NTPHydrogenFluidProperties::cp_from_v_e(Real v, Real e) const
{
  return cp_from_T(T_from_v_e(v, e));
}

void
NTPHydrogenFluidProperties::cp_from_v_e(
    Real v, Real e, Real & cp, Real & dcp_dv, Real & dcp_de) const
{
  const Real T = T_from_v_e(v, e);
  cp = cp_from_T(T);
  dcp_dv = 0;
  dcp_de = dcp_dT_from_T(T) / (cp - _R_specific);
}

Real
NTPHydrogenFluidProperties::cv_from_v_e(Real v, Real e) const
{
  return cp_from_v_e(v, e) - _R_specific;
}

void
NTPHydrogenFluidProperties::cv_from_v_e(
    Real v, Real e, Real & cv, Real & dcv_dv, Real & dcv_de) const
{
  cp_from_v_e(v, e, cv, dcv_dv, dcv_de);
  cv -= _R_specific;
}

Real
NTPHydrogenFluidProperties::mu_from_v_e(Real v, Real e) const
{
  return mu_from_p_T(0, T_from_v_e(v, e));
}

void
NTPHydrogenFluidProperties::mu_from_v_e(
    Real v, Real e, Real & mu, Real & dmu_dv, Real & dmu_de) const
{
  const Real T = T_from_v_e(v, e);
  Real dmu_dT;
  muT(T, mu, dmu_dT);
  dmu_dv = 0;
  dmu_de = dmu_dT / (cp_from_T(T) - _R_specific);
}

Real
NTPHydrogenFluidProperties::k_from_v_e(Real v, Real e) const
{
  return k_from_p_T(0, T_from_v_e(v, e));
}

void
NTPHydrogenFluidProperties::k_from_v_e(Real v, Real e, Real & k, Real & dk_dv, Real & dk_de) const
{
  const Real T = T_from_v_e(v, e);
  Real dk_dp, dk_dT;
  k_from_p_T(0, T, k, dk_dp, dk_dT);
  dk_dv = 0;
  dk_de = dk_dT / (cp_from_T(T) - _R_specific);
}

Real
NTPHydrogenFluidProperties::e_from_v_h(Real /*v*/, Real h) const
{
  return h - _R_specific * T_from_h(h);
}

void
NTPHydrogenFluidProperties::e_from_v_h(Real v, Real h, Real & e, Real & de_dv, Real & de_dh) const
{
  const Real T = T_from_h(h);
  e = e_from_v_h(v, h);
  de_dv = 0;
  de_dh = 1 - _R_specific / cp_from_T(T);
}

Real
NTPHydrogenFluidProperties::e_from_p_rho(Real p, Real rho) const
{
  return e_from_p_T(p, p / (rho * _R_specific));
}

void
NTPHydrogenFluidProperties::e_from_p_rho(
    Real p, Real rho, Real & e, Real & de_dp, Real & de_drho) const
{
  const Real T = p / (rho * _R_specific);
  const Real cv = cp_from_T(T) - _R_specific;
  e = e_from_p_T(p, T);
  de_dp = cv / (rho * _R_specific);
  de_drho = -cv * T / rho;
}

// ---------------------------------------------------------------------------- AD with derivatives
//
// The value carries the full AD derivative through the chain rule; the returned partial
// derivatives are evaluated from the second derivatives where they are needed (rho), and are
// otherwise treated as constants, as in the SinglePhaseFluidProperties defaults.

namespace
{
/// Build an AD value from the Real value and partial derivatives in (a, b)
ADReal
chain(Real value, Real d_da, Real d_db, const ADReal & a, const ADReal & b)
{
  ADReal result = value;
  result.derivatives() = a.derivatives() * d_da + b.derivatives() * d_db;
  return result;
}
}

#define NTP_H2_AD_PT(want)                                                                         \
  void NTPHydrogenFluidProperties::want##_from_p_T(                                                \
      const ADReal & p, const ADReal & T, ADReal & val, ADReal & d_dp, ADReal & d_dT) const        \
  {                                                                                                \
    Real v, dp, dT;                                                                                \
    want##_from_p_T(p.value(), T.value(), v, dp, dT);                                              \
    val = chain(v, dp, dT, p, T);                                                                  \
    d_dp = dp;                                                                                     \
    d_dT = dT;                                                                                     \
  }

NTP_H2_AD_PT(h)
NTP_H2_AD_PT(e)
NTP_H2_AD_PT(cp)
NTP_H2_AD_PT(cv)
NTP_H2_AD_PT(mu)
NTP_H2_AD_PT(k)

#undef NTP_H2_AD_PT

void
NTPHydrogenFluidProperties::rho_from_p_T(
    const ADReal & p, const ADReal & T, ADReal & rho, ADReal & drho_dp, ADReal & drho_dT) const
{
  // rho = p / (R T) is simple enough to differentiate fully
  rho = p / (_R_specific * T);
  drho_dp = 1 / (_R_specific * T);
  drho_dT = -rho / T;
}

void
NTPHydrogenFluidProperties::T_from_p_h(
    const ADReal & p, const ADReal & h, ADReal & T, ADReal & dT_dp, ADReal & dT_dh) const
{
  Real v, dp, dh;
  T_from_p_h(p.value(), h.value(), v, dp, dh);
  T = chain(v, dp, dh, p, h);
  dT_dp = dp;
  dT_dh = dh;
}
