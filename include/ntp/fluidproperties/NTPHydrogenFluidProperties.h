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

#pragma once

#include "SinglePhaseFluidProperties.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

/**
 * Hot hydrogen gas for nuclear thermal propulsion: a thermally ideal, calorically imperfect gas.
 *
 *   rho = p / (R T),  R = R_u / M,  M = 2.01588e-3 kg/mol
 *   cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4                (NASA 7-coefficient form)
 *   h/R  = a1 T + a2 T^2/2 + a3 T^3/3 + a4 T^4/4 + a5 T^5/5 + a6,  e = h - R T
 *   s/R  = a1 ln T + a2 T + a3 T^2/2 + a4 T^3/3 + a5 T^4/4 + a7 - ln(p / 1e5)
 *   mu   = mu_0 (T/T_0)^1.5 (T_0 + S) / (T + S)                 (Sutherland)
 *   k    = mu (cp + 5/4 R)                                      (Eucken)
 *
 * The default coefficients are the GRI-Mech 3.0 fit for H2 (switch at 1000 K). With them the
 * object reproduces the ideal-gas NASA hydrogen model of the reference reduced-order channel model
 * exactly. Dissociation, ortho/para conversion and real-gas effects are not included. All
 * properties have analytic derivatives; T(e) and T(h) are found by Newton iteration on the monotone
 * polynomial.
 */
class NTPHydrogenFluidProperties : public SinglePhaseFluidProperties
{
public:
  static InputParameters validParams();

  NTPHydrogenFluidProperties(const InputParameters & parameters);

  using SinglePhaseFluidProperties::c_from_p_T;
  using SinglePhaseFluidProperties::cp_from_p_T;
  using SinglePhaseFluidProperties::cv_from_p_T;
  using SinglePhaseFluidProperties::e_from_p_T;
  using SinglePhaseFluidProperties::h_from_p_T;
  using SinglePhaseFluidProperties::k_from_p_T;
  using SinglePhaseFluidProperties::mu_from_p_T;
  using SinglePhaseFluidProperties::rho_from_p_T;
  using SinglePhaseFluidProperties::T_from_p_h;

  virtual std::string fluidName() const override { return "hydrogen"; }
  virtual Real molarMass() const override { return _molar_mass; }

  // (v, e) interface, used by compressible and weakly-compressible flow objects
  virtual Real p_from_v_e(Real v, Real e) const override;
  virtual void p_from_v_e(Real v, Real e, Real & p, Real & dp_dv, Real & dp_de) const override;
  virtual Real T_from_v_e(Real v, Real e) const override;
  virtual void T_from_v_e(Real v, Real e, Real & T, Real & dT_dv, Real & dT_de) const override;
  virtual Real c_from_v_e(Real v, Real e) const override;
  virtual void c_from_v_e(Real v, Real e, Real & c, Real & dc_dv, Real & dc_de) const override;
  virtual Real cp_from_v_e(Real v, Real e) const override;
  virtual void cp_from_v_e(Real v, Real e, Real & cp, Real & dcp_dv, Real & dcp_de) const override;
  virtual Real cv_from_v_e(Real v, Real e) const override;
  virtual void cv_from_v_e(Real v, Real e, Real & cv, Real & dcv_dv, Real & dcv_de) const override;
  virtual Real mu_from_v_e(Real v, Real e) const override;
  virtual void mu_from_v_e(Real v, Real e, Real & mu, Real & dmu_dv, Real & dmu_de) const override;
  virtual Real k_from_v_e(Real v, Real e) const override;
  virtual void k_from_v_e(Real v, Real e, Real & k, Real & dk_dv, Real & dk_de) const override;
  virtual Real e_from_v_h(Real v, Real h) const override;
  virtual void e_from_v_h(Real v, Real h, Real & e, Real & de_dv, Real & de_dh) const override;
  virtual Real e_from_p_rho(Real p, Real rho) const override;
  virtual void
  e_from_p_rho(Real p, Real rho, Real & e, Real & de_dp, Real & de_drho) const override;

  // (p, T) interface, used by the incompressible / weakly-compressible NSFV objects
  virtual Real rho_from_p_T(Real p, Real T) const override;
  virtual void
  rho_from_p_T(Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const override;
  virtual Real v_from_p_T(Real p, Real T) const override;
  virtual void v_from_p_T(Real p, Real T, Real & v, Real & dv_dp, Real & dv_dT) const override;
  virtual Real e_from_p_T(Real p, Real T) const override;
  virtual void e_from_p_T(Real p, Real T, Real & e, Real & de_dp, Real & de_dT) const override;
  virtual Real h_from_p_T(Real p, Real T) const override;
  virtual void h_from_p_T(Real p, Real T, Real & h, Real & dh_dp, Real & dh_dT) const override;
  virtual Real s_from_p_T(Real p, Real T) const override;
  virtual void s_from_p_T(Real p, Real T, Real & s, Real & ds_dp, Real & ds_dT) const override;
  virtual Real cp_from_p_T(Real p, Real T) const override;
  virtual void cp_from_p_T(Real p, Real T, Real & cp, Real & dcp_dp, Real & dcp_dT) const override;
  virtual Real cv_from_p_T(Real p, Real T) const override;
  virtual void cv_from_p_T(Real p, Real T, Real & cv, Real & dcv_dp, Real & dcv_dT) const override;
  virtual Real c_from_p_T(Real p, Real T) const override;
  virtual void c_from_p_T(Real p, Real T, Real & c, Real & dc_dp, Real & dc_dT) const override;
  virtual Real mu_from_p_T(Real p, Real T) const override;
  virtual void mu_from_p_T(Real p, Real T, Real & mu, Real & dmu_dp, Real & dmu_dT) const override;
  virtual Real k_from_p_T(Real p, Real T) const override;
  virtual void k_from_p_T(Real p, Real T, Real & k, Real & dk_dp, Real & dk_dT) const override;
  virtual Real beta_from_p_T(Real p, Real T) const override;
  virtual void
  beta_from_p_T(Real p, Real T, Real & beta, Real & dbeta_dp, Real & dbeta_dT) const override;
  virtual Real gamma_from_p_T(Real p, Real T) const override;
  virtual void
  gamma_from_p_T(Real p, Real T, Real & gamma, Real & dgamma_dp, Real & dgamma_dT) const override;

  // AD overloads that also return the derivatives, used by the NSFV enthalpy objects
  virtual void rho_from_p_T(const ADReal & p,
                            const ADReal & T,
                            ADReal & rho,
                            ADReal & drho_dp,
                            ADReal & drho_dT) const override;
  virtual void h_from_p_T(const ADReal & p,
                          const ADReal & T,
                          ADReal & h,
                          ADReal & dh_dp,
                          ADReal & dh_dT) const override;
  virtual void e_from_p_T(const ADReal & p,
                          const ADReal & T,
                          ADReal & e,
                          ADReal & de_dp,
                          ADReal & de_dT) const override;
  virtual void cp_from_p_T(const ADReal & p,
                           const ADReal & T,
                           ADReal & cp,
                           ADReal & dcp_dp,
                           ADReal & dcp_dT) const override;
  virtual void cv_from_p_T(const ADReal & p,
                           const ADReal & T,
                           ADReal & cv,
                           ADReal & dcv_dp,
                           ADReal & dcv_dT) const override;
  virtual void mu_from_p_T(const ADReal & p,
                           const ADReal & T,
                           ADReal & mu,
                           ADReal & dmu_dp,
                           ADReal & dmu_dT) const override;
  virtual void k_from_p_T(const ADReal & p,
                          const ADReal & T,
                          ADReal & k,
                          ADReal & dk_dp,
                          ADReal & dk_dT) const override;
  virtual void T_from_p_h(const ADReal & p,
                          const ADReal & h,
                          ADReal & T,
                          ADReal & dT_dp,
                          ADReal & dT_dh) const override;

  // (p, h) interface, used when solving for enthalpy
  virtual Real T_from_p_h(Real p, Real h) const override;
  virtual void T_from_p_h(Real p, Real h, Real & T, Real & dT_dp, Real & dT_dh) const override;

  /// Specific enthalpy as a function of temperature only [J/kg]
  Real h_from_T(Real T) const;
  /// Isobaric specific heat as a function of temperature only [J/(kg K)]
  Real cp_from_T(Real T) const;
  /// Derivative of cp with respect to temperature [J/(kg K^2)]
  Real dcp_dT_from_T(Real T) const;
  /// Temperature from specific enthalpy [K]
  Real T_from_h(Real h) const;

protected:
  /// Polynomial coefficients for T (low or high branch)
  const std::vector<Real> & coeffs(Real T) const { return T < _T_switch ? _a_low : _a_high; }
  /// Temperature clamped to the range of validity
  Real clamp(Real T) const;
  /// Dynamic viscosity and its temperature derivative
  void muT(Real T, Real & mu, Real & dmu_dT) const;

  const Real _molar_mass;
  const Real _R_specific;
  const std::vector<Real> & _a_low;
  const std::vector<Real> & _a_high;
  const Real _T_switch;
  const Real _mu_0;
  const Real _T_mu_0;
  const Real _S_mu;
  const Real _eucken;
  const Real _T_min;
  const Real _T_max;
};

#pragma GCC diagnostic pop
