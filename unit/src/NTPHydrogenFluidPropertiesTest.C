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

#include "MooseObjectUnitTest.h"
#include "NTPHydrogenFluidProperties.h"

#include <cmath>
#include <functional>

// These tests check NTPHydrogenFluidProperties against reference values of the same ideal-gas
// NASA hydrogen model, check every analytic derivative against a central finite difference, and
// check that the inversions T(h) and T(v, e) return the temperature they started from.

class NTPHydrogenFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  NTPHydrogenFluidPropertiesTest() : MooseObjectUnitTest("OpenPronghornApp")
  {
    InputParameters pars = _factory.getValidParams("NTPHydrogenFluidProperties");
    _fe_problem->addUserObject("NTPHydrogenFluidProperties", "h2", pars);
    _fp = &_fe_problem->getUserObject<NTPHydrogenFluidProperties>("h2");
  }

protected:
  /// Compare the analytic derivatives of f(a, b) with central differences
  void checkDerivatives(const std::function<void(Real, Real, Real &, Real &, Real &)> & f,
                        const std::function<Real(Real, Real)> & value,
                        Real a,
                        Real b,
                        Real tol = 1e-6)
  {
    Real val, df_da, df_db;
    f(a, b, val, df_da, df_db);
    EXPECT_NEAR(val, value(a, b), 1e-12 * std::abs(val));

    const Real da = 1e-6 * a;
    const Real db = 1e-6 * b;
    const Real fd_a = (value(a + da, b) - value(a - da, b)) / (2 * da);
    const Real fd_b = (value(a, b + db) - value(a, b - db)) / (2 * db);
    // the scale guards derivatives that are exactly zero (e.g. dh/dp of an ideal gas)
    const Real scale_a = std::max(std::abs(fd_a), std::abs(val / a));
    const Real scale_b = std::max(std::abs(fd_b), std::abs(val / b));
    EXPECT_NEAR(df_da, fd_a, tol * scale_a);
    EXPECT_NEAR(df_db, fd_b, tol * scale_b);
  }

  const NTPHydrogenFluidProperties * _fp;
};

TEST_F(NTPHydrogenFluidPropertiesTest, matchesNtpPbrIdealNasa)
{
  const Real p = 7e6;
  // reference model, ideal-gas NASA hydrogen at 7 MPa
  const Real T[] = {30, 300, 1000, 3000};
  const Real h[] = {-3.48182769172e6, 2.64700801625e4, 1.02617883227e7, 4.40144202483e7};
  const Real cp[] = {1.05865358097e4, 1.43117571454e4, 1.49627687967e4, 1.83867744524e4};
  const Real rho[] = {5.65727481872e1, 5.65727481872, 1.69718244562, 5.65727481872e-1};
  const Real mu[] = {1.01029190752e-6, 8.76e-6, 1.84999425020e-5, 3.35448484920e-5};
  const Real k[] = {1.59041560965e-2, 1.70534080509e-1, 3.72188733374e-1, 7.89725506362e-1};

  for (unsigned int i = 0; i < 4; ++i)
  {
    EXPECT_NEAR(_fp->h_from_p_T(p, T[i]), h[i], 1e-8 * std::abs(h[i]));
    EXPECT_NEAR(_fp->cp_from_p_T(p, T[i]), cp[i], 1e-8 * cp[i]);
    EXPECT_NEAR(_fp->rho_from_p_T(p, T[i]), rho[i], 1e-9 * rho[i]);
    EXPECT_NEAR(_fp->mu_from_p_T(p, T[i]), mu[i], 1e-9 * mu[i]);
    EXPECT_NEAR(_fp->k_from_p_T(p, T[i]), k[i], 1e-9 * k[i]);
  }
}

TEST_F(NTPHydrogenFluidPropertiesTest, thermodynamicConsistency)
{
  const Real R = 8.314462618 / 2.01588e-3;
  for (const Real T : {50.0, 800.0, 2500.0, 3500.0})
  {
    const Real p = 5e6;
    // ideal gas: cp - cv = R, e = h - R T, c^2 = gamma R T
    EXPECT_NEAR(_fp->cp_from_p_T(p, T) - _fp->cv_from_p_T(p, T), R, 1e-9 * R);
    EXPECT_NEAR(_fp->e_from_p_T(p, T), _fp->h_from_p_T(p, T) - R * T, 1e-6);
    const Real c = _fp->c_from_p_T(p, T);
    EXPECT_NEAR(c * c, _fp->gamma_from_p_T(p, T) * R * T, 1e-9 * c * c);
    // inversions
    EXPECT_NEAR(_fp->T_from_p_h(p, _fp->h_from_p_T(p, T)), T, 1e-8 * T);
    const Real v = _fp->v_from_p_T(p, T);
    EXPECT_NEAR(_fp->T_from_v_e(v, _fp->e_from_p_T(p, T)), T, 1e-8 * T);
    EXPECT_NEAR(_fp->p_from_v_e(v, _fp->e_from_p_T(p, T)), p, 1e-8 * p);
  }
}

TEST_F(NTPHydrogenFluidPropertiesTest, derivatives)
{
  const auto fp = _fp;
  // (p, T) pairs on both polynomial branches
  for (const Real T : {120.0, 2800.0})
  {
    const Real p = 6e6;
#define PT_CHECK(prop)                                                                             \
  checkDerivatives([fp](Real a, Real b, Real & v, Real & da, Real & db)                            \
                   { fp->prop##_from_p_T(a, b, v, da, db); },                                      \
                   [fp](Real a, Real b) { return fp->prop##_from_p_T(a, b); },                     \
                   p,                                                                              \
                   T)
    PT_CHECK(rho);
    PT_CHECK(v);
    PT_CHECK(h);
    PT_CHECK(e);
    PT_CHECK(s);
    PT_CHECK(cp);
    PT_CHECK(cv);
    PT_CHECK(c);
    PT_CHECK(mu);
    PT_CHECK(k);
    PT_CHECK(gamma);
    PT_CHECK(beta);
#undef PT_CHECK

    const Real v = _fp->v_from_p_T(p, T);
    const Real e = _fp->e_from_p_T(p, T);
#define VE_CHECK(prop)                                                                             \
  checkDerivatives([fp](Real a, Real b, Real & val, Real & da, Real & db)                          \
                   { fp->prop##_from_v_e(a, b, val, da, db); },                                    \
                   [fp](Real a, Real b) { return fp->prop##_from_v_e(a, b); },                     \
                   v,                                                                              \
                   e)
    VE_CHECK(p);
    VE_CHECK(T);
    VE_CHECK(c);
    VE_CHECK(cp);
    VE_CHECK(cv);
    VE_CHECK(mu);
    VE_CHECK(k);
#undef VE_CHECK

    const Real h = _fp->h_from_p_T(p, T);
    checkDerivatives([fp](Real a, Real b, Real & val, Real & da, Real & db)
                     { fp->T_from_p_h(a, b, val, da, db); },
                     [fp](Real a, Real b) { return fp->T_from_p_h(a, b); },
                     p,
                     h);
    checkDerivatives([fp](Real a, Real b, Real & val, Real & da, Real & db)
                     { fp->e_from_v_h(a, b, val, da, db); },
                     [fp](Real a, Real b) { return fp->e_from_v_h(a, b); },
                     v,
                     h);
    const Real rho = _fp->rho_from_p_T(p, T);
    checkDerivatives([fp](Real a, Real b, Real & val, Real & da, Real & db)
                     { fp->e_from_p_rho(a, b, val, da, db); },
                     [fp](Real a, Real b) { return fp->e_from_p_rho(a, b); },
                     p,
                     rho);
  }
}

TEST_F(NTPHydrogenFluidPropertiesTest, adDerivativesPropagate)
{
  // The AD overloads with derivative outputs must carry dval/dT through the AD value
  ADReal p = 6e6, T = 2500;
  Moose::derivInsert(T.derivatives(), 0, 1.0);
  ADReal h, dh_dp, dh_dT;
  _fp->h_from_p_T(p, T, h, dh_dp, dh_dT);
  const Real cp = _fp->cp_from_p_T(6e6, 2500);
  EXPECT_NEAR(h.derivatives()[0], cp, 1e-12 * cp);
  EXPECT_NEAR(MetaPhysicL::raw_value(dh_dT), cp, 1e-12 * cp);
}
