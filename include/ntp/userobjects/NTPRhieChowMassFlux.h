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

#pragma once

#include "RhieChowMassFlux.h"

/**
 * Rhie-Chow face mass flux of the segregated (SIMPLE / PIMPLE) solvers for flows with steps of
 * many orders of magnitude in the momentum diagonal, as between the free channels and the packed
 * beds of a particle-bed fuel element.
 *
 * With the cell momentum balance u_C = -M_C (grad p_C + f_C), M = V / A the cell mobility and
 * f = H / V the force density of the off-diagonal and source terms, the base class interpolates
 * rho M and rho M f linearly to the faces. Across a step in M this gives the face the mobility of
 * the more permeable side: a porous layer one cell thick loses its flow resistance. This object
 * replaces the internal-face coupling by the series form of the two half segments between the
 * cell centers (NTPSeriesInterpolation),
 *
 *   m_f = -(rho M)_f [ (grad p)_f + f_f ] . n,
 *   (rho M)_f = 1 / ((1 - gC) / (rho M)_C + gC / (rho M)_N),   f_f = (1 - gC) f_C + gC f_N,
 *
 * which is the exact face flux of two uniform half segments in series and reduces to the base
 * interpolation for smooth fields. Boundary faces keep the base treatment, except that on
 * boundaries with an extrapolated (zero normal gradient) pressure condition the face coefficient
 * is made isotropic ('isotropic_zero_flux_boundaries'): the anisotropic pressure diffusion would
 * otherwise carry its tangential part through walls not aligned with the axes.
 *
 * The cell pressure gradient of the momentum equations (Green-Gauss) takes the face pressure as
 * the linear average of the two cells. Across a mobility jump that face pressure is far from the
 * one of flux continuity, which is almost the pressure of the permeable cell, so the permeable
 * cell gets a large spurious gradient and velocity. With 'series_pressure_gradient' the gradient
 * is corrected, before each momentum assembly, with the face pressure
 * p_f = (w_C p_C + w_N p_N) / (w_C + w_N), w = rho M / d the half-segment conductances.
 *
 * It runs on NONLINEAR, i.e. before every linear system assembly of the segregated solvers, and
 * updates the coupling fields of the pressure assembly after each computation of H/A (the
 * executioners compute them right before assembling the pressure equation). Declare it under the
 * name of the flow Physics' Rhie-Chow object, by default 'ins_rhie_chow_interpolator'; the Physics
 * then does not add its own.
 */
class NTPRhieChowMassFlux : public RhieChowMassFlux
{
public:
  static InputParameters validParams();

  NTPRhieChowMassFlux(const InputParameters & params);

  /// The base class clears the coupling fields here; they are rebuilt with H/A instead
  virtual void initialize() override {}
  virtual void execute() override;

protected:
  /// Series face coupling of the internal faces (face coefficient and H/A flux)
  void computeSeriesCoupling();
  /// Make the face coefficient isotropic on boundaries with an extrapolated or symmetry pressure
  /// condition
  void makeZeroFluxBoundariesIsotropic();
  /// Replace the internal-face coefficient by its normal part n . (A n)
  void makeInternalFacesNormal();
  /// Correct the cell pressure gradient with the series face pressures, once per recomputation
  void correctPressureGradient();

  /// Whether the internal-face coupling uses the series form
  const bool _series;
  /// Whether the coefficients of zero-flux pressure boundaries are made isotropic
  const bool _isotropic_walls;
  /// Whether the cell pressure gradient uses the series face pressures
  const bool _series_gradient;
  /// Whether the internal-face coefficient is reduced to its normal part
  const bool _normal_coupling;
  /// The gradient vector corrected last (the gradients alternate between two buffers)
  const NumericVector<Number> * _last_corrected_gradient = nullptr;
};
