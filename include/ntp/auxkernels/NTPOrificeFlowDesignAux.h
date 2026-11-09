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

#include "AuxKernel.h"

class MeshDivision;

/**
 * Orifice profile of a particle-bed element designed on the slice mass flows while the element is
 * solved by the segregated solvers. The cold-frit drag multiplier chi(z), one value per axial slice
 * of a mesh division, is updated at the end of every solve (every time step with PIMPLE, every
 * fixed-point iteration with SIMPLE) in two nested loops:
 *
 *  - flow matching (every update): the slice flows m_j, measured on the plenum / cold-frit
 *    interface (NTPRhieChowSliceFlow), are driven to the targets w_j sum_k m_k,
 *
 *      ln chi_j <- ln chi_j + alpha ln(m_j / (w_j sum_k m_k)),  then shifted to a mean 'level';
 *
 *    the flow through a slice falls monotonically with its drag, at the cold density of the frit,
 *    so this loop does not see the density feedback of the hot gas;
 *  - flow shares (when the flows match the targets to 'share_tolerance'): the shares w_j, which
 *    start as the slice fission power shares (a uniform outlet temperature with an equal enthalpy
 *    rise per unit mass), are corrected from the slice gas outlet temperatures,
 *
 *      w_j <- w_j (rise_j / mean rise)^gamma,  then normalized to sum 1.
 *
 * Each change of ln chi and ln w is limited to ln max_change_factor and ln max_share_change. After
 * 'max_updates' updates the profile is frozen.
 */
class NTPOrificeFlowDesignAux : public AuxKernel
{
public:
  static InputParameters validParams();

  NTPOrificeFlowDesignAux(const InputParameters & parameters);

  /**
   * One flow-matching update of log_chi (one value per slice) towards the target flows
   * shares * sum(flows); returns the largest |ln(m_j / target_j)| before the update.
   */
  static Real matchFlows(std::vector<Real> & log_chi,
                         const std::vector<Real> & flows,
                         const std::vector<Real> & shares,
                         Real alpha,
                         Real max_factor,
                         Real level);

  /**
   * One correction of the flow shares from the slice outlet temperatures; returns whether it
   * corrected (mean rise at least min_rise).
   */
  static bool correctShares(std::vector<Real> & shares,
                            const std::vector<Real> & T_slices,
                            Real T_in,
                            Real gamma,
                            Real max_factor,
                            Real min_rise);

protected:
  virtual Real computeValue() override;

  /// One design update from the current slice flows, powers and temperatures
  void update();

  /// Axial slices
  const MeshDivision & _division;
  /// Slice flows, powers and gas outlet temperatures
  const VectorPostprocessorValue & _flows;
  const VectorPostprocessorValue & _powers;
  const VectorPostprocessorValue & _T_slices;
  /// Inlet temperature
  const Real _T_in;
  /// Exponent of the flow matching
  const Real _alpha;
  /// Largest change of chi in one update, as a factor
  const Real _max_factor;
  /// Logarithm of the geometric mean of chi
  const Real _level;
  /// Exponent of the share correction (0: shares stay the power shares)
  const Real _gamma;
  /// Largest change of a share in one correction, as a factor
  const Real _max_share_factor;
  /// Flow mismatch below which the shares are corrected
  const Real _share_tolerance;
  /// Smallest mean slice temperature rise for a share correction
  const Real _min_rise;
  /// Number of updates after which the profile is frozen
  const unsigned int _max_updates;
  /// Whether to print the design history
  const bool _verbose;
  /// The profile, log chi per slice (restartable)
  std::vector<Real> & _log_chi;
  /// The flow shares (restartable; empty until the power is known)
  std::vector<Real> & _shares;
  /// (time step, fixed-point iteration) of the last update (restartable)
  std::pair<int, unsigned int> & _last_execution;
  /// Number of updates so far (restartable)
  unsigned int & _num_updates;
};
