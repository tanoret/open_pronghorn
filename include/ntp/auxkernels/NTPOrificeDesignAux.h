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
 * Orifice profile of a particle-bed element designed while the element is solved by the segregated
 * solvers, where the orifice cannot be an unknown of the (linear) systems: the cold-frit drag
 * multiplier chi(z), one value per axial slice of a mesh division, is updated at the end of every
 * solve (every time step with PIMPLE; with SIMPLE, every fixed-point iteration, so
 * fixed_point_max_its > 1 is needed for the design to iterate) towards a uniform slice gas outlet
 * temperature,
 *
 *   chi_j <- chi_j (mean rise / rise_j)^beta,   then rescaled to a geometric mean exp(level),
 *
 * rise_j being the temperature rise of slice j over the inlet. A slice hotter than the mean carries
 * too little flow and is opened; the level gives every slice room to be opened as well as
 * throttled (the same reference as NTPSliceUniformityConstraint with 'reference = mean'). Each
 * factor is limited to [1/max_change_factor, max_change_factor], and the profile is left unchanged
 * while the mean rise is below 'minimum_mean_rise'.
 */
class NTPOrificeDesignAux : public AuxKernel
{
public:
  static InputParameters validParams();

  NTPOrificeDesignAux(const InputParameters & parameters);

  /**
   * One update of the logarithm of the profile, log_chi (one value per slice), from the slice
   * temperatures. Returns whether it updated (not below the minimum mean rise).
   */
  static bool updateProfile(std::vector<Real> & log_chi,
                            const std::vector<Real> & T_slices,
                            Real T_in,
                            Real beta,
                            Real max_factor,
                            Real level,
                            Real min_rise);

protected:
  virtual Real computeValue() override;

  /// Axial slices
  const MeshDivision & _division;
  /// Slice gas outlet temperatures
  const VectorPostprocessorValue & _T_slices;
  /// Inlet temperature
  const Real _T_in;
  /// Exponent of the update (under-relaxation)
  const Real _beta;
  /// Largest change of chi in one update, as a factor
  const Real _max_factor;
  /// Logarithm of the geometric mean of chi
  const Real _level;
  /// Smallest mean temperature rise for which the profile is updated
  const Real _min_rise;
  /// The profile, log chi per slice (restartable)
  std::vector<Real> & _log_chi;
  /// (time step, fixed-point iteration) of the last update (restartable)
  std::pair<int, unsigned int> & _last_execution;
};
