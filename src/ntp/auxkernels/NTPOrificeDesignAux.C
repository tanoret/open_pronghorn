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

#include "NTPOrificeDesignAux.h"
#include "MeshDivision.h"
#include "FEProblemBase.h"
#include "Executioner.h"
#include "FixedPointSolve.h"

#include <cmath>
#include <numeric>

registerMooseObject("OpenPronghornApp", NTPOrificeDesignAux);

InputParameters
NTPOrificeDesignAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Axial orifice profile of a particle-bed element, updated after every solve towards a "
      "uniform slice gas outlet temperature, at a fixed geometric mean.");
  params.addRequiredParam<MeshDivisionName>("mesh_division", "Axial slices of the element");
  params.addRequiredParam<VectorPostprocessorName>(
      "slice_temperatures", "Vector postprocessor with the gas outlet temperature of each slice");
  params.addParam<std::string>(
      "slice_temperature_vector", "value", "Name of the vector in 'slice_temperatures'");
  params.addRequiredParam<Real>("inlet_temperature", "Gas inlet temperature [K]");
  params.addRangeCheckedParam<Real>(
      "exponent",
      0.5,
      "exponent > 0",
      "Exponent beta of the update chi <- chi (mean rise / rise)^beta");
  params.addRangeCheckedParam<Real>("max_change_factor",
                                    1.2,
                                    "max_change_factor > 1",
                                    "Largest factor by which chi may change in one update");
  params.addParam<Real>("level", 0, "Logarithm of the geometric mean of chi");
  params.addRangeCheckedParam<Real>(
      "minimum_mean_rise",
      1,
      "minimum_mean_rise > 0",
      "Smallest mean slice temperature rise [K] for an update; below it the rises are noise");
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  return params;
}

NTPOrificeDesignAux::NTPOrificeDesignAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _division(_c_fe_problem.getMeshDivision(getParam<MeshDivisionName>("mesh_division"))),
    // broadcast: every rank applies the same update to its own copy of the profile
    _T_slices(getVectorPostprocessorValue(
        "slice_temperatures", getParam<std::string>("slice_temperature_vector"), true)),
    _T_in(getParam<Real>("inlet_temperature")),
    _beta(getParam<Real>("exponent")),
    _max_factor(getParam<Real>("max_change_factor")),
    _level(getParam<Real>("level")),
    _min_rise(getParam<Real>("minimum_mean_rise")),
    _log_chi(declareRestartableData<std::vector<Real>>(
        "log_chi", std::vector<Real>(_division.getNumDivisions(), getParam<Real>("level")))),
    _last_execution(declareRestartableData<std::pair<int, unsigned int>>(
        "last_execution", std::pair<int, unsigned int>(-1, 0)))
{
  if (isNodal())
    paramError("variable", "Must be an elemental (constant monomial or FV) variable");
  // the profile is updated at TIMESTEP_END only; other flags would only return it
  for (const auto & flag : getParam<ExecFlagEnum>("execute_on"))
    if (flag != EXEC_INITIAL && flag != EXEC_TIMESTEP_END)
      paramError("execute_on", "Only 'initial' and 'timestep_end' are supported");
}

bool
NTPOrificeDesignAux::updateProfile(std::vector<Real> & log_chi,
                                   const std::vector<Real> & T_slices,
                                   const Real T_in,
                                   const Real beta,
                                   const Real max_factor,
                                   const Real level,
                                   const Real min_rise)
{
  const auto n = log_chi.size();
  mooseAssert(T_slices.size() == n, "One temperature per slice");
  Real mean_rise = 0;
  for (const auto T : T_slices)
    mean_rise += (T - T_in) / n;
  if (mean_rise < min_rise)
    return false;
  const Real max_step = std::log(max_factor);
  for (const auto j : make_range(n))
  {
    // a slice hotter than the mean carries too little flow: open it (smaller chi)
    const Real rise = std::max(T_slices[j] - T_in, 1e-3 * mean_rise);
    const Real step = -beta * std::log(rise / mean_rise);
    log_chi[j] += std::min(std::max(step, -max_step), max_step);
  }
  const Real mean = std::accumulate(log_chi.begin(), log_chi.end(), 0.0) / n;
  for (auto & l : log_chi)
    l += level - mean;
  return true;
}

Real
NTPOrificeDesignAux::computeValue()
{
  // updated once per solve: once per time step, or per fixed-point iteration when the executioner
  // iterates
  if (_c_fe_problem.getCurrentExecuteOnFlag() == EXEC_TIMESTEP_END)
  {
    auto & executioner = *_app.getExecutioner();
    const std::pair<int, unsigned int> execution(_c_fe_problem.timeStep(),
                                                 executioner.fixedPointSolve().numFixedPointIts());
    if (execution != _last_execution)
    {
      if (_T_slices.size() != _log_chi.size())
        mooseError(name(),
                   ": 'slice_temperatures' has ",
                   _T_slices.size(),
                   " values for ",
                   _log_chi.size(),
                   " slices of 'mesh_division'");
      updateProfile(_log_chi, _T_slices, _T_in, _beta, _max_factor, _level, _min_rise);
      _last_execution = execution;
    }
  }
  const auto j = _division.divisionIndex(*_current_elem);
  if (j == MooseMeshDivision::INVALID_DIVISION_INDEX || j >= _log_chi.size())
    mooseError(name(),
               ": element ",
               _current_elem->id(),
               " (centroid ",
               _current_elem->vertex_average(),
               ") is outside the slices of 'mesh_division'");
  return std::exp(_log_chi[j]);
}
