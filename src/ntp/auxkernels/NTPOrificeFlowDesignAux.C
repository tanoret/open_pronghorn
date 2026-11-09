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

#include "NTPOrificeFlowDesignAux.h"
#include "MeshDivision.h"
#include "FEProblemBase.h"
#include "Executioner.h"
#include "FixedPointSolve.h"

#include <cmath>
#include <numeric>
#include <limits>

registerMooseObject("OpenPronghornApp", NTPOrificeFlowDesignAux);

InputParameters
NTPOrificeFlowDesignAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Axial orifice profile of a particle-bed element, updated after every solve towards target "
      "slice mass flows, at a fixed geometric mean; the target flow shares start as the slice "
      "power shares and are corrected towards a uniform slice gas outlet temperature.");
  params.addRequiredParam<MeshDivisionName>("mesh_division", "Axial slices of the element");
  params.addRequiredParam<VectorPostprocessorName>(
      "slice_flows", "Vector postprocessor with the mass flow entering each slice's cold frit");
  params.addParam<std::string>("slice_flow_vector", "flow", "Name of the vector in 'slice_flows'");
  params.addRequiredParam<VectorPostprocessorName>(
      "slice_powers", "Vector postprocessor with the fission power of each slice");
  params.addParam<std::string>(
      "slice_power_vector", "value", "Name of the vector in 'slice_powers'");
  params.addRequiredParam<VectorPostprocessorName>(
      "slice_temperatures", "Vector postprocessor with the gas outlet temperature of each slice");
  params.addParam<std::string>(
      "slice_temperature_vector", "value", "Name of the vector in 'slice_temperatures'");
  params.addRequiredParam<Real>("inlet_temperature", "Gas inlet temperature [K]");
  params.addRangeCheckedParam<Real>(
      "flow_exponent",
      1.0,
      "flow_exponent > 0",
      "Exponent alpha of the flow matching ln chi <- ln chi + alpha ln(m / m_target)");
  params.addRangeCheckedParam<Real>("max_change_factor",
                                    1.5,
                                    "max_change_factor > 1",
                                    "Largest factor by which chi may change in one update");
  params.addParam<Real>("level", 0, "Logarithm of the geometric mean of chi");
  params.addRangeCheckedParam<Real>(
      "share_exponent",
      0.5,
      "share_exponent >= 0",
      "Exponent gamma of the share correction w <- w (rise / mean rise)^gamma; 0 keeps the "
      "power shares");
  params.addRangeCheckedParam<Real>("max_share_change",
                                    1.1,
                                    "max_share_change > 1",
                                    "Largest factor by which a share may change in one correction");
  params.addRangeCheckedParam<Real>(
      "share_tolerance",
      0.02,
      "share_tolerance > 0",
      "The shares are corrected only when every slice flow is within this relative distance of "
      "its target");
  params.addRangeCheckedParam<Real>(
      "minimum_mean_rise",
      1,
      "minimum_mean_rise > 0",
      "Smallest mean slice temperature rise [K] for a share correction");
  params.addParam<unsigned int>(
      "max_updates",
      std::numeric_limits<unsigned int>::max(),
      "Number of updates after which the profile is frozen, to test the designed orifice as a "
      "fixed one");
  params.addParam<bool>("verbose", false, "Print the flow mismatch and the profile per update");
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  return params;
}

NTPOrificeFlowDesignAux::NTPOrificeFlowDesignAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _division(_c_fe_problem.getMeshDivision(getParam<MeshDivisionName>("mesh_division"))),
    // broadcast: every rank applies the same update to its own copy of the profile
    _flows(getVectorPostprocessorValue(
        "slice_flows", getParam<std::string>("slice_flow_vector"), true)),
    _powers(getVectorPostprocessorValue(
        "slice_powers", getParam<std::string>("slice_power_vector"), true)),
    _T_slices(getVectorPostprocessorValue(
        "slice_temperatures", getParam<std::string>("slice_temperature_vector"), true)),
    _T_in(getParam<Real>("inlet_temperature")),
    _alpha(getParam<Real>("flow_exponent")),
    _max_factor(getParam<Real>("max_change_factor")),
    _level(getParam<Real>("level")),
    _gamma(getParam<Real>("share_exponent")),
    _max_share_factor(getParam<Real>("max_share_change")),
    _share_tolerance(getParam<Real>("share_tolerance")),
    _min_rise(getParam<Real>("minimum_mean_rise")),
    _max_updates(getParam<unsigned int>("max_updates")),
    _verbose(getParam<bool>("verbose")),
    _log_chi(declareRestartableData<std::vector<Real>>(
        "log_chi", std::vector<Real>(_division.getNumDivisions(), getParam<Real>("level")))),
    _shares(declareRestartableData<std::vector<Real>>("shares", std::vector<Real>())),
    _last_execution(declareRestartableData<std::pair<int, unsigned int>>(
        "last_execution", std::pair<int, unsigned int>(-1, 0))),
    _num_updates(declareRestartableData<unsigned int>("num_updates", 0))
{
  if (isNodal())
    paramError("variable", "Must be an elemental (constant monomial or FV) variable");
  for (const auto & flag : getParam<ExecFlagEnum>("execute_on"))
    if (flag != EXEC_INITIAL && flag != EXEC_TIMESTEP_END)
      paramError("execute_on", "Only 'initial' and 'timestep_end' are supported");
}

Real
NTPOrificeFlowDesignAux::matchFlows(std::vector<Real> & log_chi,
                                    const std::vector<Real> & flows,
                                    const std::vector<Real> & shares,
                                    const Real alpha,
                                    const Real max_factor,
                                    const Real level)
{
  const auto n = log_chi.size();
  mooseAssert(flows.size() == n && shares.size() == n, "One flow and one share per slice");
  const Real total = std::accumulate(flows.begin(), flows.end(), 0.0);
  if (!(total > 0))
    return 0;
  const Real max_step = std::log(max_factor);
  Real mismatch = 0;
  for (const auto j : make_range(n))
  {
    // a slice carrying more than its target is throttled (larger chi); a reversed or vanishing
    // flow is opened by the largest step
    const Real target = shares[j] * total;
    const Real ratio = std::log(std::max(flows[j], 1e-6 * target) / target);
    mismatch = std::max(mismatch, std::abs(ratio));
    log_chi[j] += std::min(std::max(alpha * ratio, -max_step), max_step);
  }
  const Real mean = std::accumulate(log_chi.begin(), log_chi.end(), 0.0) / n;
  for (auto & l : log_chi)
    l += level - mean;
  return mismatch;
}

bool
NTPOrificeFlowDesignAux::correctShares(std::vector<Real> & shares,
                                       const std::vector<Real> & T_slices,
                                       const Real T_in,
                                       const Real gamma,
                                       const Real max_factor,
                                       const Real min_rise)
{
  const auto n = shares.size();
  mooseAssert(T_slices.size() == n, "One temperature per slice");
  Real mean_rise = 0;
  for (const auto T : T_slices)
    mean_rise += (T - T_in) / n;
  if (mean_rise < min_rise || gamma == 0)
    return false;
  const Real max_step = std::log(max_factor);
  // a slice hotter than the mean needs a larger share of the flow
  for (const auto j : make_range(n))
  {
    const Real rise = std::max(T_slices[j] - T_in, 1e-3 * mean_rise);
    const Real step = gamma * std::log(rise / mean_rise);
    shares[j] *= std::exp(std::min(std::max(step, -max_step), max_step));
  }
  const Real sum = std::accumulate(shares.begin(), shares.end(), 0.0);
  for (auto & w : shares)
    w /= sum;
  return true;
}

void
NTPOrificeFlowDesignAux::update()
{
  const auto n = _log_chi.size();
  if (_flows.size() != n || _powers.size() != n || _T_slices.size() != n)
    mooseError(name(),
               ": the slice vectors have ",
               _flows.size(),
               " flows, ",
               _powers.size(),
               " powers and ",
               _T_slices.size(),
               " temperatures for ",
               n,
               " slices of 'mesh_division'");

  // the target shares start as the power shares
  const Real power = std::accumulate(_powers.begin(), _powers.end(), 0.0);
  if (_shares.empty())
  {
    if (!(power > 0))
      return;
    _shares.resize(n);
    for (const auto j : make_range(n))
      _shares[j] = _powers[j] / power;
  }

  const Real mismatch = matchFlows(_log_chi, _flows, _shares, _alpha, _max_factor, _level);
  bool corrected = false;
  if (mismatch < _share_tolerance)
    corrected =
        correctShares(_shares, _T_slices, _T_in, _gamma, _max_share_factor, _min_rise);

  if (_verbose)
  {
    _console << name() << ": flow mismatch " << mismatch
             << (corrected ? ", shares corrected" : "") << "\n  chi   ";
    for (const auto l : _log_chi)
      _console << " " << std::exp(l);
    _console << "\n  share ";
    for (const auto w : _shares)
      _console << " " << w;
    _console << std::endl;
  }
}

Real
NTPOrificeFlowDesignAux::computeValue()
{
  // updated once per solve: once per time step, or per fixed-point iteration when the executioner
  // iterates
  if (_c_fe_problem.getCurrentExecuteOnFlag() == EXEC_TIMESTEP_END)
  {
    auto & executioner = *_app.getExecutioner();
    const std::pair<int, unsigned int> execution(_c_fe_problem.timeStep(),
                                                 executioner.fixedPointSolve().numFixedPointIts());
    if (execution != _last_execution && _num_updates < _max_updates)
    {
      update();
      _last_execution = execution;
      ++_num_updates;
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
