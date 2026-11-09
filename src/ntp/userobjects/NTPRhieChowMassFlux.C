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

#include "NTPRhieChowMassFlux.h"
#include "NTPSeriesInterpolation.h"
#include "FEProblemBase.h"
#include "LinearSystem.h"
#include "MooseLinearVariableFV.h"
#include "PetscVectorReader.h"
#include "LinearFVExtrapolatedPressureBC.h"
#include "LinearFVPressureSymmetryBC.h"

registerMooseObject("OpenPronghornApp", NTPRhieChowMassFlux);

InputParameters
NTPRhieChowMassFlux::validParams()
{
  InputParameters params = RhieChowMassFlux::validParams();
  params.addClassDescription(
      "Rhie-Chow face mass flux of the segregated solvers with the series (half-segment "
      "resistance) face coupling, for steps of many orders of magnitude in the momentum diagonal "
      "such as packed beds next to free channels.");
  params.addParam<MooseEnum>("face_coupling",
                             MooseEnum("series base", "series"),
                             "Internal-face coupling: the series form, or the one of "
                             "RhieChowMassFlux ('pressure_diffusion_interpolation')");
  params.addParam<bool>(
      "series_pressure_gradient",
      true,
      "Correct the cell pressure gradient of the momentum equations with the face pressure of "
      "flux continuity (series weights) instead of the linear average. Across a jump of many "
      "orders of magnitude in the mobility the linear face pressure gives the permeable cell a "
      "large spurious gradient.");
  params.addParam<bool>(
      "isotropic_zero_flux_boundaries",
      true,
      "Use an isotropic face coefficient on boundaries with an extrapolated (zero normal "
      "gradient) or symmetry pressure condition, such as walls, fixed-velocity inlets and "
      "symmetry planes. With an "
      "anisotropic momentum diagonal, the tangential part of the anisotropic pressure diffusion "
      "otherwise gives a mass flux through walls whose normal is not aligned with the axes.");
  params.addParam<bool>(
      "normal_face_coupling",
      false,
      "Replace the face coefficient of internal faces by its normal part n . (A n), the same in all "
      "directions. The pressure diffusion then has no tangential (cross) term: with the "
      "anisotropic coefficient of LinearFVAnisotropicDiffusion the tangential part uses the linear "
      "average of the two cell pressure gradients, which across a mobility jump is inconsistent "
      "with the series coupling.");
  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum = {EXEC_NONLINEAR, EXEC_TIMESTEP_END};
  return params;
}

NTPRhieChowMassFlux::NTPRhieChowMassFlux(const InputParameters & params)
  : RhieChowMassFlux(params),
    _series(getParam<MooseEnum>("face_coupling") == "series"),
    _isotropic_walls(getParam<bool>("isotropic_zero_flux_boundaries")),
    _series_gradient(getParam<bool>("series_pressure_gradient")),
    _normal_coupling(getParam<bool>("normal_face_coupling"))
{
}

void
NTPRhieChowMassFlux::execute()
{
  if (_HbyA_raw.empty() || !_pressure_system)
    return;

  // before the momentum assembly: the cell pressure gradient of the momentum force; at the end of
  // a solve: the gradient and cell velocity of the last iteration, for the output
  const bool end_of_solve = _fe_problem.getCurrentExecuteOnFlag() == EXEC_TIMESTEP_END;
  if (_series_gradient && !_momentum_systems.empty() &&
      (end_of_solve || _fe_problem.currentLinearSysNum() == _momentum_systems[0]->number()))
    correctPressureGradient();
  if (end_of_solve)
    return;

  // only for the pressure assembly below
  if ((!_series && !_isotropic_walls && !_normal_coupling) ||
      _fe_problem.currentLinearSysNum() != _pressure_system->number())
    return;

  if (_isotropic_walls)
    makeZeroFluxBoundariesIsotropic();
  if (_series)
    computeSeriesCoupling();
  if (_normal_coupling)
    makeInternalFacesNormal();
}

void
NTPRhieChowMassFlux::computeSeriesCoupling()
{
  const auto time_arg = Moose::currentState();
  std::vector<PetscVectorReader> hbya_reader;
  std::vector<PetscVectorReader> ainv_reader;
  for (const auto i : index_range(_HbyA_raw))
  {
    hbya_reader.emplace_back(*_HbyA_raw[i]);
    ainv_reader.emplace_back(*_Ainv_raw[i]);
  }
  const auto sys_num = _global_momentum_system_numbers[0];

  for (const auto & fi : _fe_problem.mesh().faceInfo())
  {
    if (!hasBlocks(fi->elemPtr()->subdomain_id()) || !fi->neighborPtr() ||
        !hasBlocks(fi->neighborPtr()->subdomain_id()) || !_vel[0]->isInternalFace(*fi))
      continue;

    const auto elem_dof = fi->elemInfo()->dofIndices()[sys_num][0];
    const auto neighbor_dof = fi->neighborInfo()->dofIndices()[sys_num][0];
    const Real elem_rho = _rho(makeElemArg(fi->elemPtr()), time_arg);
    const Real neighbor_rho = _rho(makeElemArg(fi->neighborPtr()), time_arg);

    auto & Ainv = _Ainv[fi->id()];
    Real hbya_flux = 0;
    for (const auto i : index_range(_HbyA_raw))
    {
      // cell mobilities and force densities: u_C = -M_C (grad p_C + f_C)
      const Real M_elem = ainv_reader[i](elem_dof);
      const Real M_neighbor = ainv_reader[i](neighbor_dof);
      const Real f_elem = hbya_reader[i](elem_dof) / M_elem;
      const Real f_neighbor = hbya_reader[i](neighbor_dof) / M_neighbor;

      Ainv(i) =
          NTPSeriesInterpolation::coefficient(elem_rho * M_elem, neighbor_rho * M_neighbor, *fi);
      hbya_flux += fi->normal()(i) * Ainv(i) *
                   NTPSeriesInterpolation::segmentAverage(f_elem, f_neighbor, *fi);
    }
    _HbyA_flux[fi->id()] = hbya_flux;
  }
}

void
NTPRhieChowMassFlux::makeInternalFacesNormal()
{
  // a_f = n . (A_f n) in all directions: the face flux keeps its normal coefficient, and the
  // tangential term (n o A - (n . n o A) n) . grad p of the anisotropic diffusion vanishes
  for (const auto & fi : _fe_problem.mesh().faceInfo())
  {
    if (!hasBlocks(fi->elemPtr()->subdomain_id()) || !fi->neighborPtr() ||
        !hasBlocks(fi->neighborPtr()->subdomain_id()) || !_vel[0]->isInternalFace(*fi))
      continue;
    auto & Ainv = _Ainv[fi->id()];
    const auto & n = fi->normal();
    Real a = 0;
    for (const auto i : make_range(Moose::dim))
      a += n(i) * n(i) * Ainv(i);
    for (const auto i : make_range(Moose::dim))
      Ainv(i) = a;
  }
}

void
NTPRhieChowMassFlux::makeZeroFluxBoundariesIsotropic()
{
  // On a boundary with an extrapolated pressure condition the pressure diffusion has no normal
  // part; LinearFVAnisotropicDiffusion still adds the tangential part of the face tensor,
  // (n o A - (n . n o A) n) . grad p, a flux through the boundary unless A is isotropic. The normal
  // coefficient n . (n o A) is kept, so nothing else changes.
  for (const auto & fi : _fe_problem.mesh().faceInfo())
  {
    if (_vel[0]->isInternalFace(*fi))
      continue;
    const bool elem_in = hasBlocks(fi->elemPtr()->subdomain_id());
    const bool neighbor_in = fi->neighborPtr() && hasBlocks(fi->neighborPtr()->subdomain_id());
    if (!elem_in && !neighbor_in)
      continue;
    const auto & ids = fi->boundaryIDs();
    if (ids.empty())
      continue;
    const auto * const bc = _p->getBoundaryCondition(*ids.begin());
    if (!bc || (!dynamic_cast<const LinearFVExtrapolatedPressureBC *>(bc) &&
                !dynamic_cast<const LinearFVPressureSymmetryBC *>(bc)))
      continue;

    auto & Ainv = _Ainv[fi->id()];
    const auto & n = fi->normal();
    Real a = 0;
    for (const auto i : make_range(Moose::dim))
      a += n(i) * n(i) * Ainv(i);
    for (const auto i : make_range(Moose::dim))
      Ainv(i) = a;
  }
}

void
NTPRhieChowMassFlux::correctPressureGradient()
{
  // The pressure gradients are recomputed (into alternating buffers) after every pressure solve;
  // correct each new set once.
  auto & grad_container = _pressure_system->linearFVGradientContainer();
  if (grad_container.empty() || grad_container[0].get() == _last_corrected_gradient)
    return;
  _last_corrected_gradient = grad_container[0].get();

  const auto time_arg = Moose::currentState();
  const auto p_sys = _global_pressure_system_number;
  const auto m_sys = _global_momentum_system_numbers[0];
  const auto & p_solution = *_pressure_system->system().current_local_solution;

  // accumulate the corrections, then add them to the gradient vectors
  std::vector<std::vector<Real>> values(grad_container.size());
  std::vector<dof_id_type> dofs;
  {
    // the readers must be released before computeCellVelocity() below uses the vectors
    std::vector<PetscVectorReader> ainv_reader;
    for (const auto i : index_range(_Ainv_raw))
      ainv_reader.emplace_back(*_Ainv_raw[i]);
    for (const auto & fi : _fe_problem.mesh().faceInfo())
    {
      if (!hasBlocks(fi->elemPtr()->subdomain_id()) || !fi->neighborPtr() ||
          !hasBlocks(fi->neighborPtr()->subdomain_id()) || !_p->isInternalFace(*fi))
        continue;
      const auto & n = fi->normal();
      const auto elem_m = fi->elemInfo()->dofIndices()[m_sys][0];
      const auto neighbor_m = fi->neighborInfo()->dofIndices()[m_sys][0];
      // normal mobilities of the two cells, times density
      Real M_elem = 0, M_neighbor = 0;
      for (const auto i : index_range(_Ainv_raw))
      {
        M_elem += n(i) * n(i) * ainv_reader[i](elem_m);
        M_neighbor += n(i) * n(i) * ainv_reader[i](neighbor_m);
      }
      M_elem *= _rho(makeElemArg(fi->elemPtr()), time_arg);
      M_neighbor *= _rho(makeElemArg(fi->neighborPtr()), time_arg);
      if (M_elem <= 0 || M_neighbor <= 0)
        continue;

      // half-segment conductances and the face pressure of flux continuity
      const Real gc = fi->gC();
      const Real w_elem = M_elem / (1 - gc);
      const Real w_neighbor = M_neighbor / gc;
      const auto elem_p = fi->elemInfo()->dofIndices()[p_sys][0];
      const auto neighbor_p = fi->neighborInfo()->dofIndices()[p_sys][0];
      const Real p_elem = p_solution(elem_p);
      const Real p_neighbor = p_solution(neighbor_p);
      const Real p_series = (w_elem * p_elem + w_neighbor * p_neighbor) / (w_elem + w_neighbor);
      const Real p_linear = gc * p_elem + (1 - gc) * p_neighbor;
      const Real dp = p_series - p_linear;
      if (dp == 0)
        continue;

      const Real area = fi->faceArea() * fi->faceCoord();
      const Real V_elem = fi->elemInfo()->volume() * fi->elemInfo()->coordFactor();
      const Real V_neighbor = fi->neighborInfo()->volume() * fi->neighborInfo()->coordFactor();
      const bool own_elem = fi->elemPtr()->processor_id() == processor_id();
      const bool own_neighbor = fi->neighborPtr()->processor_id() == processor_id();
      for (const auto i : index_range(grad_container))
      {
        if (own_elem)
          values[i].push_back(n(i) * area * dp / V_elem);
        if (own_neighbor)
          values[i].push_back(-n(i) * area * dp / V_neighbor);
      }
      if (own_elem)
        dofs.push_back(elem_p);
      if (own_neighbor)
        dofs.push_back(neighbor_p);
    }
  }
  for (const auto i : index_range(grad_container))
  {
    auto & grad = const_cast<NumericVector<Number> &>(*grad_container[i]);
    grad.add_vector(values[i].data(), dofs);
    grad.close();
  }

  // the executioner reconstructed the cell velocity from the uncorrected gradient right after
  // the pressure solve; redo it with the corrected one
  computeCellVelocity();
}
