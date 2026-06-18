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

#include "CorrosionPlatingFlowAction.h"
#include "FEProblemBase.h"
#include "Factory.h"

#include "MoltenSaltCorrosionModel.h"
#include "CorrosionChemistry.h"

#include <cctype>
#include <cmath>

registerMooseAction("OpenPronghornApp", CorrosionPlatingFlowAction, "add_variable");
registerMooseAction("OpenPronghornApp", CorrosionPlatingFlowAction, "add_ic");
registerMooseAction("OpenPronghornApp", CorrosionPlatingFlowAction, "add_linear_fv_kernel");
registerMooseAction("OpenPronghornApp", CorrosionPlatingFlowAction, "add_linear_fv_bc");

namespace
{
std::string
canonicalElementToken(const std::string & token)
{
  std::string s = token;
  for (std::size_t i = 0; i < s.size(); ++i)
    s[i] = (i == 0) ? std::toupper(static_cast<unsigned char>(s[i]))
                    : std::tolower(static_cast<unsigned char>(s[i]));
  return s;
}
}

InputParameters
CorrosionPlatingFlowAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Sets up salt-side molten salt corrosion (linear finite-volume passive "
                             "scalars plus a Butler-Volmer wall reaction) for a flowing MSR.");

  params.addParam<DataFileName>("database",
                                "corrosion_database.json",
                                "The JSON corrosion database. Defaults to data/corrosion_database.json.");

  MultiMooseEnum elements("Cr Fe Ni Mo Nb Tc Ru Ag Sb Te", "Cr Fe Ni");
  params.addParam<MultiMooseEnum>("elements", elements, "The metal elements to track.");
  params.addParam<std::vector<std::string>>(
      "release_variables",
      {},
      "Existing variable each element's corrosion is released into (one per element), instead of "
      "creating a new c_<El>. Use this to feed the wall corrosion product into a radiolysis-tracked "
      "cation (e.g. 'Cr_II'), coupling corrosion and radiolysis; the named variable's transport and "
      "chemistry are then owned by the radiolysis action.");

  params.addRequiredParam<MooseFunctorName>("temperature", "Temperature functor [K].");
  params.addParam<Real>("reference_temperature",
                        923.15,
                        "Scalar temperature [K] for the kinetics seed and the Arrhenius anchor.");
  params.addParam<bool>("temperature_dependent_kinetics",
                        true,
                        "Make the exchange current Arrhenius in the local temperature field.");

  params.addParam<std::string>("material_class", "hastelloy_n", "Alloy class.");
  params.addParam<std::string>("salt_class", "fluoride_fuel", "Salt class.");
  params.addParam<std::string>("redox_class", "oxidizing_fef2", "Redox/overpotential class.");
  params.addParam<std::string>("position_class", "nominal", "Loop position class.");
  params.addParam<Real>("flow_factor", 1.0, "Circulation/mass-transfer factor.");
  params.addParam<Real>("delta_T_C", 0.0, "Loop thermal gradient [C].");
  params.addParam<Real>("applied_overpotential", 0.1, "Imposed metal-minus-salt potential [V].");
  params.addParam<Real>("reference_concentration", 1.0, "Reference concentration c_ref [mol/m^3].");

  params.addRequiredParam<std::vector<BoundaryName>>("reaction_boundary",
                                                     "The corroding wall boundary(ies).");

  params.addParam<UserObjectName>("rhie_chow_user_object",
                                  "The Rhie-Chow user object for advection by the segregated flow.");
  params.addParam<RealVectorValue>("velocity", "A prescribed constant advection velocity.");
  MooseEnum interp("average upwind", "upwind");
  params.addParam<MooseEnum>("advected_interp_method", interp, "Advected interpolation method.");
  params.addParam<MooseFunctorName>(
      "diffusivity", "Effective (molecular + turbulent) diffusivity functor [m^2/s].");

  params.addParam<std::vector<BoundaryName>>("inlet_boundary", {}, "Inlet boundary(ies).");
  params.addParam<Real>("inlet_concentration", 0.0, "Inlet cation concentration [mol/m^3].");
  params.addParam<std::vector<BoundaryName>>("outlet_boundary", {}, "Outlet boundary(ies).");

  params.addParam<Real>("initial_concentration", 0.0, "Initial salt-ion concentration [mol/m^3].");

  params.addParam<std::vector<SubdomainName>>("block", {}, "Blocks (default: whole mesh).");
  params.addParam<bool>("time_derivative", true, "Add a time derivative to each species.");
  params.addParam<bool>("verbose", false, "Print the generated species and kinetics.");
  return params;
}

CorrosionPlatingFlowAction::CorrosionPlatingFlowAction(const InputParameters & parameters)
  : Action(parameters),
    _temperature(getParam<MooseFunctorName>("temperature")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _temperature_dependent(getParam<bool>("temperature_dependent_kinetics")),
    _transient(getParam<bool>("time_derivative"))
{
  buildPlan();
}

void
CorrosionPlatingFlowAction::buildPlan()
{
  Corrosion::MoltenSaltCorrosionDatabase db(getParam<DataFileName>("database"));
  Corrosion::MoltenSaltCorrosionModel model(db);

  Corrosion::CorrosionFeatures features;
  features.material_class = getParam<std::string>("material_class");
  features.salt_class = getParam<std::string>("salt_class");
  features.redox_class = getParam<std::string>("redox_class");
  features.position_class = getParam<std::string>("position_class");
  features.flow_factor = getParam<Real>("flow_factor");
  features.delta_T_C = getParam<Real>("delta_T_C");
  features.temperature_K = _reference_temperature;

  const Real rate_um_y = model.corrosionRateUmY(features);
  const Real density = db.density(features.material_class);
  const Real eta = getParam<Real>("applied_overpotential");
  const Real c_ref = getParam<Real>("reference_concentration");
  const Real initial = getParam<Real>("initial_concentration");
  const Real f = Corrosion::faraday / (Corrosion::R_gas * _reference_temperature);

  const auto & elements = getParam<MultiMooseEnum>("elements");
  const auto & release = getParam<std::vector<std::string>>("release_variables");
  if (!release.empty() && release.size() != elements.size())
    paramError("release_variables",
               "Must list one variable per tracked element (or leave empty to create c_<El>).");

  unsigned int index = 0;
  for (const auto & el : elements)
  {
    const Corrosion::ElementProperties props = db.element(el.name());
    ElementPlan plan;
    plan.name = canonicalElementToken(el.name());
    plan.owns_variable = release.empty() || release[index].empty();
    plan.salt_var = plan.owns_variable ? "c_" + plan.name : release[index];
    plan.valence = props.valence;
    plan.molar_mass = props.molar_mass_g_mol;
    plan.alpha_a = props.alpha_a;
    plan.alpha_c = props.alpha_c;
    plan.E0 = props.E0_V;
    plan.c_ref = c_ref;
    plan.initial_salt = initial;

    // Seed the exchange current so that, at the applied overpotential and c = c_ref, the
    // Butler-Volmer current Faradaically equals the calibrated dissolution rate.
    const Real i_ref_a_m2 =
        Corrosion::umYToCorrosionCurrent(rate_um_y, plan.valence, plan.molar_mass, density) * 1.0e4;
    const Real bracket = std::exp(plan.alpha_a * plan.valence * f * eta) -
                         std::exp(-plan.alpha_c * plan.valence * f * eta);
    plan.i0 = (std::abs(bracket) > 1.0e-30) ? i_ref_a_m2 / bracket : i_ref_a_m2;

    _elements.push_back(plan);
    ++index;
  }

  if (getParam<bool>("verbose"))
  {
    _console << "[CorrosionPlatingFlow] reference corrosion rate " << rate_um_y << " um/y; "
             << _elements.size() << " elements\n";
    for (const auto & e : _elements)
      _console << "  " << e.name << " -> " << e.salt_var << "  i0=" << e.i0 << " A/m^2\n";
    _console << std::flush;
  }
}

void
CorrosionPlatingFlowAction::maybeAssignBlocks(InputParameters & params) const
{
  const auto & blocks = getParam<std::vector<SubdomainName>>("block");
  if (!blocks.empty())
    params.set<std::vector<SubdomainName>>("block") = blocks;
}

void
CorrosionPlatingFlowAction::act()
{
  if (_current_task == "add_variable")
    addVariables();
  else if (_current_task == "add_ic")
    addInitialConditions();
  else if (_current_task == "add_linear_fv_kernel")
    addKernels();
  else if (_current_task == "add_linear_fv_bc")
    addBoundaryConditions();
}

void
CorrosionPlatingFlowAction::addVariables()
{
  _problem->needFV();
  for (const auto & e : _elements)
  {
    if (!e.owns_variable)
      continue; // the variable (e.g. a radiolysis cation) is created elsewhere
    auto params = _factory.getValidParams("MooseLinearVariableFVReal");
    params.set<SolverSystemName>("solver_sys") = systemName(e.salt_var);
    maybeAssignBlocks(params);
    _problem->addVariable("MooseLinearVariableFVReal", e.salt_var, params);
  }
}

void
CorrosionPlatingFlowAction::addInitialConditions()
{
  for (const auto & e : _elements)
  {
    if (!e.owns_variable)
      continue;
    auto params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = e.salt_var;
    params.set<Real>("value") = e.initial_salt;
    maybeAssignBlocks(params);
    _problem->addInitialCondition("ConstantIC", e.salt_var + "_ic", params);
  }
}

void
CorrosionPlatingFlowAction::addKernels()
{
  const bool has_rhie_chow = isParamValid("rhie_chow_user_object");
  const bool has_velocity = isParamValid("velocity");
  const bool has_diffusion = isParamValid("diffusivity");

  for (const auto & e : _elements)
  {
    if (!e.owns_variable)
      continue; // transport and chemistry of the released variable are owned elsewhere

    if (_transient)
    {
      auto params = _factory.getValidParams("LinearFVTimeDerivative");
      params.set<LinearVariableName>("variable") = e.salt_var;
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("LinearFVTimeDerivative", "corr_" + e.salt_var + "_time", params);
    }

    if (has_rhie_chow)
    {
      auto params = _factory.getValidParams("LinearFVScalarAdvection");
      params.set<LinearVariableName>("variable") = e.salt_var;
      params.set<UserObjectName>("rhie_chow_user_object") =
          getParam<UserObjectName>("rhie_chow_user_object");
      params.set<MooseEnum>("advected_interp_method") = getParam<MooseEnum>("advected_interp_method");
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel(
          "LinearFVScalarAdvection", "corr_" + e.salt_var + "_advection", params);
    }
    else if (has_velocity)
    {
      auto params = _factory.getValidParams("LinearFVAdvection");
      params.set<LinearVariableName>("variable") = e.salt_var;
      params.set<RealVectorValue>("velocity") = getParam<RealVectorValue>("velocity");
      params.set<MooseEnum>("advected_interp_method") = getParam<MooseEnum>("advected_interp_method");
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("LinearFVAdvection", "corr_" + e.salt_var + "_advection", params);
    }

    if (has_diffusion)
    {
      auto params = _factory.getValidParams("LinearFVDiffusion");
      params.set<LinearVariableName>("variable") = e.salt_var;
      params.set<MooseFunctorName>("diffusion_coeff") = getParam<MooseFunctorName>("diffusivity");
      params.set<bool>("use_nonorthogonal_correction") = false;
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("LinearFVDiffusion", "corr_" + e.salt_var + "_diffusion", params);
    }
  }
}

void
CorrosionPlatingFlowAction::addBoundaryConditions()
{
  const auto & walls = getParam<std::vector<BoundaryName>>("reaction_boundary");
  const auto & inlets = getParam<std::vector<BoundaryName>>("inlet_boundary");
  const auto & outlets = getParam<std::vector<BoundaryName>>("outlet_boundary");
  const Real eta = getParam<Real>("applied_overpotential");
  const Real inlet_c = getParam<Real>("inlet_concentration");
  const bool use_temperature = _temperature_dependent;

  // The activation energy for the Arrhenius exchange current comes from the calibrated correlation.
  Corrosion::MoltenSaltCorrosionDatabase db(getParam<DataFileName>("database"));
  const Real Ea = use_temperature ? db.parameter("Ea_corr_kJ_mol") : 0.0;

  for (const auto & e : _elements)
  {
    {
      auto params = _factory.getValidParams("CorrosionLinearFVButlerVolmerBC");
      params.set<LinearVariableName>("variable") = e.salt_var;
      params.set<std::vector<BoundaryName>>("boundary") = walls;
      params.set<Real>("valence") = e.valence;
      params.set<Real>("exchange_current_density") = e.i0;
      params.set<Real>("E0") = e.E0;
      params.set<Real>("alpha_a") = e.alpha_a;
      params.set<Real>("alpha_c") = e.alpha_c;
      params.set<Real>("c_ref") = e.c_ref;
      params.set<Real>("reference_temperature") = _reference_temperature;
      params.set<Real>("activation_energy") = Ea;
      params.set<MooseFunctorName>("temperature") = _temperature;
      params.set<MooseFunctorName>("metal_potential") = std::to_string(eta);
      _problem->addLinearFVBC(
          "CorrosionLinearFVButlerVolmerBC", "corr_" + e.salt_var + "_wall", params);
    }

    if (!inlets.empty() && e.owns_variable)
    {
      auto params = _factory.getValidParams("LinearFVAdvectionDiffusionFunctorDirichletBC");
      params.set<LinearVariableName>("variable") = e.salt_var;
      params.set<std::vector<BoundaryName>>("boundary") = inlets;
      params.set<MooseFunctorName>("functor") = std::to_string(inlet_c);
      _problem->addLinearFVBC(
          "LinearFVAdvectionDiffusionFunctorDirichletBC", "corr_" + e.salt_var + "_inlet", params);
    }
    if (!outlets.empty() && e.owns_variable)
    {
      auto params = _factory.getValidParams("LinearFVAdvectionDiffusionOutflowBC");
      params.set<LinearVariableName>("variable") = e.salt_var;
      params.set<std::vector<BoundaryName>>("boundary") = outlets;
      params.set<bool>("use_two_term_expansion") = false;
      _problem->addLinearFVBC(
          "LinearFVAdvectionDiffusionOutflowBC", "corr_" + e.salt_var + "_outlet", params);
    }
  }
}
