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

#include "CorrosionPlatingAction.h"
#include "FEProblemBase.h"
#include "Factory.h"

#include "MoltenSaltCorrosionModel.h"
#include "CorrosionChemistry.h"

#include <cctype>
#include <cmath>
#include <map>
#include <sstream>

registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_variable");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_aux_variable");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_ic");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_kernel");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_interface_kernel");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_bc");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_aux_kernel");

namespace
{
// Format a Real with full precision so it can be used as a constant functor name.
std::string
numberToFunctor(Real value)
{
  std::ostringstream oss;
  oss.precision(17);
  oss << value;
  return oss.str();
}

// Normalize an element token (the MultiMooseEnum upper-cases its entries) to the canonical chemical
// symbol case, e.g. "CR" -> "Cr", so the generated variables are c_Cr, cs_Cr.
std::string
canonicalElement(const std::string & token)
{
  std::string s = token;
  for (std::size_t i = 0; i < s.size(); ++i)
    s[i] = (i == 0) ? std::toupper(static_cast<unsigned char>(s[i]))
                    : std::tolower(static_cast<unsigned char>(s[i]));
  return s;
}
}

InputParameters
CorrosionPlatingAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Creates the variables, transport kernels, potential equation and "
                             "Butler-Volmer electrode kinetics of a molten salt corrosion and "
                             "plating problem from the built-in corrosion database.");

  params.addParam<DataFileName>("database",
                                "corrosion_database.json",
                                "The JSON corrosion database (elements, material tables, calibrated "
                                "correlation parameters). Defaults to data/corrosion_database.json.");

  MultiMooseEnum elements("Cr Fe Ni Mo Nb Tc Ru Ag Sb Te", "Cr Fe Ni");
  params.addParam<MultiMooseEnum>("elements", elements, "The metal elements to track.");
  MooseEnum recession_element("Cr Fe Ni Mo Nb Tc Ru Ag Sb Te", "Cr");
  params.addParam<MooseEnum>(
      "recession_element", recession_element, "The element controlling the recession/plating output.");

  MooseEnum topology("two_block salt_only solid_only");
  params.addRequiredParam<MooseEnum>("topology", topology, "Which phases are modeled as mesh blocks.");
  params.addParam<std::vector<SubdomainName>>("salt_block", {}, "The salt (electrolyte) block(s).");
  params.addParam<std::vector<SubdomainName>>("solid_block", {}, "The solid (metal) block(s).");
  params.addParam<BoundaryName>("interface_boundary",
                                "The metal/salt interface sideset (two_block topology).");
  params.addParam<BoundaryName>("reaction_boundary",
                                "The external electrode boundary (salt_only/solid_only topology).");

  params.addRequiredParam<MooseFunctorName>("temperature", "The temperature functor [K].");
  params.addParam<Real>("reference_temperature",
                        923.15,
                        "Scalar temperature [K] for the kinetics seed and the Butler-Volmer "
                        "exponent. Set this to the case operating temperature so the interfacial "
                        "current reproduces the calibrated rate; defaults to 650 C.");

  // Case feature selectors driving the ported effective-correlation kinetics.
  params.addParam<std::string>("material_class", "hastelloy_n", "Alloy class.");
  params.addParam<std::string>("salt_class", "fluoride_fuel", "Salt class.");
  params.addParam<std::string>("redox_class", "purified_baseline", "Redox/overpotential class.");
  params.addParam<std::string>("position_class", "nominal", "Loop position class.");
  params.addParam<std::string>("surface_class", "metal", "Deposition surface class.");
  params.addParam<Real>("flow_factor", 1.0, "Circulation/mass-transfer factor.");
  params.addParam<Real>("delta_T_C", 0.0, "Loop thermal gradient [C].");

  // Electrode kinetics and operating point.
  params.addParam<Real>("applied_overpotential",
                        0.1,
                        "Imposed metal-minus-salt potential difference [V]. Also used to seed the "
                        "exchange current so the interfacial current reproduces the calibrated "
                        "dissolution rate at this overpotential.");
  params.addParam<MooseFunctorName>(
      "applied_potential",
      "Counter-phase potential functor [V] for single-domain topologies. Defaults to "
      "'applied_overpotential' when not supplied.");
  params.addParam<MooseFunctorName>(
      "counter_concentration",
      "Fixed salt-ion concentration functor [mol/m^3] of the absent salt phase (solid_only).");

  // Electric potential (electrophoresis) solve.
  params.addParam<bool>("solve_potential",
                        false,
                        "Solve the salt electric potential (current continuity) and include "
                        "electromigration. Requires 'pin_potential_boundary'.");
  params.addParam<bool>("supporting_electrolyte",
                        true,
                        "Treat the salt as a supporting electrolyte: the Butler-Volmer current is "
                        "not fed back into the potential equation (the dilute tracers only migrate "
                        "in the field). Set false to couple the interfacial current to the salt "
                        "charge balance.");
  params.addParam<BoundaryName>("pin_potential_boundary",
                                "Boundary at which the salt potential is pinned (fixes the "
                                "otherwise-undetermined constant mode).");
  params.addParam<Real>("pin_potential_value", 0.0, "Pinned salt potential value [V].");

  // Advection by an external (one-way) velocity field.
  params.addParam<MooseFunctorName>("velocity_x", "x-velocity functor [m/s] for salt-ion advection.");
  params.addParam<MooseFunctorName>("velocity_y", "y-velocity functor [m/s] for salt-ion advection.");
  params.addParam<MooseFunctorName>("velocity_z", "z-velocity functor [m/s] for salt-ion advection.");

  params.addParam<Real>("salt_diffusivity",
                        "Override the salt-ion diffusivity [m^2/s] for all elements (e.g. an "
                        "effective dispersion coefficient under flow). Defaults to the per-element "
                        "database value.");
  params.addParam<Real>("solid_diffusivity",
                        "Override the solid-state diffusivity [m^2/s] of the tracked metals in the "
                        "alloy. Defaults to the calibrated chromium diffusivity D_s(T).");

  // Initial concentrations.
  params.addParam<Real>("default_salt_concentration", 1.0, "Default salt-ion IC [mol/m^3].");
  params.addParam<Real>("default_solid_concentration", 1.0e4, "Default solid-metal IC [mol/m^3].");
  params.addParam<Real>("reference_concentration",
                        "Reference concentration c_ref [mol/m^3] at which the exchange current is "
                        "defined (the Butler-Volmer cathodic branch carries c/c_ref). Defaults to "
                        "the initial salt concentration; set it explicitly to start the salt below "
                        "c_ref (fresh salt) while keeping the calibrated rate at c_ref.");
  params.addParam<std::vector<std::string>>(
      "initial_condition_variables", {}, "Generated variable names with a specified initial value.");
  params.addParam<std::vector<Real>>(
      "initial_condition_values", {}, "Initial values matching 'initial_condition_variables'.");

  params.addParam<bool>("transient", true, "Add a time derivative to every transported species.");
  params.addParam<bool>("verbose", false, "Print the generated variables and kinetics.");

  return params;
}

CorrosionPlatingAction::CorrosionPlatingAction(const InputParameters & parameters)
  : Action(parameters),
    _topology(getParam<MooseEnum>("topology") == "two_block"
                  ? Topology::TwoBlock
                  : (getParam<MooseEnum>("topology") == "salt_only" ? Topology::SaltOnly
                                                                    : Topology::SolidOnly)),
    _temperature(getParam<MooseFunctorName>("temperature")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _solve_potential(getParam<bool>("solve_potential")),
    _supporting_electrolyte(getParam<bool>("supporting_electrolyte")),
    _transient(getParam<bool>("transient"))
{
  // Topology-dependent requirements.
  if (_topology == Topology::TwoBlock)
  {
    if (getParam<std::vector<SubdomainName>>("salt_block").empty() ||
        getParam<std::vector<SubdomainName>>("solid_block").empty())
      mooseError("CorrosionPlating: 'two_block' topology requires both 'salt_block' and "
                 "'solid_block'.");
    if (!isParamValid("interface_boundary"))
      mooseError("CorrosionPlating: 'two_block' topology requires 'interface_boundary'.");
  }
  else if (!isParamValid("reaction_boundary"))
    mooseError("CorrosionPlating: single-domain topology requires 'reaction_boundary'.");

  if (_solve_potential && !isParamValid("pin_potential_boundary"))
    mooseError("CorrosionPlating: 'solve_potential' requires 'pin_potential_boundary' to fix the "
               "constant mode of the potential.");

  buildPlan();
}

void
CorrosionPlatingAction::buildPlan()
{
  Corrosion::MoltenSaltCorrosionDatabase db(getParam<DataFileName>("database"));
  Corrosion::MoltenSaltCorrosionModel model(db);

  // Build the case features from the class selectors.
  Corrosion::CorrosionFeatures features;
  features.material_class = getParam<std::string>("material_class");
  features.salt_class = getParam<std::string>("salt_class");
  features.redox_class = getParam<std::string>("redox_class");
  features.position_class = getParam<std::string>("position_class");
  features.surface_class = getParam<std::string>("surface_class");
  features.flow_factor = getParam<Real>("flow_factor");
  features.delta_T_C = getParam<Real>("delta_T_C");
  // Evaluate the calibrated rate at the scalar reference temperature (the field temperature functor
  // drives the spatial transport coefficients separately).
  features.temperature_K = _reference_temperature;

  const Real rate_um_y = model.corrosionRateUmY(features);
  const Real density = db.density(features.material_class);
  const Real solid_diffusivity = isParamValid("solid_diffusivity")
                                     ? getParam<Real>("solid_diffusivity")
                                     : model.crDiffusivityM2S(features);
  const Real eta = getParam<Real>("applied_overpotential");

  // Initial-condition overrides keyed by generated variable name.
  const auto & ic_vars = getParam<std::vector<std::string>>("initial_condition_variables");
  const auto & ic_vals = getParam<std::vector<Real>>("initial_condition_values");
  if (ic_vars.size() != ic_vals.size())
    paramError("initial_condition_values",
               "'initial_condition_variables' and 'initial_condition_values' must match in length.");
  std::map<std::string, Real> ic;
  for (const auto i : index_range(ic_vars))
    ic[ic_vars[i]] = ic_vals[i];
  const Real default_salt = getParam<Real>("default_salt_concentration");
  const Real default_solid = getParam<Real>("default_solid_concentration");

  const std::string recession_name = canonicalElement(getParam<MooseEnum>("recession_element"));

  for (const auto & el : getParam<MultiMooseEnum>("elements"))
  {
    const Corrosion::ElementProperties props = db.element(el.name());
    ElementPlan plan;
    plan.name = canonicalElement(el.name());
    plan.salt_var = "c_" + plan.name;
    plan.solid_var = "cs_" + plan.name;
    plan.valence = props.valence;
    plan.molar_mass = props.molar_mass_g_mol;
    plan.diffusivity =
        isParamValid("salt_diffusivity") ? getParam<Real>("salt_diffusivity") : props.diffusivity_m2_s;
    plan.solid_diffusivity = solid_diffusivity;
    plan.alpha_a = props.alpha_a;
    plan.alpha_c = props.alpha_c;
    plan.E0 = props.E0_V;
    plan.initial_salt = ic.count(plan.salt_var) ? ic.at(plan.salt_var) : default_salt;
    plan.initial_solid = ic.count(plan.solid_var) ? ic.at(plan.solid_var) : default_solid;
    // The exchange current is defined at c_ref; default it to the initial salt concentration so the
    // reproduction bridge is exact at the initial state, or take an explicit reference.
    if (isParamValid("reference_concentration"))
      plan.c_ref = getParam<Real>("reference_concentration");
    else
      plan.c_ref = plan.initial_salt > 0.0 ? plan.initial_salt : 1.0;

    // Seed the exchange current density so the Butler-Volmer current at the applied overpotential
    // Faradaically equals the calibrated dissolution rate (i0 = i_ref / bracket(eta)).
    const Real i_ref_a_cm2 =
        Corrosion::umYToCorrosionCurrent(rate_um_y, plan.valence, plan.molar_mass, density);
    const Real i_ref_a_m2 = i_ref_a_cm2 * 1.0e4;
    const Real f = Corrosion::faraday / (Corrosion::R_gas * _reference_temperature);
    const Real bracket = std::exp(plan.alpha_a * plan.valence * f * eta) -
                         std::exp(-plan.alpha_c * plan.valence * f * eta);
    plan.i0 = (std::abs(bracket) > 1.0e-30) ? i_ref_a_m2 / bracket : i_ref_a_m2;

    if (plan.name == recession_name)
      _recession_index = _elements.size();
    _elements.push_back(plan);
  }

  if (getParam<bool>("verbose"))
  {
    _console << "[CorrosionPlating] reference corrosion rate " << rate_um_y << " um/y; "
             << _elements.size() << " elements:\n";
    for (const auto & e : _elements)
      _console << "  " << e.name << " z=" << e.valence << " D=" << e.diffusivity
               << " i0=" << e.i0 << " A/m^2\n";
    _console << std::flush;
  }
}

void
CorrosionPlatingAction::assignBlocks(InputParameters & params, bool solid) const
{
  const auto & blocks = solid ? getParam<std::vector<SubdomainName>>("solid_block")
                              : getParam<std::vector<SubdomainName>>("salt_block");
  if (!blocks.empty())
    params.set<std::vector<SubdomainName>>("block") = blocks;
}

std::string
CorrosionPlatingAction::saltPotential() const
{
  return "phi_salt";
}

std::string
CorrosionPlatingAction::solidPotential() const
{
  return "phi_solid";
}

void
CorrosionPlatingAction::act()
{
  if (_current_task == "add_variable")
    addVariables();
  else if (_current_task == "add_aux_variable")
    addAuxVariables();
  else if (_current_task == "add_ic")
    addInitialConditions();
  else if (_current_task == "add_kernel")
    addKernels();
  else if (_current_task == "add_interface_kernel")
    addInterfaceKernels();
  else if (_current_task == "add_bc")
    addBoundaryConditions();
  else if (_current_task == "add_aux_kernel")
    addAuxKernels();
}

void
CorrosionPlatingAction::addVariables()
{
  for (const auto & e : _elements)
  {
    if (modelSalt())
    {
      auto params = _factory.getValidParams("MooseVariable");
      assignBlocks(params, false);
      _problem->addVariable("MooseVariable", e.salt_var, params);
    }
    if (modelSolid())
    {
      auto params = _factory.getValidParams("MooseVariable");
      assignBlocks(params, true);
      _problem->addVariable("MooseVariable", e.solid_var, params);
    }
  }

  if (_solve_potential && modelSalt())
  {
    auto params = _factory.getValidParams("MooseVariable");
    assignBlocks(params, false);
    _problem->addVariable("MooseVariable", saltPotential(), params);
  }
}

void
CorrosionPlatingAction::addAuxVariables()
{
  // Recession/plating-thickness tracker and instantaneous penetration-rate readout on the
  // recession-controlling element.
  for (const std::string & name : {std::string("recession_um"), std::string("corrosion_rate_um_y")})
  {
    auto params = _factory.getValidParams("MooseVariable");
    params.set<MooseEnum>("family") = "LAGRANGE";
    params.set<MooseEnum>("order") = "FIRST";
    _problem->addAuxVariable("MooseVariable", name, params);
  }
}

void
CorrosionPlatingAction::addInitialConditions()
{
  auto add_ic = [&](const std::string & var, Real value, bool solid)
  {
    auto params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = var;
    params.set<Real>("value") = value;
    assignBlocks(params, solid);
    _problem->addInitialCondition("ConstantIC", var + "_ic", params);
  };

  for (const auto & e : _elements)
  {
    if (modelSalt())
      add_ic(e.salt_var, e.initial_salt, false);
    if (modelSolid())
      add_ic(e.solid_var, e.initial_solid, true);
  }

  if (_solve_potential && modelSalt())
  {
    auto params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = saltPotential();
    params.set<Real>("value") = getParam<Real>("pin_potential_value");
    assignBlocks(params, false);
    _problem->addInitialCondition("ConstantIC", saltPotential() + "_ic", params);
  }
}

void
CorrosionPlatingAction::addKernels()
{
  const bool has_velocity = isParamValid("velocity_x") || isParamValid("velocity_y") ||
                            isParamValid("velocity_z");

  for (const auto & e : _elements)
  {
    if (modelSalt())
    {
      if (_transient)
      {
        auto params = _factory.getValidParams("ADTimeDerivative");
        params.set<NonlinearVariableName>("variable") = e.salt_var;
        assignBlocks(params, false);
        _problem->addKernel("ADTimeDerivative", e.salt_var + "_time", params);
      }
      {
        auto params = _factory.getValidParams("CorrosionNernstPlanckFlux");
        params.set<NonlinearVariableName>("variable") = e.salt_var;
        params.set<Real>("diffusivity") = e.diffusivity;
        params.set<Real>("valence") = e.valence;
        params.set<MooseFunctorName>("temperature") = _temperature;
        if (_solve_potential)
          params.set<std::vector<VariableName>>("potential") = {saltPotential()};
        assignBlocks(params, false);
        _problem->addKernel("CorrosionNernstPlanckFlux", e.salt_var + "_transport", params);
      }
      if (has_velocity)
      {
        auto params = _factory.getValidParams("CorrosionAdvection");
        params.set<NonlinearVariableName>("variable") = e.salt_var;
        if (isParamValid("velocity_x"))
          params.set<MooseFunctorName>("vel_x") = getParam<MooseFunctorName>("velocity_x");
        if (isParamValid("velocity_y"))
          params.set<MooseFunctorName>("vel_y") = getParam<MooseFunctorName>("velocity_y");
        if (isParamValid("velocity_z"))
          params.set<MooseFunctorName>("vel_z") = getParam<MooseFunctorName>("velocity_z");
        assignBlocks(params, false);
        _problem->addKernel("CorrosionAdvection", e.salt_var + "_advection", params);
      }
    }

    if (modelSolid())
    {
      if (_transient)
      {
        auto params = _factory.getValidParams("ADTimeDerivative");
        params.set<NonlinearVariableName>("variable") = e.solid_var;
        assignBlocks(params, true);
        _problem->addKernel("ADTimeDerivative", e.solid_var + "_time", params);
      }
      {
        auto params = _factory.getValidParams("CorrosionSolidDiffusion");
        params.set<NonlinearVariableName>("variable") = e.solid_var;
        params.set<Real>("diffusivity") = e.solid_diffusivity;
        assignBlocks(params, true);
        _problem->addKernel("CorrosionSolidDiffusion", e.solid_var + "_diffusion", params);
      }
    }
  }

  // Current-continuity equation for the salt potential.
  if (_solve_potential && modelSalt())
  {
    std::vector<VariableName> ions;
    std::vector<Real> valences;
    std::vector<Real> diffusivities;
    for (const auto & e : _elements)
    {
      ions.push_back(e.salt_var);
      valences.push_back(e.valence);
      diffusivities.push_back(e.diffusivity);
    }
    auto params = _factory.getValidParams("CorrosionCurrentContinuity");
    params.set<NonlinearVariableName>("variable") = saltPotential();
    params.set<std::vector<VariableName>>("concentrations") = ions;
    params.set<std::vector<Real>>("valences") = valences;
    params.set<std::vector<Real>>("diffusivities") = diffusivities;
    params.set<MooseFunctorName>("temperature") = _temperature;
    assignBlocks(params, false);
    _problem->addKernel("CorrosionCurrentContinuity", "phi_salt_current", params);
  }
}

void
CorrosionPlatingAction::addInterfaceKernels()
{
  if (_topology != Topology::TwoBlock)
    return;

  const BoundaryName interface = getParam<BoundaryName>("interface_boundary");
  const Real eta = getParam<Real>("applied_overpotential");
  const Real temperature = _reference_temperature;

  for (const auto & e : _elements)
  {
    auto params = _factory.getValidParams("ButlerVolmerInterface");
    params.set<NonlinearVariableName>("variable") = e.salt_var;
    params.set<std::vector<VariableName>>("neighbor_var") = {e.solid_var};
    params.set<std::vector<BoundaryName>>("boundary") = {interface};
    params.set<Real>("valence") = e.valence;
    params.set<Real>("exchange_current_density") = e.i0;
    params.set<Real>("E0") = e.E0;
    params.set<Real>("alpha_a") = e.alpha_a;
    params.set<Real>("alpha_c") = e.alpha_c;
    params.set<Real>("temperature") = temperature;
    params.set<Real>("c_ref") = e.c_ref;
    // The metal is treated as equipotential at the applied overpotential; the salt potential is
    // solved (coupled) when requested.
    params.set<Real>("metal_potential") = eta;
    params.set<Real>("salt_potential") = 0.0;
    if (_solve_potential)
      params.set<std::vector<VariableName>>("phi_salt") = {saltPotential()};
    _problem->addInterfaceKernel("ButlerVolmerInterface", "bv_" + e.name, params);
  }
}

void
CorrosionPlatingAction::addBoundaryConditions()
{
  const Real eta = getParam<Real>("applied_overpotential");
  const Real temperature = _reference_temperature;
  const MooseFunctorName applied = isParamValid("applied_potential")
                                       ? getParam<MooseFunctorName>("applied_potential")
                                       : MooseFunctorName(numberToFunctor(eta));

  // Pin the salt potential constant mode.
  if (_solve_potential && modelSalt())
  {
    auto params = _factory.getValidParams("ADDirichletBC");
    params.set<NonlinearVariableName>("variable") = saltPotential();
    params.set<std::vector<BoundaryName>>("boundary") = {
        getParam<BoundaryName>("pin_potential_boundary")};
    params.set<Real>("value") = getParam<Real>("pin_potential_value");
    _problem->addBoundaryCondition("ADDirichletBC", "phi_salt_pin", params);
  }

  if (_topology == Topology::TwoBlock)
  {
    // The interfacial current is fed into the salt charge balance unless the salt is treated as a
    // supporting electrolyte.
    if (_solve_potential && !_supporting_electrolyte)
    {
      const BoundaryName interface = getParam<BoundaryName>("interface_boundary");
      for (const auto & e : _elements)
      {
        auto params = _factory.getValidParams("ButlerVolmerBC");
        params.set<NonlinearVariableName>("variable") = saltPotential();
        params.set<std::vector<BoundaryName>>("boundary") = {interface};
        params.set<MooseEnum>("flux_type") = "charge";
        params.set<std::vector<VariableName>>("concentration") = {e.salt_var};
        params.set<Real>("valence") = e.valence;
        params.set<Real>("exchange_current_density") = e.i0;
        params.set<Real>("E0") = e.E0;
        params.set<Real>("alpha_a") = e.alpha_a;
        params.set<Real>("alpha_c") = e.alpha_c;
        params.set<Real>("temperature") = temperature;
        params.set<Real>("c_ref") = e.c_ref;
        params.set<MooseFunctorName>("applied_potential") = numberToFunctor(eta);
        _problem->addBoundaryCondition("ButlerVolmerBC", "bv_charge_" + e.name, params);
      }
    }
    return;
  }

  // Single-domain topologies: Butler-Volmer species flux on the external electrode boundary.
  const BoundaryName boundary = getParam<BoundaryName>("reaction_boundary");
  const bool metal_domain = (_topology == Topology::SolidOnly);
  for (const auto & e : _elements)
  {
    auto params = _factory.getValidParams("ButlerVolmerBC");
    params.set<NonlinearVariableName>("variable") = metal_domain ? e.solid_var : e.salt_var;
    params.set<std::vector<BoundaryName>>("boundary") = {boundary};
    params.set<MooseEnum>("flux_type") = "species";
    params.set<bool>("metal_domain") = metal_domain;
    params.set<Real>("valence") = e.valence;
    params.set<Real>("exchange_current_density") = e.i0;
    params.set<Real>("E0") = e.E0;
    params.set<Real>("alpha_a") = e.alpha_a;
    params.set<Real>("alpha_c") = e.alpha_c;
    params.set<Real>("temperature") = temperature;
    params.set<Real>("c_ref") = e.c_ref;
    params.set<MooseFunctorName>("applied_potential") = applied;
    if (metal_domain && isParamValid("counter_concentration"))
      params.set<MooseFunctorName>("counter_concentration") =
          getParam<MooseFunctorName>("counter_concentration");
    if (_solve_potential && !metal_domain)
      params.set<std::vector<VariableName>>("potential") = {saltPotential()};
    _problem->addBoundaryCondition("ButlerVolmerBC", "bv_" + e.name, params);
  }
}

void
CorrosionPlatingAction::addAuxKernels()
{
  // Recession / plating-thickness output on the controlling element's electrode boundary.
  const ElementPlan & e = _elements[_recession_index];
  Corrosion::MoltenSaltCorrosionDatabase db(getParam<DataFileName>("database"));
  const Real density = db.density(getParam<std::string>("material_class"));
  const Real eta = getParam<Real>("applied_overpotential");

  const bool metal_domain = (_topology == Topology::SolidOnly);
  const VariableName conc = metal_domain ? e.solid_var : e.salt_var;
  const BoundaryName boundary = (_topology == Topology::TwoBlock)
                                    ? getParam<BoundaryName>("interface_boundary")
                                    : getParam<BoundaryName>("reaction_boundary");

  auto make_aux = [&](const std::string & var,
                      const std::string & mode,
                      const std::string & name,
                      const std::string & execute_on)
  {
    auto params = _factory.getValidParams("CorrosionRateAux");
    params.set<AuxVariableName>("variable") = var;
    params.set<std::vector<BoundaryName>>("boundary") = {boundary};
    params.set<MooseEnum>("mode") = mode;
    params.set<ExecFlagEnum>("execute_on") = execute_on;
    params.set<std::vector<VariableName>>("concentration") = {conc};
    params.set<Real>("valence") = e.valence;
    params.set<Real>("molar_mass") = e.molar_mass;
    params.set<Real>("density") = density;
    params.set<Real>("exchange_current_density") = e.i0;
    params.set<Real>("E0") = e.E0;
    params.set<Real>("alpha_a") = e.alpha_a;
    params.set<Real>("alpha_c") = e.alpha_c;
    params.set<Real>("temperature") = _reference_temperature;
    params.set<Real>("c_ref") = e.c_ref;
    params.set<Real>("metal_potential_value") = eta;
    if (_solve_potential && modelSalt())
      params.set<std::vector<VariableName>>("salt_potential") = {saltPotential()};
    else
      params.set<Real>("salt_potential_value") = 0.0;
    _problem->addAuxKernel("CorrosionRateAux", name, params);
  };

  // The recession integrates over time (steps only); the rate is a snapshot (also at the initial
  // state, where it equals the calibrated reference rate).
  make_aux("recession_um", "recession", "recession_aux", "TIMESTEP_END");
  make_aux("corrosion_rate_um_y", "penetration_rate", "corrosion_rate_aux", "INITIAL TIMESTEP_END");
}
