//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FissionProductInventoryAction.h"
#include "FEProblemBase.h"
#include "Factory.h"
#include "LinearFVAdvectionUtils.h"

#include <algorithm>
#include <sstream>

registerMooseAction("OpenPronghornApp", FissionProductInventoryAction, "add_variable");
registerMooseAction("OpenPronghornApp",
                    FissionProductInventoryAction,
                    "add_interpolation_method_physics");
registerMooseAction("OpenPronghornApp", FissionProductInventoryAction, "add_ic");
registerMooseAction("OpenPronghornApp", FissionProductInventoryAction, "add_linear_fv_kernel");
registerMooseAction("OpenPronghornApp", FissionProductInventoryAction, "add_linear_fv_bc");

namespace
{
// Format a Real with full precision so it can be used as a constant functor name.
std::string
fpNumberToFunctor(Real value)
{
  std::ostringstream oss;
  oss.precision(17);
  oss << value;
  return oss.str();
}

const std::vector<std::string> fp_destinations = {"offgas", "plated", "processed"};
}

InputParameters
FissionProductInventoryAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Creates the variables and kernels of a fission-product inventory (Bateman depletion) "
      "problem in a circulating fuel salt, with MSR removal pathways (off-gas sparging and "
      "noble-metal plating), from the built-in nuclide database.");

  params.addParam<DataFileName>(
      "database",
      "fission_product_database.json",
      "The JSON nuclide database (half-lives, decay branchings, capture cross sections, fission "
      "yields, chemical groups). Defaults to the bundled data/fission_product_database.json "
      "(ENDF/B-VIII.0).");
  params.addRequiredParam<std::vector<std::string>>(
      "nuclides",
      "The nuclides to track (e.g. 'I135 Xe135'). Decay and capture links between "
      "tracked nuclides are added automatically.");

  params.addParam<std::string>(
      "yield_set", "U235_thermal", "The fission yield set of the database (e.g. U235_fast).");
  MooseEnum yield_type("effective cumulative independent", "effective");
  params.addParam<MooseEnum>(
      "yield_type",
      yield_type,
      "Yield used for the direct fission source. 'effective' (default) is the cumulative yield "
      "minus the part fed by tracked parents, so that no fission production is double counted and "
      "untracked short-lived precursors are lumped into their first tracked descendant. "
      "'cumulative' ignores tracked parents (use for a single nuclide in saturation) and "
      "'independent' uses the direct yield only.");

  params.addParam<MooseFunctorName>("fission_rate_density",
                                    "Volumetric fission rate density [fissions/(m^3 s)].");
  params.addParam<MooseFunctorName>(
      "power_density",
      "Volumetric fission power density [W/m^3]; converted to a fission rate density with "
      "'energy_per_fission'. Alternative to 'fission_rate_density'.");
  params.addRangeCheckedParam<Real>(
      "energy_per_fission", 200.0, "energy_per_fission>0", "Recoverable energy per fission [MeV].");
  params.addParam<MooseFunctorName>(
      "neutron_flux", "0", "One-group neutron flux [1/(m^2 s)] driving neutron capture.");

  params.addParam<MultiMooseEnum>(
      "removal_groups",
      MultiMooseEnum("noble_gas noble_metal halogen salt_seeking"),
      "Chemical groups with a first-order removal pathway (one entry per pathway; a group may be "
      "listed several times).");
  params.addParam<std::vector<MooseFunctorName>>(
      "removal_rates", {}, "First-order removal rate [1/s] of each pathway in 'removal_groups'.");
  params.addParam<MultiMooseEnum>(
      "removal_destinations",
      MultiMooseEnum("offgas plated processed"),
      "Where each removal pathway sends the inventory: 'offgas' (sparging / cover gas), 'plated' "
      "(metal surfaces) or 'processed' (chemical clean-up, not tracked).");
  params.addParam<bool>(
      "track_removed_inventories",
      true,
      "Create <nuclide>_offgas and <nuclide>_plated variables holding the removed inventory per "
      "unit salt volume [mol/m^3]; they keep decaying, so the decay of plated or stripped parents "
      "feeds the corresponding daughters.");

  params.addParam<std::vector<BoundaryName>>(
      "plating_boundaries",
      {},
      "Boundaries on which the noble metals plate with the mass-transfer flux J = k_m N "
      "(requires 'diffusivity').");
  params.addParam<MooseFunctorName>("plating_mass_transfer_coefficient",
                                    "Wall mass-transfer coefficient k_m [m/s] for "
                                    "'plating_boundaries'.");

  params.addParam<MooseFunctorName>(
      "diffusivity",
      "Diffusivity of the dissolved/entrained nuclides [m^2/s] (turbulent "
      "diffusivity in a flowing salt). Diffusion is added only when supplied.");
  params.addParam<UserObjectName>(
      "rhie_chow_user_object",
      "The Rhie-Chow user object providing the mass flux; when supplied the salt inventories are "
      "advected by the segregated flow solution.");
  params.addParam<RealVectorValue>(
      "velocity",
      "A prescribed constant advection velocity, used when 'rhie_chow_user_object' is not given.");
  OpenPronghorn::addAdvectedInterpolationParams(params);

  params.addParam<std::vector<std::string>>(
      "initial_condition_nuclides", {}, "Nuclides with a nonzero initial salt inventory.");
  params.addParam<std::vector<Real>>("initial_condition_values",
                                     {},
                                     "Initial inventories [mol/m^3] matching "
                                     "'initial_condition_nuclides'.");

  params.addParam<std::vector<SubdomainName>>(
      "block", {}, "The blocks the inventory problem is defined on (default: whole mesh).");
  params.addParam<bool>(
      "time_derivative", true, "Add a time-derivative kernel to every inventory.");
  params.addParam<bool>(
      "verbose", false, "Print the generated inventories, effective yields and system names.");
  return params;
}

FissionProductInventoryAction::FissionProductInventoryAction(const InputParameters & params)
  : Action(params),
    _db(getParam<DataFileName>("database")),
    _nuclides(getParam<std::vector<std::string>>("nuclides"))
{
  if (isParamValid("fission_rate_density") == isParamValid("power_density"))
    paramError("fission_rate_density",
               "Exactly one of 'fission_rate_density' and 'power_density' must be supplied.");

  const auto & y = _db.yieldSets();
  if (std::find(y.begin(), y.end(), getParam<std::string>("yield_set")) == y.end())
    paramError("yield_set",
               "Yield set '" + getParam<std::string>("yield_set") +
                   "' is not in the database; available: " + Moose::stringify(y));

  const auto n_groups = getParam<MultiMooseEnum>("removal_groups").size();
  if (getParam<std::vector<MooseFunctorName>>("removal_rates").size() != n_groups ||
      getParam<MultiMooseEnum>("removal_destinations").size() != n_groups)
    paramError("removal_rates",
               "'removal_groups', 'removal_rates' and 'removal_destinations' must have equal "
               "length.");
  {
    // Each (group, destination) pair carries one rate; a second one would be silently ignored.
    const auto & groups = getParam<MultiMooseEnum>("removal_groups");
    const auto & destinations = getParam<MultiMooseEnum>("removal_destinations");
    std::set<std::pair<std::string, std::string>> pathways;
    for (const auto i : make_range(n_groups))
      if (!pathways.emplace(MooseUtils::toLower(groups[i]), MooseUtils::toLower(destinations[i]))
               .second)
        paramError("removal_groups",
                   "The group '",
                   groups[i],
                   "' is removed to '",
                   destinations[i],
                   "' more than once. Combine the rates of such pathways into one functor.");
  }

  if (!getParam<std::vector<BoundaryName>>("plating_boundaries").empty())
  {
    if (!isParamValid("plating_mass_transfer_coefficient"))
      paramError("plating_mass_transfer_coefficient",
                 "Required when 'plating_boundaries' is supplied.");
    if (!isParamValid("diffusivity"))
      paramError(
          "plating_boundaries",
          "Wall plating is a boundary flux of the diffusion operator; supply 'diffusivity'.");
  }

  std::set<std::string> seen;
  for (const auto & n : _nuclides)
  {
    if (!_db.hasNuclide(n))
      paramError("nuclides",
                 "Nuclide '" + n + "' is not in the database '" + _db.filename() +
                     "'. Available: " + Moose::stringify(_db.nuclideNames()));
    if (!seen.insert(n).second)
      paramError("nuclides", "Nuclide '" + n + "' is listed more than once.");
  }

  buildPlan();
  orderPlan();
}

MooseFunctorName
FissionProductInventoryAction::groupRemoval(const std::string & group,
                                            const std::string & location) const
{
  const auto & g = getParam<MultiMooseEnum>("removal_groups");
  const auto & r = getParam<std::vector<MooseFunctorName>>("removal_rates");
  const auto & d = getParam<MultiMooseEnum>("removal_destinations");
  for (const auto i : make_range(g.size()))
    if (MooseUtils::toLower(g[i]) == group && MooseUtils::toLower(d[i]) == location)
      return r[i];
  return "";
}

void
FissionProductInventoryAction::buildPlan()
{
  const auto & ic_nuc = getParam<std::vector<std::string>>("initial_condition_nuclides");
  const auto & ic_val = getParam<std::vector<Real>>("initial_condition_values");
  if (ic_nuc.size() != ic_val.size())
    paramError("initial_condition_values",
               "'initial_condition_nuclides' and 'initial_condition_values' must have equal "
               "length.");
  std::map<std::string, Real> ic;
  for (const auto i : index_range(ic_nuc))
  {
    if (std::find(_nuclides.begin(), _nuclides.end(), ic_nuc[i]) == _nuclides.end())
      paramError("initial_condition_nuclides", "Nuclide '" + ic_nuc[i] + "' is not in 'nuclides'.");
    ic[ic_nuc[i]] = ic_val[i];
  }

  const std::set<std::string> tracked(_nuclides.begin(), _nuclides.end());
  const auto & yield_set = getParam<std::string>("yield_set");
  const std::string yield_type = getParam<MooseEnum>("yield_type");

  // Tracked parents feeding each tracked nuclide (decay branches and capture products)
  std::map<std::string, std::vector<ParentFeed>> feeds;
  for (const auto & p : _nuclides)
  {
    const auto & pd = _db.nuclide(p);
    for (const auto & branch : pd.decay)
      if (tracked.count(branch.daughter))
        feeds[branch.daughter].push_back({p, pd.decayConstant() * branch.branching, 0.0});
    if (tracked.count(pd.capture_product) && pd.capture_cross_section > 0.0)
      feeds[pd.capture_product].push_back({p, 0.0, pd.capture_cross_section});
  }

  // Salt inventories
  for (const auto & n : _nuclides)
  {
    const auto & d = _db.nuclide(n);
    InventoryPlan plan;
    plan.name = n;
    plan.nuclide = n;
    plan.location = "salt";
    plan.initial_condition = ic.count(n) ? ic.at(n) : 0.0;
    plan.lambda = d.decayConstant();
    plan.sigma = d.capture_cross_section;
    plan.feeds = feeds[n];
    plan.transported = true;

    const auto it = d.yields.find(yield_set);
    const Real independent = it != d.yields.end() ? it->second.first : 0.0;
    const Real cumulative = it != d.yields.end() ? it->second.second : 0.0;
    if (yield_type == "independent")
      plan.yield = independent;
    else if (yield_type == "cumulative")
      plan.yield = cumulative;
    else
    {
      // Cumulative yield minus the fission production that reaches this nuclide through tracked
      // parents by decay: CY_i = IY_i + sum_p b_pi CY_p, restricted to the tracked p.
      Real fed = 0.0;
      for (const auto & p : _nuclides)
        for (const auto & branch : _db.nuclide(p).decay)
          if (branch.daughter == n)
          {
            const auto pit = _db.nuclide(p).yields.find(yield_set);
            if (pit != _db.nuclide(p).yields.end())
              fed += branch.branching * pit->second.second;
          }
      plan.yield = std::max(cumulative - fed, independent);
    }

    for (const auto & dest : fp_destinations)
    {
      const auto r = groupRemoval(d.group, dest);
      if (!r.empty())
        plan.removal.push_back(r);
    }
    _plan.push_back(plan);
  }

  if (!getParam<bool>("track_removed_inventories"))
    return;

  // Removed inventories: every nuclide removed to a location, plus (recursively) its tracked
  // daughters, so that decay on the surface / in the off-gas system conserves the chain.
  for (const std::string location : {"offgas", "plated"})
  {
    std::set<std::string> members;
    for (const auto & n : _nuclides)
      if (!groupRemoval(_db.nuclide(n).group, location).empty())
        members.insert(n);
    bool grew = true;
    while (grew)
    {
      grew = false;
      for (const auto & m : std::set<std::string>(members))
        for (const auto & branch : _db.nuclide(m).decay)
          if (tracked.count(branch.daughter) && members.insert(branch.daughter).second)
            grew = true;
    }

    for (const auto & n : _nuclides)
    {
      if (!members.count(n))
        continue;
      const auto & d = _db.nuclide(n);
      InventoryPlan plan;
      plan.name = n + "_" + location;
      plan.nuclide = n;
      plan.location = location;
      plan.lambda = d.decayConstant();
      plan.transfer_rate = groupRemoval(d.group, location);
      // Removed inventories are out of the neutron flux: decay feeds only.
      for (const auto & f : feeds[n])
        if (f.decay_coefficient > 0.0 && members.count(f.parent))
          plan.feeds.push_back({f.parent + "_" + location, f.decay_coefficient, 0.0});
      _plan.push_back(plan);
    }
  }
}

void
FissionProductInventoryAction::orderPlan()
{
  // Topological sort, parent before daughter, stable with respect to the input order.
  std::vector<InventoryPlan> ordered;
  std::set<std::string> placed;
  auto remaining = _plan;
  while (!remaining.empty())
  {
    bool progress = false;
    for (auto it = remaining.begin(); it != remaining.end();)
    {
      bool ready = true;
      for (const auto & f : it->feeds)
        if (!placed.count(f.parent) && f.parent != it->name)
          ready = false;
      if (ready)
      {
        placed.insert(it->name);
        ordered.push_back(*it);
        it = remaining.erase(it);
        progress = true;
      }
      else
        ++it;
    }
    if (!progress)
    {
      // A cycle (e.g. through capture); keep the remaining input order. The fixed-point
      // iteration then converges the coupling.
      mooseWarning("FissionProductInventory: the nuclide links contain a cycle; the systems "
                   "cannot be ordered parent-before-daughter and need several fixed-point sweeps.");
      for (const auto & p : remaining)
        ordered.push_back(p);
      break;
    }
  }
  _plan = ordered;
}

void
FissionProductInventoryAction::maybeAssignBlocks(InputParameters & params) const
{
  const auto & blocks = getParam<std::vector<SubdomainName>>("block");
  if (!blocks.empty())
    params.set<std::vector<SubdomainName>>("block") = blocks;
}

void
FissionProductInventoryAction::act()
{
  if (_current_task == "add_variable")
    addVariables();
  else if (_current_task == "add_ic")
    addInitialConditions();
  else if (_current_task == "add_interpolation_method_physics")
  {
    if (isParamValid("rhie_chow_user_object") || isParamValid("velocity"))
      OpenPronghorn::addAdvectedInterpolationMethod(parameters(), *_problem, _factory);
  }
  else if (_current_task == "add_linear_fv_kernel")
    addKernels();
  else if (_current_task == "add_linear_fv_bc")
    addBoundaryConditions();
}

void
FissionProductInventoryAction::addVariables()
{
  if (getParam<bool>("verbose"))
  {
    _console << "[FissionProductInventory] " << _plan.size()
             << " inventories, ordered parent before daughter (list the systems in this order "
                "in [Problem] linear_sys_names and the executioner):\n";
    for (const auto & p : _plan)
    {
      _console << "  " << std::setw(14) << std::left << p.name << " -> " << systemName(p.name);
      if (p.location == "salt")
        _console << "  yield=" << p.yield << " lambda=" << p.lambda << " 1/s sigma=" << p.sigma
                 << " b";
      _console << "\n";
    }
    _console << "  linear_sys_names = '";
    for (const auto i : index_range(_plan))
      _console << (i ? " " : "") << systemName(_plan[i].name);
    _console << "'\n" << std::flush;
  }

  _problem->needFV();
  for (const auto & p : _plan)
  {
    auto params = _factory.getValidParams("MooseLinearVariableFVReal");
    params.set<SolverSystemName>("solver_sys") = systemName(p.name);
    maybeAssignBlocks(params);
    _problem->addVariable("MooseLinearVariableFVReal", p.name, params);
  }
}

void
FissionProductInventoryAction::addInitialConditions()
{
  for (const auto & p : _plan)
  {
    auto params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = p.name;
    params.set<Real>("value") = p.initial_condition;
    maybeAssignBlocks(params);
    _problem->addInitialCondition("ConstantIC", "fp_" + p.name + "_ic", params);
  }
}

void
FissionProductInventoryAction::addKernels()
{
  const bool has_rhie_chow = isParamValid("rhie_chow_user_object");
  const bool has_velocity = isParamValid("velocity");
  const bool has_diffusion = isParamValid("diffusivity");
  const auto & flux = getParam<MooseFunctorName>("neutron_flux");

  // The fission source density and the factor converting it to [mol/(m^3 s)] per unit yield.
  const bool by_power = isParamValid("power_density");
  const MooseFunctorName source = by_power ? getParam<MooseFunctorName>("power_density")
                                           : getParam<MooseFunctorName>("fission_rate_density");
  const Real per_unit_yield =
      by_power ? 1.0 / (FissionProducts::avogadro * getParam<Real>("energy_per_fission") *
                        FissionProducts::MeV_to_J)
               : 1.0 / FissionProducts::avogadro;

  for (const auto & p : _plan)
  {
    const std::string tag = "fp_" + p.name;

    if (getParam<bool>("time_derivative"))
    {
      auto params = _factory.getValidParams("LinearFVTimeDerivative");
      params.set<LinearVariableName>("variable") = p.name;
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("LinearFVTimeDerivative", tag + "_time", params);
    }

    if (p.transported)
    {
      if (has_rhie_chow)
      {
        auto params = _factory.getValidParams("LinearFVScalarAdvection");
        params.set<LinearVariableName>("variable") = p.name;
        params.set<UserObjectName>("rhie_chow_user_object") =
            getParam<UserObjectName>("rhie_chow_user_object");
        params.set<InterpolationMethodName>("advected_interp_method_name") =
            OpenPronghorn::advectedInterpolationMethodName(parameters());
        maybeAssignBlocks(params);
        _problem->addLinearFVKernel("LinearFVScalarAdvection", tag + "_advection", params);
      }
      else if (has_velocity)
      {
        auto params = _factory.getValidParams("LinearFVAdvection");
        params.set<LinearVariableName>("variable") = p.name;
        params.set<RealVectorValue>("velocity") = getParam<RealVectorValue>("velocity");
        params.set<InterpolationMethodName>("advected_interp_method_name") =
            OpenPronghorn::advectedInterpolationMethodName(parameters());
        maybeAssignBlocks(params);
        _problem->addLinearFVKernel("LinearFVAdvection", tag + "_advection", params);
      }
      if (has_diffusion)
      {
        auto params = _factory.getValidParams("LinearFVDiffusion");
        params.set<LinearVariableName>("variable") = p.name;
        params.set<MooseFunctorName>("diffusion_coeff") = getParam<MooseFunctorName>("diffusivity");
        params.set<bool>("use_nonorthogonal_correction") = false;
        maybeAssignBlocks(params);
        _problem->addLinearFVKernel("LinearFVDiffusion", tag + "_diffusion", params);
      }
    }

    // Direct fission production
    if (p.location == "salt" && p.yield > 0.0)
    {
      auto params = _factory.getValidParams("LinearFVSource");
      params.set<LinearVariableName>("variable") = p.name;
      params.set<MooseFunctorName>("source_density") = source;
      params.set<MooseFunctorName>("scaling_factor") = fpNumberToFunctor(p.yield * per_unit_yield);
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("LinearFVSource", tag + "_fission", params);
    }

    // Production from tracked parents (decay and capture)
    for (const auto & f : p.feeds)
    {
      auto params = _factory.getValidParams("FissionProductFeed");
      params.set<LinearVariableName>("variable") = p.name;
      params.set<MooseFunctorName>("parent") = f.parent;
      params.set<Real>("decay_coefficient") = f.decay_coefficient;
      params.set<Real>("capture_cross_section") = f.capture_cross_section;
      params.set<MooseFunctorName>("neutron_flux") = flux;
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("FissionProductFeed", tag + "_from_" + f.parent, params);
    }

    // Transfer from the salt into a removed inventory
    if (!p.transfer_rate.empty())
    {
      auto params = _factory.getValidParams("FissionProductFeed");
      params.set<LinearVariableName>("variable") = p.name;
      params.set<MooseFunctorName>("parent") = p.nuclide;
      params.set<MooseFunctorName>("transfer_rate") = p.transfer_rate;
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("FissionProductFeed", tag + "_removal_from_salt", params);
    }

    // Loss by decay, capture and removal
    if (p.lambda > 0.0 || (p.location == "salt" && (p.sigma > 0.0 || !p.removal.empty())))
    {
      auto params = _factory.getValidParams("FissionProductLoss");
      params.set<LinearVariableName>("variable") = p.name;
      params.set<Real>("decay_constant") = p.lambda;
      if (p.location == "salt")
      {
        params.set<Real>("capture_cross_section") = p.sigma;
        params.set<MooseFunctorName>("neutron_flux") = flux;
        params.set<std::vector<MooseFunctorName>>("removal_rates") = p.removal;
      }
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("FissionProductLoss", tag + "_loss", params);
    }
  }
}

void
FissionProductInventoryAction::addBoundaryConditions()
{
  const auto & boundaries = getParam<std::vector<BoundaryName>>("plating_boundaries");
  if (boundaries.empty())
    return;
  for (const auto & p : _plan)
  {
    if (p.location != "salt" || _db.nuclide(p.nuclide).group != "noble_metal")
      continue;
    // Plating flux J = k_m N into the wall (a perfect sink): the degassing BC with no back
    // pressure.
    auto params = _factory.getValidParams("MSRDegassingBC");
    params.set<LinearVariableName>("variable") = p.name;
    params.set<std::vector<BoundaryName>>("boundary") = boundaries;
    params.set<MooseFunctorName>("mass_transfer_coefficient") =
        getParam<MooseFunctorName>("plating_mass_transfer_coefficient");
    _problem->addLinearFVBC("MSRDegassingBC", "fp_" + p.name + "_wall_plating", params);
  }
}
