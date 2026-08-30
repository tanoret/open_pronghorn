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
#include "CorrosionPlatingActionUtils.h"
#include "FEProblemBase.h"
#include "Factory.h"

#include "MoltenSaltCorrosionModel.h"
#include "CorrosionChemistry.h"
#include "AdvancedCorrosionModelData.h"
#include "MSTDBTCData.h"
#include "MSTDBTCStandardStateCorrosionModel.h"

#include <array>
#include <cctype>
#include <cmath>
#include <set>

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
  params.addClassDescription(
      "Sets up salt-side molten salt corrosion (linear finite-volume passive scalars plus a "
      "Butler-Volmer wall reaction) for a flowing MSR.");

  params.addParam<DataFileName>("database",
                                "corrosion_database.json",
                                "The JSON corrosion database. Defaults to "
                                "data/corrosion_database.json.");

  MooseEnum kinetics_model("reduced_empirical mstdb_tc_standard_state", "reduced_empirical");
  params.addParam<MooseEnum>(
      "kinetics_model",
      kinetics_model,
      "Static model used to seed the wall Butler-Volmer exchange currents. In MSTDB-TC mode the "
      "total dissolution-front rate is apportioned among Cr, Fe, and Ni.");
  params.addParam<DataFileName>(
      "advanced_database",
      "Advanced-model closure parameters and mandatory MSTDB-TC edition/SHA-256 binding. "
      "Required when kinetics_model=mstdb_tc_standard_state.");
  params.addParam<FileName>(
      "fluoride_database",
      "Authorized external MSTDB-TC V3.1 fluoride *_No_Func.dat file. Required in MSTDB mode.");
  params.addParam<FileName>(
      "chloride_database",
      "Authorized external MSTDB-TC V3.1 chloride *_No_Func.dat file. Required in MSTDB mode.");
  params.addParam<bool>(
      "allow_extrapolation",
      false,
      "Permit Gibbs-function evaluation beyond the last MSTDB-TC interval. This never disables "
      "the version or SHA-256 checks.");

  MultiMooseEnum elements("Cr Fe Ni Mo Nb Tc Ru Ag Sb Te", "Cr Fe Ni");
  params.addParam<MultiMooseEnum>("elements", elements, "The metal elements to track.");
  params.addParam<std::vector<std::string>>(
      "release_variables",
      {},
      "One entry per tracked element: the existing variable that element's corrosion is released "
      "into (e.g. the radiolysis-tracked 'Cr_II'), or 'none' for an element that gets its own "
      "c_<El> variable created and transported here. Leave empty to create c_<El> for every "
      "element. Releasing into a radiolysis cation couples corrosion and radiolysis; that "
      "variable's transport and chemistry are then owned by the radiolysis action. Example for "
      "'Cr Fe Ni': "
      "'Cr_II none none' feeds chromium into the radiolysis Cr(II) and tracks iron and nickel as "
      "their own c_Fe / c_Ni.");

  params.addRequiredParam<MooseFunctorName>("temperature", "Temperature functor [K].");
  params.addParam<Real>("reference_temperature",
                        923.15,
                        "Scalar temperature [K] for the kinetics seed and the Arrhenius anchor.");
  params.addParam<Real>(
      "reference_cold_temperature",
      "Reference cold-side salt temperature [K]. Required in MSTDB mode and must not exceed "
      "reference_temperature.");
  params.addParam<Real>(
      "reference_exposure_time",
      "Reference physical exposure time [s]. Required and nonnegative in MSTDB mode.");
  params.addParam<Real>(
      "area_to_salt_mass",
      "Explicit wetted-area to salt-mass ratio [cm^2/g]. Required and positive in MSTDB mode.");
  params.addParam<Real>(
      "inventory_coupling_factor",
      "Explicit fraction applied only to dissolved inventory in the modeled salt [0,1]. It does "
      "not scale static cold-capture mass gain. Required in MSTDB mode.");
  params.addParam<bool>("temperature_dependent_kinetics",
                        true,
                        "Make the exchange current Arrhenius in the local temperature field. In "
                        "MSTDB mode this is an explicitly reduced post-seed scaling using the base "
                        "corrosion database activation energy; MSTDB is not reevaluated locally.");

  params.addParam<std::string>("material_class", "hastelloy_n", "Alloy class.");
  params.addParam<std::string>("salt_class", "fluoride_fuel", "Salt class.");
  params.addParam<std::string>("redox_class", "oxidizing_fef2", "Redox/overpotential class.");
  params.addParam<std::string>("position_class", "nominal", "Loop position class.");
  params.addParam<Real>("flow_factor", 1.0, "Circulation/mass-transfer factor.");
  params.addParam<Real>("delta_T_C", 0.0, "Loop thermal gradient [C].");
  params.addParam<Real>(
      "applied_overpotential",
      0.1,
      "Imposed metal-minus-salt potential difference [V]. The element E0 is subtracted from this "
      "value both when seeding and when evaluating the Butler-Volmer wall law.");
  params.addParam<Real>(
      "reference_concentration",
      1.0,
      "Explicit runtime reference concentration c_ref [mol/m^3], replacing the provenance-bound "
      "base database default without disabling its provenance check.");

  params.addRequiredParam<std::vector<BoundaryName>>("reaction_boundary",
                                                     "The corroding wall boundary(ies).");

  params.addParam<UserObjectName>(
      "rhie_chow_user_object", "The Rhie-Chow user object for advection by the segregated flow.");
  params.addParam<RealVectorValue>("velocity", "A prescribed constant advection velocity.");
  MooseEnum interp("average upwind", "upwind");
  params.addParam<MooseEnum>("advected_interp_method", interp, "Advected interpolation method.");
  params.addParam<MooseFunctorName>(
      "diffusivity",
      "Explicit effective (molecular + turbulent) diffusivity functor [m^2/s], replacing the "
      "provenance-bound base database default without disabling its provenance check.");

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
    _use_mstdb_tc(getParam<MooseEnum>("kinetics_model") == "mstdb_tc_standard_state"),
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

  Real rate_um_y = 0.0;
  std::array<Real, 3> source_fractions{{1.0, 1.0, 1.0}};
  std::array<Real, 3> affinities{{0.0, 0.0, 0.0}};
  std::string calibration_id;
  std::string mstdb_version;
  std::string fluoride_sha256;
  std::string chloride_sha256;

  const auto & elements = getParam<MultiMooseEnum>("elements");
  if (!_use_mstdb_tc)
    // Keep the legacy path and its one-total-rate-per-selected-element seeding unchanged.
    rate_um_y = model.corrosionRateUmY(features);
  else
  {
    const std::set<std::string> required_elements{"Cr", "Fe", "Ni"};
    std::set<std::string> selected_elements;
    unsigned int selected_count = 0;
    for (const auto & element : elements)
    {
      selected_elements.insert(canonicalElementToken(element.name()));
      ++selected_count;
    }
    if (selected_count != required_elements.size() || selected_elements != required_elements)
      paramError("elements",
                 "kinetics_model=mstdb_tc_standard_state requires exactly 'Cr Fe Ni'. The "
                 "calibrated source fractions are normalized over those three elements.");

    for (const auto * parameter : {"advanced_database",
                                   "fluoride_database",
                                   "chloride_database",
                                   "reference_cold_temperature",
                                   "reference_exposure_time",
                                   "area_to_salt_mass",
                                   "inventory_coupling_factor"})
      if (!isParamValid(parameter))
        paramError(parameter,
                   "This parameter is required when "
                   "kinetics_model=mstdb_tc_standard_state.");

    const Real cold_temperature = getParam<Real>("reference_cold_temperature");
    const Real exposure_time = getParam<Real>("reference_exposure_time");
    const Real area_to_salt_mass = getParam<Real>("area_to_salt_mass");
    const Real inventory_coupling = getParam<Real>("inventory_coupling_factor");
    if (!std::isfinite(_reference_temperature) || _reference_temperature <= 0.0)
      paramError("reference_temperature", "Must be finite and positive in MSTDB mode.");
    if (!std::isfinite(cold_temperature) || cold_temperature <= 0.0)
      paramError("reference_cold_temperature", "Must be finite and positive.");
    if (cold_temperature > _reference_temperature)
      paramError("reference_cold_temperature", "Must not exceed reference_temperature.");
    if (!std::isfinite(exposure_time) || exposure_time < 0.0)
      paramError("reference_exposure_time", "Must be finite and nonnegative.");
    if (!std::isfinite(area_to_salt_mass) || area_to_salt_mass <= 0.0)
      paramError("area_to_salt_mass", "Must be finite and positive.");
    if (!std::isfinite(inventory_coupling) || inventory_coupling < 0.0 || inventory_coupling > 1.0)
      paramError("inventory_coupling_factor", "Must be finite and in [0,1].");

    Corrosion::AdvancedCorrosionModelDatabase advanced(getParam<DataFileName>("advanced_database"));
    advanced.validateBaseModel(db);
    if (advanced.expectedMSTDBVersion() != "3.1")
      paramError("advanced_database",
                 "This action release requires an MSTDB-TC V3.1 calibration binding; found '",
                 advanced.expectedMSTDBVersion(),
                 "'. Recalibration and revalidation are required before changing editions.");
    Corrosion::MSTDBTCPair thermodynamics(getParam<FileName>("fluoride_database"),
                                          getParam<FileName>("chloride_database"),
                                          advanced.expectedMSTDBVersion(),
                                          advanced.expectedFluorideSHA256(),
                                          advanced.expectedChlorideSHA256(),
                                          false,
                                          getParam<bool>("allow_extrapolation"));
    Corrosion::MSTDBTCStandardStateCorrosionModel mstdb_model(db, advanced, thermodynamics);
    Corrosion::MSTDBTCCorrosionFeatures mstdb_features;
    mstdb_features.hot_temperature_K = _reference_temperature;
    mstdb_features.cold_temperature_K = cold_temperature;
    mstdb_features.exposure_s = exposure_time;
    mstdb_features.flow_factor = features.flow_factor;
    mstdb_features.area_to_salt_mass_cm2_g = area_to_salt_mass;
    mstdb_features.inventory_coupling_factor = inventory_coupling;
    mstdb_features.material_class = features.material_class;
    mstdb_features.salt_class = features.salt_class;
    mstdb_features.redox_class = features.redox_class;

    const auto endpoint = mstdb_model.evaluate(mstdb_features);
    if (!std::isfinite(endpoint.front_rate_um_y) || endpoint.front_rate_um_y < 0.0)
      paramError("advanced_database",
                 "The MSTDB-TC endpoint produced a nonfinite or negative front rate. Check the "
                 "calibrated parameters and thermodynamic input range.");
    Real fraction_sum = 0.0;
    for (const auto fraction : endpoint.source_fraction)
    {
      if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0)
        paramError("advanced_database",
                   "The MSTDB-TC endpoint produced a source fraction outside [0,1]. Check the "
                   "calibrated parameters and thermodynamic input range.");
      fraction_sum += fraction;
    }
    if (std::abs(fraction_sum - 1.0) > 1.0e-10)
      paramError("advanced_database",
                 "The MSTDB-TC Cr/Fe/Ni source fractions sum to ",
                 fraction_sum,
                 " instead of one.");
    rate_um_y = endpoint.front_rate_um_y;
    source_fractions = endpoint.source_fraction;
    affinities = endpoint.affinity_log_k_over_q;
    calibration_id = advanced.calibrationId();
    mstdb_version = thermodynamics.version();
    fluoride_sha256 = thermodynamics.fluoride().sha256();
    chloride_sha256 = thermodynamics.chloride().sha256();
  }

  const Real density = db.density(features.material_class);
  const Real eta = getParam<Real>("applied_overpotential");
  const Real c_ref = getParam<Real>("reference_concentration");
  const Real initial = getParam<Real>("initial_concentration");

  const auto & release = getParam<std::vector<std::string>>("release_variables");
  if (!release.empty() && release.size() != elements.size())
    paramError("release_variables",
               "Must list one entry per tracked element: the existing variable each element's "
               "corrosion is released into, or 'none' for an element that gets its own c_<El> "
               "variable. Leave the whole parameter empty to create c_<El> for every element.");

  unsigned int index = 0;
  for (const auto & el : elements)
  {
    const Corrosion::ElementProperties props = db.element(el.name());
    ElementPlan plan;
    plan.name = canonicalElementToken(el.name());
    // An element owns its own variable when no release target is given (empty list, empty entry, or
    // the explicit sentinel 'none'); otherwise it is released into the named existing variable.
    plan.owns_variable = release.empty() || release[index].empty() || release[index] == "none";
    plan.salt_var = plan.owns_variable ? "c_" + plan.name : release[index];
    plan.valence = props.valence;
    plan.molar_mass = props.molar_mass_g_mol;
    plan.alpha_a = props.alpha_a;
    plan.alpha_c = props.alpha_c;
    plan.E0 = props.E0_V;
    plan.c_ref = c_ref;
    plan.initial_salt = initial;

    if (!_use_mstdb_tc)
      plan.planned_rate_um_y = rate_um_y;
    else if (plan.name == "Cr")
      plan.planned_rate_um_y = rate_um_y * source_fractions[0];
    else if (plan.name == "Fe")
      plan.planned_rate_um_y = rate_um_y * source_fractions[1];
    else if (plan.name == "Ni")
      plan.planned_rate_um_y = rate_um_y * source_fractions[2];
    else
      paramError("elements", "Internal error: unsupported MSTDB-TC element '", plan.name, "'.");

    // The input is the metal-minus-salt potential difference.  Match the boundary object's E0
    // subtraction and clipped exponents so custom nonzero-E0 data reproduce the reference current.
    const Real i_ref_a_m2 =
        Corrosion::umYToCorrosionCurrent(plan.planned_rate_um_y,
                                         plan.valence,
                                         plan.molar_mass,
                                         density) *
        1.0e4;
    const Real bracket = Corrosion::ActionKinetics::butlerVolmerSeedBracket(
        eta, plan.E0, plan.valence, plan.alpha_a, plan.alpha_c, _reference_temperature);
    if (i_ref_a_m2 > 0.0 && (!std::isfinite(bracket) || bracket <= 1.0e-30))
      paramError("applied_overpotential",
                 "The metal-minus-salt potential difference must produce a positive anodic "
                 "Butler-Volmer bracket for ",
                 plan.name,
                 " (applied_overpotential=",
                 eta,
                 " V, E0=",
                 plan.E0,
                 " V). A positive dissolution-rate seed cannot be represented otherwise.");
    plan.i0 = i_ref_a_m2 == 0.0 ? 0.0 : i_ref_a_m2 / bracket;

    _elements.push_back(plan);
    ++index;
  }

  if (getParam<bool>("verbose"))
  {
    if (!_use_mstdb_tc)
    {
      _console << "[CorrosionPlatingFlow] reference corrosion rate " << rate_um_y << " um/y; "
               << _elements.size() << " elements\n";
      for (const auto & e : _elements)
        _console << "  " << e.name << " -> " << e.salt_var << "  i0=" << e.i0 << " A/m^2\n";
    }
    else
    {
      _console << "[CorrosionPlatingFlow] kinetics_model=mstdb_tc_standard_state"
               << " calibration=" << calibration_id << " MSTDB-TC=" << mstdb_version << '\n'
               << "  fluoride_sha256=" << fluoride_sha256 << '\n'
               << "  chloride_sha256=" << chloride_sha256 << '\n'
               << "  total front rate=" << rate_um_y << " um/y; hot="
               << _reference_temperature << " K, cold="
               << getParam<Real>("reference_cold_temperature") << " K, exposure="
               << getParam<Real>("reference_exposure_time") << " s, area/salt="
               << getParam<Real>("area_to_salt_mass") << " cm^2/g, inventory coupling="
               << getParam<Real>("inventory_coupling_factor") << '\n';
      for (const auto & element : _elements)
      {
        const unsigned int i = element.name == "Cr" ? 0 : (element.name == "Fe" ? 1 : 2);
        _console << "  " << element.name << " -> " << element.salt_var
                 << " source_fraction=" << source_fractions[i]
                 << " affinity_log_K_over_Q=" << affinities[i]
                 << " planned_rate=" << element.planned_rate_um_y << " um/y"
                 << " i0=" << element.i0 << " A/m^2\n";
      }
    }
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
      params.set<MooseEnum>("advected_interp_method") =
          getParam<MooseEnum>("advected_interp_method");
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel(
          "LinearFVScalarAdvection", "corr_" + e.salt_var + "_advection", params);
    }
    else if (has_velocity)
    {
      auto params = _factory.getValidParams("LinearFVAdvection");
      params.set<LinearVariableName>("variable") = e.salt_var;
      params.set<RealVectorValue>("velocity") = getParam<RealVectorValue>("velocity");
      params.set<MooseEnum>("advected_interp_method") =
          getParam<MooseEnum>("advected_interp_method");
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
      params.set<MooseFunctorName>("metal_potential") =
          Corrosion::ActionKinetics::realFunctorName(eta);
      _problem->addLinearFVBC(
          "CorrosionLinearFVButlerVolmerBC", "corr_" + e.salt_var + "_wall", params);
    }

    if (!inlets.empty() && e.owns_variable)
    {
      auto params = _factory.getValidParams("LinearFVAdvectionDiffusionFunctorDirichletBC");
      params.set<LinearVariableName>("variable") = e.salt_var;
      params.set<std::vector<BoundaryName>>("boundary") = inlets;
      params.set<MooseFunctorName>("functor") = Corrosion::ActionKinetics::realFunctorName(inlet_c);
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
