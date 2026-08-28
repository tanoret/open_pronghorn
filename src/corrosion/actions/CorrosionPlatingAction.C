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
#include <map>
#include <set>

registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_variable");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_aux_variable");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_ic");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_kernel");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_interface_kernel");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_bc");
registerMooseAction("OpenPronghornApp", CorrosionPlatingAction, "add_aux_kernel");

namespace
{
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
  params.addClassDescription(
      "Creates the variables, transport kernels, potential equation and Butler-Volmer electrode "
      "kinetics of a molten salt corrosion and plating problem using either the reduced empirical "
      "correlation or a provenance-bound MSTDB-TC standard-state seed.");

  params.addParam<DataFileName>("database",
                                "corrosion_database.json",
                                "The JSON corrosion database (elements, material tables, calibrated "
                                "correlation parameters). Defaults to data/corrosion_database.json.");

  MooseEnum kinetics_model("reduced_empirical mstdb_tc_standard_state", "reduced_empirical");
  params.addParam<MooseEnum>(
      "kinetics_model",
      kinetics_model,
      "Static model used to seed the Butler-Volmer exchange currents. The default preserves the "
      "validated reduced-empirical correlation. 'mstdb_tc_standard_state' evaluates the calibrated "
      "MSTDB-TC standard-state/Nernst engineering model at the reference endpoint and apportions "
      "its total dissolution-front rate among Cr, Fe, and Ni.");
  params.addParam<DataFileName>(
      "advanced_database",
      "Advanced-model closure parameters and the mandatory MSTDB-TC edition/SHA-256 binding. "
      "Required when kinetics_model=mstdb_tc_standard_state.");
  params.addParam<FileName>(
      "fluoride_database",
      "Authorized external MSTDB-TC V3.1 fluoride *_No_Func.dat file. Required in MSTDB mode; "
      "the file hash must match advanced_database.");
  params.addParam<FileName>(
      "chloride_database",
      "Authorized external MSTDB-TC V3.1 chloride *_No_Func.dat file. Required in MSTDB mode; "
      "the file hash must match advanced_database.");
  params.addParam<bool>(
      "allow_extrapolation",
      false,
      "Permit Gibbs-function evaluation beyond the last MSTDB-TC interval. This never disables "
      "the version or SHA-256 checks.");

  MultiMooseEnum elements("Cr Fe Ni Mo Nb Tc Ru Ag Sb Te", "Cr Fe Ni");
  params.addParam<MultiMooseEnum>("elements", elements, "The metal elements to track.");
  MooseEnum recession_element("Cr Fe Ni Mo Nb Tc Ru Ag Sb Te", "Cr");
  params.addParam<MooseEnum>(
      "recession_element",
      recession_element,
      "The element controlling the legacy reduced-model recession/plating output. MSTDB-TC mode "
      "instead reports the summed Cr/Fe/Ni front rate and recession.");

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
                        "Finite imposed metal-minus-salt potential difference [V]. The element E0 is "
                        "subtracted from this value when seeding the exchange current, exactly as "
                        "in the generated Butler-Volmer objects; a positive dissolution seed "
                        "therefore requires a positive anodic bracket.");
  params.addParam<MooseFunctorName>(
      "applied_potential",
      "Absolute potential functor [V] of the absent counter phase in a single-domain topology. "
      "When omitted, salt_only uses applied_overpotential plus the pinned salt-potential value "
      "when that potential is solved (or just applied_overpotential when it is not solved); "
      "solid_only uses -applied_overpotential for the external salt. A supplied functor is an "
      "explicit absolute potential and is never shifted automatically.");
  params.addParam<MooseFunctorName>(
      "counter_concentration",
      "Salt-ion concentration functor [mol/m^3] of the absent salt phase. Required for "
      "solid_only and used identically by the Butler-Volmer boundary and rate diagnostics.");

  // Electric potential (electrophoresis) solve.
  params.addParam<bool>("solve_potential",
                        false,
                        "Solve the salt electric potential (current continuity) and include "
                        "electromigration. Requires 'pin_potential_boundary' and a modeled salt "
                        "phase, so it is unavailable for solid_only.");
  params.addParam<bool>("supporting_electrolyte",
                        true,
                        "Treat the salt as a supporting electrolyte: the Butler-Volmer current is "
                        "not fed back into the potential equation (the dilute tracers only migrate "
                        "in the field). Set false to couple the interfacial current to the salt "
                        "charge balance.");
  params.addParam<BoundaryName>("pin_potential_boundary",
                                "Boundary at which the salt potential is pinned (fixes the "
                                "otherwise-undetermined constant mode).");
  params.addParam<Real>(
      "pin_potential_value",
      0.0,
      "Finite pinned salt potential value [V]. The generated default metal potential carries this "
      "same absolute offset so metal minus salt remains applied_overpotential; their sum must also "
      "be representable as a finite Real.");

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
    _use_mstdb_tc(getParam<MooseEnum>("kinetics_model") == "mstdb_tc_standard_state"),
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

  if (_topology == Topology::SolidOnly && !isParamValid("counter_concentration"))
    paramError("counter_concentration",
               "The solid_only topology requires the absent salt-ion concentration as a functor.");

  if (_topology == Topology::SolidOnly && _solve_potential)
    paramError("solve_potential",
               "The solid_only topology has no modeled salt phase or salt-potential variable. "
               "Set solve_potential=false or use salt_only/two_block.");

  if (_solve_potential && !isParamValid("pin_potential_boundary"))
    mooseError("CorrosionPlating: 'solve_potential' requires 'pin_potential_boundary' to fix the "
               "constant mode of the potential.");

  const Real applied_overpotential = getParam<Real>("applied_overpotential");
  if (!std::isfinite(applied_overpotential))
    paramError("applied_overpotential", "Must be finite.");
  if (_solve_potential)
  {
    const Real salt_pin = getParam<Real>("pin_potential_value");
    if (!std::isfinite(salt_pin))
      paramError("pin_potential_value",
                 "Must be finite when solve_potential=true because it sets the absolute gauge of "
                 "every generated salt and metal potential.");
    if (!std::isfinite(Corrosion::ActionKinetics::defaultMetalPotential(
            applied_overpotential, true, salt_pin)))
      paramError("pin_potential_value",
                 "Its sum with applied_overpotential must be finite so the generated default "
                 "metal potential is representable.");
  }

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

  Real rate_um_y = 0.0;
  std::array<Real, 3> source_fractions{{1.0, 1.0, 1.0}};
  std::array<Real, 3> affinities{{0.0, 0.0, 0.0}};
  std::string calibration_id;
  std::string mstdb_version;
  std::string fluoride_sha256;
  std::string chloride_sha256;

  if (!_use_mstdb_tc)
    // Keep the legacy path and its one-total-rate-per-selected-element seeding unchanged.
    rate_um_y = model.corrosionRateUmY(features);
  else
  {
    const std::set<std::string> required_elements{"Cr", "Fe", "Ni"};
    std::set<std::string> selected_elements;
    unsigned int selected_count = 0;
    for (const auto & element : getParam<MultiMooseEnum>("elements"))
    {
      selected_elements.insert(canonicalElement(element.name()));
      ++selected_count;
    }
    if (selected_count != required_elements.size() || selected_elements != required_elements)
      paramError("elements",
                 "kinetics_model=mstdb_tc_standard_state requires exactly 'Cr Fe Ni'. The "
                 "calibrated source fractions are normalized over those three elements, so "
                 "omitting one or adding another would discard or duplicate front recession.");

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
    if (!std::isfinite(inventory_coupling) || inventory_coupling < 0.0 ||
        inventory_coupling > 1.0)
      paramError("inventory_coupling_factor", "Must be finite and in [0,1].");

    Corrosion::AdvancedCorrosionModelDatabase advanced(
        getParam<DataFileName>("advanced_database"));
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

    // The input is the metal-minus-salt potential difference.  The boundary objects subtract E0,
    // so use the same effective overpotential and exponent clipping here; otherwise a custom
    // nonzero-E0 database would not reproduce the requested reference current.
    const Real i_ref_a_cm2 =
        Corrosion::umYToCorrosionCurrent(plan.planned_rate_um_y,
                                         plan.valence,
                                         plan.molar_mass,
                                         density);
    const Real i_ref_a_m2 = i_ref_a_cm2 * 1.0e4;
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

    if (plan.name == recession_name)
      _recession_index = _elements.size();
    _elements.push_back(plan);
  }

  if (getParam<bool>("verbose"))
  {
    if (!_use_mstdb_tc)
    {
      _console << "[CorrosionPlating] reference corrosion rate " << rate_um_y << " um/y; "
               << _elements.size() << " elements:\n";
      for (const auto & e : _elements)
        _console << "  " << e.name << " z=" << e.valence << " D=" << e.diffusivity
                 << " i0=" << e.i0 << " A/m^2\n";
    }
    else
    {
      _console << "[CorrosionPlating] kinetics_model=mstdb_tc_standard_state"
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
        _console << "  " << element.name << " source_fraction=" << source_fractions[i]
                 << " affinity_log_K_over_Q=" << affinities[i]
                 << " planned_rate=" << element.planned_rate_um_y << " um/y"
                 << " z=" << element.valence << " D=" << element.diffusivity
                 << " i0=" << element.i0 << " A/m^2\n";
      }
    }
    if (_solve_potential)
    {
      const Real salt_pin = getParam<Real>("pin_potential_value");
      _console << "  salt_potential_pin=" << salt_pin << " V; ";
      if (_topology == Topology::SaltOnly && isParamValid("applied_potential"))
        _console << "external_metal_potential=custom_absolute (not shifted by the Action); ";
      else
        _console << "default_metal_potential="
                 << Corrosion::ActionKinetics::defaultMetalPotential(eta, true, salt_pin) << " V; ";
      _console << "seed_metal_minus_salt=" << eta << " V\n";
    }
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
  const Real default_metal_potential = Corrosion::ActionKinetics::defaultMetalPotential(
      eta, _solve_potential, getParam<Real>("pin_potential_value"));
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
    // The metal carries the pinned salt's absolute offset so their reference difference is eta.
    params.set<Real>("metal_potential") = default_metal_potential;
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
  const bool metal_domain = (_topology == Topology::SolidOnly);
  const Real default_metal_potential = Corrosion::ActionKinetics::defaultMetalPotential(
      eta, _solve_potential, getParam<Real>("pin_potential_value"));
  // ButlerVolmerBC interprets applied_potential as the absent counter-phase absolute potential.
  // salt_only carries the modeled salt pin into the external metal; solid_only retains salt=-eta
  // against the modeled metal's zero fallback. Both defaults therefore give metal-salt=eta.
  const Real default_counter_potential = metal_domain ? -eta : default_metal_potential;
  const MooseFunctorName applied =
      isParamValid("applied_potential")
          ? getParam<MooseFunctorName>("applied_potential")
          : MooseFunctorName(Corrosion::ActionKinetics::realFunctorName(default_counter_potential));

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
        params.set<MooseFunctorName>("applied_potential") =
            Corrosion::ActionKinetics::realFunctorName(default_metal_potential);
        _problem->addBoundaryCondition("ButlerVolmerBC", "bv_charge_" + e.name, params);
      }
    }
    return;
  }

  // Single-domain topologies: Butler-Volmer species flux on the external electrode boundary.
  const BoundaryName boundary = getParam<BoundaryName>("reaction_boundary");
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
  // Recession / plating-thickness output on the electrode boundary. The reduced model retains its
  // legacy controlling-element diagnostic; MSTDB mode reports the sum of all source-resolved rates.
  const ElementPlan & e = _elements[_recession_index];
  Corrosion::MoltenSaltCorrosionDatabase db(getParam<DataFileName>("database"));
  const Real density = db.density(getParam<std::string>("material_class"));
  const Real eta = getParam<Real>("applied_overpotential");

  const bool metal_domain = (_topology == Topology::SolidOnly);
  const Real default_metal_potential = Corrosion::ActionKinetics::defaultMetalPotential(
      eta, _solve_potential, getParam<Real>("pin_potential_value"));
  const Real default_counter_potential = metal_domain ? -eta : default_metal_potential;
  const MooseFunctorName applied =
      isParamValid("applied_potential")
          ? getParam<MooseFunctorName>("applied_potential")
          : MooseFunctorName(Corrosion::ActionKinetics::realFunctorName(default_counter_potential));
  const BoundaryName boundary = (_topology == Topology::TwoBlock)
                                    ? getParam<BoundaryName>("interface_boundary")
                                    : getParam<BoundaryName>("reaction_boundary");

  if (_use_mstdb_tc)
  {
    std::vector<VariableName> concentrations;
    std::vector<Real> valences;
    std::vector<Real> molar_masses;
    std::vector<Real> exchange_current_densities;
    std::vector<Real> E0_values;
    std::vector<Real> alpha_a_values;
    std::vector<Real> alpha_c_values;
    std::vector<Real> c_ref_values;
    for (const auto & element : _elements)
    {
      concentrations.push_back(element.salt_var);
      valences.push_back(element.valence);
      molar_masses.push_back(element.molar_mass);
      exchange_current_densities.push_back(element.i0);
      E0_values.push_back(element.E0);
      alpha_a_values.push_back(element.alpha_a);
      alpha_c_values.push_back(element.alpha_c);
      c_ref_values.push_back(element.c_ref);
    }

    auto make_total_aux = [&](const std::string & var,
                              const std::string & mode,
                              const std::string & name,
                              const std::string & execute_on)
    {
      auto params = _factory.getValidParams("CorrosionTotalRateAux");
      params.set<AuxVariableName>("variable") = var;
      params.set<std::vector<BoundaryName>>("boundary") = {boundary};
      params.set<MooseEnum>("mode") = mode;
      params.set<ExecFlagEnum>("execute_on") = execute_on;
      if (metal_domain)
        params.set<MooseFunctorName>("concentration_functor") =
            getParam<MooseFunctorName>("counter_concentration");
      else
        params.set<std::vector<VariableName>>("concentrations") = concentrations;
      params.set<std::vector<Real>>("valences") = valences;
      params.set<std::vector<Real>>("molar_masses") = molar_masses;
      params.set<std::vector<Real>>("exchange_current_densities") =
          exchange_current_densities;
      params.set<std::vector<Real>>("E0_values") = E0_values;
      params.set<std::vector<Real>>("alpha_a_values") = alpha_a_values;
      params.set<std::vector<Real>>("alpha_c_values") = alpha_c_values;
      params.set<std::vector<Real>>("c_ref_values") = c_ref_values;
      params.set<Real>("density") = density;
      params.set<Real>("temperature") = _reference_temperature;
      if (_topology == Topology::TwoBlock)
        params.set<Real>("metal_potential_value") = default_metal_potential;
      else if (metal_domain)
      {
        params.set<MooseFunctorName>("salt_potential_functor") = applied;
        params.set<Real>("metal_potential_value") = 0.0;
      }
      else
        params.set<MooseFunctorName>("metal_potential_functor") = applied;
      if (_solve_potential && modelSalt())
        params.set<std::vector<VariableName>>("salt_potential") = {saltPotential()};
      else if (!metal_domain)
        params.set<Real>("salt_potential_value") = 0.0;
      _problem->addAuxKernel("CorrosionTotalRateAux", name, params);
    };

    make_total_aux("recession_um", "recession", "recession_aux", "TIMESTEP_END");
    make_total_aux("corrosion_rate_um_y",
                   "penetration_rate",
                   "corrosion_rate_aux",
                   "INITIAL TIMESTEP_END");
    return;
  }

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
    if (metal_domain)
      params.set<MooseFunctorName>("concentration_functor") =
          getParam<MooseFunctorName>("counter_concentration");
    else
      params.set<std::vector<VariableName>>("concentration") = {e.salt_var};
    params.set<Real>("valence") = e.valence;
    params.set<Real>("molar_mass") = e.molar_mass;
    params.set<Real>("density") = density;
    params.set<Real>("exchange_current_density") = e.i0;
    params.set<Real>("E0") = e.E0;
    params.set<Real>("alpha_a") = e.alpha_a;
    params.set<Real>("alpha_c") = e.alpha_c;
    params.set<Real>("temperature") = _reference_temperature;
    params.set<Real>("c_ref") = e.c_ref;
    if (_topology == Topology::TwoBlock)
      params.set<Real>("metal_potential_value") = default_metal_potential;
    else if (metal_domain)
    {
      params.set<MooseFunctorName>("salt_potential_functor") = applied;
      params.set<Real>("metal_potential_value") = 0.0;
    }
    else
      params.set<MooseFunctorName>("metal_potential_functor") = applied;
    if (_solve_potential && modelSalt())
      params.set<std::vector<VariableName>>("salt_potential") = {saltPotential()};
    else if (!metal_domain)
      params.set<Real>("salt_potential_value") = 0.0;
    _problem->addAuxKernel("CorrosionRateAux", name, params);
  };

  // The recession integrates over time (steps only); the rate is a snapshot (also at the initial
  // state, where it equals the calibrated reference rate).
  make_aux("recession_um", "recession", "recession_aux", "TIMESTEP_END");
  make_aux("corrosion_rate_um_y", "penetration_rate", "corrosion_rate_aux", "INITIAL TIMESTEP_END");
}
