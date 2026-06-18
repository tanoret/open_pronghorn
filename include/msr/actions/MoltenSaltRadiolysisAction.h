//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"
#include "MoltenSaltRadiolysisData.h"

#include <map>
#include <string>
#include <vector>

/**
 * Action that assembles a molten salt radiolysis problem from the built-in chemistry database.
 *
 * Given a salt kernel (chloride or fluoride), an optional set of metals (Zn, Cr, U) and optional
 * gas-phase diatomics (Cl2, F2), the action creates one linear finite-volume passive-scalar
 * variable per species and adds the kernels that integrate the radiolysis problem:
 *
 *  - LinearFVTimeDerivative for every species,
 *  - LinearFVScalarAdvection (when a Rhie-Chow user object is supplied) and LinearFVDiffusion for
 *    the dissolved species (flow coupling and molecular diffusion),
 *  - LinearFVSource for the radiolytic source S = G * dose_rate / (100 eV * N_A),
 *  - MSRReaction for each species' contribution to every mass-action reaction it participates in,
 *  - MSRGasExchange for each gas-phase diatomic (dissolved <-> headspace).
 *
 * Each species is placed in its own linear solver system named "<species>_sys"; these system names
 * must be listed in [Problem] linear_sys_names and in the executioner so that the segregated
 * fixed-point solve iterates the (lagged) reaction coupling. The action prints the generated system
 * names when 'verbose = true'.
 */
class MoltenSaltRadiolysisAction : public Action
{
public:
  static InputParameters validParams();

  MoltenSaltRadiolysisAction(const InputParameters & params);

  virtual void act() override;

protected:
  /// Description of one species to be created
  struct SpeciesPlan
  {
    std::string name;
    bool is_gas;
    Real initial_condition;
    Real g_value; // radiolytic yield G [molecules/100eV] (zero if none)
  };

  /// A dissolved/headspace diatomic pair with its Henry coefficient
  struct GasPair
  {
    std::string liquid;
    std::string vapor;
    Real kH;
  };

  /// Build the species list, reaction network and source terms from the database and parameters
  void buildPlan();

  /// Create one passive-scalar variable per species
  void addVariables();

  /// Create the constant initial conditions for the species
  void addInitialConditions();

  /// Add the transport, source, reaction and gas-exchange kernels
  void addKernels();

  /// Add buoyant advection and dispersion for a gas-phase species
  void addGasTransport(const std::string & name,
                       bool has_rhie_chow,
                       bool has_velocity,
                       bool has_buoyancy,
                       bool has_gas_dispersion);

  /// Add the MSRReaction kernels for a single reaction
  void addReactionKernels(const MSR::ReactionData & reaction, unsigned int index);

  /// Convenience: linear solver system name for a species
  std::string systemName(const std::string & species) const { return species + "_sys"; }

  /// Set the block restriction on a set of child parameters, when one was provided
  void maybeAssignBlocks(InputParameters & params) const;

  /// Salt kernel ("chloride" or "fluoride")
  const std::string _salt;
  /// Temperature functor name [K]
  const MooseFunctorName _temperature;

  /// The species to create
  std::vector<SpeciesPlan> _species;
  /// Fast lookup of whether a species is tracked (and gas phase)
  std::map<std::string, bool> _is_gas;
  /// The full reaction network (core plus metal templates)
  std::vector<MSR::ReactionData> _reactions;
  /// The diatomic gas pairs with an active headspace
  std::vector<GasPair> _gas_pairs;
};
