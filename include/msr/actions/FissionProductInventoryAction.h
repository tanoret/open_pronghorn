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
#include "FissionProductData.h"

#include <map>
#include <set>
#include <string>
#include <vector>

/**
 * Action that assembles a fission-product inventory (depletion) problem in a circulating fuel salt
 * from the built-in nuclide database.
 *
 * Given a list of nuclides, the action creates one linear finite-volume variable per nuclide
 * (molar inventory [mol/m^3]) and the kernels of the Bateman equations with MSR removal terms:
 *
 *   dN_i/dt + div(u N_i) - div(D grad N_i) =
 *       Y_i F + sum_p (lambda_p b_pi + sigma_p phi delta_pi) N_p - (lambda_i + sigma_i phi + r_i)
 * N_i
 *
 * where F is the fission rate density, Y_i the effective yield of nuclide i, b_pi the decay
 * branching of a tracked parent p to i, sigma the capture cross section, phi the one-group flux and
 * r_i the sum of the first-order removal rates of i's chemical group (noble-gas sparging,
 * noble-metal plating and noble-metal entrainment in the off-gas). Optionally the removed
 * inventories are tracked as <nuclide>_offgas and <nuclide>_plated variables that keep decaying
 * (so plated Mo-99 becomes plated Tc-99). Noble metals may also plate on boundaries through a
 * mass-transfer flux J = k_m N.
 *
 * Each nuclide is placed in its own linear system "<nuclide>_sys"; the action orders them parent
 * before daughter (printed with 'verbose'), so that a single Gauss-Seidel sweep per time step
 * solves the lower-triangular chain exactly when the systems are listed in that order.
 */
class FissionProductInventoryAction : public Action
{
public:
  static InputParameters validParams();

  FissionProductInventoryAction(const InputParameters & params);

  virtual void act() override;

protected:
  /// A source of nuclide inventory from a tracked parent inventory
  struct ParentFeed
  {
    std::string parent;
    Real decay_coefficient = 0.0;
    Real capture_cross_section = 0.0;
  };

  /// Description of one variable to be created
  struct InventoryPlan
  {
    /// Variable name
    std::string name;
    /// Nuclide this inventory belongs to
    std::string nuclide;
    /// "salt", "offgas" or "plated"
    std::string location;
    Real initial_condition = 0.0;
    /// Effective fission yield (salt inventories only)
    Real yield = 0.0;
    /// Decay constant [1/s]
    Real lambda = 0.0;
    /// Capture cross section [barn] (salt inventories only)
    Real sigma = 0.0;
    /// Removal-rate functors (salt inventories only)
    std::vector<MooseFunctorName> removal;
    /// Parent feeds by decay or capture
    std::vector<ParentFeed> feeds;
    /// Transfer from the salt inventory (removed inventories only)
    MooseFunctorName transfer_rate;
    /// Whether this inventory is transported with the salt
    bool transported = false;
  };

  void buildPlan();
  void orderPlan();
  void addVariables();
  void addInitialConditions();
  void addKernels();
  void addBoundaryConditions();
  void maybeAssignBlocks(InputParameters & params) const;
  std::string systemName(const std::string & variable) const { return variable + "_sys"; }

  /// Name of the removal-rate functor for a chemical group and pathway, empty if none
  MooseFunctorName groupRemoval(const std::string & group, const std::string & location) const;

  /// The nuclide database
  const FissionProducts::FissionProductDatabase _db;
  /// Tracked nuclides
  std::vector<std::string> _nuclides;
  /// Variables to create, ordered parent before daughter
  std::vector<InventoryPlan> _plan;
};
