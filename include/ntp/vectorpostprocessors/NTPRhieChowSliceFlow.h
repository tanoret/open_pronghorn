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

#pragma once

#include "GeneralVectorPostprocessor.h"
#include "BlockRestrictable.h"

class RhieChowMassFlux;
class MeshDivision;

/**
 * Mass flow rate through the interface between two blocks, one value per division of a mesh
 * division, from the face mass fluxes of the segregated solvers' Rhie-Chow object: the sum over
 * the interface faces of m_f A_f, positive from 'from_block' into 'to_block'. Each face is
 * attributed to the division of its 'to_block' cell. These are the fluxes the discrete mass and
 * energy balances see, so the flows are exact; between the plenum and the cold frit of a
 * particle-bed element they are the slice flows of the orifice design (NTPOrificeFlowDesignAux).
 */
class NTPRhieChowSliceFlow : public GeneralVectorPostprocessor
{
public:
  static InputParameters validParams();

  NTPRhieChowSliceFlow(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;

protected:
  /// The divisions (slices)
  const MeshDivision & _division;
  /// Subdomains the flow leaves and enters
  std::set<SubdomainID> _from_ids;
  std::set<SubdomainID> _to_ids;
  /// The Rhie-Chow object holding the face fluxes
  const RhieChowMassFlux * _rc;
  /// Flow per division
  VectorPostprocessorValue & _flow;
  /// Index of each division
  VectorPostprocessorValue & _division_index;
};
