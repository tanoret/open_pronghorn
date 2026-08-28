//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "AdvancedCorrosionModelUserObject.h"
#include "AdvancedCorrosionModelData.h"
#include "MSTDBTCData.h"
#include "MSTDBTCStandardStateCorrosionModel.h"
#include "MoltenSaltCorrosionData.h"

/**
 * Static MSTDB-TC standard-state corrosion endpoint.
 *
 * Raw MSTDB-TC files are external inputs. Their edition and SHA-256 digests are checked against the
 * advanced-model parameter database before the endpoint is made available to other MOOSE objects.
 */
class MSTDBTCCorrosionUserObject : public AdvancedCorrosionModelUserObject
{
public:
  static InputParameters validParams();

  MSTDBTCCorrosionUserObject(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual void finalize() override {}

  virtual Real scalarValue(const std::string & quantity) const override;

  const Corrosion::MSTDBTCCorrosionFeatures & features() const { return _features; }
  const Corrosion::MSTDBTCCorrosionResult & result() const { return _result; }
  const Corrosion::MSTDBTCPair & thermodynamics() const { return _thermodynamics; }
  const Corrosion::AdvancedCorrosionModelDatabase & parameterDatabase() const
  {
    return _advanced_database;
  }

protected:
  const Corrosion::MoltenSaltCorrosionDatabase _base_database;
  const Corrosion::AdvancedCorrosionModelDatabase _advanced_database;
  const Corrosion::MSTDBTCPair _thermodynamics;
  const Corrosion::MSTDBTCStandardStateCorrosionModel _model;
  const Corrosion::MSTDBTCCorrosionFeatures _features;
  const Corrosion::MSTDBTCCorrosionResult _result;
};

