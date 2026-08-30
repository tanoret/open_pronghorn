//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "AdvancedCorrosionModelUserObject.h"
#include "AdvancedCorrosionModelData.h"
#include "DRIDNModel.h"

/**
 * Restartable time integrator for the lumped DRIDN corrosion model.
 *
 * The complete physical Context is supplied explicitly. The wrapper never selects geometry,
 * inventory, or closure behavior from experiment identifiers or requested output quantities.
 */
class DRIDNCorrosionUserObject : public AdvancedCorrosionModelUserObject
{
public:
  static InputParameters validParams();

  DRIDNCorrosionUserObject(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}

  virtual Real scalarValue(const std::string & quantity) const override;

  const Corrosion::DRIDNModel::Context & context() const { return _context; }
  const Corrosion::DRIDNModel::State & state() const { return _state; }
  Corrosion::DRIDNModel::Outputs outputs() const;
  Real elapsedTimeYears() const { return _elapsed_time_y; }

protected:
  const Corrosion::AdvancedCorrosionModelDatabase _advanced_database;
  const Corrosion::DRIDNModel::Parameters _parameters;
  const Corrosion::DRIDNModel::Context _context;
  const Corrosion::DRIDNModel::IntegrationOptions _integration_options;
  const Corrosion::DRIDNModel::ModelOptions _model_options;
  const Corrosion::DRIDNModel _model;

  Corrosion::DRIDNModel::State & _state;
  Real & _elapsed_time_y;
  unsigned int & _accepted_steps;
  unsigned int & _rejected_steps;
};
