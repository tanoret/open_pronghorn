//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "GeneralPostprocessor.h"

#include <string>

class AdvancedCorrosionModelUserObject;

/** Exposes one selected scalar from an advanced corrosion-model UserObject. */
class AdvancedCorrosionModelPostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  AdvancedCorrosionModelPostprocessor(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}
  virtual Real getValue() const override { return _value; }

protected:
  const AdvancedCorrosionModelUserObject & _model;
  const std::string _quantity;
  Real _value = 0.0;
};
