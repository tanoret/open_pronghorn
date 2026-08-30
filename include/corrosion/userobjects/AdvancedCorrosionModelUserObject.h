//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "GeneralUserObject.h"

#include <string>

/**
 * Common scalar-output interface for the advanced lumped corrosion models.
 *
 * The interface lets one postprocessor expose both the static MSTDB-TC endpoint and the transient
 * DRIDN state without teaching the postprocessor either model's implementation details.
 */
class AdvancedCorrosionModelUserObject : public GeneralUserObject
{
public:
  static InputParameters validParams();

  AdvancedCorrosionModelUserObject(const InputParameters & parameters);

  /// Return one named scalar result. Implementations must reject unsupported names.
  virtual Real scalarValue(const std::string & quantity) const = 0;
};
