//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "AdvancedCorrosionModelUserObject.h"

InputParameters
AdvancedCorrosionModelUserObject::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Common scalar-output interface for advanced lumped molten-salt corrosion models.");
  return params;
}

AdvancedCorrosionModelUserObject::AdvancedCorrosionModelUserObject(
    const InputParameters & parameters)
  : GeneralUserObject(parameters)
{
}

