//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "AuxKernel.h"
#include "MoltenSaltCorrosionData.h"
#include "MoltenSaltCorrosionModel.h"

/**
 * Evaluates IGC depth from an already accumulated uniform corrosion-depth variable.
 */
class MoltenSaltIGCDepthAux : public AuxKernel
{
public:
  static InputParameters validParams();

  MoltenSaltIGCDepthAux(const InputParameters & parameters);

protected:
  Real computeValue() override;

  const VariableValue & _corrosion_depth;
  Corrosion::MoltenSaltCorrosionDatabase _database;
  Corrosion::MoltenSaltCorrosionModel _model;
  Corrosion::CorrosionFeatures _features;
};
