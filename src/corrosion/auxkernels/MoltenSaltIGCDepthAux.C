//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MoltenSaltIGCDepthAux.h"

registerMooseObject("OpenPronghornApp", MoltenSaltIGCDepthAux);

InputParameters
MoltenSaltIGCDepthAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Evaluates MoltenSaltCorrosionModel::igcDepthFromCorrosionDepthUm() "
                             "from an accumulated uniform corrosion-depth variable.");
  params.addRequiredCoupledVar("corrosion_depth",
                               "Accumulated uniform corrosion-depth variable [um].");
  params.addParam<DataFileName>(
      "database", "corrosion_database.json", "JSON molten-salt corrosion database.");
  params.addParam<std::string>("salt_class", "generic_salt", "Salt class.");
  params.addParam<std::string>(
      "redox_class", "purified_baseline", "Redox/effective-overpotential class.");
  return params;
}

MoltenSaltIGCDepthAux::MoltenSaltIGCDepthAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _corrosion_depth(coupledValue("corrosion_depth")),
    _database(getParam<DataFileName>("database")),
    _model(_database)
{
  if (isNodal())
    paramError("variable", "MoltenSaltIGCDepthAux requires an elemental auxiliary variable.");

  _features.salt_class = getParam<std::string>("salt_class");
  _features.redox_class = getParam<std::string>("redox_class");
}

Real
MoltenSaltIGCDepthAux::computeValue()
{
  return _model.igcDepthFromCorrosionDepthUm(_features, _corrosion_depth[_qp]);
}
