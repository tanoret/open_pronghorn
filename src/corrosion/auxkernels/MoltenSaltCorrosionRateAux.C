//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MoltenSaltCorrosionRateAux.h"

registerMooseObject("OpenPronghornApp", MoltenSaltCorrosionRateAux);

InputParameters
MoltenSaltCorrosionRateAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Evaluates MoltenSaltCorrosionModel::corrosionRateUmY() with a local temperature functor.");
  params.addRequiredParam<MooseFunctorName>("temperature", "Local temperature functor [K].");
  params.addParam<DataFileName>(
      "database", "corrosion_database.json", "JSON molten-salt corrosion database.");
  params.addParam<std::string>("material_class", "generic_metal", "Alloy class.");
  params.addParam<std::string>("salt_class", "generic_salt", "Salt class.");
  params.addParam<std::string>(
      "redox_class", "purified_baseline", "Redox/effective-overpotential class.");
  params.addParam<std::string>("position_class", "nominal", "Loop position class.");
  params.addParam<Real>("flow_factor", 0.75, "Circulation/mass-transfer factor.");
  params.addParam<Real>("delta_T_C", 0.0, "Loop thermal gradient [C].");
  return params;
}

MoltenSaltCorrosionRateAux::MoltenSaltCorrosionRateAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _temperature(getFunctor<Real>("temperature")),
    _database(getParam<DataFileName>("database")),
    _model(_database)
{
  if (isNodal())
    paramError("variable", "MoltenSaltCorrosionRateAux requires an elemental auxiliary variable.");

  _features.material_class = getParam<std::string>("material_class");
  _features.salt_class = getParam<std::string>("salt_class");
  _features.redox_class = getParam<std::string>("redox_class");
  _features.position_class = getParam<std::string>("position_class");
  _features.flow_factor = getParam<Real>("flow_factor");
  _features.delta_T_C = getParam<Real>("delta_T_C");
}

Real
MoltenSaltCorrosionRateAux::computeValue()
{
  const Moose::ElemQpArg qp_arg = {_current_elem, _qp, _qrule, _q_point[_qp]};
  _features.temperature_K = _temperature(qp_arg, determineState());
  return _model.corrosionRateUmY(_features);
}
