//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* OpenPronghorn is powered by the MOOSE Framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Copyright 2025, Battelle Energy Alliance, LLC
//* ALL RIGHTS RESERVED
//*

#include "NTPDarcyMassFlux.h"

registerMooseObject("OpenPronghornApp", NTPDarcyMassFlux);

InputParameters
NTPDarcyMassFlux::validParams()
{
  InputParameters params = NTPDarcyFluxInterface::validParams();
  params.addClassDescription("Gas mass conservation with the Darcy-Forchheimer face mass flux, "
                             "including the momentum flux of free channels.");
  return params;
}

NTPDarcyMassFlux::NTPDarcyMassFlux(const InputParameters & params) : NTPDarcyFluxInterface(params)
{
}

ADReal
NTPDarcyMassFlux::computeQpResidual()
{
  return faceMassFlux(determineState());
}
