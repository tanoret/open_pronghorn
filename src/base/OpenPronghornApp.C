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

#include "OpenPronghornApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

#include "OpenPronghornRevision.h"

InputParameters
OpenPronghornApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

OpenPronghornApp::OpenPronghornApp(const InputParameters & parameters) : MooseApp(parameters)
{
  OpenPronghornApp::registerAll(_factory, _action_factory, _syntax);
}

OpenPronghornApp::~OpenPronghornApp() {}

void
OpenPronghornApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<OpenPronghornApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"OpenPronghornApp"});
  Registry::registerActionsTo(af, {"OpenPronghornApp"});

  // Register the application data directory so that DataFileName parameters (e.g. the molten salt
  // radiolysis chemistry database) resolve against open_pronghorn/data.
  registerAppDataFilePath("open_pronghorn");

  // Molten salt radiolysis action: builds species variables and chemistry kernels
  registerSyntax("MoltenSaltRadiolysisAction", "MoltenSaltRadiolysis");

  // Molten salt corrosion and plating action: builds species variables, transport, the electric
  // potential equation and the Butler-Volmer electrode kinetics (finite-element / Newton)
  registerSyntax("CorrosionPlatingAction", "CorrosionPlating");

  // Flow-coupled corrosion action: salt-side corrosion products as linear finite-volume passive
  // scalars with a Butler-Volmer wall reaction, for coupling into a flowing MSR (SIMPLE/PIMPLE)
  registerSyntax("CorrosionPlatingFlowAction", "CorrosionPlatingFlow");
}

void
OpenPronghornApp::registerApps()
{
  registerApp(OpenPronghornApp);
  ModulesApp::registerApps();
}

std::string
OpenPronghornApp::getInstallableInputs() const
{
  return OPENPRONGHORN_INSTALLABLE_DIRS;
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
OpenPronghornApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  OpenPronghornApp::registerAll(f, af, s);
}
extern "C" void
OpenPronghornApp__registerApps()
{
  OpenPronghornApp::registerApps();
}
