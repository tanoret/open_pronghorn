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

#pragma once

#include "ADKernelGrad.h"

/**
 * Nernst-Planck diffusion and migration flux of one dissolved cation M^{z+}.
 *
 * The molar flux is N = -D grad c - (z F / R T) D c grad(phi), so the weak (grad-test) residual of
 * div(N) contributes
 *   (grad test, D grad c) + (grad test, (z F / R T) D c grad phi).
 * Advection is added separately by CorrosionAdvection, and the Butler-Volmer flux enters as the
 * natural boundary term handled by the interface kernel or boundary condition.
 */
class CorrosionNernstPlanckFlux : public ADKernelGrad
{
public:
  static InputParameters validParams();

  CorrosionNernstPlanckFlux(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  /// Molecular diffusivity of the cation [m^2/s].
  const Real _diffusivity;
  /// Charge number z of the cation.
  const Real _valence;
  /// Temperature functor [K].
  const Moose::Functor<Real> & _temperature;
  /// Whether the migration term (electric-field coupling) is active.
  const bool _has_potential;
  /// Gradient of the electric potential [V/m].
  const ADVariableGradient & _grad_potential;
};
