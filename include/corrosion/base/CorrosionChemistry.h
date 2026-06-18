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

#include "MooseTypes.h"

/**
 * Chemical constants and Faradaic unit conversions for the molten salt corrosion and plating
 * framework. These mirror the validated reference model in msr_corrosion_bv/chemistry.py so that the
 * mechanistic Butler-Volmer objects and the ported effective correlation share one definition of the
 * electrochemistry.
 *
 * The Faradaic conversions take the valence, molar mass and density explicitly (rather than looking
 * them up by material class) so that they remain pure and independent of the database; the class
 * lookups live in MoltenSaltCorrosionData.
 */
namespace Corrosion
{

/// Universal gas constant [J/(mol K)] (matches chemistry.py R_GAS).
const Real R_gas = 8.31446261815324;
/// Faraday constant [C/mol] (matches chemistry.py FARADAY).
const Real faraday = 96485.33212;
/// Seconds in a Julian year [s] (matches chemistry.py SEC_PER_YEAR).
const Real seconds_per_year = 365.25 * 24.0 * 3600.0;
/// Hours in a Julian year [h] (matches chemistry.py HOURS_PER_YEAR).
const Real hours_per_year = 365.25 * 24.0;
/// Reference temperature for the calibrated correlation [K], 650 C (matches chemistry.py T_REF_K).
const Real T_ref_K = 650.0 + 273.15;

/// Exponential with the argument clipped to [lo, hi], matching model.py exp_clip. The clip keeps the
/// effective Butler-Volmer rate finite when the driving force is large (Newton robustness).
Real expClip(Real x, Real lo = -60.0, Real hi = 60.0);

/// Natural logarithm of x floored at a small positive value, matching chemistry.py safe_log.
Real safeLog(Real x, Real floor = 1.0e-30);

/// Convert an areal mass change [mg/cm^2] to an equivalent dense-layer thickness [um].
Real mgCm2ToUm(Real value_mg_cm2, Real density_g_cm3);

/// Convert a penetration depth [um] to an equivalent areal mass change [mg/cm^2].
Real umToMgCm2(Real value_um, Real density_g_cm3);

/**
 * Convert an anodic current density [A/cm^2] to an equivalent penetration rate [um/y].
 *
 * The conversion is Faradaic and assumes the selected species controls metal loss:
 *   cm_per_s = current * molar_mass / (valence * F * density)
 * (port of chemistry.py::corrosion_current_to_um_y).
 */
Real corrosionCurrentToUmY(Real current_a_cm2, Real valence, Real molar_mass_g_mol, Real density_g_cm3);

/// Inverse Faradaic conversion: penetration rate [um/y] to anodic current density [A/cm^2]
/// (port of chemistry.py::um_y_to_corrosion_current).
Real umYToCorrosionCurrent(Real rate_um_y, Real valence, Real molar_mass_g_mol, Real density_g_cm3);

} // namespace Corrosion
