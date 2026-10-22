//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVAdvectionDiffusionBC.h"
#include "LinearFVCHTBCInterface.h"
#include "NS.h"

/**
 * Conjugate-heat-transfer condition for a turbulent wall resolved with wall functions, used on
 * both sides of the interface (parameter `side`).
 *
 * The wall heat-transfer coefficient comes from the thermal law of the wall (Jayatilleke 1969):
 *
 *   T+ = Pr_t (ln(E y+)/kappa + P(Pr/Pr_t)),
 *   P  = 9.24 ((Pr/Pr_t)^0.75 - 1) (1 + 0.28 exp(-0.007 Pr/Pr_t))
 *
 * h_wf = rho cp u_tau / T+ in the log layer (y+ >= 30), k / y in the viscous sublayer
 * (y+ <= 5) and a linear blend in between. The P-function makes the coefficient correct for
 * Prandtl numbers far from one: for a molten fluoride salt (Pr ~ 15-20) the near-wall cell
 * conduction k_eff / y overpredicts h by more than an order of magnitude.
 *
 * Between the fluid cell P (temperature T_P) and the solid cell S (T_S) there are two
 * resistances in series: the wall function 1/h_wf and the half solid cell d_S/k_S.
 *   side = fluid: outward flux h_wf (T_P - T_w), with T_w the solid surface temperature
 *                 (interface_temperature_solid_<bnd>); hands T_P to the solid.
 *   side = solid: outward flux U (T_S - T_P), U = 1/(1/h_wf + d_S/k_S), with T_P from the fluid
 *                 (interface_temperature_fluid_<bnd>); hands T_w = T_S + U d_S/k_S (T_P - T_S)
 *                 to the fluid.
 * Each side is implicit in its own temperature and lags the other, so neither sub-problem is
 * left without a temperature reference (a solid cooled only by the fluid stays well posed), no
 * stabilizing coefficient is needed, and the two fluxes are identical at convergence.
 */
class LinearFVTurbulentWallCHTBC : public LinearFVAdvectionDiffusionBC,
                                   public LinearFVCHTBCInterface
{
public:
  static InputParameters validParams();

  LinearFVTurbulentWallCHTBC(const InputParameters & parameters);

  virtual Real computeBoundaryValue() const override;
  virtual Real computeBoundaryNormalGradient() const override;
  virtual Real computeBoundaryValueMatrixContribution() const override;
  virtual Real computeBoundaryValueRHSContribution() const override;
  virtual Real computeBoundaryGradientMatrixContribution() const override;
  virtual Real computeBoundaryGradientRHSContribution() const override;
  virtual bool includesMaterialPropertyMultiplier() const override { return true; }

  /// Jayatilleke P-function
  static Real jayatillekeP(Real Pr, Real Pr_t);

protected:
  /// Wall-function heat-transfer coefficient (W/m^2/K) of the fluid cell next to the face
  Real wallFunctionHTC() const;
  /// Coefficient between this side's cell and the other side's reference temperature
  Real couplingCoefficient() const;
  /// Temperature of the cell on this side of the face
  Real cellTemperature() const;
  /// Reference temperature from the other side (T_w for the fluid, T_P for the solid)
  Real otherTemperature() const;

  /// Which side of the interface this condition is applied to
  const bool _is_fluid_side;
  /// Velocity components
  const Moose::Functor<Real> & _u;
  const Moose::Functor<Real> * const _v;
  const Moose::Functor<Real> * const _w;
  /// Fluid properties
  const Moose::Functor<Real> & _rho;
  const Moose::Functor<Real> & _mu;
  const Moose::Functor<Real> & _cp;
  /// Fluid molecular conductivity (sets the Prandtl number and the sublayer coefficient)
  const Moose::Functor<Real> & _k_fluid;
  /// Conductivity used by this side's conduction kernel on the face
  const Moose::Functor<Real> & _k_kernel;
  /// Turbulent kinetic energy (non-equilibrium treatment only)
  const Moose::Functor<Real> * const _tke;
  /// Turbulent Prandtl number
  const Real _Pr_t;
  const Real _C_mu;
  const NS::WallTreatmentEnum _wall_treatment;
  /// Temperature handed over by the other side through the CHT handler
  const Moose::Functor<Real> & _other_temperature;
};
