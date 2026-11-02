//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "MooseTypes.h"

#include <array>
#include <cstddef>
#include <limits>

namespace Corrosion
{

/**
 * Pure C++ Dynamic Redox-Inventory-Depletion Network (DRIDN).
 *
 * This class owns no MOOSE variables.  It evaluates and advances the 19 scalar states of the
 * well-mixed DRIDN model, leaving a later MOOSE wrapper to map those states to scalar variables.
 * The right-hand side uses years, micrometres, ppm, and mg/cm2, exactly as documented below. A
 * MOOSE residual whose time is seconds must therefore divide rhs() by seconds_per_year.
 *
 * The default model is deliberately stricter than the original Python research implementation:
 *  - geometry and inventory coupling are explicit and never inferred from response metadata;
 *  - the Nernst activity floor is separate from the actual dissolved inventory;
 *  - sink rates use the nonnegative physical inventory, not a fictitious concentration floor;
 *  - an affinity gate makes the dissolution branch vanish at and below equilibrium; and
 *  - advance() projects the surface/redox box constraints while rejecting negative inventories.
 *
 * legacyCompatibility() is provided only for algebraic parity tests at admissible states.  It
 * restores the Python model's always-positive charge-transfer branch but does not restore its
 * nonconservative negative-inventory projection.  The safer affinity-gated default has not been
 * refitted or validated against the legacy data set: a production validation claim requires a new,
 * mode-specific fit.  A legacy fitted vector and a corrected-mode result must never be combined
 * into one regression golden.
 */
class DRIDNModel
{
public:
  static constexpr std::size_t n_elements = 3;
  static constexpr std::size_t n_states = 19;
  static constexpr Real seconds_per_year = 31557600.0;

  enum class Element : std::size_t
  {
    Cr = 0,
    Fe = 1,
    Ni = 2
  };

  enum class InventoryScale
  {
    Explicit,
    MSRE,
    Loop
  };

  enum class DepositionClosure
  {
    FuelLike,
    FLiNaK
  };

  enum class ChargeTransferMode
  {
    /// Corrected dissolution-only law: zero at and below the weighted equilibrium affinity.
    AffinityGatedDissolution,
    /// Original Python law: an always-positive exponential penetration rate.
    LegacyIrreversible
  };

  using Triplet = std::array<Real, n_elements>;
  using State = std::array<Real, n_states>;

  /**
   * Fixed state ordering shared by rhs(), advance(), and every MOOSE wrapper:
   *
   *  - 0..2: Cr/Fe/Ni surface availability [-];
   *  - 3..5: dissolved Cr/Fe/Ni inventory [ppm];
   *  - 6..8: cumulative Cr/Fe/Ni source [ppm];
   *  - 9..11: coupon Cr/Fe/Ni deposit [mg/cm2];
   *  - 12: reaction-front depth [um];
   *  - 13: mass-equivalent recession [um];
   *  - 14: grain-boundary diffusion depth squared [um2];
   *  - 15: dimensionless log redox shift; and
   *  - 16..18: bulk-captured Cr/Fe/Ni inventory [ppm].
   *
   * rhs() returns each quantity per year where applicable.
   */
  static constexpr std::size_t surfaceIndex(const Element element)
  {
    return static_cast<std::size_t>(element);
  }
  static constexpr std::size_t dissolvedIndex(const Element element)
  {
    return 3 + static_cast<std::size_t>(element);
  }
  static constexpr std::size_t cumulativeSourceIndex(const Element element)
  {
    return 6 + static_cast<std::size_t>(element);
  }
  static constexpr std::size_t couponDepositIndex(const Element element)
  {
    return 9 + static_cast<std::size_t>(element);
  }
  static constexpr std::size_t front_depth = 12;
  static constexpr std::size_t mass_recession = 13;
  static constexpr std::size_t grain_boundary_depth_squared = 14;
  static constexpr std::size_t redox_shift = 15;
  static constexpr std::size_t bulkCaptureIndex(const Element element)
  {
    return 16 + static_cast<std::size_t>(element);
  }

  /**
   * The 26 DRIDN closure parameters, stored in their native/raw form.
   *
   * Member defaults are documented physical priors, not a validated calibration.  Production
   * wrappers must replace all 26 values from AdvancedCorrosionModelDatabase and record the
   * calibration identifier; partial overlay of a fitted vector is intentionally unsupported.
   */
  struct Parameters
  {
    Real log_rate_scale = 0.0;
    Real affinity_feedback_scale = 0.035;
    Real log_surface_reservoir_um = 3.6888794541139363;        // log(40 um)
    Real log_surface_replenishment_y_inv = 0.6931471805599453; // log(2/y)
    Real surface_availability_exponent = 0.8;
    Real surface_reservoir_cr_exponent = 0.8;
    Real log_dynamic_cr_exchange_bias = 0.0;
    Real log_dynamic_fe_capture_bias = 0.30;
    Real inventory_inhibition_scale = 0.75;
    Real log_redox_relaxation_y_inv = 0.0;
    Real redox_buffer_retention = 0.75;
    Real redox_consumption_per_um = 0.005;
    Real log_stress_interfacial_factor = 0.0;
    Real log_fluoride_impurity_interfacial_factor = 0.0;
    Real log_deposition_rate_y_inv_fuel = -0.6931471805599453;  // log(0.5/y)
    Real log_deposition_rate_y_inv_flinak = 0.6931471805599453; // log(2/y)
    Real bulk_capture_area_multiplier_fuel = 1.0;
    Real bulk_capture_area_multiplier_flinak = 0.5;
    Real log_bulk_precipitation_rate_scale = 0.0;
    Real log_inventory_scale_msre = 0.0;
    Real log_inventory_scale_loop = 0.0;
    Real log_deposit_area_scale_fuel = 0.0;
    Real log_deposit_area_scale_flinak = 0.0;
    Real log_mass_loss_scale = 0.0;
    Real gb_dynamic_scale = 1.0;
    Real damage_affinity_scale = 0.03;

    void validate() const;
  };

  /**
   * Previously hard-coded numerical/closure constants.  They are explicit here so sensitivity
   * studies cannot silently change the calibrated model definition.
   */
  struct ClosureConstants
  {
    Real surface_availability_floor = 0.015;
    Triplet replenishment_relative{{1.0, 3.0, 1.8}};
    Real minimum_capture_flow_factor = 0.05;
    Real bulk_precipitation_prefactor_y_inv = 0.35;
    Real maximum_bulk_precipitation_rate_y_inv = 200.0;
    Real minimum_mass_loss_fraction = 0.01;
    Real log_rate_clip = 60.0;
    Real relative_affinity_clip = 80.0;

    void validate() const;
  };

  /**
   * All time-independent terms for one well-mixed exposure.
   *
   * Upstream thermochemistry supplies affinity_baseline and log_charge_base_no_redox.  For the
   * corrected model it must aggregate repeated reaction species (the Fe self-buffer reaction has
   * zero standard-state Gibbs energy) and must evaluate the actual Be/BeF2 branch when present. No
   * response kind, experiment identifier, or source-name heuristic is accepted here.
   */
  struct Context
  {
    Triplet mass_fractions{{0.07, 0.05, 0.71}};      // alloy mass fractions [-]
    Triplet log_exchange_offsets{{0.0, 0.0, 0.0}};  // log rate offsets [-]
    Triplet affinity_baseline{{0.0, 0.0, 0.0}};     // dimensionless -Delta_r G/(RT)
    Triplet cold_capture_fraction{{0.0, 0.0, 0.0}}; // fractions [0,1]
    Triplet initial_dissolved_ppm{{0.0, 0.0, 0.0}}; // physical inventory [ppm]

    Real cr_fraction_ratio = 1.0;          // alloy Cr fraction / reference Cr fraction [-]
    Real density_g_cm3 = 8.3;              // alloy density [g/cm3]
    Real flow_factor = 1.0;                // explicit positive circulation factor [-]
    Real selectivity_scale = 0.0;          // affinity selectivity coefficient [-]
    Real product_activity_floor_ppm = 1.0; // Nernst-only activity floor [ppm]
    Real redox_shift_initial = 0.0;         // dimensionless log redox shift
    Real log_charge_base_no_redox = 0.0;   // log of charge-transfer rate [um/y]
    Real mass_transfer_rate_um_y = 1.0;    // limiting front rate [um/y]
    Real inventory_capacity_ppm = 1.0;     // inhibition/precipitation scale [ppm]
    Real area_to_salt_mass_cm2_g = 1.0;    // wetted area / salt mass [cm2/g]
    Real explicit_inventory_scale = 1.0;   // explicit geometry multiplier [-]
    Real inventory_coupling_factor = 1.0;  // source coupled to modeled inventory [0,1]
    Real deposit_area_factor = 1.0;        // coupon deposit-area mapping [-]
    Real mass_loss_fraction = 1.0;         // recession reported as mass loss [0,1]
    Real cr_diffusion_cm2_s = 0.0;         // solid Cr diffusivity [cm2/s]
    Real front_damage_multiplier = 1.0;    // reaction-front IGC multiplier [-]
    Real gb_length_multiplier = 1.0;       // grain-boundary length multiplier [-]

    InventoryScale inventory_scale = InventoryScale::Explicit;
    DepositionClosure deposition_closure = DepositionClosure::FuelLike;
    bool transient_redox = false;
    bool stress_interfacial_activation = false;
    bool fluoride_impurity_interfacial_activation = false;
    bool chloride_salt = false;

    void validate(const ClosureConstants & closures) const;
  };

  struct ModelOptions
  {
    ChargeTransferMode charge_transfer = ChargeTransferMode::AffinityGatedDissolution;

    void validate() const;

    static ModelOptions legacyCompatibility()
    {
      ModelOptions options;
      options.charge_transfer = ChargeTransferMode::LegacyIrreversible;
      return options;
    }
  };

  struct Rates
  {
    Triplet affinity{{0.0, 0.0, 0.0}};
    Triplet species_fraction{{0.0, 0.0, 0.0}};
    Triplet source_rate_ppm_y{{0.0, 0.0, 0.0}};
    Triplet deposition_rate_y_inv{{0.0, 0.0, 0.0}};
    Triplet bulk_capture_rate_y_inv{{0.0, 0.0, 0.0}};

    Real product_feedback = 0.0;
    Real surface_feedback = 0.0;
    Real effective_affinity = 0.0;
    Real charge_transfer_rate_um_y = 0.0;
    Real mass_transfer_rate_um_y = 0.0;
    Real front_rate_um_y = 0.0;
    Real bulk_precipitation_rate_y_inv = 0.0;
    Real deposit_conversion_mg_cm2_per_ppm = 0.0;
    Real loss_fraction = 0.0;
  };

  struct IntegrationOptions
  {
    Real relative_tolerance = 2.0e-7;
    Real absolute_tolerance = 1.0e-10;
    Real maximum_step_y = std::numeric_limits<Real>::infinity();
    Real minimum_step_y = 1.0e-14;
    unsigned int nominal_maximum_step_intervals = 120;
    unsigned int maximum_step_attempts = 200000;

    void validate() const;
  };

  struct IntegrationStats
  {
    unsigned int accepted_steps = 0;
    unsigned int rejected_steps = 0;
    Real final_step_y = 0.0;
  };

  struct Outputs
  {
    Real front_depth_um = 0.0;
    Real mass_recession_um = 0.0;
    Real mass_loss_mg_cm2 = 0.0;
    Real mass_gain_mg_cm2 = 0.0;
    Real igc_depth_um = 0.0;
    Real average_corrosion_rate_um_y = 0.0;
    Real instantaneous_front_rate_um_y = 0.0;
    Triplet dissolved_ppm{{0.0, 0.0, 0.0}};
    Triplet cumulative_source_ppm{{0.0, 0.0, 0.0}};
    Triplet coupon_deposit_mg_cm2{{0.0, 0.0, 0.0}};
    Triplet bulk_captured_ppm{{0.0, 0.0, 0.0}};
    Triplet surface_availability{{0.0, 0.0, 0.0}};
    Real redox_shift = 0.0;
    Real mass_balance_relative_error = 0.0;
  };

  DRIDNModel();
  explicit DRIDNModel(const Parameters & parameters);
  DRIDNModel(const Parameters & parameters, const ClosureConstants & closures);
  DRIDNModel(const Parameters & parameters,
             const ClosureConstants & closures,
             const ModelOptions & options);

  const Parameters & parameters() const { return _parameters; }
  const ClosureConstants & closures() const { return _closures; }
  const ModelOptions & options() const { return _options; }

  State initialState(const Context & context) const;
  Rates rates(const Context & context, const State & state) const;
  State rhs(const Context & context, const State & state) const;

  /**
   * Advance state by an explicit duration in years using adaptive RK4 step doubling.
   *
   * Nonphysical stages are rejected.  Element inventories are conserved across every accepted
   * step.  No default exposure duration is inferred by this class.
   */
  IntegrationStats advance(const Context & context, State & state, Real duration_y) const;
  IntegrationStats advance(const Context & context,
                           State & state,
                           Real duration_y,
                           const IntegrationOptions & options) const;

  Outputs outputs(const Context & context, const State & state, Real elapsed_time_y) const;

private:
  Rates ratesImpl(const Context & context, const State & state) const;
  State rhsImpl(const Context & context, const State & state) const;
  bool admissible(const Context & context, const State & state, Real tolerance) const;
  bool rk4Step(const Context & context,
               const State & state,
               Real step_y,
               Real admissibility_tolerance,
               State & result) const;
  void restoreConservation(const Context & context,
                           const Triplet & inventory_invariant,
                           State & state,
                           Real tolerance) const;
  Real inventoryScale(const Context & context) const;
  Real depositAreaScale(const Context & context) const;
  Real depositionBaseRate(const Context & context) const;
  Real bulkCaptureAreaMultiplier(const Context & context) const;

  Parameters _parameters;
  ClosureConstants _closures;
  ModelOptions _options;
};

} // namespace Corrosion
