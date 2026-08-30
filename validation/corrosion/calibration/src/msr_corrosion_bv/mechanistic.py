"""Reduced-mechanistic corrosion, mass-transfer, and deposition model for MSRs.

This module complements :mod:`msr_corrosion_bv.model`.  The effective Butler-Volmer
model remains the electrochemical calibration layer; this model reuses its fitted
redox log-driving-force and chromium solid-state diffusion coefficient while making
several transport and inventory mechanisms explicit:

* chromium activity dependence of alloy dissolution,
* charge-transfer and salt-side mass-transfer resistances in series,
* non-congruent mass loss for selective chromium removal,
* diffusion-length-controlled intergranular penetration,
* hot-to-cold supersaturation-driven deposition, and
* chromium inventory closure through a surface-area/salt-mass ratio.

It is deliberately a reduced mechanistic model, not a first-principles chemistry or
CFD solver.  Salt chemical-potential shifts, geometry ratios, capture efficiency,
and grain-boundary enhancement factors remain calibrated because the legacy
validation data generally do not contain complete salt compositions, velocities,
surface areas, or salt inventories.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .chemistry import R_GAS, SEC_PER_YEAR, T_REF_K, cr_weight_fraction, density
from .model import MoltenSaltBVModel


MECHANISTIC_ACTIVE_ROLES = {"direct", "range", "upper", "lower"}
MECHANISTIC_DIRECT_ROLES = {"direct", "range"}
MECHANISTIC_RESPONSE_KINDS = {
    "corrosion_rate_um_y",
    "corrosion_depth_um",
    "mass_loss_mg_cm2",
    "mass_gain_mg_cm2",
    "salt_cr_ppm",
    "salt_fe_decrease_ppm",
    "cr_diffusion_cm2_s",
    "redox_acceleration_ratio",
    "redox_acceleration_qualitative",
    "igc_depth_um",
}


@dataclass(frozen=True)
class MechanisticParameterSpec:
    name: str
    initial: float
    lower: float
    upper: float
    prior: float
    prior_sigma: float
    description: str


# These parameters are physically dimensioned where practical.  The reference
# electrochemical activation energy and redox offsets are intentionally NOT fitted
# here; they are inherited from the effective BV model so the comparison isolates
# the mechanistic transport/material layer.
MECHANISTIC_PARAMETER_SPECS: tuple[MechanisticParameterSpec, ...] = (
    MechanisticParameterSpec("log_front_rate0_um_y", math.log(1.0), math.log(0.05), math.log(20.0), math.log(1.0), 1.0, "Reference Cr-controlled dissolution-front rate at 650 C, log(um/y)."),
    MechanisticParameterSpec("cr_activity_exponent", 1.6, 0.2, 4.0, 1.6, 1.0, "Exponent mapping alloy Cr mass fraction to an effective Cr activity ratio."),
    MechanisticParameterSpec("delta_mu_flinak_kJ_mol", -8.0, -25.0, 20.0, -8.0, 10.0, "Effective FLiNaK chemical-potential shift relative to fuel fluoride; negative values increase dissolution drive under this sign convention."),
    MechanisticParameterSpec("delta_mu_flibe_kJ_mol", 2.0, -15.0, 20.0, 2.0, 10.0, "Effective FLiBe chemical-potential shift relative to fuel fluoride."),
    MechanisticParameterSpec("delta_mu_fluoroborate_kJ_mol", -10.0, -30.0, 20.0, -10.0, 10.0, "Effective fluoroborate chemical-potential shift relative to fuel fluoride."),
    MechanisticParameterSpec("delta_mu_chloride_kJ_mol", -5.0, -30.0, 20.0, -5.0, 10.0, "Effective chloride chemical-potential shift relative to fuel fluoride."),
    MechanisticParameterSpec("log_mass_transfer_cap_um_y", math.log(2000.0), math.log(30.0), math.log(1.0e5), math.log(2000.0), 1.5, "Reference salt-side mass-transfer capacity, log(um/y)."),
    MechanisticParameterSpec("Ea_mass_transfer_kJ_mol", 25.0, 0.0, 100.0, 25.0, 30.0, "Apparent activation energy for salt-side mass transfer."),
    MechanisticParameterSpec("flow_mass_transfer_exponent", 0.6, 0.1, 1.5, 0.6, 0.4, "Flow-factor exponent for salt-side mass transfer."),
    MechanisticParameterSpec("gb_diffusion_multiplier", 0.65, 0.1, 2.0, 0.65, 0.4, "Multiplier on the Cr solid-state diffusion length for ordinary IGC/void penetration."),
    MechanisticParameterSpec("multi_alloy_gb_multiplier", 9.0, 1.0, 30.0, 9.0, 6.0, "Grain-boundary/galvanic penetration multiplier in multi-alloy tests."),
    MechanisticParameterSpec("x750_gb_multiplier", 12.0, 1.0, 30.0, 12.0, 8.0, "Additional X-750 susceptibility multiplier in multi-alloy tests."),
    MechanisticParameterSpec("chloride_reaction_multiplier", 6.0, 1.0, 20.0, 6.0, 5.0, "Reaction-front multiplier for impure chloride IGC."),
    MechanisticParameterSpec("chloride_gb_multiplier", 4.0, 1.0, 20.0, 4.0, 5.0, "Diffusion-length multiplier for impure chloride IGC."),
    MechanisticParameterSpec("tellurium_gb_multiplier", 8.0, 1.0, 30.0, 8.0, 8.0, "Diffusion-length multiplier for tellurium-assisted intergranular attack."),
    MechanisticParameterSpec("cold_leg_capture_efficiency", 0.75, 0.2, 1.0, 0.75, 0.2, "Fraction of hot-side transported inventory captured by cold-side supersaturation."),
    MechanisticParameterSpec("effective_solution_enthalpy_fuel_kJ_mol", 70.0, 10.0, 250.0, 70.0, 50.0, "Effective enthalpy controlling temperature-dependent dissolved-metal capacity in fuel fluoride."),
    MechanisticParameterSpec("effective_solution_enthalpy_flinak_kJ_mol", 140.0, 10.0, 300.0, 140.0, 80.0, "Effective enthalpy controlling temperature-dependent dissolved-metal capacity in FLiNaK."),
    MechanisticParameterSpec("log_area_to_salt_mass_msre_cm2_g", math.log(0.13), math.log(0.01), math.log(2.0), math.log(0.13), 1.0, "Effective exposed-area/salt-mass ratio for MSRE-scale inventory closure, log(cm2/g)."),
    MechanisticParameterSpec("log_area_to_salt_mass_loop_cm2_g", math.log(3.0), math.log(0.1), math.log(30.0), math.log(3.0), 1.2, "Effective exposed-area/salt-mass ratio for small-loop inventory closure, log(cm2/g)."),
    MechanisticParameterSpec("fe_to_cr_inventory_ratio", 0.20, 0.02, 1.0, 0.20, 0.25, "Effective magnitude ratio between salt Fe decrease and Cr increase."),
    MechanisticParameterSpec("mass_loss_fraction_logit", math.log(0.72 / 0.28), -4.0, 4.0, math.log(0.72 / 0.28), 1.0, "Reference fraction of the affected alloy layer removed as net mass."),
    MechanisticParameterSpec("mass_loss_cr_exponent", -0.25, -2.0, 1.0, -0.25, 0.5, "Cr-content exponent for non-congruent mass loss."),
    MechanisticParameterSpec("flinak_high_cr_selectivity_exponent", -0.8, -2.5, 0.5, -0.8, 0.7, "Additional non-congruent selective-leaching exponent for high-Cr alloys in FLiNaK."),
    MechanisticParameterSpec("log_dissolved_cr_capacity_ppm", math.log(1800.0), math.log(300.0), math.log(10000.0), math.log(1800.0), 0.8, "Effective dissolved-Cr inventory capacity used to saturate long-duration loop inventories, log(ppm)."),
)

MECHANISTIC_PARAMETER_NAMES = tuple(spec.name for spec in MECHANISTIC_PARAMETER_SPECS)


def mechanistic_initial_vector() -> np.ndarray:
    return np.asarray([spec.initial for spec in MECHANISTIC_PARAMETER_SPECS], dtype=float)


def mechanistic_lower_bounds() -> np.ndarray:
    return np.asarray([spec.lower for spec in MECHANISTIC_PARAMETER_SPECS], dtype=float)


def mechanistic_upper_bounds() -> np.ndarray:
    return np.asarray([spec.upper for spec in MECHANISTIC_PARAMETER_SPECS], dtype=float)


def mechanistic_vector_to_params(vector: Iterable[float]) -> dict[str, float]:
    return {name: float(value) for name, value in zip(MECHANISTIC_PARAMETER_NAMES, vector)}


def _is_positive_finite(value: Any) -> bool:
    try:
        return bool(np.isfinite(float(value)) and float(value) > 0.0)
    except Exception:
        return False


def _row_value(row: pd.Series | Mapping[str, Any], name: str, default: Any = None) -> Any:
    value = row.get(name, default)
    if value is None:
        return default
    try:
        if bool(pd.isna(value)):
            return default
    except Exception:
        pass
    return value


def mechanistic_scope_reason(row: pd.Series | Mapping[str, Any]) -> str:
    """Explain whether a validation row belongs to the reduced-mechanistic scope."""
    role = str(_row_value(row, "fit_role", ""))
    if role not in MECHANISTIC_ACTIVE_ROLES:
        return "not an active numerical/inequality response constraint"
    kind = str(_row_value(row, "response_kind", ""))
    if kind not in MECHANISTIC_RESPONSE_KINDS:
        return f"response kind '{kind}' remains in the effective/auxiliary model"
    text = " ".join(
        str(_row_value(row, key, "")).lower()
        for key in ("model_mapping", "experiment_class", "observable", "location_or_condition")
    )
    if "system metal-transfer" in text or "metal transfer rate inferred" in text:
        return "system-level inferred transfer proxy; no coupon-scale mechanistic geometry is available"
    if "external oxidation" in text or "under heaters in air" in text:
        return "external air oxidation is outside the molten-salt mechanistic model"
    return "supported"


def mechanistic_is_supported(row: pd.Series | Mapping[str, Any]) -> bool:
    return mechanistic_scope_reason(row) == "supported"


class MechanisticCorrosionModel:
    """Reduced mechanistic layer sharing electrochemical calibration with a BV model."""

    def __init__(self, base_model: MoltenSaltBVModel, params: Mapping[str, float] | None = None) -> None:
        self.base_model = base_model
        self.params = mechanistic_vector_to_params(mechanistic_initial_vector()) if params is None else dict(params)

    @classmethod
    def from_vector(cls, base_model: MoltenSaltBVModel, vector: Iterable[float]) -> "MechanisticCorrosionModel":
        return cls(base_model, mechanistic_vector_to_params(vector))

    def _param(self, name: str) -> float:
        return float(self.params[name])

    def salt_delta_mu_kJ_mol(self, salt_class: str) -> float:
        mapping = {
            "flinak": "delta_mu_flinak_kJ_mol",
            "flibe": "delta_mu_flibe_kJ_mol",
            "fluoroborate": "delta_mu_fluoroborate_kJ_mol",
            "chloride": "delta_mu_chloride_kJ_mol",
        }
        key = mapping.get(salt_class)
        return 0.0 if key is None else self._param(key)

    def redox_log_drive(self, redox_class: str) -> float:
        # The effective model defines its redox term directly in log-current space.
        # Subtracting the purified reference preserves that calibrated electrochemical
        # boundary condition without introducing a second unidentifiable redox fit.
        return self.base_model.redox_offset(redox_class) - self.base_model.redox_offset("purified_baseline")

    def cr_diffusion_cm2_s(self, row: pd.Series | Mapping[str, Any]) -> float:
        return self.base_model.cr_diffusion_cm2_s(row)

    def diffusion_length_um(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = float(_row_value(row, "time_years", 0.0) or 0.0)
        if time_y <= 0.0:
            return 0.0
        D = max(self.cr_diffusion_cm2_s(row), 0.0)
        return 2.0 * math.sqrt(D * time_y * SEC_PER_YEAR / math.pi) * 1.0e4

    def dissolution_front_rate_um_y(
        self,
        row: pd.Series | Mapping[str, Any],
        *,
        redox_override: str | None = None,
        temperature_K: float | None = None,
    ) -> float:
        p = self.params
        T = float(temperature_K if temperature_K is not None else (_row_value(row, "temperature_K", T_REF_K) or T_REF_K))
        material = str(_row_value(row, "material_class", "generic_metal"))
        salt = str(_row_value(row, "salt_class", "generic_salt"))
        redox = redox_override or str(_row_value(row, "redox_class", "purified_baseline"))
        flow = max(float(_row_value(row, "flow_factor", 0.75) or 0.75), 1.0e-3)
        if salt == "no_salt" or flow < 0.05:
            return 1.0e-9

        cr_ratio = max(cr_weight_fraction(material) / 0.07, 1.0e-6)
        Ea_corr = float(self.base_model.params["Ea_corr_kJ_mol"])
        thermal = Ea_corr * 1000.0 / R_GAS * (1.0 / T_REF_K - 1.0 / T)
        salt_drive = -self.salt_delta_mu_kJ_mol(salt) * 1000.0 / (R_GAS * T)
        log_charge_transfer = (
            p["log_front_rate0_um_y"]
            + thermal
            + p["cr_activity_exponent"] * math.log(cr_ratio)
            + salt_drive
            + self.redox_log_drive(redox)
        )
        rate_charge_transfer = float(np.exp(np.clip(log_charge_transfer, -60.0, 60.0)))

        thermal_mt = p["Ea_mass_transfer_kJ_mol"] * 1000.0 / R_GAS * (1.0 / T_REF_K - 1.0 / T)
        log_mass_transfer = (
            p["log_mass_transfer_cap_um_y"]
            + thermal_mt
            + p["flow_mass_transfer_exponent"] * math.log(flow)
        )
        rate_mass_transfer = float(np.exp(np.clip(log_mass_transfer, -60.0, 60.0)))
        return 1.0 / (1.0 / max(rate_charge_transfer, 1.0e-30) + 1.0 / max(rate_mass_transfer, 1.0e-30))

    def mass_loss_fraction(self, row: pd.Series | Mapping[str, Any]) -> float:
        material = str(_row_value(row, "material_class", "generic_metal"))
        salt = str(_row_value(row, "salt_class", "generic_salt"))
        cr_ratio = max(cr_weight_fraction(material) / 0.07, 1.0e-6)
        logit = self._param("mass_loss_fraction_logit")
        reference_fraction = 1.0 / (1.0 + math.exp(-logit))
        value = reference_fraction * cr_ratio ** self._param("mass_loss_cr_exponent")
        if salt == "flinak" and cr_ratio > 1.0:
            value *= cr_ratio ** self._param("flinak_high_cr_selectivity_exponent")
        return float(np.clip(value, 0.05, 1.0))

    def corrosion_rate_um_y(self, row: pd.Series | Mapping[str, Any], redox_override: str | None = None) -> float:
        return self.dissolution_front_rate_um_y(row, redox_override=redox_override) * self.mass_loss_fraction(row)

    def corrosion_depth_um(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = float(_row_value(row, "time_years", np.nan) or np.nan)
        if not np.isfinite(time_y):
            return np.nan
        # This is the selective depletion/reaction-front depth, not dense-layer
        # equivalent recession.  Net mass recession is handled separately.
        return self.dissolution_front_rate_um_y(row) * max(time_y, 0.0)

    def mass_loss_mg_cm2(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = float(_row_value(row, "time_years", np.nan) or np.nan)
        if not np.isfinite(time_y):
            return np.nan
        material = str(_row_value(row, "material_class", "generic_metal"))
        depth_equiv = self.corrosion_rate_um_y(row) * max(time_y, 0.0)
        return depth_equiv * density(material) * 0.1

    def _cold_temperature_K(self, row: pd.Series | Mapping[str, Any], hot_temperature_K: float) -> float:
        cold_C = _row_value(row, "temperature_cold_C", None)
        if cold_C is not None:
            try:
                if np.isfinite(float(cold_C)):
                    return float(cold_C) + 273.15
            except Exception:
                pass
        dT = max(float(_row_value(row, "delta_T_C", 0.0) or 0.0), 0.0)
        return max(hot_temperature_K - dT, 1.0)

    def _supersaturation_fraction(self, salt_class: str, hot_K: float, cold_K: float) -> float:
        if hot_K <= cold_K:
            return 0.0
        if salt_class == "flinak":
            H_kJ_mol = self._param("effective_solution_enthalpy_flinak_kJ_mol")
        else:
            H_kJ_mol = self._param("effective_solution_enthalpy_fuel_kJ_mol")
        exponent = -H_kJ_mol * 1000.0 / R_GAS * max(1.0 / cold_K - 1.0 / hot_K, 0.0)
        return float(np.clip(1.0 - math.exp(exponent), 0.0, 1.0))

    def mass_gain_mg_cm2(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = float(_row_value(row, "time_years", np.nan) or np.nan)
        if not np.isfinite(time_y):
            return np.nan
        dT = max(float(_row_value(row, "delta_T_C", 0.0) or 0.0), 0.0)
        if dT <= 0.0:
            return 0.0
        T_row = float(_row_value(row, "temperature_K", T_REF_K) or T_REF_K)
        position = str(_row_value(row, "position_class", "nominal"))
        hot_K = T_row + dT if position == "cold_leg" else T_row
        cold_K = T_row if position == "cold_leg" else self._cold_temperature_K(row, hot_K)
        salt = str(_row_value(row, "salt_class", "generic_salt"))
        material = str(_row_value(row, "material_class", "generic_metal"))

        donor_front_rate = self.dissolution_front_rate_um_y(row, temperature_K=hot_K)
        donor_mass_mg_cm2 = donor_front_rate * max(time_y, 0.0) * density(material) * 0.1
        capture = self._param("cold_leg_capture_efficiency")
        supersaturation = self._supersaturation_fraction(salt, hot_K, cold_K)
        return donor_mass_mg_cm2 * capture * supersaturation

    def _inventory_source_material(self, row: pd.Series | Mapping[str, Any]) -> str:
        explicit = _row_value(row, "inventory_source_material", None)
        if explicit:
            return str(explicit)
        material = str(_row_value(row, "material_class", "generic_metal"))
        salt = str(_row_value(row, "salt_class", "generic_salt"))
        # The legacy inventory rows store the measured species (Cr) rather than the
        # donor structural alloy.  HN is the documented donor in the fluoride-fuel
        # loop/MSRE inventory cases represented by the current validation table.
        if material == "generic_metal" and salt == "fluoride_fuel":
            return "hastelloy_n"
        return material

    def salt_cr_ppm(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = float(_row_value(row, "time_years", np.nan) or np.nan)
        if not np.isfinite(time_y):
            return np.nan
        donor_material = self._inventory_source_material(row)
        donor_row = dict(row) if not isinstance(row, pd.Series) else row.to_dict()
        donor_row["material_class"] = donor_material
        T_hot = float(_row_value(row, "temperature_K", T_REF_K) or T_REF_K)
        front_depth_um = self.dissolution_front_rate_um_y(donor_row) * max(time_y, 0.0)
        cr_mass_g_cm2 = front_depth_um * 1.0e-4 * density(donor_material) * cr_weight_fraction(donor_material)

        family = str(_row_value(row, "experiment_family", ""))
        source = str(_row_value(row, "source_id", ""))
        if "MSRE" in family or source.startswith("ORNL-TM-3"):
            area_to_mass = math.exp(self._param("log_area_to_salt_mass_msre_cm2_g"))
        else:
            area_to_mass = math.exp(self._param("log_area_to_salt_mass_loop_cm2_g"))

        retained = 1.0
        dT = max(float(_row_value(row, "delta_T_C", 0.0) or 0.0), 0.0)
        if dT > 0.0:
            cold_K = self._cold_temperature_K(row, T_hot)
            supersaturation = self._supersaturation_fraction(str(_row_value(row, "salt_class", "generic_salt")), T_hot, cold_K)
            retained = max(0.0, 1.0 - self._param("cold_leg_capture_efficiency") * supersaturation)

        raw_ppm = cr_mass_g_cm2 * area_to_mass * retained * 1.0e6
        capacity_ppm = math.exp(self._param("log_dissolved_cr_capacity_ppm"))
        return capacity_ppm * (1.0 - math.exp(-raw_ppm / max(capacity_ppm, 1.0e-30)))

    def salt_fe_decrease_ppm(self, row: pd.Series | Mapping[str, Any]) -> float:
        return self.salt_cr_ppm(row) * self._param("fe_to_cr_inventory_ratio")

    def redox_acceleration_ratio(self, row: pd.Series | Mapping[str, Any]) -> float:
        oxidized = self.corrosion_rate_um_y(row, redox_override="oxidizing_fef2")
        baseline = self.corrosion_rate_um_y(row, redox_override="purified_baseline")
        return oxidized / max(baseline, 1.0e-30)

    def igc_depth_um(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = float(_row_value(row, "time_years", np.nan) or np.nan)
        if not np.isfinite(time_y):
            return np.nan
        front_depth = self.dissolution_front_rate_um_y(row) * max(time_y, 0.0)
        diffusion_length = self.diffusion_length_um(row)
        redox = str(_row_value(row, "redox_class", "purified_baseline"))
        salt = str(_row_value(row, "salt_class", "generic_salt"))
        material = str(_row_value(row, "material_class", "generic_metal"))

        if redox == "multi_alloy":
            multiplier = self._param("multi_alloy_gb_multiplier")
            if material == "alloy_x750":
                multiplier *= self._param("x750_gb_multiplier")
            return front_depth + multiplier * diffusion_length
        if salt == "chloride":
            return (
                self._param("chloride_reaction_multiplier") * front_depth
                + self._param("chloride_gb_multiplier") * diffusion_length
            )
        if redox == "tellurium":
            return front_depth + self._param("tellurium_gb_multiplier") * diffusion_length
        return front_depth + self._param("gb_diffusion_multiplier") * diffusion_length

    def predict_response(self, row: pd.Series | Mapping[str, Any]) -> float:
        kind = str(_row_value(row, "response_kind", ""))
        if kind == "corrosion_rate_um_y":
            return self.corrosion_rate_um_y(row)
        if kind == "corrosion_depth_um":
            return self.corrosion_depth_um(row)
        if kind == "mass_loss_mg_cm2":
            return self.mass_loss_mg_cm2(row)
        if kind == "mass_gain_mg_cm2":
            return self.mass_gain_mg_cm2(row)
        if kind == "salt_cr_ppm":
            return self.salt_cr_ppm(row)
        if kind == "salt_fe_decrease_ppm":
            return self.salt_fe_decrease_ppm(row)
        if kind == "cr_diffusion_cm2_s":
            return self.cr_diffusion_cm2_s(row)
        if kind in {"redox_acceleration_ratio", "redox_acceleration_qualitative"}:
            return self.redox_acceleration_ratio(row)
        if kind == "igc_depth_um":
            return self.igc_depth_um(row)
        return np.nan

    def parameter_table(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "parameter": spec.name,
                    "value": self.params[spec.name],
                    "prior": spec.prior,
                    "prior_sigma": spec.prior_sigma,
                    "lower": spec.lower,
                    "upper": spec.upper,
                    "description": spec.description,
                }
                for spec in MECHANISTIC_PARAMETER_SPECS
            ]
        )


def mechanistic_prior_residuals(params: Mapping[str, float], strength: float = 0.25) -> np.ndarray:
    values = []
    for spec in MECHANISTIC_PARAMETER_SPECS:
        if spec.prior_sigma > 0.0:
            values.append(strength * (float(params[spec.name]) - spec.prior) / spec.prior_sigma)
    return np.asarray(values, dtype=float)


def mechanistic_residuals_for_targets(
    model: MechanisticCorrosionModel,
    targets: pd.DataFrame,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    residuals: list[float] = []
    details: list[dict[str, Any]] = []
    for _, row in targets.iterrows():
        if not mechanistic_is_supported(row):
            continue
        pred = model.predict_response(row)
        if not _is_positive_finite(pred):
            continue
        pred = float(pred)
        role = str(row.get("fit_role", ""))
        sigma = max(float(row.get("default_sigma_ln", 0.75) or 0.75), 0.15)
        quality = float(row.get("quality_weight", 1.0) or 1.0)
        if not np.isfinite(quality) or quality <= 0.0:
            quality = 0.5
        weight = math.sqrt(quality)
        target_low = row.get("target_low")
        target_mid = row.get("target_mid")
        target_high = row.get("target_high")
        res: float | None = None

        if role == "direct" and _is_positive_finite(target_mid):
            res = math.log(pred / float(target_mid)) / sigma * weight
        elif (
            role == "range"
            and _is_positive_finite(target_low)
            and _is_positive_finite(target_high)
        ):
            lower, upper = float(target_low), float(target_high)
            if lower <= pred <= upper:
                res = 0.0
            else:
                reference = lower if pred < lower else upper
                res = math.log(pred / reference) / sigma * weight
        elif role == "upper":
            upper = target_high if _is_positive_finite(target_high) else target_mid
            if _is_positive_finite(upper):
                res = max(0.0, math.log(pred / float(upper))) / sigma * weight
        elif role == "lower":
            lower = target_low if _is_positive_finite(target_low) else target_mid
            if _is_positive_finite(lower):
                res = max(0.0, math.log(float(lower) / pred)) / sigma * weight
        if res is None:
            continue
        residuals.append(float(res))
        details.append(
            {
                "measurement_id": row.get("measurement_id"),
                "case_id": row.get("case_id"),
                "response_kind": row.get("response_kind"),
                "fit_role": role,
                "prediction": pred,
                "target_low": target_low,
                "target_mid": target_mid,
                "target_high": target_high,
                "residual_ln_weighted": float(res),
            }
        )
    return np.asarray(residuals, dtype=float), details


def _objective_vector(
    vector: np.ndarray,
    base_model: MoltenSaltBVModel,
    targets: pd.DataFrame,
    prior_strength: float,
) -> np.ndarray:
    model = MechanisticCorrosionModel.from_vector(base_model, vector)
    data_residuals, _ = mechanistic_residuals_for_targets(model, targets)
    prior_residuals = mechanistic_prior_residuals(model.params, strength=prior_strength)
    return np.concatenate([data_residuals, prior_residuals])


def fit_mechanistic_model(
    targets: pd.DataFrame,
    base_model: MoltenSaltBVModel,
    *,
    prior_strength: float = 0.25,
    max_nfev: int = 5000,
) -> dict[str, Any]:
    result = least_squares(
        _objective_vector,
        mechanistic_initial_vector(),
        args=(base_model, targets, prior_strength),
        bounds=(mechanistic_lower_bounds(), mechanistic_upper_bounds()),
        loss="soft_l1",
        f_scale=1.0,
        x_scale="jac",
        max_nfev=max_nfev,
        verbose=0,
    )
    model = MechanisticCorrosionModel.from_vector(base_model, result.x)
    residuals, details = mechanistic_residuals_for_targets(model, targets)
    return {
        "model": model,
        "optimizer_result": result,
        "residuals": residuals,
        "residual_details": pd.DataFrame(details),
        "parameter_table": model.parameter_table(),
        "prior_strength": float(prior_strength),
    }


def _constraint_pass(row: pd.Series, prediction: float) -> bool:
    if not _is_positive_finite(prediction):
        return False
    role = str(row.get("fit_role", ""))
    low = row.get("target_low")
    mid = row.get("target_mid")
    high = row.get("target_high")
    pred = float(prediction)
    if role == "direct" and _is_positive_finite(mid):
        return abs(math.log(pred / float(mid))) <= math.log(2.0)
    if role == "range" and _is_positive_finite(low) and _is_positive_finite(high):
        return float(low) <= pred <= float(high)
    if role == "upper":
        limit = high if _is_positive_finite(high) else mid
        return bool(_is_positive_finite(limit) and pred <= float(limit))
    if role == "lower":
        limit = low if _is_positive_finite(low) else mid
        return bool(_is_positive_finite(limit) and pred >= float(limit))
    return False


def compare_mechanistic_to_effective(
    targets: pd.DataFrame,
    mechanistic_model: MechanisticCorrosionModel,
    effective_model: MoltenSaltBVModel,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for _, row in targets.iterrows():
        reason = mechanistic_scope_reason(row)
        if reason != "supported":
            continue
        mech = float(mechanistic_model.predict_response(row))
        effective = float(effective_model.predict_response(row))
        item = row.to_dict()
        item["mechanistic_scope_reason"] = reason
        item["effective_prediction"] = effective
        item["mechanistic_prediction"] = mech
        item["effective_constraint_pass"] = _constraint_pass(row, effective)
        item["mechanistic_constraint_pass"] = _constraint_pass(row, mech)
        mid = row.get("target_mid")
        if _is_positive_finite(mid):
            target = float(mid)
            item["effective_factor_error"] = max(effective / target, target / effective)
            item["mechanistic_factor_error"] = max(mech / target, target / mech)
        else:
            item["effective_factor_error"] = np.nan
            item["mechanistic_factor_error"] = np.nan
        rows.append(item)
    comparison = pd.DataFrame(rows)

    metrics: dict[str, Any] = {
        "n_mechanistic_constraints": int(len(comparison)),
        "n_mechanistic_parameters": int(len(MECHANISTIC_PARAMETER_SPECS)),
        "shared_effective_inputs": [
            "Ea_corr_kJ_mol",
            "redox class log-driving-force offsets",
            "log_Dcr_ref_cm2_s",
            "Ea_Dcr_kJ_mol",
        ],
    }
    direct = comparison[comparison["fit_role"].isin(MECHANISTIC_DIRECT_ROLES)].copy()
    direct = direct[
        pd.to_numeric(direct["target_mid"], errors="coerce").gt(0)
        & pd.to_numeric(direct["effective_prediction"], errors="coerce").gt(0)
        & pd.to_numeric(direct["mechanistic_prediction"], errors="coerce").gt(0)
    ]
    metrics["n_direct_or_range_targets"] = int(len(direct))
    for label, column in (("effective", "effective_prediction"), ("mechanistic", "mechanistic_prediction")):
        if len(direct):
            ln_error = np.log(direct[column].astype(float).to_numpy() / direct["target_mid"].astype(float).to_numpy())
            metrics[f"{label}_log_rmse_direct"] = float(np.sqrt(np.mean(ln_error**2)))
            metrics[f"{label}_median_factor_error_direct"] = float(np.exp(np.median(np.abs(ln_error))))
            metrics[f"{label}_within_factor_2_direct"] = float(np.mean(np.abs(ln_error) <= math.log(2.0)))
            metrics[f"{label}_within_factor_5_direct"] = float(np.mean(np.abs(ln_error) <= math.log(5.0)))
        else:
            metrics[f"{label}_log_rmse_direct"] = None
            metrics[f"{label}_median_factor_error_direct"] = None
            metrics[f"{label}_within_factor_2_direct"] = None
            metrics[f"{label}_within_factor_5_direct"] = None
    metrics["effective_constraint_pass_fraction"] = float(comparison["effective_constraint_pass"].mean()) if len(comparison) else None
    metrics["mechanistic_constraint_pass_fraction"] = float(comparison["mechanistic_constraint_pass"].mean()) if len(comparison) else None
    if len(direct):
        metrics["mechanistic_better_factor_error_count"] = int(
            (direct["mechanistic_factor_error"] < direct["effective_factor_error"]).sum()
        )
    else:
        metrics["mechanistic_better_factor_error_count"] = 0
    if metrics.get("effective_log_rmse_direct") and metrics.get("mechanistic_log_rmse_direct") is not None:
        metrics["log_rmse_reduction_fraction"] = float(
            1.0 - metrics["mechanistic_log_rmse_direct"] / metrics["effective_log_rmse_direct"]
        )
    return comparison, metrics


def save_mechanistic_outputs(outputs: Mapping[str, Any], output_dir: str | Path) -> None:
    output_dir = Path(output_dir)
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    outputs["comparison"].to_csv(results_dir / "mechanistic_validation_comparison.csv", index=False)
    outputs["fit"]["parameter_table"].to_csv(results_dir / "mechanistic_parameters.csv", index=False)
    with open(results_dir / "mechanistic_parameters.json", "w", encoding="utf-8") as handle:
        json.dump(outputs["fit"]["model"].params, handle, indent=2)
    with open(results_dir / "mechanistic_comparison_metrics.json", "w", encoding="utf-8") as handle:
        json.dump(outputs["metrics"], handle, indent=2)
