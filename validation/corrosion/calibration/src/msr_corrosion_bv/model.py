"""Effective Butler-Volmer/transport model for molten-salt corrosion and plating.

The sparse validation workbook generally lacks measured electrochemical potentials,
activities, velocities, geometries, and loop volumes. This model therefore uses an
"effective Butler-Volmer" form: unknown redox potentials and activity terms enter as
calibrated overpotential offsets, while temperature, redox, salt chemistry, alloy
class, flow/circulation, and thermal-gradient factors remain explicit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable, Mapping
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .chemistry import (
    FARADAY,
    R_GAS,
    T_REF_K,
    cr_weight_fraction,
    density,
    safe_log,
    um_to_mg_cm2,
)


def exp_clip(x: float | np.ndarray, lo: float = -60.0, hi: float = 60.0) -> float | np.ndarray:
    return np.exp(np.clip(x, lo, hi))

@lru_cache(maxsize=1)
def solid_diffusivities() -> dict[str, Any]:
    """Load source-specific solid diffusivities from the production database."""
    db_path = Path(__file__).resolve().parents[5] / "data" / "corrosion_database.json"
    with db_path.open(encoding="utf-8") as stream:
        database = json.load(stream)
    return database["solid_diffusivities"]

@lru_cache(maxsize=1)
def solid_diffusivity_fallbacks() -> dict[str, str]:
    """Load explicit solid-diffusivity material fallbacks from the production database."""
    db_path = Path(__file__).resolve().parents[5] / "data" / "corrosion_database.json"
    with db_path.open(encoding="utf-8") as stream:
        database = json.load(stream)
    return database.get("solid_diffusivity_fallbacks", {})

def _case_insensitive_get(mapping: Mapping[str, Any], key: str) -> Any:
    """Return a mapping value using case-insensitive string-key lookup."""
    key_lower = key.lower()
    for candidate, value in mapping.items():
        if candidate.lower() == key_lower:
            return value
    return None

@dataclass(frozen=True)
class ParameterSpec:
    """One fitted coefficient and the bounds/prior that define its calibration."""
    name: str
    initial: float
    lower: float
    upper: float
    prior: float
    prior_sigma: float
    description: str


# Parameters are expressed in model space. Activation energies are in kJ/mol.
# Every entry has a finite prior width and therefore contributes a residual to
# the objective. This regularization is necessary because 59 coefficients are
# inferred from only 37 active calibration rows.
PARAMETER_SPECS: tuple[ParameterSpec, ...] = (
    ParameterSpec("log_rate0_um_y", math.log(4.0), math.log(0.02), math.log(200.0), math.log(4.0), 1.5, "Reference anodic dissolution rate for Hastelloy N in fuel fluoride at 650 C, um/y."),
    ParameterSpec("Ea_corr_kJ_mol", 60.0, 20.0, 180.0, 60.0, 30.0, "Apparent activation energy for anodic dissolution."),
    ParameterSpec("gamma_flow_corr", 0.35, 0.15, 2.0, 0.35, 0.45, "Power-law exponent for circulation/mass-transfer effect on corrosion."),
    ParameterSpec("theta_dT_corr", 0.25, 0.05, 3.0, 0.25, 0.60, "Thermal-gradient enhancement of corrosion/mass transfer."),
    ParameterSpec("hot_leg_bonus", 0.15, 0.0, 1.8, 0.15, 0.55, "Additional log-rate offset for hot-leg/hottest specimens."),
    ParameterSpec("cold_leg_corrosion_penalty", -0.25, -2.5, 1.0, -0.25, 0.70, "Log-rate offset for cold-leg corrosion before deposition."),
    ParameterSpec("log_mt_cap_um_y", math.log(1000.0), math.log(10.0), math.log(1.0e5), math.log(1000.0), 1.8, "Reference transport-capacity rate; kinetic and transport rates are combined harmonically."),
    ParameterSpec("Ea_mt_kJ_mol", 25.0, 0.0, 120.0, 25.0, 35.0, "Apparent activation energy for transport-capacity branch."),
    ParameterSpec("gamma_flow_mt", 0.75, -0.5, 2.5, 0.75, 0.70, "Power-law flow exponent for transport capacity."),
    ParameterSpec("theta_dT_mt", 0.60, -1.0, 3.0, 0.60, 0.90, "Thermal-gradient enhancement of transport capacity."),

    # Material offsets relative to Hastelloy N.
    ParameterSpec("mat_modified_hastelloy_n", -0.20, -2.5, 2.0, -0.20, 0.80, "Modified Hastelloy N log-rate offset."),
    ParameterSpec("mat_gh3535", 0.05, -2.0, 2.0, 0.05, 0.80, "GH3535 log-rate offset."),
    ParameterSpec("mat_stainless_304", 1.20, -1.5, 4.0, 1.20, 0.85, "304/304L stainless log-rate offset."),
    ParameterSpec("mat_stainless_316", 1.05, -1.5, 4.0, 1.05, 0.85, "316/316H/316L stainless log-rate offset."),
    ParameterSpec("mat_alloy_x750", 2.40, -1.0, 7.0, 2.40, 1.20, "Alloy X750 log-rate offset."),
    ParameterSpec("mat_in625", 0.70, -2.0, 3.5, 0.70, 1.00, "IN625 log-rate offset."),
    ParameterSpec("mat_ni_alloy", 0.40, -2.0, 3.5, 0.40, 1.00, "Generic high-Ni alloy log-rate offset."),
    ParameterSpec("mat_generic_metal", 0.70, -2.0, 3.5, 0.70, 1.20, "Generic metal log-rate offset."),

    # Salt offsets relative to fuel fluoride.
    ParameterSpec("salt_flinak", 0.70, -2.0, 4.0, 0.70, 0.90, "FLiNaK log-rate offset."),
    ParameterSpec("salt_flibe", -0.35, -3.0, 3.0, -0.35, 0.90, "FLiBe log-rate offset."),
    ParameterSpec("salt_fluoroborate", 0.95, -2.0, 4.5, 0.95, 1.00, "Fluoroborate coolant log-rate offset."),
    ParameterSpec("salt_chloride", 1.65, -1.0, 5.5, 1.65, 1.05, "Chloride salt log-rate offset."),
    ParameterSpec("salt_generic_salt", 0.40, -2.0, 4.0, 0.40, 1.25, "Generic/unspecified salt log-rate offset."),
    ParameterSpec("salt_no_salt", -8.0, -15.0, -2.0, -8.0, 1.50, "No-salt control offset; keeps molten-salt corrosion near zero."),

    # Redox/effective overpotential offsets. Additive in log-rate is equivalent
    # to alpha*n*F*eta/(RT) in the anodic branch of Butler-Volmer.
    ParameterSpec("redox_purified_baseline", 0.0, -1.0, 1.0, 0.0, 0.35, "Purified baseline redox offset."),
    ParameterSpec("redox_msre_or_fuel_baseline", -0.35, -3.0, 2.0, -0.35, 0.80, "MSRE/fuel baseline redox offset."),
    ParameterSpec("redox_oxidizing_fef2", 1.80, -0.5, 5.0, 1.80, 0.80, "FeF2 oxidizing perturbation offset."),
    ParameterSpec("redox_reducing_be", -2.00, -6.0, 0.5, -2.00, 0.90, "Be-reduced salt offset."),
    ParameterSpec("redox_impure_moisture", 1.75, -0.5, 5.0, 1.75, 0.90, "Moisture/impure salt oxidizing offset."),
    ParameterSpec("redox_chloride_unspecified", 0.65, -1.5, 3.5, 0.65, 1.00, "Unspecified chloride redox/impurity offset."),
    ParameterSpec("redox_tellurium", 0.40, -2.0, 4.0, 0.40, 1.10, "Tellurium-bearing salt offset."),
    ParameterSpec("redox_stressed", 1.40, -0.5, 5.0, 1.40, 1.00, "Stress-assisted corrosion offset."),
    ParameterSpec("redox_multi_alloy", 1.25, -0.5, 6.0, 1.25, 1.10, "Multi-alloy/plating/galvanic interaction offset."),
    ParameterSpec("redox_fission_product", -0.20, -3.0, 3.0, -0.20, 1.00, "Fission-product/noble-metal environment offset."),
    ParameterSpec("redox_gas_control", -6.0, -12.0, -1.0, -6.0, 1.40, "Gas control offset."),

    # Damage/dealloying morphology factors.
    ParameterSpec("damage_base_log", 0.10, -1.5, 2.5, 0.10, 0.70, "Base multiplier for void/IGC damage depth beyond uniform recession."),
    ParameterSpec("damage_chloride_log", 0.70, -1.0, 4.0, 0.70, 1.00, "Additional chloride IGC/void damage multiplier."),
    ParameterSpec("damage_tellurium_log", 2.00, 0.0, 5.5, 2.00, 1.00, "Tellurium IGC multiplier."),
    ParameterSpec("damage_stress_log", 1.10, -0.5, 4.5, 1.10, 1.00, "Stress-assisted cracking multiplier."),
    ParameterSpec("damage_multi_alloy_log", 1.00, -0.5, 4.5, 1.00, 1.00, "Multi-alloy galvanic/plating damage multiplier."),

    # Deposition/plating branch.
    ParameterSpec("log_dep0_um_y", math.log(5.0), math.log(0.01), math.log(400.0), math.log(5.0), 1.7, "Reference deposition/plating rate, um/y."),
    ParameterSpec("Ea_dep_kJ_mol", 20.0, -40.0, 120.0, 20.0, 45.0, "Apparent temperature coefficient for deposition branch."),
    ParameterSpec("gamma_flow_dep", 0.60, 0.15, 2.5, 0.60, 0.70, "Flow exponent for plating/deposition."),
    ParameterSpec("theta_dT_dep", 0.95, 0.0, 4.0, 0.95, 0.85, "Thermal-gradient/cold-leg deposition enhancement."),
    ParameterSpec("dep_cold_bonus", 0.80, -1.0, 4.0, 0.80, 0.90, "Cold-leg deposition bonus."),
    ParameterSpec("dep_graphite_offset", -0.60, -3.5, 1.5, -0.60, 0.90, "Graphite surface deposition offset relative to metal."),
    ParameterSpec("dep_turbulent_bonus", 0.45, -1.0, 2.5, 0.45, 0.70, "Turbulent metal deposition offset."),
    ParameterSpec("dep_laminar_bonus", 0.00, -1.5, 1.5, 0.00, 0.60, "Laminar metal deposition offset."),
    ParameterSpec("dep_multi_alloy_bonus", 0.90, -1.0, 4.0, 0.90, 1.00, "Multi-alloy donor/plating source enhancement."),
    ParameterSpec("dep_salt_flinak", 0.40, -3.0, 4.0, 0.40, 1.10, "Deposition-only FLiNaK offset relative to fuel fluoride."),
    ParameterSpec("dep_salt_flibe", -0.20, -3.0, 4.0, -0.20, 1.10, "Deposition-only FLiBe offset relative to fuel fluoride."),
    ParameterSpec("dep_salt_fluoroborate", 0.20, -3.0, 4.0, 0.20, 1.10, "Deposition-only fluoroborate offset relative to fuel fluoride."),
    ParameterSpec("dep_salt_chloride", 0.30, -3.0, 4.0, 0.30, 1.20, "Deposition-only chloride offset relative to fuel fluoride."),

    # Inventory, diffusion, and fission-product submodels.
    ParameterSpec("log_ppm_scale_msre", math.log(9.0), math.log(0.01), math.log(500.0), math.log(9.0), 1.5, "Cr ppm per effective um source depth for MSRE-scale inventory closure."),
    ParameterSpec("log_ppm_scale_loop", math.log(120.0), math.log(0.1), math.log(5000.0), math.log(120.0), 1.8, "Cr ppm per effective um source depth for small loop inventory closure."),
    ParameterSpec("log_fe_to_cr_ppm_ratio", math.log(0.25), math.log(0.01), math.log(2.0), math.log(0.25), 1.2, "Magnitude ratio of Fe decrease to Cr increase in NCL inventory closure."),
    ParameterSpec("logit_offgas_fraction", math.log(0.001 / (1.0 - 0.001)), math.log(1.0e-6 / (1.0 - 1.0e-6)), math.log(0.20 / 0.80), math.log(0.001 / (1.0 - 0.001)), 2.0, "Noble-metal off-gas escape fraction as a logit."),
    ParameterSpec("log_te_soluble_ppm", math.log(2.0), math.log(0.01), math.log(100.0), math.log(2.0), 1.5, "Effective soluble tellurium upper-bound concentration."),
    ParameterSpec("log_te_threshold_ratio", math.log(150.0), math.log(5.0), math.log(2000.0), math.log(150.0), 0.90, "Critical U(IV)/U(III) ratio for telluride formation."),
)

PARAMETER_NAMES = tuple(p.name for p in PARAMETER_SPECS)


def initial_parameter_vector() -> np.ndarray:
    return np.array([p.initial for p in PARAMETER_SPECS], dtype=float)


def lower_bounds() -> np.ndarray:
    return np.array([p.lower for p in PARAMETER_SPECS], dtype=float)


def upper_bounds() -> np.ndarray:
    return np.array([p.upper for p in PARAMETER_SPECS], dtype=float)


def vector_to_params(vector: Iterable[float]) -> dict[str, float]:
    return {name: float(v) for name, v in zip(PARAMETER_NAMES, vector)}


def params_to_vector(params: Mapping[str, float]) -> np.ndarray:
    return np.array([params[name] for name in PARAMETER_NAMES], dtype=float)


def priors_as_residuals(params: Mapping[str, float], scale: float = 1.0) -> np.ndarray:
    values = []
    for spec in PARAMETER_SPECS:
        if spec.prior_sigma > 0:
            values.append(scale * (params[spec.name] - spec.prior) / spec.prior_sigma)
    return np.asarray(values, dtype=float)


class MoltenSaltBVModel:
    """Effective source-term correlation with corrosion and plating branches.

    Class offsets absorb unresolved electrochemical potentials, activities,
    geometry, and transport details.  Predictions are therefore conditional on
    the workbook's categorical mapping and should not be presented as a
    first-principles electrode model.
    """

    def __init__(self, params: Mapping[str, float] | None = None) -> None:
        self.params = vector_to_params(initial_parameter_vector()) if params is None else dict(params)

    @classmethod
    def from_vector(cls, vector: Iterable[float]) -> "MoltenSaltBVModel":
        return cls(vector_to_params(vector))

    def to_vector(self) -> np.ndarray:
        return params_to_vector(self.params)

    def update(self, params: Mapping[str, float]) -> None:
        self.params.update(params)

    def _param(self, name: str) -> float:
        return float(self.params[name])

    @staticmethod
    def _feature(row: pd.Series | Mapping[str, Any], name: str, default: Any = None) -> Any:
        if isinstance(row, pd.Series):
            value = row.get(name, default)
        else:
            value = row.get(name, default)
        if isinstance(value, float) and math.isnan(value):
            return default
        return default if value is None else value

    def material_offset(self, material_class: str) -> float:
        mapping = {
            "modified_hastelloy_n": "mat_modified_hastelloy_n",
            "gh3535": "mat_gh3535",
            "stainless_304": "mat_stainless_304",
            "stainless_304l": "mat_stainless_304",
            "stainless_316": "mat_stainless_316",
            "stainless_316h": "mat_stainless_316",
            "stainless_316l": "mat_stainless_316",
            "alloy_x750": "mat_alloy_x750",
            "in625": "mat_in625",
            "ni_alloy": "mat_ni_alloy",
            "generic_metal": "mat_generic_metal",
            "graphite": "mat_generic_metal",
        }
        key = mapping.get(material_class, "mat_generic_metal")
        return 0.0 if material_class == "hastelloy_n" else self._param(key)

    def salt_offset(self, salt_class: str) -> float:
        mapping = {
            "flinak": "salt_flinak",
            "flibe": "salt_flibe",
            "fluoroborate": "salt_fluoroborate",
            "chloride": "salt_chloride",
            "generic_salt": "salt_generic_salt",
            "no_salt": "salt_no_salt",
        }
        key = mapping.get(salt_class, None)
        return 0.0 if key is None else self._param(key)

    def redox_offset(self, redox_class: str) -> float:
        mapping = {
            "purified_baseline": "redox_purified_baseline",
            "msre_or_fuel_baseline": "redox_msre_or_fuel_baseline",
            "oxidizing_fef2": "redox_oxidizing_fef2",
            "reducing_be": "redox_reducing_be",
            "impure_moisture": "redox_impure_moisture",
            "chloride_unspecified": "redox_chloride_unspecified",
            "tellurium": "redox_tellurium",
            "stressed": "redox_stressed",
            "multi_alloy": "redox_multi_alloy",
            "fission_product": "redox_fission_product",
            "gas_control": "redox_gas_control",
        }
        key = mapping.get(redox_class, "redox_purified_baseline")
        return self._param(key)

    def position_offset(self, position_class: str) -> float:
        if position_class == "hot_leg":
            return self._param("hot_leg_bonus")
        if position_class == "cold_leg":
            return self._param("cold_leg_corrosion_penalty")
        return 0.0

    def thermal_term(self, temperature_K: float, Ea_kJ_mol: float) -> float:
        if not np.isfinite(temperature_K) or temperature_K <= 0:
            temperature_K = T_REF_K
        return (Ea_kJ_mol * 1000.0 / R_GAS) * (1.0 / T_REF_K - 1.0 / temperature_K)

    def bv_overpotential_equivalent_V(self, redox_class: str, temperature_K: float, alpha_n: float = 1.0) -> float:
        """Return the effective overpotential implied by the fitted redox offset.

        Since the corrosion workbook lacks measured potentials, the redox offset is
        the identifiable quantity. This helper converts that log-current offset to
        an equivalent eta using offset = alpha*n*F*eta/(R*T).
        """
        return self.redox_offset(redox_class) * R_GAS * temperature_K / (alpha_n * FARADAY)

    def corrosion_rate_um_y(self, row: pd.Series | Mapping[str, Any], redox_override: str | None = None) -> float:
        p = self.params
        T = float(self._feature(row, "temperature_K", T_REF_K) or T_REF_K)
        material = str(self._feature(row, "material_class", "generic_metal"))
        salt = str(self._feature(row, "salt_class", "generic_salt"))
        redox = redox_override or str(self._feature(row, "redox_class", "purified_baseline"))
        position = str(self._feature(row, "position_class", "nominal"))
        flow = max(float(self._feature(row, "flow_factor", 0.75) or 0.75), 1.0e-3)
        dT = max(float(self._feature(row, "delta_T_C", 0.0) or 0.0), 0.0)

        if salt == "no_salt" or flow < 0.05:
            # Keep finite for log plots while indicating negligible molten-salt corrosion.
            return 1.0e-6

        log_kin = (
            p["log_rate0_um_y"]
            + self.thermal_term(T, p["Ea_corr_kJ_mol"])
            + self.material_offset(material)
            + self.salt_offset(salt)
            + self.redox_offset(redox)
            + self.position_offset(position)
            + p["gamma_flow_corr"] * math.log(flow)
            + p["theta_dT_corr"] * math.log1p(dT / 100.0)
        )
        rate_kin = float(exp_clip(log_kin))

        log_mt = (
            p["log_mt_cap_um_y"]
            + self.thermal_term(T, p["Ea_mt_kJ_mol"])
            + 0.35 * self.salt_offset(salt)
            + p["gamma_flow_mt"] * math.log(flow)
            + p["theta_dT_mt"] * math.log1p(dT / 100.0)
        )
        rate_mt = float(exp_clip(log_mt))
        # Harmonic mean: kinetic-limited at low current, transport-limited at high driving force.
        return 1.0 / (1.0 / max(rate_kin, 1.0e-12) + 1.0 / max(rate_mt, 1.0e-12))

    def damage_multiplier(self, row: pd.Series | Mapping[str, Any]) -> float:
        redox = str(self._feature(row, "redox_class", "purified_baseline"))
        salt = str(self._feature(row, "salt_class", "generic_salt"))
        value = self._param("damage_base_log")
        if salt == "chloride" or redox in {"impure_moisture", "chloride_unspecified"}:
            value += self._param("damage_chloride_log")
        if redox == "tellurium":
            value += self._param("damage_tellurium_log")
        if redox == "stressed":
            value += self._param("damage_stress_log")
        if redox == "multi_alloy":
            value += self._param("damage_multi_alloy_log")
        return float(exp_clip(value, -10.0, 10.0))

    def corrosion_depth_um(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = float(self._feature(row, "time_years", np.nan) or np.nan)
        if not np.isfinite(time_y):
            return np.nan
        return self.corrosion_rate_um_y(row) * max(time_y, 0.0)

    def igc_depth_um(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = float(self._feature(row, "time_years", np.nan) or np.nan)
        if not np.isfinite(time_y):
            return np.nan
        uniform = self.corrosion_rate_um_y(row) * max(time_y, 0.0)
        # Mixed linear/parabolic morphology. The square-root term allows deep IGC/void
        # penetration with limited mass loss while remaining tied to BV corrosion drive.
        return uniform * self.damage_multiplier(row) + math.sqrt(max(uniform, 0.0)) * self.damage_multiplier(row)

    def mass_loss_mg_cm2(self, row: pd.Series | Mapping[str, Any]) -> float:
        material = str(self._feature(row, "material_class", "generic_metal"))
        return um_to_mg_cm2(self.corrosion_depth_um(row), material)

    def deposition_salt_offset(self, salt_class: str) -> float:
        mapping = {
            "flinak": "dep_salt_flinak",
            "flibe": "dep_salt_flibe",
            "fluoroborate": "dep_salt_fluoroborate",
            "chloride": "dep_salt_chloride",
        }
        key = mapping.get(salt_class)
        return 0.0 if key is None else self._param(key)

    def deposition_rate_um_y(self, row: pd.Series | Mapping[str, Any], surface_override: str | None = None) -> float:
        p = self.params
        T = float(self._feature(row, "temperature_K", T_REF_K) or T_REF_K)
        salt = str(self._feature(row, "salt_class", "generic_salt"))
        redox = str(self._feature(row, "redox_class", "purified_baseline"))
        position = str(self._feature(row, "position_class", "nominal"))
        surface = surface_override or str(self._feature(row, "surface_class", "metal"))
        flow = max(float(self._feature(row, "flow_factor", 0.75) or 0.75), 1.0e-3)
        dT = max(float(self._feature(row, "delta_T_C", 0.0) or 0.0), 0.0)
        if salt == "no_salt":
            return 1.0e-9

        surface_offset = 0.0
        if surface == "graphite":
            surface_offset += p["dep_graphite_offset"]
        elif surface == "turbulent_metal":
            surface_offset += p["dep_turbulent_bonus"]
        elif surface == "laminar_metal":
            surface_offset += p["dep_laminar_bonus"]
        if redox == "multi_alloy":
            surface_offset += p["dep_multi_alloy_bonus"]
        # Cathodic/plating branch: redox drive enters with a reduced sign because
        # noble metals and corrosion products deposit as cathodic/precipitation processes.
        redox_cathodic = -0.30 * self.redox_offset(redox)
        cold_bonus = p["dep_cold_bonus"] if position == "cold_leg" else 0.0
        log_dep = (
            p["log_dep0_um_y"]
            + self.thermal_term(T, p["Ea_dep_kJ_mol"])
            + 0.20 * self.salt_offset(salt)
            + self.deposition_salt_offset(salt)
            + redox_cathodic
            + p["gamma_flow_dep"] * math.log(flow)
            + p["theta_dT_dep"] * math.log1p(dT / 100.0)
            + cold_bonus
            + surface_offset
        )
        return float(exp_clip(log_dep))

    def deposition_depth_um(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = float(self._feature(row, "time_years", np.nan) or np.nan)
        if not np.isfinite(time_y):
            return np.nan
        return self.deposition_rate_um_y(row) * max(time_y, 0.0)

    def mass_gain_mg_cm2(self, row: pd.Series | Mapping[str, Any]) -> float:
        material = str(self._feature(row, "material_class", "generic_metal"))
        return um_to_mg_cm2(self.deposition_depth_um(row), material)

    def redox_acceleration_ratio(self, row: pd.Series | Mapping[str, Any]) -> float:
        oxidized = self.corrosion_rate_um_y(row, redox_override="oxidizing_fef2")
        baseline = self.corrosion_rate_um_y(row, redox_override="purified_baseline")
        return oxidized / max(baseline, 1.0e-12)

    def cr_diffusion_cm2_s(self, row: pd.Series | Mapping[str, Any]) -> float:
        T = float(self._feature(row, "temperature_K", T_REF_K) or T_REF_K)
        material = str(self._feature(row, "material_class", "generic_metal"))

        diffusivities = solid_diffusivities()
        fallbacks = solid_diffusivity_fallbacks()

        material_data = _case_insensitive_get(diffusivities, material)

        if material_data is None or "Cr" not in material_data:
            fallback_material = _case_insensitive_get(fallbacks, material)
            if fallback_material is not None:
                material_data = _case_insensitive_get(diffusivities, fallback_material)

        if material_data is None or "Cr" not in material_data:
            material_data = diffusivities["generic_metal"]

        props = material_data["Cr"]
        D0_cm2_s = float(props["D0_cm2_s"])
        Q_kJ_mol = float(props["Q_kJ_mol"])

        logD = math.log(D0_cm2_s) - Q_kJ_mol * 1000.0 / (R_GAS * T)
        return float(exp_clip(logD, -80.0, -20.0))

    def salt_cr_ppm(self, row: pd.Series | Mapping[str, Any]) -> float:
        family = str(self._feature(row, "experiment_family", ""))
        source = str(self._feature(row, "source_id", ""))
        depth = max(self.corrosion_depth_um(row), 1.0e-9)
        material = str(self._feature(row, "material_class", "generic_metal"))
        cr_factor = cr_weight_fraction(material) / 0.07
        if "MSRE" in family or source.startswith("ORNL-TM-3"):
            scale = float(exp_clip(self._param("log_ppm_scale_msre")))
        else:
            scale = float(exp_clip(self._param("log_ppm_scale_loop")))
        return scale * depth * cr_factor

    def salt_fe_decrease_ppm(self, row: pd.Series | Mapping[str, Any]) -> float:
        return self.salt_cr_ppm(row) * float(exp_clip(self._param("log_fe_to_cr_ppm_ratio")))

    def offgas_fraction_percent(self, row: pd.Series | Mapping[str, Any] | None = None) -> float:
        logit = self._param("logit_offgas_fraction")
        frac = 1.0 / (1.0 + math.exp(-logit))
        return 100.0 * frac

    def te_soluble_ppm(self, row: pd.Series | Mapping[str, Any] | None = None) -> float:
        return float(exp_clip(self._param("log_te_soluble_ppm")))

    def te_redox_threshold_ratio(self, row: pd.Series | Mapping[str, Any] | None = None) -> float:
        return float(exp_clip(self._param("log_te_threshold_ratio")))

    def predict_response(self, row: pd.Series | Mapping[str, Any]) -> float:
        kind = str(self._feature(row, "response_kind", "corrosion_rate_um_y"))
        if kind == "corrosion_rate_um_y":
            return self.corrosion_rate_um_y(row)
        if kind == "corrosion_depth_um":
            return self.corrosion_depth_um(row)
        if kind == "igc_depth_um":
            return self.igc_depth_um(row)
        if kind == "mass_loss_mg_cm2":
            return self.mass_loss_mg_cm2(row)
        if kind == "mass_gain_mg_cm2":
            return self.mass_gain_mg_cm2(row)
        if kind == "redox_acceleration_ratio" or kind == "redox_acceleration_qualitative":
            return self.redox_acceleration_ratio(row)
        if kind == "salt_cr_ppm":
            return self.salt_cr_ppm(row)
        if kind == "salt_fe_decrease_ppm":
            return self.salt_fe_decrease_ppm(row)
        if kind == "cr_diffusion_cm2_s":
            return self.cr_diffusion_cm2_s(row)
        if kind == "offgas_fraction_percent":
            return self.offgas_fraction_percent(row)
        if kind == "te_soluble_ppm":
            return self.te_soluble_ppm(row)
        if kind == "te_redox_threshold_ratio":
            return self.te_redox_threshold_ratio(row)
        if kind == "noble_metal_deposition_ranking":
            # Return turbulent/graphite ratio as a scalar diagnostic.
            turbulent = self.deposition_rate_um_y(row, surface_override="turbulent_metal")
            graphite = self.deposition_rate_um_y(row, surface_override="graphite")
            return turbulent / max(graphite, 1.0e-12)
        return np.nan

    def predict_targets(self, targets: pd.DataFrame) -> pd.DataFrame:
        out = targets.copy()
        preds = []
        etas = []
        rates = []
        dep_rates = []
        for _, row in out.iterrows():
            preds.append(self.predict_response(row))
            T = row.get("temperature_K", T_REF_K)
            if not np.isfinite(T):
                T = T_REF_K
            etas.append(self.bv_overpotential_equivalent_V(str(row.get("redox_class", "purified_baseline")), float(T)))
            rates.append(self.corrosion_rate_um_y(row))
            dep_rates.append(self.deposition_rate_um_y(row))
        out["prediction"] = preds
        out["pred_corrosion_rate_um_y"] = rates
        out["pred_deposition_rate_um_y"] = dep_rates
        out["effective_eta_V"] = etas
        out["prediction_units"] = out["target_units_model"]
        return out

    def predict_cases(self, cases: pd.DataFrame) -> pd.DataFrame:
        out = cases.copy()
        rates = []
        depths = []
        dep_rates = []
        dep_depths = []
        ppms = []
        etas = []
        for _, row in out.iterrows():
            rates.append(self.corrosion_rate_um_y(row))
            depths.append(self.corrosion_depth_um(row))
            dep_rates.append(self.deposition_rate_um_y(row))
            dep_depths.append(self.deposition_depth_um(row))
            ppms.append(self.salt_cr_ppm(row) if row.get("salt_class") != "no_salt" and np.isfinite(row.get("time_years", np.nan)) else np.nan)
            T = row.get("temperature_K", T_REF_K)
            if not np.isfinite(T):
                T = T_REF_K
            etas.append(self.bv_overpotential_equivalent_V(str(row.get("redox_class", "purified_baseline")), float(T)))
        out["pred_corrosion_rate_um_y"] = rates
        out["pred_corrosion_depth_um"] = depths
        out["pred_deposition_rate_um_y"] = dep_rates
        out["pred_deposition_depth_um"] = dep_depths
        out["pred_salt_cr_ppm"] = ppms
        out["effective_eta_V"] = etas
        return out

    def deposition_ranking(self, base_row: pd.Series | Mapping[str, Any]) -> dict[str, float]:
        return {
            "graphite": self.deposition_rate_um_y(base_row, "graphite"),
            "laminar_metal": self.deposition_rate_um_y(base_row, "laminar_metal"),
            "metal": self.deposition_rate_um_y(base_row, "metal"),
            "turbulent_metal": self.deposition_rate_um_y(base_row, "turbulent_metal"),
        }

    def simulate_loop(
        self,
        segments: pd.DataFrame,
        duration_h: float,
        dt_h: float = 24.0,
        salt_volume_cm3: float = 1000.0,
        initial_cr_ppm: float = 0.0,
    ) -> pd.DataFrame:
        """Simulate a simplified multi-segment loop.

        Parameters
        ----------
        segments:
            DataFrame with at least material_class, salt_class, temperature_C or
            temperature_K, flow_factor, delta_T_C, redox_class, surface_area_cm2,
            and position_class. Missing geometry defaults to unit area.
        duration_h:
            Total exposure time.
        dt_h:
            Time step.
        salt_volume_cm3:
            Effective well-mixed salt volume. Used for normalized inventory closure.
        initial_cr_ppm:
            Initial chromium concentration.

        Returns
        -------
        DataFrame with time, salt Cr ppm, cumulative dissolution and deposition
        by segment. This is a compact engineering ODE intended for validation and
        sensitivity studies, not CFD.
        """
        if duration_h <= 0:
            raise ValueError("duration_h must be positive")
        if dt_h <= 0:
            raise ValueError("dt_h must be positive")
        segs = segments.copy().reset_index(drop=True)
        if "temperature_K" not in segs:
            segs["temperature_K"] = segs.get("temperature_C", 650.0) + 273.15
        if "surface_area_cm2" not in segs:
            segs["surface_area_cm2"] = 1.0
        for col, default in {
            "flow_factor": 1.0,
            "delta_T_C": 0.0,
            "redox_class": "purified_baseline",
            "position_class": "nominal",
            "surface_class": "metal",
        }.items():
            if col not in segs:
                segs[col] = default

        n_steps = int(math.ceil(duration_h / dt_h))
        cr_ppm = float(initial_cr_ppm)
        cum_diss_mg = np.zeros(len(segs))
        cum_dep_mg = np.zeros(len(segs))
        rows: list[dict[str, float | int | str]] = []
        for step in range(n_steps + 1):
            time_h = min(step * dt_h, duration_h)
            rows.append(
                {
                    "time_h": time_h,
                    "salt_cr_ppm": cr_ppm,
                    "total_cr_dissolved_mg": float(cum_diss_mg.sum()),
                    "total_deposit_mg": float(cum_dep_mg.sum()),
                }
            )
            if step == n_steps:
                break
            actual_dt_h = min(dt_h, duration_h - time_h)
            dt_y = actual_dt_h / (365.25 * 24.0)
            for i, seg in segs.iterrows():
                area = float(seg.get("surface_area_cm2", 1.0) or 1.0)
                mat = str(seg.get("material_class", "generic_metal"))
                corr_depth = self.corrosion_rate_um_y(seg) * dt_y
                total_mg = um_to_mg_cm2(corr_depth, mat) * area
                cr_mg = total_mg * cr_weight_fraction(mat)
                cum_diss_mg[i] += cr_mg
                # Convert mg Cr to ppm by mass using salt density approximated as 2 g/cm3.
                salt_mass_mg = salt_volume_cm3 * 2.0 * 1000.0
                cr_ppm += cr_mg / max(salt_mass_mg, 1.0) * 1.0e6
                dep_depth = self.deposition_rate_um_y(seg) * dt_y
                dep_mg = um_to_mg_cm2(dep_depth, mat) * area
                # Deposition cannot exceed a relaxed available inventory in this simplified closure.
                available_mg = max(cr_ppm, 0.0) / 1.0e6 * salt_mass_mg
                dep_mg = min(dep_mg, 0.25 * available_mg)
                cum_dep_mg[i] += dep_mg
                cr_ppm -= dep_mg / max(salt_mass_mg, 1.0) * 1.0e6
                cr_ppm = max(cr_ppm, 0.0)
        return pd.DataFrame(rows)

    def parameter_table(self) -> pd.DataFrame:
        rows = []
        for spec in PARAMETER_SPECS:
            rows.append(
                {
                    "parameter": spec.name,
                    "value": self.params[spec.name],
                    "prior": spec.prior,
                    "prior_sigma": spec.prior_sigma,
                    "lower": spec.lower,
                    "upper": spec.upper,
                    "description": spec.description,
                }
            )
        return pd.DataFrame(rows)
