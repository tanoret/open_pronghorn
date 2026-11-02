"""Mass loss, diffusion, deposition, inventories, and response mappings."""

from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np
import pandas as pd

from .chemistry import R_GAS, SEC_PER_YEAR, density
from .mechanistic import _row_value
from .thermochemical_equilibrium import ADVANCED_SUPPORTED_SALTS
from .thermochemical_data import (
    ALLOY_ELEMENT_MASS_FRACTIONS,
    ELEMENT_MOLAR_MASS_G_MOL,
    THERMOCHEMICAL_PARAMETER_SPECS,
)


class ThermochemicalTransportMixin:
    _ADVANCED_POSITIONS = frozenset(
        {"cold_leg", "hot_leg", "core_or_reactor", "gas_or_offgas", "nominal"}
    )

    @staticmethod
    def _exposure_time_years(row: pd.Series | Mapping[str, Any]) -> float:
        raw = _row_value(row, "time_years", None)
        if raw is None:
            raise ValueError("time_years is required for advanced endpoint physics")
        value = float(raw)
        if not np.isfinite(value) or value < 0.0:
            raise ValueError("time_years must be finite and nonnegative")
        return value

    def mass_loss_fraction(self, row: pd.Series | Mapping[str, Any]) -> float:
        material, salt, _ = self._validate_advanced_labels(row)
        cr_fraction = ALLOY_ELEMENT_MASS_FRACTIONS.get(material, ALLOY_ELEMENT_MASS_FRACTIONS["generic_metal"])["Cr"]
        cr_ratio = max(cr_fraction / 0.07, 1.0e-6)
        logit = self._param("mass_loss_fraction_logit")
        reference_fraction = 1.0 / (1.0 + math.exp(-logit))
        value = reference_fraction * cr_ratio ** self._param("mass_loss_cr_exponent")
        if salt == "flinak" and cr_ratio > 1.0:
            value *= cr_ratio ** self._param("flinak_high_cr_selectivity_exponent")
        return float(np.clip(value, 0.03, 1.0))

    def corrosion_rate_um_y(self, row: pd.Series | Mapping[str, Any], redox_override: str | None = None) -> float:
        return self.dissolution_front_rate_um_y(row, redox_override=redox_override) * self.mass_loss_fraction(row)

    def corrosion_depth_um(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = self._exposure_time_years(row)
        return self.dissolution_front_rate_um_y(row) * time_y

    def mass_loss_mg_cm2(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = self._exposure_time_years(row)
        material = str(_row_value(row, "material_class", "generic_metal"))
        return self.corrosion_rate_um_y(row) * time_y * density(material) * 0.1

    def cr_diffusion_cm2_s(self, row: pd.Series | Mapping[str, Any]) -> float:
        self._validate_advanced_labels(row)
        return self.base_model.cr_diffusion_cm2_s(row)

    def diffusion_length_um(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = self._exposure_time_years(row)
        if time_y == 0.0:
            return 0.0
        return 2.0 * math.sqrt(max(self.cr_diffusion_cm2_s(row), 0.0) * time_y * SEC_PER_YEAR / math.pi) * 1.0e4

    def _hot_cold_temperatures_K(
        self, row: pd.Series | Mapping[str, Any]
    ) -> tuple[float, float]:
        """Normalize a row temperature to the explicit hot/cold salt pair."""
        raw_temperature = _row_value(row, "temperature_K", None)
        raw_delta = _row_value(row, "delta_T_C", None)
        if raw_temperature is None:
            raise ValueError("temperature_K is required for advanced thermochemical physics")
        if raw_delta is None:
            raise ValueError("delta_T_C is required for advanced thermochemical physics")
        row_temperature_K = float(raw_temperature)
        delta_T_C = float(raw_delta)
        position = str(_row_value(row, "position_class", "")).strip()
        if not np.isfinite(row_temperature_K) or row_temperature_K <= 0.0:
            raise ValueError("temperature_K must be finite and positive")
        if not np.isfinite(delta_T_C) or delta_T_C < 0.0:
            raise ValueError("delta_T_C must be finite and nonnegative")
        if position not in self._ADVANCED_POSITIONS:
            raise ValueError(
                f"unsupported advanced position_class: {position!r}; "
                f"expected one of {sorted(self._ADVANCED_POSITIONS)}"
            )
        explicit_cold_C = _row_value(row, "temperature_cold_C", None)
        explicit_cold_K: float | None = None
        if explicit_cold_C is not None:
            explicit_cold_K = float(explicit_cold_C) + 273.15
            if not np.isfinite(explicit_cold_K) or explicit_cold_K <= 0.0:
                raise ValueError("temperature_cold_C must be finite and physically positive")
        if position == "cold_leg":
            hot_K = row_temperature_K + delta_T_C
            expected_cold_K = row_temperature_K
        else:
            hot_K = row_temperature_K
            expected_cold_K = row_temperature_K - delta_T_C
        if explicit_cold_K is not None and not math.isclose(
            explicit_cold_K, expected_cold_K, rel_tol=0.0, abs_tol=1.0e-6
        ):
            raise ValueError(
                "temperature_cold_C is inconsistent with temperature_K, delta_T_C, "
                "and position_class"
            )
        cold_K = expected_cold_K if explicit_cold_K is None else explicit_cold_K
        if not np.isfinite(hot_K) or hot_K <= 0.0:
            raise ValueError("hot temperature must be finite and positive")
        if not np.isfinite(cold_K) or cold_K <= 0.0:
            raise ValueError("cold temperature must be finite and positive")
        if cold_K > hot_K:
            raise ValueError("cold temperature must not exceed hot temperature")
        return float(hot_K), float(cold_K)

    def saturation_activity(self, element: str, salt_class: str, temperature_K: float) -> float:
        if salt_class not in ADVANCED_SUPPORTED_SALTS:
            raise ValueError(f"unsupported advanced salt_class: {salt_class!r}")
        if not np.isfinite(temperature_K) or temperature_K <= 0.0:
            raise ValueError("temperature_K must be finite and positive")
        database, species_map = self._database_and_species(salt_class)
        species = species_map[element]
        delta_g = database.standard_gibbs_J_mol(species["solid"], temperature_K) - database.standard_gibbs_J_mol(
            species["dissolved"], temperature_K
        )
        return float(np.clip(math.exp(np.clip(delta_g / (R_GAS * temperature_K), -80.0, 0.0)), 1.0e-30, 1.0))

    def cold_capture_fraction(self, element: str, salt_class: str, hot_K: float, cold_K: float) -> float:
        if not np.isfinite(hot_K) or hot_K <= 0.0:
            raise ValueError("hot_K must be finite and positive")
        if not np.isfinite(cold_K) or cold_K <= 0.0:
            raise ValueError("cold_K must be finite and positive")
        if cold_K > hot_K:
            raise ValueError("cold_K must not exceed hot_K")
        if hot_K == cold_K:
            return 0.0
        hot = self.saturation_activity(element, salt_class, hot_K)
        cold = self.saturation_activity(element, salt_class, cold_K)
        return float(np.clip(1.0 - cold / max(hot, 1.0e-30), 0.0, 1.0))

    def deposition_species_fractions(self, row: pd.Series | Mapping[str, Any]) -> dict[str, float]:
        """Return the normalized Cr/Fe/Ni composition of the captured cold-leg deposit."""
        hot_K, cold_K = self._hot_cold_temperatures_K(row)
        if hot_K == cold_K:
            return {element: 0.0 for element in ("Cr", "Fe", "Ni")}
        salt_class = str(_row_value(row, "salt_class", "generic_salt"))
        source = self.species_flux_fractions(row, hot_K)
        weights = {
            element: source[element] * self.cold_capture_fraction(element, salt_class, hot_K, cold_K)
            for element in source
        }
        total = sum(weights.values())
        if total <= 0.0:
            return {element: 0.0 for element in ("Cr", "Fe", "Ni")}
        return {element: weights[element] / total for element in ("Cr", "Fe", "Ni")}

    def mass_gain_mg_cm2(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = self._exposure_time_years(row)
        hot_K, cold_K = self._hot_cold_temperatures_K(row)
        if hot_K == cold_K:
            return 0.0
        salt_class = str(_row_value(row, "salt_class", "generic_salt"))
        material = str(_row_value(row, "material_class", "generic_metal"))
        donor_row = row.to_dict() if isinstance(row, pd.Series) else dict(row)
        donor_row["temperature_K"] = hot_K
        donor_rate = self.dissolution_front_rate_um_y(donor_row)
        donor_mass = donor_rate * time_y * density(material) * 0.1
        fractions = self.species_flux_fractions(row, hot_K)
        weighted_capture = sum(
            fractions[element] * self.cold_capture_fraction(element, salt_class, hot_K, cold_K)
            for element in fractions
        )
        return donor_mass * (self._param("deposition_capture_area_factor_flinak") if salt_class == "flinak" else self._param("deposition_capture_area_factor_fuel")) * weighted_capture

    def _inventory_source_material(self, row: pd.Series | Mapping[str, Any]) -> str:
        explicit = _row_value(row, "inventory_source_material", None)
        if not explicit:
            raise ValueError("inventory_source_material is required for advanced inventory physics")
        return str(explicit)

    def _area_to_salt_mass(self, row: pd.Series | Mapping[str, Any]) -> float:
        explicit = _row_value(row, "area_to_salt_mass_cm2_g", None)
        if explicit is not None:
            value = float(explicit)
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError("area_to_salt_mass_cm2_g must be finite and positive")
            return value
        raise ValueError("area_to_salt_mass_cm2_g is required for advanced inventory physics")

    def _inventory_coupling_factor(self, row: pd.Series | Mapping[str, Any]) -> float:
        explicit = _row_value(row, "inventory_coupling_factor", None)
        if explicit is None:
            raise ValueError("inventory_coupling_factor is required for advanced inventory physics")
        value = float(explicit)
        if not np.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError("inventory_coupling_factor must be finite and in [0, 1]")
        return value

    def species_inventory_ppm(self, row: pd.Series | Mapping[str, Any]) -> dict[str, float]:
        time_y = self._exposure_time_years(row)
        source_material = self._inventory_source_material(row)
        donor_row = row.to_dict() if isinstance(row, pd.Series) else dict(row)
        donor_row["material_class"] = source_material
        T_hot, cold_K = self._hot_cold_temperatures_K(row)
        donor_row["temperature_K"] = T_hot
        net_depth_um = self.dissolution_front_rate_um_y(donor_row) * time_y
        total_mass_g_cm2 = net_depth_um * 1.0e-4 * density(source_material)
        fractions = self.species_flux_fractions(donor_row, T_hot)
        salt_class = str(_row_value(row, "salt_class", "generic_salt"))
        capture_factor = (self._param("deposition_capture_area_factor_flinak") if salt_class == "flinak" else self._param("deposition_capture_area_factor_fuel"))
        area_to_mass = self._area_to_salt_mass(row)
        coupling_factor = self._inventory_coupling_factor(row)
        raw_inventory: dict[str, float] = {}
        for element, fraction in fractions.items():
            retained = 1.0
            if T_hot > cold_K:
                retained = max(
                    0.0,
                    1.0 - min(1.0, capture_factor * self.cold_capture_fraction(element, salt_class, T_hot, cold_K)),
                )
            raw_inventory[element] = (
                total_mass_g_cm2
                * fraction
                * area_to_mass
                * coupling_factor
                * retained
                * 1.0e6
            )

        # The capacity is a total dissolved-corrosion-product capacity, not an
        # independent cap for every species.  Preserve the mechanistic source
        # fractions while smoothly saturating the combined Cr+Fe+Ni inventory.
        raw_total = sum(raw_inventory.values())
        if raw_total <= 0.0:
            return {element: 0.0 for element in ("Cr", "Fe", "Ni")}
        capacity = math.exp(self._param("log_inventory_capacity_ppm"))
        dissolved_total = capacity * (1.0 - math.exp(-raw_total / max(capacity, 1.0e-30)))
        scale = dissolved_total / raw_total
        return {element: raw_inventory[element] * scale for element in ("Cr", "Fe", "Ni")}

    def salt_cr_ppm(self, row: pd.Series | Mapping[str, Any]) -> float:
        return self.species_inventory_ppm(row)["Cr"]

    def salt_fe_decrease_ppm(self, row: pd.Series | Mapping[str, Any]) -> float:
        inventory = self.species_inventory_ppm(row)
        cr_moles = max(inventory["Cr"], 0.0) / ELEMENT_MOLAR_MASS_G_MOL["Cr"]
        ni_moles = max(inventory["Ni"], 0.0) / ELEMENT_MOLAR_MASS_G_MOL["Ni"]
        fe_produced_moles = max(inventory["Fe"], 0.0) / ELEMENT_MOLAR_MASS_G_MOL["Fe"]
        cementation_decrease = max(0.0, cr_moles + ni_moles - fe_produced_moles) * ELEMENT_MOLAR_MASS_G_MOL["Fe"]
        return min(math.exp(self._param("log_initial_fe2_ppm")), cementation_decrease)

    def redox_acceleration_ratio(self, row: pd.Series | Mapping[str, Any]) -> float:
        oxidized = self.corrosion_rate_um_y(row, redox_override="oxidizing_fef2")
        baseline = self.corrosion_rate_um_y(row, redox_override="purified_baseline")
        return oxidized / max(baseline, 1.0e-30)

    def igc_depth_um(self, row: pd.Series | Mapping[str, Any]) -> float:
        time_y = self._exposure_time_years(row)
        front = self.dissolution_front_rate_um_y(row) * time_y
        diffusion = self.diffusion_length_um(row)
        redox = str(_row_value(row, "redox_class", "purified_baseline"))
        salt = str(_row_value(row, "salt_class", "generic_salt"))
        material = str(_row_value(row, "material_class", "generic_metal"))
        if redox == "multi_alloy":
            multiplier = self._param("multi_alloy_gb_multiplier")
            if material == "alloy_x750":
                multiplier *= self._param("x750_gb_multiplier")
            return front + multiplier * diffusion
        if salt == "chloride":
            return self._param("chloride_reaction_multiplier") * front + self._param("chloride_gb_multiplier") * diffusion
        if redox == "tellurium":
            return front + self._param("tellurium_gb_multiplier") * diffusion
        return front + self._param("gb_diffusion_multiplier") * diffusion

    def predict_response(self, row: pd.Series | Mapping[str, Any]) -> float:
        # A cold-leg measurement reports the cold temperature in
        # ``temperature_K``.  The governing dissolution, diffusion, and redox
        # terms are nevertheless evaluated at the explicit hot temperature;
        # only deposition/capture uses the preserved cold temperature.  Put
        # every response through the same canonical hot-leg representation so
        # nested endpoint methods cannot accidentally reinterpret a cold-leg
        # row at its raw measurement temperature.
        hot_K, cold_K = self._hot_cold_temperatures_K(row)
        physics_row = row.to_dict() if isinstance(row, pd.Series) else dict(row)
        physics_row["temperature_K"] = hot_K
        physics_row["temperature_cold_C"] = cold_K - 273.15
        physics_row["delta_T_C"] = hot_K - cold_K
        physics_row["position_class"] = "hot_leg"

        kind = str(_row_value(physics_row, "response_kind", ""))
        if kind == "corrosion_rate_um_y":
            return self.corrosion_rate_um_y(physics_row)
        if kind == "corrosion_depth_um":
            return self.corrosion_depth_um(physics_row)
        if kind == "mass_loss_mg_cm2":
            return self.mass_loss_mg_cm2(physics_row)
        if kind == "mass_gain_mg_cm2":
            return self.mass_gain_mg_cm2(physics_row)
        if kind == "salt_cr_ppm":
            return self.salt_cr_ppm(physics_row)
        if kind == "salt_fe_decrease_ppm":
            return self.salt_fe_decrease_ppm(physics_row)
        if kind == "cr_diffusion_cm2_s":
            return self.cr_diffusion_cm2_s(physics_row)
        if kind in {"redox_acceleration_ratio", "redox_acceleration_qualitative"}:
            return self.redox_acceleration_ratio(physics_row)
        if kind == "igc_depth_um":
            return self.igc_depth_um(physics_row)
        return np.nan

    def diagnostic_row(self, row: pd.Series | Mapping[str, Any]) -> dict[str, Any]:
        hot_K, cold_K = self._hot_cold_temperatures_K(row)
        physics_row = row.to_dict() if isinstance(row, pd.Series) else dict(row)
        physics_row["temperature_K"] = hot_K
        physics_row["temperature_cold_C"] = cold_K - 273.15
        physics_row["position_class"] = "hot_leg"
        T = hot_K
        salt = str(_row_value(row, "salt_class", "generic_salt"))
        fractions = self.species_flux_fractions(physics_row, T)
        inventory = self.species_inventory_ppm(physics_row)
        deposit_fractions = self.deposition_species_fractions(physics_row)
        result: dict[str, Any] = {
            "measurement_id": _row_value(row, "measurement_id", None),
            "case_id": _row_value(row, "case_id", None),
            "salt_class": salt,
            "redox_class": _row_value(row, "redox_class", None),
            "temperature_K": T,
            "pred_corrosion_rate_um_y": self.corrosion_rate_um_y(physics_row),
            "pred_cr_ppm": inventory["Cr"],
            "pred_fe_ppm": inventory["Fe"],
            "pred_ni_ppm": inventory["Ni"],
            "pred_total_corrosion_product_ppm": sum(inventory.values()),
            "pred_corrosion_depth_um": self.corrosion_depth_um(physics_row),
            "pred_igc_depth_um": self.igc_depth_um(physics_row),
            "pred_mass_gain_mg_cm2": self.mass_gain_mg_cm2(physics_row),
        }
        for element in ("Cr", "Fe", "Ni"):
            result[f"ln_K_over_Q_{element}"] = self.reaction_log_K_over_Q(element, physics_row, T)
            result[f"dissolution_fraction_{element}"] = fractions[element]
            result[f"saturation_activity_{element}"] = self.saturation_activity(element, salt, T)
            result[f"deposition_fraction_{element}"] = deposit_fractions[element]
        return result

    def predict_cases(self, cases: pd.DataFrame) -> pd.DataFrame:
        """Return corrosion, damage, deposition, and Cr/Fe/Ni diagnostics for cases."""
        return pd.DataFrame([self.diagnostic_row(row) for _, row in cases.iterrows()])

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
                for spec in THERMOCHEMICAL_PARAMETER_SPECS
            ]
        )
