"""MSTDB/Nernst reaction driving forces and selective-dissolution kinetics."""

from __future__ import annotations

import math
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from .chemistry import R_GAS, T_REF_K
from .mechanistic import _row_value
from .model import MoltenSaltBVModel
from .mstdb import MSTDBPair
from .thermochemical_data import (
    ALLOY_ELEMENT_MASS_FRACTIONS,
    CHLORIDE_SPECIES,
    ELEMENT_MOLAR_MASS_G_MOL,
    FLUORIDE_SPECIES,
    HALIDE_MOLAR_MASS_G_MOL,
    _alloy_atomic_activities,
    _logsumexp,
    _ppm_to_activity,
    thermochemical_initial_vector,
    thermochemical_vector_to_params,
)


ADVANCED_SUPPORTED_MATERIALS = frozenset(
    {
        "hastelloy_n",
        "modified_hastelloy_n",
        "gh3535",
        "stainless_304",
        "stainless_304l",
        "stainless_316",
        "stainless_316h",
        "stainless_316l",
        "alloy_x750",
        "in625",
        "ni_alloy",
        "generic_metal",
    }
)
ADVANCED_SUPPORTED_SALTS = frozenset(
    {"fluoride_fuel", "flibe", "flinak", "fluoroborate", "chloride"}
)
ADVANCED_SUPPORTED_REDOX = frozenset(
    {
        "purified_baseline",
        "msre_or_fuel_baseline",
        "oxidizing_fef2",
        "reducing_be",
        "impure_moisture",
        "chloride_unspecified",
        "tellurium",
        "stressed",
        "multi_alloy",
        "fission_product",
        "gas_control",
    }
)

class ThermochemicalEquilibriumMixin:
    """Mechanistic corrosion model driven by MSTDB-TC reaction affinities."""

    _ADVANCED_SUPPORTED_MATERIALS = ADVANCED_SUPPORTED_MATERIALS
    _ADVANCED_SUPPORTED_SALTS = ADVANCED_SUPPORTED_SALTS
    _ADVANCED_SUPPORTED_REDOX = ADVANCED_SUPPORTED_REDOX

    def __init__(
        self,
        base_model: MoltenSaltBVModel,
        mstdb: MSTDBPair,
        params: Mapping[str, float] | None = None,
    ) -> None:
        self.base_model = base_model
        self.mstdb = mstdb
        self.params = thermochemical_vector_to_params(thermochemical_initial_vector()) if params is None else dict(params)

    @classmethod
    def from_vector(
        cls,
        base_model: MoltenSaltBVModel,
        mstdb: MSTDBPair,
        vector: Iterable[float],
    ) -> "MSTDBThermochemicalCorrosionModel":
        return cls(base_model, mstdb, thermochemical_vector_to_params(vector))

    def _param(self, name: str) -> float:
        return float(self.params[name])

    @staticmethod
    def _validate_advanced_labels(
        row: pd.Series | Mapping[str, Any], *, redox_override: str | None = None
    ) -> tuple[str, str, str]:
        material = str(_row_value(row, "material_class", "")).strip()
        salt = str(_row_value(row, "salt_class", "")).strip()
        redox = str(
            redox_override
            if redox_override is not None
            else _row_value(row, "redox_class", "")
        ).strip()
        if material not in ADVANCED_SUPPORTED_MATERIALS:
            raise ValueError(f"unsupported advanced material_class: {material!r}")
        if salt not in ADVANCED_SUPPORTED_SALTS:
            raise ValueError(f"unsupported advanced salt_class: {salt!r}")
        if redox not in ADVANCED_SUPPORTED_REDOX:
            raise ValueError(f"unsupported advanced redox_class: {redox!r}")
        return material, salt, redox

    @staticmethod
    def _salt_family(salt_class: str) -> str:
        return "chloride" if salt_class == "chloride" else "fluoride"

    def _database_and_species(self, salt_class: str) -> tuple[Any, Mapping[str, Mapping[str, str]]]:
        if salt_class == "chloride":
            return self.mstdb.chloride, CHLORIDE_SPECIES
        return self.mstdb.fluoride, FLUORIDE_SPECIES

    def redox_log_shift(self, redox_class: str) -> float:
        return self.base_model.redox_offset(redox_class) - self.base_model.redox_offset("purified_baseline")

    def _product_floor_ppm(self) -> float:
        return math.exp(self._param("log_product_floor_ppm"))

    def _fe2_buffer_ppm(self, salt_class: str, redox_class: str) -> float:
        baseline = {
            "flinak": 100.0,
            "flibe": 30.0,
            "fluoroborate": 100.0,
            "chloride": 250.0,
            "generic_salt": 100.0,
        }.get(salt_class, 100.0)
        if redox_class == "oxidizing_fef2":
            return 500.0
        shift = self.redox_log_shift(redox_class)
        return float(np.clip(baseline * math.exp(shift), 1.0e-3, 1.0e6))

    def _uranium_ratio(self, redox_class: str) -> float:
        # The effective redox offset is used only to infer an explicit Nernst
        # activity ratio when no measured potential is present.  For M(II), a
        # ratio multiplier exp(delta/2) produces a reaction-affinity shift delta.
        return 100.0 * math.exp(self.redox_log_shift(redox_class) / 2.0)

    def reaction_log_K_over_Q(
        self,
        element: str,
        row: pd.Series | Mapping[str, Any],
        temperature_K: float,
        *,
        redox_override: str | None = None,
        product_ppm: float | None = None,
    ) -> float:
        if not np.isfinite(temperature_K) or temperature_K <= 0.0:
            raise ValueError("temperature_K must be finite and positive")
        _, salt_class, redox_class = self._validate_advanced_labels(
            row, redox_override=redox_override
        )
        database, species_map = self._database_and_species(salt_class)
        species = species_map[element]
        metal_activity = _alloy_atomic_activities(str(_row_value(row, "material_class", "generic_metal")))[element]
        family = self._salt_family(salt_class)
        product_activity = _ppm_to_activity(
            self._product_floor_ppm() if product_ppm is None else product_ppm,
            HALIDE_MOLAR_MASS_G_MOL[family][element],
            salt_class,
        )

        if salt_class == "fluoride_fuel":
            reaction = {
                species["dissolved"]: 1.0,
                "UF3_L1(liq)": 2.0,
                species["metal"]: -1.0,
                "UF4_L1(liq)": -2.0,
            }
            lnK = database.equilibrium_log_constant(reaction, temperature_K)
            ratio = max(self._uranium_ratio(redox_class), 1.0e-30)
            return lnK - math.log(product_activity) + math.log(metal_activity) + 2.0 * math.log(ratio)

        if salt_class == "flibe" and redox_class == "reducing_be":
            reaction = {
                species["dissolved"]: 1.0,
                "Be_S1(s)": 1.0,
                species["metal"]: -1.0,
                "BeF2_L1(liq)": -1.0,
            }
            lnK = database.equilibrium_log_constant(reaction, temperature_K)
            # Be metal and BeF2 are treated as unit-activity buffer phases.
            return lnK - math.log(product_activity) + math.log(metal_activity)

        # FeX2/Fe is an explicit impurity/redox buffer for non-fuel salts.
        fe_species = species_map["Fe"]
        if element == "Fe":
            # The nominal buffer reaction is
            #
            #     FeX2 + Fe <=> Fe + FeX2,
            #
            # an exact identity whose net stoichiometric vector and reaction
            # Gibbs energy are zero. Constructing it with a dict silently
            # collapses duplicate keys and creates a fictitious Fe affinity.
            return 0.0
        reaction = {
            species["dissolved"]: 1.0,
            fe_species["metal"]: 1.0,
            species["metal"]: -1.0,
            fe_species["dissolved"]: -1.0,
        }
        lnK = database.equilibrium_log_constant(reaction, temperature_K)
        fe_metal_activity = _alloy_atomic_activities(str(_row_value(row, "material_class", "generic_metal")))["Fe"]
        fe2_activity = _ppm_to_activity(
            self._fe2_buffer_ppm(salt_class, redox_class),
            HALIDE_MOLAR_MASS_G_MOL[family]["Fe"],
            salt_class,
        )
        return (
            lnK
            - math.log(product_activity)
            - math.log(fe_metal_activity)
            + math.log(metal_activity)
            + math.log(fe2_activity)
        )

    def _salt_activity_correction(self, salt_class: str) -> float:
        mapping = {
            "flinak": "log_gamma_cr_flinak",
            "flibe": "log_gamma_cr_flibe",
            "fluoroborate": "log_gamma_cr_fluoroborate",
            "chloride": "log_gamma_cr_chloride",
        }
        key = mapping.get(salt_class)
        return 0.0 if key is None else self._param(key)

    def thermochemical_salt_log_drive(
        self,
        row: pd.Series | Mapping[str, Any],
        temperature_K: float,
    ) -> float:
        """MSTDB standard-state and activity contribution relative to fuel fluoride.

        The four fitted quantities are logarithms of corrosion-product activity
        coefficients, not arbitrary energy shifts.  A native Thermochimica SUBQ
        calculation can replace them directly when complete salt compositions are
        supplied.
        """
        salt = str(_row_value(row, "salt_class", "generic_salt"))
        if salt in {"fluoride_fuel", "no_salt"}:
            return 0.0
        material = str(_row_value(row, "material_class", "generic_metal"))
        target_row = {"material_class": material, "salt_class": salt, "redox_class": "purified_baseline"}
        reference_row = {"material_class": material, "salt_class": "fluoride_fuel", "redox_class": "purified_baseline"}
        raw_difference = self.reaction_log_K_over_Q("Cr", target_row, temperature_K) - self.reaction_log_K_over_Q(
            "Cr", reference_row, temperature_K
        )
        return raw_difference - self._salt_activity_correction(salt)

    def _species_log_rates(
        self,
        row: pd.Series | Mapping[str, Any],
        temperature_K: float,
        *,
        redox_override: str | None = None,
    ) -> dict[str, float]:
        material = str(_row_value(row, "material_class", "generic_metal"))
        mass_fractions = ALLOY_ELEMENT_MASS_FRACTIONS.get(material, ALLOY_ELEMENT_MASS_FRACTIONS["generic_metal"])
        affinity = {
            element: self.reaction_log_K_over_Q(element, row, temperature_K, redox_override=redox_override)
            for element in ("Cr", "Fe", "Ni")
        }
        reference = affinity["Cr"]
        scale = self._param("selectivity_affinity_scale")
        offsets = {"Cr": 0.0, "Fe": self._param("log_exchange_fe_relative"), "Ni": self._param("log_exchange_ni_relative")}
        return {
            element: math.log(max(mass_fractions.get(element, 0.0), 1.0e-16))
            + offsets[element]
            + scale * float(np.clip(affinity[element] - reference, -80.0, 80.0))
            for element in ("Cr", "Fe", "Ni")
        }

    def species_flux_fractions(
        self,
        row: pd.Series | Mapping[str, Any],
        temperature_K: float | None = None,
        *,
        redox_override: str | None = None,
    ) -> dict[str, float]:
        raw_temperature = temperature_K if temperature_K is not None else _row_value(row, "temperature_K", None)
        if raw_temperature is None:
            raise ValueError("temperature_K is required for advanced thermochemical physics")
        T = float(raw_temperature)
        if not np.isfinite(T) or T <= 0.0:
            raise ValueError("temperature_K must be finite and positive")
        logs = self._species_log_rates(row, T, redox_override=redox_override)
        norm = _logsumexp(logs.values())
        return {element: math.exp(value - norm) for element, value in logs.items()}

    def dissolution_front_rate_um_y(
        self,
        row: pd.Series | Mapping[str, Any],
        *,
        redox_override: str | None = None,
        temperature_K: float | None = None,
    ) -> float:
        raw_temperature = temperature_K if temperature_K is not None else _row_value(row, "temperature_K", None)
        if raw_temperature is None:
            raise ValueError("temperature_K is required for advanced thermochemical physics")
        T = float(raw_temperature)
        if not np.isfinite(T) or T <= 0.0:
            raise ValueError("temperature_K must be finite and positive")
        material, salt_class, redox = self._validate_advanced_labels(
            row, redox_override=redox_override
        )
        raw_flow = _row_value(row, "flow_factor", None)
        if raw_flow is None:
            raise ValueError("flow_factor is required for advanced thermochemical physics")
        flow = float(raw_flow)
        if not np.isfinite(flow) or flow <= 0.0:
            raise ValueError("flow_factor must be finite and positive")
        cr_fraction = ALLOY_ELEMENT_MASS_FRACTIONS.get(material, ALLOY_ELEMENT_MASS_FRACTIONS["generic_metal"])["Cr"]
        cr_ratio = max(cr_fraction / 0.07, 1.0e-6)
        Ea_corr = float(self.base_model.params["Ea_corr_kJ_mol"])
        thermal = Ea_corr * 1000.0 / R_GAS * (1.0 / T_REF_K - 1.0 / T)
        log_charge = (
            self._param("log_front_rate0_um_y")
            + thermal
            + self._param("cr_activity_exponent") * math.log(cr_ratio)
            + self.thermochemical_salt_log_drive(row, T)
            + self.redox_log_shift(redox)
        )
        rate_charge = float(np.exp(np.clip(log_charge, -60.0, 60.0)))
        thermal_mt = self._param("Ea_mass_transfer_kJ_mol") * 1000.0 / R_GAS * (1.0 / T_REF_K - 1.0 / T)
        log_mt = self._param("log_mass_transfer_cap_um_y") + thermal_mt + self._param("flow_mass_transfer_exponent") * math.log(flow)
        rate_mt = float(np.exp(np.clip(log_mt, -60.0, 60.0)))
        return 1.0 / (1.0 / max(rate_charge, 1.0e-30) + 1.0 / max(rate_mt, 1.0e-30))
