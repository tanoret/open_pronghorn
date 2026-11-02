"""State-evolving redox, inventory, and surface-depletion corrosion model.

The Dynamic Redox-Inventory-Depletion Network (DRIDN) is a third model layer
for molten-salt corrosion.  It retains the MSTDB-TC standard-state reaction
thermodynamics from :mod:`msr_corrosion_bv.thermochemical`, but advances the
following state variables through the exposure:

* Cr/Fe/Ni surface availability in the alloy,
* dissolved Cr/Fe/Ni inventories in a well-mixed salt compartment,
* a relaxing/consumable redox-buffer coordinate,
* cumulative reaction-front and equivalent mass-recession depths,
* species-resolved cold-leg deposition, and
* a parabolic grain-boundary penetration coordinate.

The model is intentionally reduced order.  It is a transparent dynamic network,
not CFD and not a spatially resolved CALPHAD diffusion calculation.  Its purpose
is to test whether feedback among surface depletion, Nernst product buildup,
redox-buffer consumption, and cold-leg capture improves endpoint prediction and
provides physically interpretable transients for the validation experiments.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .chemistry import R_GAS, SEC_PER_YEAR, T_REF_K, density
from .mechanistic import _row_value
from .thermochemical_data import (
    ALLOY_ELEMENT_MASS_FRACTIONS,
    ELEMENT_MOLAR_MASS_G_MOL,
)
from .thermochemical_model import MSTDBThermochemicalCorrosionModel

ELEMENTS: tuple[str, ...] = ("Cr", "Fe", "Ni")
ELEMENT_INDEX = {element: index for index, element in enumerate(ELEMENTS)}


def _explicit_bool(value: Any, name: str) -> bool:
    """Parse a required physical switch without truth-testing arbitrary text."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1"}:
            return True
        if normalized in {"false", "0"}:
            return False
    raise ValueError(f"{name} must be an explicit boolean")


@dataclass(frozen=True)
class DynamicParameterSpec:
    name: str
    initial: float
    lower: float
    upper: float
    prior: float
    prior_sigma: float
    description: str


DYNAMIC_PARAMETER_SPECS: tuple[DynamicParameterSpec, ...] = (
    DynamicParameterSpec(
        "log_rate_scale",
        0.0,
        -1.5,
        1.5,
        0.0,
        0.50,
        "Log multiplier on the thermochemical charge-transfer rate.",
    ),
    DynamicParameterSpec(
        "affinity_feedback_scale",
        0.035,
        0.0,
        0.30,
        0.035,
        0.060,
        "Scale converting Nernst product buildup into dynamic rate suppression.",
    ),
    DynamicParameterSpec(
        "log_surface_reservoir_um",
        math.log(40.0),
        math.log(0.5),
        math.log(500.0),
        math.log(40.0),
        1.0,
        "Effective near-surface alloy reservoir thickness, log(um).",
    ),
    DynamicParameterSpec(
        "log_surface_replenishment_y_inv",
        math.log(2.0),
        math.log(0.005),
        math.log(200.0),
        math.log(2.0),
        1.3,
        "Reference solid-state replenishment rate of the surface reservoir, log(1/y).",
    ),
    DynamicParameterSpec(
        "surface_availability_exponent",
        0.8,
        0.0,
        3.0,
        0.8,
        0.8,
        "Exponent coupling depleted surface availability to interfacial current.",
    ),
    DynamicParameterSpec(
        "surface_reservoir_cr_exponent",
        0.8,
        0.0,
        2.5,
        0.8,
        0.7,
        "Scaling of effective surface-reservoir thickness with bulk alloy Cr fraction.",
    ),
    DynamicParameterSpec(
        "log_dynamic_cr_exchange_bias",
        0.007,
        -0.5,
        0.5,
        0.0,
        0.15,
        "Small dynamic interfacial-transfer correction for selective Cr dissolution.",
    ),
    DynamicParameterSpec(
        "log_dynamic_fe_capture_bias",
        0.30,
        -0.5,
        0.5,
        0.30,
        0.20,
        "Fe-rich cold-leg capture bias constrained by the 316H FLiNaK deposit observation.",
    ),
    DynamicParameterSpec(
        "inventory_inhibition_scale",
        0.75,
        0.0,
        6.0,
        0.75,
        1.2,
        "Strength of dissolved-product inhibition of salt-side mass transfer.",
    ),
    DynamicParameterSpec(
        "log_redox_relaxation_y_inv",
        math.log(1.0),
        math.log(0.005),
        math.log(200.0),
        math.log(1.0),
        1.4,
        "Relaxation rate of transient FeF2, Be, or moisture redox perturbations, log(1/y).",
    ),
    DynamicParameterSpec(
        "redox_buffer_retention",
        0.75,
        0.0,
        1.0,
        0.75,
        0.25,
        "Long-time fraction of the imposed transient redox perturbation retained by buffering.",
    ),
    DynamicParameterSpec(
        "redox_consumption_per_um",
        0.005,
        0.0,
        0.25,
        0.005,
        0.025,
        "Dimensionless redox-buffer consumption per micrometre of reaction-front advance.",
    ),
    DynamicParameterSpec(
        "log_stress_interfacial_factor",
        0.0,
        0.0,
        math.log(12.0),
        0.0,
        0.80,
        "Log multiplier representing stress-opened active area and interfacial transport.",
    ),
    DynamicParameterSpec(
        "log_fluoride_impurity_interfacial_factor",
        0.0,
        0.0,
        math.log(12.0),
        0.0,
        0.80,
        "Log multiplier representing impurity-activated fluoride interfacial transport.",
    ),
    DynamicParameterSpec(
        "log_deposition_rate_y_inv_fuel",
        math.log(0.5),
        math.log(1.0e-3),
        math.log(100.0),
        math.log(0.5),
        1.5,
        "First-order cold-leg capture rate for fuel/fluoroborate salts, log(1/y).",
    ),
    DynamicParameterSpec(
        "log_deposition_rate_y_inv_flinak",
        math.log(2.0),
        math.log(1.0e-3),
        math.log(200.0),
        math.log(2.0),
        1.5,
        "First-order cold-leg capture rate for FLiNaK, log(1/y).",
    ),
    DynamicParameterSpec(
        "bulk_capture_area_multiplier_fuel",
        1.0,
        0.0,
        12.0,
        1.0,
        2.0,
        "Uninstrumented-to-coupon cold-surface capture multiplier for fuel/fluoroborate loops.",
    ),
    DynamicParameterSpec(
        "bulk_capture_area_multiplier_flinak",
        0.5,
        0.0,
        12.0,
        0.5,
        2.0,
        "Uninstrumented-to-coupon cold-surface capture multiplier for FLiNaK loops.",
    ),
    DynamicParameterSpec(
        "log_bulk_precipitation_rate_scale",
        0.0,
        math.log(0.05),
        math.log(50.0),
        0.0,
        1.0,
        "Log multiplier on nonlinear bulk precipitation under dissolved-metal supersaturation.",
    ),
    DynamicParameterSpec(
        "log_inventory_scale_msre",
        0.0,
        -2.0,
        2.0,
        0.0,
        0.75,
        "Log multiplier on MSRE-scale exposed-area-to-salt-mass source terms.",
    ),
    DynamicParameterSpec(
        "log_inventory_scale_loop",
        0.0,
        -2.0,
        2.0,
        0.0,
        0.75,
        "Log multiplier on small-loop exposed-area-to-salt-mass source terms.",
    ),
    DynamicParameterSpec(
        "log_deposit_area_scale_fuel",
        0.0,
        -math.log(10.0),
        math.log(10.0),
        0.0,
        0.80,
        "Log correction to the observed cold-leg area mapping for fuel/fluoroborate salts (within one decade).",
    ),
    DynamicParameterSpec(
        "log_deposit_area_scale_flinak",
        0.0,
        -math.log(10.0),
        math.log(10.0),
        0.0,
        0.80,
        "Log correction to the observed cold-leg area mapping for FLiNaK (within one decade).",
    ),
    DynamicParameterSpec(
        "log_mass_loss_scale",
        0.0,
        -1.0,
        1.0,
        0.0,
        0.45,
        "Log multiplier on non-congruent equivalent mass recession.",
    ),
    DynamicParameterSpec(
        "gb_dynamic_scale",
        1.0,
        0.2,
        3.0,
        1.0,
        0.55,
        "Multiplier on inherited grain-boundary diffusion lengths.",
    ),
    DynamicParameterSpec(
        "damage_affinity_scale",
        0.03,
        0.0,
        0.40,
        0.03,
        0.08,
        "Positive-redox enhancement of parabolic grain-boundary penetration.",
    ),
)

DYNAMIC_PARAMETER_NAMES = tuple(spec.name for spec in DYNAMIC_PARAMETER_SPECS)
DYNAMIC_FIXED_PARAMETER_NAMES: tuple[str, ...] = (
    # Weakly identifiable feedback terms are fixed at physically motivated
    # priors so the calibrated dynamic subset remains smaller than the number
    # of direct/range validation targets.
    "affinity_feedback_scale",
    "surface_reservoir_cr_exponent",
    "log_dynamic_cr_exchange_bias",
    "log_dynamic_fe_capture_bias",
    "inventory_inhibition_scale",
    "redox_consumption_per_um",
    "bulk_capture_area_multiplier_fuel",
    "bulk_capture_area_multiplier_flinak",
    "log_bulk_precipitation_rate_scale",
    "damage_affinity_scale",
)
DYNAMIC_FIT_PARAMETER_NAMES = tuple(
    name for name in DYNAMIC_PARAMETER_NAMES if name not in DYNAMIC_FIXED_PARAMETER_NAMES
)
_DYNAMIC_SPEC_BY_NAME = {spec.name: spec for spec in DYNAMIC_PARAMETER_SPECS}


def dynamic_initial_vector() -> np.ndarray:
    return np.asarray([spec.initial for spec in DYNAMIC_PARAMETER_SPECS], dtype=float)


def dynamic_lower_bounds() -> np.ndarray:
    return np.asarray([spec.lower for spec in DYNAMIC_PARAMETER_SPECS], dtype=float)


def dynamic_upper_bounds() -> np.ndarray:
    return np.asarray([spec.upper for spec in DYNAMIC_PARAMETER_SPECS], dtype=float)


def dynamic_vector_to_params(vector: Iterable[float]) -> dict[str, float]:
    return {name: float(value) for name, value in zip(DYNAMIC_PARAMETER_NAMES, vector)}


def dynamic_params_to_vector(params: Mapping[str, float]) -> np.ndarray:
    return np.asarray([float(params[name]) for name in DYNAMIC_PARAMETER_NAMES], dtype=float)


def dynamic_fit_initial_vector() -> np.ndarray:
    return np.asarray([_DYNAMIC_SPEC_BY_NAME[name].initial for name in DYNAMIC_FIT_PARAMETER_NAMES], dtype=float)


def dynamic_fit_lower_bounds() -> np.ndarray:
    return np.asarray([_DYNAMIC_SPEC_BY_NAME[name].lower for name in DYNAMIC_FIT_PARAMETER_NAMES], dtype=float)


def dynamic_fit_upper_bounds() -> np.ndarray:
    return np.asarray([_DYNAMIC_SPEC_BY_NAME[name].upper for name in DYNAMIC_FIT_PARAMETER_NAMES], dtype=float)


def dynamic_fit_params_to_vector(params: Mapping[str, float]) -> np.ndarray:
    return np.asarray([float(params[name]) for name in DYNAMIC_FIT_PARAMETER_NAMES], dtype=float)


def dynamic_fit_vector_to_params(
    vector: Iterable[float],
    base_params: Mapping[str, float] | None = None,
) -> dict[str, float]:
    params = (
        dynamic_vector_to_params(dynamic_initial_vector())
        if base_params is None
        else {name: float(base_params[name]) for name in DYNAMIC_PARAMETER_NAMES}
    )
    params.update({name: float(value) for name, value in zip(DYNAMIC_FIT_PARAMETER_NAMES, vector)})
    return params


@dataclass(frozen=True)
class DynamicCaseContext:
    """Thermochemical and geometric terms that do not change during fitting."""

    row: Mapping[str, Any]
    measurement_id: str
    material: str
    salt_class: str
    redox_class: str
    position_class: str
    temperature_K: float
    cold_temperature_K: float
    delta_T_C: float
    flow_factor: float
    time_years: float
    density_g_cm3: float
    mass_fractions: np.ndarray
    cr_fraction_ratio: float
    inventory_scale: str
    explicit_inventory_scale: float
    deposition_closure: str
    log_exchange_offsets: np.ndarray
    selectivity_scale: float
    product_floor_ppm: float
    initial_dissolved_ppm: np.ndarray
    affinity_baseline: np.ndarray
    redox_shift_initial: float
    log_charge_base_no_redox: float
    mass_transfer_rate_um_y: float
    inventory_capacity_ppm: float
    area_to_salt_mass_cm2_g: float
    inventory_coupling_factor: float
    cold_capture_fraction: np.ndarray
    deposit_area_factor: float
    mass_loss_fraction: float
    cr_diffusion_cm2_s: float
    front_damage_multiplier: float
    gb_length_multiplier: float
    transient_redox: bool
    stress_interfacial_activation: bool
    fluoride_impurity_interfacial_activation: bool
    chloride_salt: bool


@dataclass
class DynamicSimulationResult:
    """Endpoint and optional trajectory from one validation-case simulation."""

    front_depth_um: float
    mass_recession_um: float
    mass_loss_mg_cm2: float
    mass_gain_mg_cm2: float
    igc_depth_um: float
    corrosion_rate_um_y: float
    dissolved_ppm: dict[str, float]
    cumulative_source_ppm: dict[str, float]
    deposit_mg_cm2: dict[str, float]
    bulk_captured_ppm: dict[str, float]
    surface_availability: dict[str, float]
    redox_shift: float
    mass_balance_relative_error: float
    trajectory: pd.DataFrame | None = None

    @property
    def deposit_fractions(self) -> dict[str, float]:
        total = sum(self.deposit_mg_cm2.values())
        if total <= 0.0:
            return {element: 0.0 for element in ELEMENTS}
        return {element: self.deposit_mg_cm2[element] / total for element in ELEMENTS}

    @property
    def cumulative_source_fractions(self) -> dict[str, float]:
        total = sum(self.cumulative_source_ppm.values())
        if total <= 0.0:
            return {element: 0.0 for element in ELEMENTS}
        return {element: self.cumulative_source_ppm[element] / total for element in ELEMENTS}


class DynamicRedoxInventoryDepletionModel:
    """State-evolving dynamic network built on an MSTDB thermochemical model."""

    def __init__(
        self,
        thermochemical_model: MSTDBThermochemicalCorrosionModel,
        params: Mapping[str, float] | None = None,
        *,
        integration_steps: int = 120,
    ) -> None:
        self.thermochemical_model = thermochemical_model
        self.params = dynamic_vector_to_params(dynamic_initial_vector())
        if params is not None:
            self.params.update({str(name): float(value) for name, value in params.items()})
        self.integration_steps = max(int(integration_steps), 12)

    @classmethod
    def from_vector(
        cls,
        thermochemical_model: MSTDBThermochemicalCorrosionModel,
        vector: Iterable[float],
        *,
        integration_steps: int = 120,
    ) -> "DynamicRedoxInventoryDepletionModel":
        return cls(
            thermochemical_model,
            dynamic_vector_to_params(vector),
            integration_steps=integration_steps,
        )

    def _param(self, name: str) -> float:
        return float(self.params[name])

    def build_context(
        self,
        row: pd.Series | Mapping[str, Any],
        *,
        redox_override: str | None = None,
    ) -> DynamicCaseContext:
        thermo = self.thermochemical_model
        original = row.to_dict() if isinstance(row, pd.Series) else dict(row)
        required_context = (
            "measurement_id",
            "salt_class",
            "redox_class",
            "position_class",
            "temperature_K",
            "delta_T_C",
            "flow_factor",
            "time_years",
            "inventory_source_material",
            "inventory_coupling_factor",
            "inventory_scale",
            "area_to_salt_mass_cm2_g",
            "deposition_closure",
            "deposit_area_factor",
            "initial_dissolved_Cr_ppm",
            "initial_dissolved_Fe_ppm",
            "initial_dissolved_Ni_ppm",
            "explicit_inventory_scale",
            "transient_redox",
            "stress_interfacial_activation",
            "fluoride_impurity_interfacial_activation",
            "chloride_salt",
        )
        missing = [name for name in required_context if _row_value(original, name, None) is None]
        if missing:
            measurement_id = str(_row_value(original, "measurement_id", "<unknown>"))
            raise ValueError(
                f"Advanced case {measurement_id} is missing explicit context fields: {missing}"
            )
        measurement_id = str(_row_value(original, "measurement_id")).strip()
        salt_class = str(_row_value(original, "salt_class")).strip()
        material = str(_row_value(original, "inventory_source_material")).strip()
        if not measurement_id or not salt_class:
            raise ValueError("measurement_id and salt_class must be nonempty")
        if material not in thermo._ADVANCED_SUPPORTED_MATERIALS:
            raise ValueError(f"inventory_source_material is not recognized: {material!r}")
        inventory_coupling_factor = float(_row_value(original, "inventory_coupling_factor"))
        if not np.isfinite(inventory_coupling_factor) or not 0.0 <= inventory_coupling_factor <= 1.0:
            raise ValueError("inventory_coupling_factor must be in [0, 1]")
        inventory_scale = str(_row_value(original, "inventory_scale")).lower()
        if inventory_scale not in {"msre", "loop", "explicit"}:
            raise ValueError("inventory_scale must be 'msre', 'loop', or 'explicit'")
        explicit_inventory_scale = float(_row_value(original, "explicit_inventory_scale"))
        if not np.isfinite(explicit_inventory_scale) or explicit_inventory_scale <= 0.0:
            raise ValueError("explicit_inventory_scale must be finite and positive")
        area_to_mass = float(_row_value(original, "area_to_salt_mass_cm2_g"))
        if not np.isfinite(area_to_mass) or area_to_mass <= 0.0:
            raise ValueError("area_to_salt_mass_cm2_g must be finite and positive")
        deposition_closure = str(_row_value(original, "deposition_closure")).lower()
        if deposition_closure not in {"fuel", "flinak"}:
            raise ValueError("deposition_closure must be 'fuel' or 'flinak'")
        area_factor = float(_row_value(original, "deposit_area_factor"))
        if not np.isfinite(area_factor) or area_factor <= 0.0:
            raise ValueError("deposit_area_factor must be finite and positive")
        initial_dissolved = np.asarray(
            [
                float(_row_value(original, f"initial_dissolved_{element}_ppm"))
                for element in ELEMENTS
            ],
            dtype=float,
        )
        if not np.all(np.isfinite(initial_dissolved)) or np.any(initial_dissolved < 0.0):
            raise ValueError("initial_dissolved_Cr/Fe/Ni_ppm must be finite and nonnegative")
        transient_redox = _explicit_bool(
            _row_value(original, "transient_redox"), "transient_redox"
        )
        stress_interfacial_activation = _explicit_bool(
            _row_value(original, "stress_interfacial_activation"),
            "stress_interfacial_activation",
        )
        fluoride_impurity_interfacial_activation = _explicit_bool(
            _row_value(original, "fluoride_impurity_interfacial_activation"),
            "fluoride_impurity_interfacial_activation",
        )
        chloride_salt = _explicit_bool(_row_value(original, "chloride_salt"), "chloride_salt")
        sim_row = dict(original)
        sim_row["material_class"] = material
        redox_class = redox_override or str(_row_value(original, "redox_class")).strip()
        if not redox_class:
            raise ValueError("redox_class must be nonempty")
        sim_row["redox_class"] = redox_class
        position = str(_row_value(original, "position_class")).strip()
        if position not in thermo._ADVANCED_POSITIONS:
            raise ValueError(f"unsupported advanced position_class: {position!r}")
        thermo._validate_advanced_labels(sim_row)
        temperature_K = float(_row_value(original, "temperature_K"))
        if not np.isfinite(temperature_K) or temperature_K <= 0.0:
            raise ValueError("temperature_K must be finite and positive")
        delta_T_C = float(_row_value(original, "delta_T_C"))
        if not np.isfinite(delta_T_C) or delta_T_C < 0.0:
            raise ValueError("delta_T_C must be finite and nonnegative")
        hot_temperature_K, cold_temperature_K = thermo._hot_cold_temperatures_K(original)
        flow = float(_row_value(original, "flow_factor"))
        if not np.isfinite(flow) or flow <= 0.0:
            raise ValueError("flow_factor must be finite and positive")
        time_years = float(_row_value(original, "time_years"))
        if not np.isfinite(time_years) or time_years < 0.0:
            raise ValueError("time_years must be finite and nonnegative")

        mass_fraction_mapping = ALLOY_ELEMENT_MASS_FRACTIONS.get(
            material,
            ALLOY_ELEMENT_MASS_FRACTIONS["generic_metal"],
        )
        mass_fractions = np.asarray([mass_fraction_mapping[element] for element in ELEMENTS], dtype=float)
        exchange_offsets = np.asarray(
            [
                0.0,
                float(thermo.params["log_exchange_fe_relative"]),
                float(thermo.params["log_exchange_ni_relative"]),
            ],
            dtype=float,
        )
        floor = max(math.exp(float(thermo.params["log_product_floor_ppm"])), 1.0e-8)
        baseline_affinity = np.asarray(
            [
                thermo.reaction_log_K_over_Q(
                    element,
                    sim_row,
                    hot_temperature_K,
                    redox_override="purified_baseline",
                    product_ppm=floor,
                )
                for element in ELEMENTS
            ],
            dtype=float,
        )
        redox_shift = thermo.redox_log_shift(redox_class)

        cr_fraction = max(mass_fraction_mapping["Cr"], 1.0e-12)
        cr_ratio = max(cr_fraction / 0.07, 1.0e-8)
        ea_corr = float(thermo.base_model.params["Ea_corr_kJ_mol"])
        thermal = ea_corr * 1000.0 / R_GAS * (1.0 / T_REF_K - 1.0 / hot_temperature_K)
        log_charge_no_redox = (
            float(thermo.params["log_front_rate0_um_y"])
            + thermal
            + float(thermo.params["cr_activity_exponent"]) * math.log(cr_ratio)
            + thermo.thermochemical_salt_log_drive(sim_row, hot_temperature_K)
        )
        thermal_mt = float(thermo.params["Ea_mass_transfer_kJ_mol"]) * 1000.0 / R_GAS * (
            1.0 / T_REF_K - 1.0 / hot_temperature_K
        )
        mass_transfer_rate = math.exp(
            float(thermo.params["log_mass_transfer_cap_um_y"])
            + thermal_mt
            + float(thermo.params["flow_mass_transfer_exponent"]) * math.log(flow)
        )
        capture = np.asarray(
            [
                thermo.cold_capture_fraction(element, salt_class, hot_temperature_K, cold_temperature_K)
                if cold_temperature_K < hot_temperature_K
                else 0.0
                for element in ELEMENTS
            ],
            dtype=float,
        )
        capacity = math.exp(float(thermo.params["log_inventory_capacity_ppm"]))

        front_multiplier = 1.0
        gb_multiplier = float(thermo.params["gb_diffusion_multiplier"])
        if redox_class == "multi_alloy":
            gb_multiplier = float(thermo.params["multi_alloy_gb_multiplier"])
            if material == "alloy_x750":
                gb_multiplier *= float(thermo.params["x750_gb_multiplier"])
        elif salt_class == "chloride":
            front_multiplier = float(thermo.params["chloride_reaction_multiplier"])
            gb_multiplier = float(thermo.params["chloride_gb_multiplier"])
        elif redox_class == "tellurium":
            gb_multiplier = float(thermo.params["tellurium_gb_multiplier"])

        return DynamicCaseContext(
            row=original,
            measurement_id=measurement_id,
            material=material,
            salt_class=salt_class,
            redox_class=redox_class,
            position_class=position,
            temperature_K=hot_temperature_K,
            cold_temperature_K=cold_temperature_K,
            delta_T_C=delta_T_C,
            flow_factor=flow,
            time_years=time_years,
            density_g_cm3=density(material),
            mass_fractions=mass_fractions,
            cr_fraction_ratio=cr_ratio,
            inventory_scale=inventory_scale,
            explicit_inventory_scale=explicit_inventory_scale,
            deposition_closure=deposition_closure,
            log_exchange_offsets=exchange_offsets,
            selectivity_scale=float(thermo.params["selectivity_affinity_scale"]),
            product_floor_ppm=floor,
            initial_dissolved_ppm=initial_dissolved,
            affinity_baseline=baseline_affinity,
            redox_shift_initial=redox_shift,
            log_charge_base_no_redox=log_charge_no_redox,
            mass_transfer_rate_um_y=mass_transfer_rate,
            inventory_capacity_ppm=capacity,
            area_to_salt_mass_cm2_g=area_to_mass,
            inventory_coupling_factor=inventory_coupling_factor,
            cold_capture_fraction=capture,
            deposit_area_factor=area_factor,
            mass_loss_fraction=thermo.mass_loss_fraction(sim_row),
            cr_diffusion_cm2_s=thermo.cr_diffusion_cm2_s(sim_row),
            front_damage_multiplier=front_multiplier,
            gb_length_multiplier=gb_multiplier,
            transient_redox=transient_redox,
            stress_interfacial_activation=stress_interfacial_activation,
            fluoride_impurity_interfacial_activation=fluoride_impurity_interfacial_activation,
            chloride_salt=chloride_salt,
        )

    def _species_fractions(
        self,
        context: DynamicCaseContext,
        dissolved_ppm: np.ndarray,
        surface_availability: np.ndarray,
        redox_shift: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        product_ratio = np.maximum(dissolved_ppm / context.product_floor_ppm, 1.0)
        affinity = context.affinity_baseline + redox_shift - np.log(product_ratio)
        relative = affinity - affinity[ELEMENT_INDEX["Cr"]]
        if not np.all(np.isfinite(context.mass_fractions)) or np.any(context.mass_fractions < 0.0):
            raise ValueError("mass_fractions must be finite and nonnegative")
        positive_mass = context.mass_fractions > 0.0
        if not np.any(positive_mass):
            raise ValueError("at least one Cr/Fe/Ni mass fraction must be positive")
        log_mass_fractions = np.full(len(ELEMENTS), -np.inf, dtype=float)
        log_mass_fractions[positive_mass] = np.log(context.mass_fractions[positive_mass])
        log_weights = (
            log_mass_fractions
            + context.log_exchange_offsets
            + context.selectivity_scale * np.clip(relative, -80.0, 80.0)
            + self._param("surface_availability_exponent")
            * np.log(np.maximum(surface_availability, 1.0e-8))
        )
        log_weights[ELEMENT_INDEX["Cr"]] += self._param("log_dynamic_cr_exchange_bias")
        maximum = float(np.max(log_weights))
        weights = np.exp(log_weights - maximum)
        fractions = weights / max(float(weights.sum()), 1.0e-30)
        return fractions, affinity

    def simulate_context(
        self,
        context: DynamicCaseContext,
        *,
        return_trajectory: bool = False,
    ) -> DynamicSimulationResult:
        time_total = context.time_years
        reservoir_um = (
            math.exp(self._param("log_surface_reservoir_um"))
            * context.cr_fraction_ratio ** self._param("surface_reservoir_cr_exponent")
        )
        replenishment = math.exp(self._param("log_surface_replenishment_y_inv"))
        replenishment_relative = np.asarray([1.0, 3.0, 1.8], dtype=float)
        if context.inventory_scale == "msre":
            inventory_scale = math.exp(self._param("log_inventory_scale_msre"))
        elif context.inventory_scale == "loop":
            inventory_scale = math.exp(self._param("log_inventory_scale_loop"))
        else:
            inventory_scale = context.explicit_inventory_scale
        deposit_area_scale = math.exp(
            self._param(
                "log_deposit_area_scale_flinak"
                if context.deposition_closure == "flinak"
                else "log_deposit_area_scale_fuel"
            )
        )
        mass_loss_scale = math.exp(self._param("log_mass_loss_scale"))
        redox_relaxation = math.exp(self._param("log_redox_relaxation_y_inv"))
        redox_target = (
            self._param("redox_buffer_retention") * context.redox_shift_initial
            if context.transient_redox
            else context.redox_shift_initial
        )
        deposition_rate_constant = math.exp(
            self._param(
                "log_deposition_rate_y_inv_flinak"
                if context.deposition_closure == "flinak"
                else "log_deposition_rate_y_inv_fuel"
            )
        )
        deposition_rate_constants = (
            deposition_rate_constant
            * context.cold_capture_fraction
            * math.sqrt(max(context.flow_factor, 0.05))
        )
        deposition_rate_constants[ELEMENT_INDEX["Fe"]] *= math.exp(
            self._param("log_dynamic_fe_capture_bias")
        )
        bulk_capture_area_multiplier = self._param(
            "bulk_capture_area_multiplier_flinak"
            if context.deposition_closure == "flinak"
            else "bulk_capture_area_multiplier_fuel"
        )
        bulk_capture_rate_constants = deposition_rate_constants * bulk_capture_area_multiplier
        area_to_mass = context.area_to_salt_mass_cm2_g * inventory_scale
        diffusion_um2_y = max(context.cr_diffusion_cm2_s, 0.0) * 1.0e8 * SEC_PER_YEAR
        loss_fraction = float(np.clip(context.mass_loss_fraction * mass_loss_scale, 0.01, 1.0))
        deposit_conversion = (
            1.0e-3
            / max(area_to_mass, 1.0e-30)
            * context.deposit_area_factor
            * deposit_area_scale
        )
        redox_lower = min(0.0, context.redox_shift_initial)
        redox_upper = max(0.0, context.redox_shift_initial)

        # State order: surface[3], dissolved[3], cumulative source[3],
        # coupon deposit mass[3], front depth, mass recession, gb depth squared,
        # redox, and cumulative uninstrumented/bulk capture[3].
        y0 = np.zeros(19, dtype=float)
        y0[0:3] = 1.0
        y0[3:6] = context.initial_dissolved_ppm
        y0[15] = context.redox_shift_initial

        def rates_from_state(state: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, float]:
            surface = np.clip(state[0:3], 0.015, 1.0)
            dissolved = np.maximum(state[3:6], 0.0)
            redox_shift = float(np.clip(state[15], redox_lower, redox_upper))
            fractions, _ = self._species_fractions(context, dissolved, surface, redox_shift)
            product_feedback = -float(
                np.dot(fractions, np.log(np.maximum(dissolved / context.product_floor_ppm, 1.0)))
            )
            surface_feedback = float(np.dot(fractions, np.log(np.maximum(surface, 1.0e-8))))
            interfacial_log_factor = 0.0
            if context.stress_interfacial_activation:
                interfacial_log_factor += self._param("log_stress_interfacial_factor")
            if context.fluoride_impurity_interfacial_activation and not context.chloride_salt:
                interfacial_log_factor += self._param(
                    "log_fluoride_impurity_interfacial_factor"
                )
            log_charge_rate = (
                context.log_charge_base_no_redox
                + redox_shift
                + self._param("log_rate_scale")
                + interfacial_log_factor
                + self._param("affinity_feedback_scale") * product_feedback
                + self._param("surface_availability_exponent") * surface_feedback
            )
            charge_rate = math.exp(float(np.clip(log_charge_rate, -60.0, 60.0)))
            inventory_ratio = max(float(dissolved.sum()), 0.0) / max(context.inventory_capacity_ppm, 1.0e-30)
            mass_transfer_rate = (
                context.mass_transfer_rate_um_y
                * math.exp(interfacial_log_factor)
                / (1.0 + self._param("inventory_inhibition_scale") * inventory_ratio)
            )
            total_rate = 1.0 / (
                1.0 / max(charge_rate, 1.0e-30)
                + 1.0 / max(mass_transfer_rate, 1.0e-30)
            )
            source_rate_ppm_y = (
                total_rate
                * context.density_g_cm3
                * 100.0
                * area_to_mass
                * context.inventory_coupling_factor
                * fractions
            )
            bulk_capacity_rate = min(
                200.0,
                0.35
                * math.exp(self._param("log_bulk_precipitation_rate_scale"))
                * (
                    float(dissolved.sum())
                    / max(context.inventory_capacity_ppm, 1.0e-30)
                )
                ** 2,
            )
            return total_rate, fractions, source_rate_ppm_y, bulk_capacity_rate

        if time_total == 0.0:
            zeros = {element: 0.0 for element in ELEMENTS}
            initial = {
                element: float(context.initial_dissolved_ppm[index])
                for index, element in enumerate(ELEMENTS)
            }
            surface = {element: 1.0 for element in ELEMENTS}
            trajectory = None
            if return_trajectory:
                instantaneous_rate, _, _, _ = rates_from_state(y0)
                trajectory_row: dict[str, float] = {
                    "time_years": 0.0,
                    "front_depth_um": 0.0,
                    "mass_recession_um": 0.0,
                    "mass_loss_mg_cm2": 0.0,
                    "mass_gain_mg_cm2": 0.0,
                    "igc_depth_um": 0.0,
                    "instantaneous_corrosion_rate_um_y": float(instantaneous_rate),
                    "redox_log_shift": context.redox_shift_initial,
                }
                for element in ELEMENTS:
                    trajectory_row[f"dissolved_{element}_ppm"] = initial[element]
                    trajectory_row[f"surface_{element}_availability"] = 1.0
                    trajectory_row[f"deposit_{element}_mg_cm2"] = 0.0
                    trajectory_row[f"cumulative_source_{element}_ppm"] = 0.0
                    trajectory_row[f"bulk_captured_{element}_ppm"] = 0.0
                trajectory = pd.DataFrame([trajectory_row])
            return DynamicSimulationResult(
                front_depth_um=0.0,
                mass_recession_um=0.0,
                mass_loss_mg_cm2=0.0,
                mass_gain_mg_cm2=0.0,
                igc_depth_um=0.0,
                corrosion_rate_um_y=0.0,
                dissolved_ppm=initial,
                cumulative_source_ppm=dict(zeros),
                deposit_mg_cm2=dict(zeros),
                bulk_captured_ppm=dict(zeros),
                surface_availability=surface,
                redox_shift=context.redox_shift_initial,
                mass_balance_relative_error=0.0,
                trajectory=trajectory,
            )

        def rhs(_time: float, state: np.ndarray) -> np.ndarray:
            total_rate, fractions, source_rate_ppm_y, bulk_capacity_rate = rates_from_state(state)
            surface = np.clip(state[0:3], 0.015, 1.0)
            dissolved = np.maximum(state[3:6], 0.0)
            redox_shift = float(np.clip(state[15], redox_lower, redox_upper))
            derivative = np.zeros_like(state)
            derivative[0:3] = (
                replenishment * replenishment_relative * (1.0 - surface)
                - total_rate * fractions / max(reservoir_um, 1.0e-30)
            )
            for index in range(3):
                if state[index] <= 0.015 and derivative[index] < 0.0:
                    derivative[index] = 0.0
                if state[index] >= 1.0 and derivative[index] > 0.0:
                    derivative[index] = 0.0
            total_sink_rate = (
                deposition_rate_constants
                + bulk_capture_rate_constants
                + bulk_capacity_rate
            )
            derivative[3:6] = source_rate_ppm_y - total_sink_rate * dissolved
            derivative[6:9] = source_rate_ppm_y
            derivative[9:12] = deposition_rate_constants * dissolved * deposit_conversion
            derivative[16:19] = (bulk_capture_rate_constants + bulk_capacity_rate) * dissolved
            derivative[12] = total_rate
            derivative[13] = total_rate * loss_fraction
            gb_multiplier = context.gb_length_multiplier * self._param("gb_dynamic_scale")
            affinity_damage = math.exp(
                self._param("damage_affinity_scale") * max(redox_shift, 0.0)
            )
            derivative[14] = (
                4.0
                * diffusion_um2_y
                / math.pi
                * gb_multiplier**2
                * affinity_damage
            )
            if context.transient_redox:
                consumption = (
                    -math.copysign(
                        self._param("redox_consumption_per_um") * total_rate,
                        redox_shift,
                    )
                    if abs(redox_shift) > 1.0e-12
                    else 0.0
                )
                derivative[15] = redox_relaxation * (redox_target - redox_shift) + consumption
                if state[15] <= redox_lower and derivative[15] < 0.0:
                    derivative[15] = 0.0
                if state[15] >= redox_upper and derivative[15] > 0.0:
                    derivative[15] = 0.0
            return derivative

        t_eval = (
            np.linspace(0.0, time_total, self.integration_steps + 1)
            if return_trajectory
            else None
        )
        solution = solve_ivp(
            rhs,
            (0.0, time_total),
            y0,
            method="LSODA",
            t_eval=t_eval,
            rtol=2.0e-5,
            atol=1.0e-8,
            max_step=max(time_total / self.integration_steps, 1.0e-8),
        )
        if not solution.success:
            raise RuntimeError(f"DRIDN integration failed for {context.measurement_id}: {solution.message}")
        final = solution.y[:, -1]
        surface = np.clip(final[0:3], 0.015, 1.0)
        dissolved = np.maximum(final[3:6], 0.0)
        cumulative_source = np.maximum(final[6:9], 0.0)
        deposit_mg_cm2 = np.maximum(final[9:12], 0.0)
        front_depth = max(float(final[12]), 0.0)
        mass_recession = max(float(final[13]), 0.0)
        gb_depth_squared = max(float(final[14]), 0.0)
        redox_shift = float(np.clip(final[15], redox_lower, redox_upper))
        bulk_captured_ppm = np.maximum(final[16:19], 0.0)
        igc_depth = context.front_damage_multiplier * front_depth + math.sqrt(gb_depth_squared)
        local_captured_ppm = deposit_mg_cm2 / max(deposit_conversion, 1.0e-30)
        inventory_expected = context.initial_dissolved_ppm + cumulative_source
        inventory_accounted = dissolved + local_captured_ppm + bulk_captured_ppm
        mass_balance_relative_error = float(
            np.max(
                np.abs(inventory_expected - inventory_accounted)
                / np.maximum(inventory_expected, 1.0e-12)
            )
        )

        trajectory = None
        if return_trajectory:
            rows: list[dict[str, float]] = []
            for column, time_y in enumerate(solution.t):
                state = solution.y[:, column]
                rate, _, _, _ = rates_from_state(state)
                state_surface = np.clip(state[0:3], 0.015, 1.0)
                state_dissolved = np.maximum(state[3:6], 0.0)
                state_deposit = np.maximum(state[9:12], 0.0)
                state_front = max(float(state[12]), 0.0)
                state_mass = max(float(state[13]), 0.0)
                state_gb = max(float(state[14]), 0.0)
                row = {
                    "time_years": float(time_y),
                    "front_depth_um": state_front,
                    "mass_recession_um": state_mass,
                    "mass_loss_mg_cm2": state_mass * context.density_g_cm3 * 0.1,
                    "mass_gain_mg_cm2": float(state_deposit.sum()),
                    "igc_depth_um": context.front_damage_multiplier * state_front + math.sqrt(state_gb),
                    "instantaneous_corrosion_rate_um_y": float(rate),
                    "redox_log_shift": float(np.clip(state[15], redox_lower, redox_upper)),
                }
                for index, element in enumerate(ELEMENTS):
                    row[f"dissolved_{element}_ppm"] = float(state_dissolved[index])
                    row[f"surface_{element}_availability"] = float(state_surface[index])
                    row[f"deposit_{element}_mg_cm2"] = float(state_deposit[index])
                    row[f"cumulative_source_{element}_ppm"] = float(max(state[6 + index], 0.0))
                    row[f"bulk_captured_{element}_ppm"] = float(max(state[16 + index], 0.0))
                rows.append(row)
            trajectory = pd.DataFrame(rows)

        return DynamicSimulationResult(
            front_depth_um=front_depth,
            mass_recession_um=mass_recession,
            mass_loss_mg_cm2=mass_recession * context.density_g_cm3 * 0.1,
            mass_gain_mg_cm2=float(deposit_mg_cm2.sum()),
            igc_depth_um=float(igc_depth),
            corrosion_rate_um_y=mass_recession / max(time_total, 1.0e-30),
            dissolved_ppm={element: float(dissolved[index]) for index, element in enumerate(ELEMENTS)},
            cumulative_source_ppm={
                element: float(cumulative_source[index]) for index, element in enumerate(ELEMENTS)
            },
            deposit_mg_cm2={element: float(deposit_mg_cm2[index]) for index, element in enumerate(ELEMENTS)},
            bulk_captured_ppm={element: float(bulk_captured_ppm[index]) for index, element in enumerate(ELEMENTS)},
            surface_availability={element: float(surface[index]) for index, element in enumerate(ELEMENTS)},
            redox_shift=redox_shift,
            mass_balance_relative_error=mass_balance_relative_error,
            trajectory=trajectory,
        )

    def simulate(
        self,
        row: pd.Series | Mapping[str, Any],
        *,
        redox_override: str | None = None,
        return_trajectory: bool = False,
    ) -> DynamicSimulationResult:
        return self.simulate_context(
            self.build_context(row, redox_override=redox_override),
            return_trajectory=return_trajectory,
        )

    def predict_response(self, row: pd.Series | Mapping[str, Any]) -> float:
        kind = str(_row_value(row, "response_kind", ""))
        if kind == "cr_diffusion_cm2_s":
            return self.thermochemical_model.cr_diffusion_cm2_s(row)
        if kind in {"redox_acceleration_ratio", "redox_acceleration_qualitative"}:
            oxidized = self.simulate(row, redox_override="oxidizing_fef2")
            baseline = self.simulate(row, redox_override="purified_baseline")
            return oxidized.corrosion_rate_um_y / max(baseline.corrosion_rate_um_y, 1.0e-30)
        result = self.simulate(row)
        if kind == "corrosion_rate_um_y":
            return result.corrosion_rate_um_y
        if kind == "corrosion_depth_um":
            return result.front_depth_um
        if kind == "mass_loss_mg_cm2":
            return result.mass_loss_mg_cm2
        if kind == "mass_gain_mg_cm2":
            return result.mass_gain_mg_cm2
        if kind == "salt_cr_ppm":
            return result.dissolved_ppm["Cr"]
        if kind == "salt_fe_decrease_ppm":
            source = result.cumulative_source_ppm
            cr_moles = max(source["Cr"], 0.0) / ELEMENT_MOLAR_MASS_G_MOL["Cr"]
            ni_moles = max(source["Ni"], 0.0) / ELEMENT_MOLAR_MASS_G_MOL["Ni"]
            fe_moles = max(source["Fe"], 0.0) / ELEMENT_MOLAR_MASS_G_MOL["Fe"]
            decrease = max(0.0, cr_moles + ni_moles - fe_moles) * ELEMENT_MOLAR_MASS_G_MOL["Fe"]
            limit = math.exp(float(self.thermochemical_model.params["log_initial_fe2_ppm"]))
            return min(limit, decrease)
        if kind == "igc_depth_um":
            return result.igc_depth_um
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
                for spec in DYNAMIC_PARAMETER_SPECS
            ]
        )
