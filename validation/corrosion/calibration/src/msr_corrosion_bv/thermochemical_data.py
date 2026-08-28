"""Shared constants and regularized parameters for the thermochemical model."""

from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from typing import Iterable, Mapping

import numpy as np

ELEMENT_MOLAR_MASS_G_MOL: Mapping[str, float] = {
    "Cr": 51.9961,
    "Fe": 55.845,
    "Ni": 58.6934,
}
HALIDE_MOLAR_MASS_G_MOL: Mapping[str, Mapping[str, float]] = {
    "fluoride": {"Cr": 89.9929, "Fe": 93.8418, "Ni": 96.6902},
    "chloride": {"Cr": 122.9021, "Fe": 126.7510, "Ni": 129.5994},
}
MEAN_SALT_MOLAR_MASS_G_MOL: Mapping[str, float] = {
    "fluoride_fuel": 42.0,
    "flibe": 33.1,
    "flinak": 41.3,
    "fluoroborate": 103.0,
    "chloride": 71.1,
    "generic_salt": 55.0,
}

# Approximate bulk alloy mass fractions used to close the Cr/Fe/Ni source term.
# Unlisted alloying elements remain in the matrix and affect density but are not
# assigned an electrochemical dissolution branch in this reduced model.
ALLOY_ELEMENT_MASS_FRACTIONS: Mapping[str, Mapping[str, float]] = {
    "hastelloy_n": {"Cr": 0.070, "Fe": 0.050, "Ni": 0.710},
    "modified_hastelloy_n": {"Cr": 0.070, "Fe": 0.050, "Ni": 0.710},
    "gh3535": {"Cr": 0.070, "Fe": 0.050, "Ni": 0.700},
    "stainless_304": {"Cr": 0.190, "Fe": 0.690, "Ni": 0.100},
    "stainless_304l": {"Cr": 0.190, "Fe": 0.690, "Ni": 0.100},
    "stainless_316": {"Cr": 0.170, "Fe": 0.660, "Ni": 0.120},
    "stainless_316h": {"Cr": 0.170, "Fe": 0.660, "Ni": 0.120},
    "stainless_316l": {"Cr": 0.170, "Fe": 0.660, "Ni": 0.120},
    "alloy_x750": {"Cr": 0.155, "Fe": 0.070, "Ni": 0.700},
    "in625": {"Cr": 0.215, "Fe": 0.050, "Ni": 0.610},
    "ni_alloy": {"Cr": 0.150, "Fe": 0.050, "Ni": 0.700},
    "generic_metal": {"Cr": 0.120, "Fe": 0.500, "Ni": 0.300},
    "graphite": {"Cr": 0.0, "Fe": 0.0, "Ni": 0.0},
}


FLUORIDE_SPECIES: Mapping[str, Mapping[str, str]] = {
    "Cr": {"metal": "Cr_S1(s)", "dissolved": "CrF2_L1(liq)", "solid": "CrF2_P21/c_No.14(s)"},
    "Fe": {"metal": "Fe_bcc(s)", "dissolved": "FeF2_L1(liq)", "solid": "FeF2_P42/mnm_No.136(s)"},
    "Ni": {"metal": "Ni_fcc(s)", "dissolved": "NiF2_L1(liq)", "solid": "NiF2_P42/mnm_No.136(s)"},
}
CHLORIDE_SPECIES: Mapping[str, Mapping[str, str]] = {
    "Cr": {"metal": "Cr_S1(s)", "dissolved": "CrCl2_L1(liq)", "solid": "CrCl2_Pnnm_No.58(s)"},
    "Fe": {"metal": "Fe_bcc(s)", "dissolved": "FeCl2_L1(liq)", "solid": "FeCl2_R3M_No166(s)"},
    "Ni": {"metal": "Ni_Solid_FCC(s)", "dissolved": "NiCl2_L1(liq)", "solid": "NiCl2_R3M_No.166(s)"},
}


@dataclass(frozen=True)
class ThermochemicalParameterSpec:
    name: str
    initial: float
    lower: float
    upper: float
    prior: float
    prior_sigma: float
    description: str


THERMOCHEMICAL_PARAMETER_SPECS: tuple[ThermochemicalParameterSpec, ...] = (
    ThermochemicalParameterSpec("log_front_rate0_um_y", 0.0575, math.log(0.01), math.log(100.0), 0.0575, 1.0, "Reference selective-dissolution front rate at 650 C, log(um/y)."),
    ThermochemicalParameterSpec("cr_activity_exponent", 1.86, 0.2, 4.0, 1.86, 0.9, "Exponent mapping bulk alloy Cr fraction to reaction-front kinetics."),
    ThermochemicalParameterSpec("log_gamma_cr_flinak", -5.90, -12.0, 4.0, -5.90, 2.0, "Inferred Cr(II) activity-coefficient correction in FLiNaK; replaceable by native Thermochimica output."),
    ThermochemicalParameterSpec("log_gamma_cr_flibe", -5.42, -12.0, 4.0, -5.42, 2.0, "Inferred Cr(II) activity-coefficient correction in FLiBe; replaceable by native Thermochimica output."),
    ThermochemicalParameterSpec("log_gamma_cr_fluoroborate", -5.45, -12.0, 4.0, -5.45, 2.0, "Inferred Cr(II) activity-coefficient correction in fluoroborate salt."),
    ThermochemicalParameterSpec("log_gamma_cr_chloride", -6.03, -14.0, 4.0, -6.03, 2.5, "Inferred Cr(II) activity-coefficient correction in chloride salt."),
    ThermochemicalParameterSpec("selectivity_affinity_scale", 0.10, 0.0, 1.0, 0.10, 0.20, "Scale converting relative MSTDB Cr/Fe/Ni affinities into selective-dissolution fractions."),
    ThermochemicalParameterSpec("log_exchange_fe_relative", -1.0, -5.0, 2.0, -1.0, 1.0, "Fe exchange-current offset relative to Cr for species partitioning."),
    ThermochemicalParameterSpec("log_exchange_ni_relative", -2.5, -7.0, 1.0, -2.5, 1.2, "Ni exchange-current offset relative to Cr for species partitioning."),
    ThermochemicalParameterSpec("log_product_floor_ppm", math.log(1.0), math.log(1.0e-4), math.log(300.0), math.log(1.0), 2.0, "Initial dissolved corrosion-product floor in the Nernst quotient, log(ppm)."),
    ThermochemicalParameterSpec("log_mass_transfer_cap_um_y", 6.81, math.log(10.0), math.log(1.0e5), 6.81, 1.4, "Reference salt-side mass-transfer capacity, log(um/y)."),
    ThermochemicalParameterSpec("Ea_mass_transfer_kJ_mol", 23.4, 0.0, 120.0, 23.4, 35.0, "Apparent activation energy for salt-side transport."),
    ThermochemicalParameterSpec("flow_mass_transfer_exponent", 0.75, 0.05, 2.0, 0.75, 0.55, "Flow-factor exponent for salt-side transport."),
    ThermochemicalParameterSpec("gb_diffusion_multiplier", 0.63, 0.05, 3.0, 0.63, 0.6, "Multiplier on Cr solid-state diffusion length for ordinary IGC."),
    ThermochemicalParameterSpec("multi_alloy_gb_multiplier", 6.42, 1.0, 40.0, 6.42, 8.0, "Galvanic/grain-boundary penetration multiplier in multi-alloy tests."),
    ThermochemicalParameterSpec("x750_gb_multiplier", 12.62, 1.0, 40.0, 12.62, 10.0, "Additional X-750 susceptibility multiplier."),
    ThermochemicalParameterSpec("chloride_reaction_multiplier", 5.62, 0.5, 30.0, 5.62, 7.0, "Reaction-front multiplier for impure chloride IGC."),
    ThermochemicalParameterSpec("chloride_gb_multiplier", 3.74, 0.5, 30.0, 3.74, 7.0, "Diffusion-length multiplier for impure chloride IGC."),
    ThermochemicalParameterSpec("tellurium_gb_multiplier", 8.0, 1.0, 40.0, 8.0, 10.0, "Diffusion-length multiplier for tellurium-assisted IGC."),
    ThermochemicalParameterSpec("deposition_capture_area_factor_fuel", 1.6, 0.02, 8.0, 1.6, 1.5, "Fuel-salt cold-leg capture efficiency times donor/cold surface-area ratio."),
    ThermochemicalParameterSpec("deposition_capture_area_factor_flinak", 3.0, 0.02, 10.0, 3.0, 2.0, "FLiNaK cold-leg capture efficiency times donor/cold surface-area ratio."),
    ThermochemicalParameterSpec("log_area_to_salt_mass_msre_cm2_g", -1.98, math.log(0.002), math.log(5.0), -1.98, 1.2, "Effective exposed-area/salt-mass ratio for MSRE-scale inventory closure."),
    ThermochemicalParameterSpec("log_area_to_salt_mass_loop_cm2_g", 2.05, math.log(0.03), math.log(80.0), 2.05, 1.4, "Effective exposed-area/salt-mass ratio for small-loop inventory closure."),
    ThermochemicalParameterSpec("log_initial_fe2_ppm", math.log(100.0), math.log(1.0), math.log(3000.0), math.log(100.0), 1.2, "Initial dissolved Fe(II) inventory available for Cr/Fe cementation, log(ppm)."),
    ThermochemicalParameterSpec("mass_loss_fraction_logit", 1.47, -5.0, 5.0, 1.47, 1.2, "Reference fraction of the selectively affected layer removed as net mass."),
    ThermochemicalParameterSpec("mass_loss_cr_exponent", 0.23, -2.5, 1.5, 0.23, 0.7, "Cr-content exponent for non-congruent mass loss."),
    ThermochemicalParameterSpec("flinak_high_cr_selectivity_exponent", -1.55, -3.5, 0.5, -1.55, 0.9, "Additional selective-leaching exponent for high-Cr alloys in FLiNaK."),
    ThermochemicalParameterSpec("log_inventory_capacity_ppm", 6.19, math.log(100.0), math.log(2.0e4), 6.19, 1.0, "Effective mixed-salt corrosion-product inventory capacity, log(ppm)."),
)
THERMOCHEMICAL_PARAMETER_NAMES = tuple(spec.name for spec in THERMOCHEMICAL_PARAMETER_SPECS)


def thermochemical_initial_vector() -> np.ndarray:
    return np.asarray([spec.initial for spec in THERMOCHEMICAL_PARAMETER_SPECS], dtype=float)


def thermochemical_lower_bounds() -> np.ndarray:
    return np.asarray([spec.lower for spec in THERMOCHEMICAL_PARAMETER_SPECS], dtype=float)


def thermochemical_upper_bounds() -> np.ndarray:
    return np.asarray([spec.upper for spec in THERMOCHEMICAL_PARAMETER_SPECS], dtype=float)


def thermochemical_vector_to_params(vector: Iterable[float]) -> dict[str, float]:
    return {name: float(value) for name, value in zip(THERMOCHEMICAL_PARAMETER_NAMES, vector)}


def _logsumexp(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=float)
    if not len(array):
        return -math.inf
    maximum = float(np.max(array))
    return maximum + math.log(float(np.exp(array - maximum).sum()))


@functools.lru_cache(maxsize=None)
def _alloy_atomic_activities(material: str) -> dict[str, float]:
    fractions = ALLOY_ELEMENT_MASS_FRACTIONS.get(material, ALLOY_ELEMENT_MASS_FRACTIONS["generic_metal"])
    tracked_moles = {element: fractions.get(element, 0.0) / ELEMENT_MOLAR_MASS_G_MOL[element] for element in ELEMENT_MOLAR_MASS_G_MOL}
    remainder = max(0.0, 1.0 - sum(fractions.values())) / 60.0
    denominator = max(sum(tracked_moles.values()) + remainder, 1.0e-30)
    return {element: max(moles / denominator, 1.0e-16) for element, moles in tracked_moles.items()}


def _ppm_to_activity(ppm: float, species_molar_mass: float, salt_class: str) -> float:
    mean_salt_molar_mass = MEAN_SALT_MOLAR_MASS_G_MOL.get(salt_class, 55.0)
    return max(float(ppm) * 1.0e-6 * mean_salt_molar_mass / species_molar_mass, 1.0e-30)


def _activity_to_ppm(activity: float, species_molar_mass: float, salt_class: str) -> float:
    mean_salt_molar_mass = MEAN_SALT_MOLAR_MASS_G_MOL.get(salt_class, 55.0)
    return max(float(activity), 0.0) * species_molar_mass / mean_salt_molar_mass * 1.0e6
