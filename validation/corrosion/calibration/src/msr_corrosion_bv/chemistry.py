"""Chemical constants, material properties, and unit conversions for the MSR BV model."""

from __future__ import annotations

import math
import re
from typing import Mapping

R_GAS = 8.31446261815324  # J mol-1 K-1
FARADAY = 96485.33212  # C mol-1
SEC_PER_YEAR = 365.25 * 24.0 * 3600.0
HOURS_PER_YEAR = 365.25 * 24.0
T_REF_K = 650.0 + 273.15

# Densities are engineering values intended for converting sparse validation data
# to comparable units. They are not substitutes for alloy-specific measurements.
DENSITY_G_CM3: Mapping[str, float] = {
    "hastelloy_n": 8.89,
    "modified_hastelloy_n": 8.89,
    "gh3535": 8.90,
    "stainless_304": 8.00,
    "stainless_304l": 8.00,
    "stainless_316": 8.00,
    "stainless_316h": 8.00,
    "stainless_316l": 8.00,
    "alloy_x750": 8.28,
    "in625": 8.44,
    "ni_alloy": 8.50,
    "generic_metal": 8.30,
    "graphite": 1.80,
}

# Approximate chromium mass fractions; used only by the optional inventory module.
ALLOY_CR_WT_FRAC: Mapping[str, float] = {
    "hastelloy_n": 0.07,
    "modified_hastelloy_n": 0.07,
    "gh3535": 0.07,
    "stainless_304": 0.19,
    "stainless_304l": 0.19,
    "stainless_316": 0.17,
    "stainless_316h": 0.17,
    "stainless_316l": 0.17,
    "alloy_x750": 0.155,
    "in625": 0.215,
    "ni_alloy": 0.15,
    "generic_metal": 0.12,
}

MOLAR_MASS_G_MOL: Mapping[str, float] = {
    "Cr": 51.9961,
    "Fe": 55.845,
    "Ni": 58.6934,
    "Mo": 95.95,
    "Nb": 92.9064,
    "Tc": 98.0,
    "Ru": 101.07,
    "Ag": 107.8682,
    "Sb": 121.760,
    "Te": 127.60,
    "generic_metal": 55.0,
}

VALENCE: Mapping[str, int] = {
    "Cr": 2,
    "Fe": 2,
    "Ni": 2,
    "Mo": 3,
    "Nb": 4,
    "Tc": 4,
    "Ru": 3,
    "Ag": 1,
    "Sb": 3,
    "Te": 2,
    "generic_metal": 2,
}


def canonical_material(text: object) -> str:
    """Map diverse material names from legacy reports into model classes."""
    s = "" if text is None else str(text).lower()
    if "graphite" in s and not any(k in s for k in ["hastelloy", "steel", "alloy"]):
        return "graphite"
    if "modified" in s and "hastelloy" in s:
        return "modified_hastelloy_n"
    if "gh3535" in s:
        return "gh3535"
    if "hastelloy n" in s or "hn" == s.strip():
        return "hastelloy_n"
    if "x750" in s or "x-750" in s:
        return "alloy_x750"
    if "316h" in s:
        return "stainless_316h"
    if "316l" in s:
        return "stainless_316l"
    if "316" in s:
        return "stainless_316"
    if "304l" in s:
        return "stainless_304l"
    if "304" in s:
        return "stainless_304"
    if "in625" in s or "inconel 625" in s or "alloy 625" in s:
        return "in625"
    if "nickel" in s or "ni-" in s or "ni " in s:
        return "ni_alloy"
    if "alloy" in s or "steel" in s or "metal" in s:
        return "generic_metal"
    return "generic_metal"


def canonical_salt(text: object) -> str:
    """Map salt descriptions to model classes."""
    s = "" if text is None else str(text).lower()
    if "no salt" in s or "n/a" == s.strip():
        return "no_salt"
    if "nabf4" in s or "fluoroborate" in s or "bf4" in s:
        return "fluoroborate"
    if "cl" in s and ("nacl" in s or "mgcl" in s or "kcl" in s or "chloride" in s):
        return "chloride"
    if "flinak" in s or ("lif" in s and "naf" in s and "kf" in s):
        return "flinak"
    if "flibe" in s or ("lif" in s and "bef2" in s and "uf" not in s and "zr" not in s and "th" not in s):
        return "flibe"
    if any(k in s for k in ["uf4", "zrf4", "thf4", "fuel"]):
        return "fluoride_fuel"
    if "lif" in s or "bef2" in s or "fluoride" in s:
        return "fluoride_fuel"
    if "molten salt" in s or "various" in s:
        return "generic_salt"
    return "generic_salt"


def density(material_class: str) -> float:
    return DENSITY_G_CM3.get(material_class, DENSITY_G_CM3["generic_metal"])


def cr_weight_fraction(material_class: str) -> float:
    return ALLOY_CR_WT_FRAC.get(material_class, ALLOY_CR_WT_FRAC["generic_metal"])


def mil_to_um(value: float) -> float:
    return 25.4 * value


def um_to_mil(value: float) -> float:
    return value / 25.4


def mg_cm2_to_um(value_mg_cm2: float, material_class: str) -> float:
    """Convert areal mass change to an equivalent dense-layer thickness in micrometers."""
    return value_mg_cm2 / (density(material_class) * 0.1)


def um_to_mg_cm2(value_um: float, material_class: str) -> float:
    return value_um * density(material_class) * 0.1


def mg_cm2_h_to_um_y(value_mg_cm2_h: float, material_class: str) -> float:
    return mg_cm2_to_um(value_mg_cm2_h, material_class) * HOURS_PER_YEAR


def corrosion_current_to_um_y(
    current_a_cm2: float,
    material_class: str,
    species: str = "Cr",
) -> float:
    """Convert an anodic current density to equivalent penetration rate.

    The conversion is Faradaic and assumes the selected species controls metal loss.
    The calibration absorbs deviations from this idealization.
    """
    z = VALENCE.get(species, VALENCE["generic_metal"])
    molar_mass = MOLAR_MASS_G_MOL.get(species, MOLAR_MASS_G_MOL["generic_metal"])
    rho = density(material_class)
    cm_per_s = current_a_cm2 * molar_mass / (z * FARADAY * rho)
    return cm_per_s * 1.0e4 * SEC_PER_YEAR


def um_y_to_corrosion_current(
    rate_um_y: float,
    material_class: str,
    species: str = "Cr",
) -> float:
    z = VALENCE.get(species, VALENCE["generic_metal"])
    molar_mass = MOLAR_MASS_G_MOL.get(species, MOLAR_MASS_G_MOL["generic_metal"])
    rho = density(material_class)
    cm_per_s = rate_um_y / (1.0e4 * SEC_PER_YEAR)
    return cm_per_s * z * FARADAY * rho / molar_mass


def safe_log(x: float, floor: float = 1.0e-30) -> float:
    return math.log(max(float(x), floor))


def first_float(text: object, default: float | None = None) -> float | None:
    if text is None:
        return default
    if isinstance(text, (int, float)) and math.isfinite(float(text)):
        return float(text)
    s = str(text).replace(",", "")
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    if m:
        return float(m.group(0))
    return default
