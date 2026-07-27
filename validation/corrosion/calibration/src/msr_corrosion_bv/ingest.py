"""Convert the source workbook into model-ready, traceable calibration rows.

The workbook remains the authoritative source.  This module makes explicit the
otherwise easy-to-miss scientific judgments needed to map heterogeneous
literature statements into material, salt, flow, redox, position, response, and
constraint classes.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .chemistry import (
    HOURS_PER_YEAR,
    canonical_material,
    canonical_salt,
    mg_cm2_h_to_um_y,
    mg_cm2_to_um,
    mil_to_um,
)


REQUIRED_WORKBOOK_SHEETS = {
    "Validation_Cases",
    "Detailed_Measurements",
    "Source_Index",
}


@dataclass(frozen=True)
class ParsedNumber:
    low: float | None
    mid: float | None
    high: float | None
    relation: str
    raw: str


def _clean_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value).strip()


def parse_numeric_interval(value: Any) -> ParsedNumber:
    """Parse sparse values such as '40 to 100', '~5', '<5', '+500'.

    Returns low/mid/high and a relation code: exact, range, upper, lower, qualitative, none.
    """
    raw = _clean_text(value)
    if raw == "":
        return ParsedNumber(None, None, None, "none", raw)
    if isinstance(value, (int, float)) and not (isinstance(value, float) and math.isnan(value)):
        x = float(value)
        return ParsedNumber(x, x, x, "exact", raw)

    s = raw.lower().replace(",", "")
    s = s.replace("~", " about ")
    s = s.replace("approximately", "about")
    s = s.replace("approx.", "about")
    s = s.replace("+", " +")

    # Common qualitative phrases.
    qualitative_phrases = [
        "no corrosion observed",
        "below detection",
        "minimal",
        "no detectable",
        "extensive",
        "metal > graphite",
        "much below",
    ]
    if any(p in s for p in qualitative_phrases):
        nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)]
        if "below" in s and nums:
            x = nums[0]
            return ParsedNumber(None, x / 2.0, x, "upper", raw)
        if "up to" in s and nums:
            x = nums[0]
            return ParsedNumber(None, x / 2.0, x, "upper", raw)
        return ParsedNumber(None, None, None, "qualitative", raw)

    if "up to" in s:
        nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)]
        if nums:
            x = nums[0]
            return ParsedNumber(None, x / 2.0, x, "upper", raw)

    if s.startswith("<") or "less than" in s or "below" in s:
        nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)]
        if nums:
            x = nums[0]
            return ParsedNumber(None, x / 2.0, x, "upper", raw)
    if s.startswith(">") or "greater than" in s or "more than" in s:
        nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)]
        if nums:
            x = nums[0]
            return ParsedNumber(x, x * 1.5, None, "lower", raw)

    # 'a to b', 'a-b' ranges. Avoid interpreting negative signs as ranges.
    range_match = re.search(
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*(?:to|-)\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)",
        s,
    )
    if range_match:
        a = float(range_match.group(1))
        b = float(range_match.group(2))
        low, high = min(a, b), max(a, b)
        if low > 0 and high > 0:
            mid = math.sqrt(low * high)
        else:
            mid = 0.5 * (low + high)
        return ParsedNumber(low, mid, high, "range", raw)

    nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)]
    if nums:
        x = nums[0]
        if "about" in s:
            return ParsedNumber(0.8 * x, x, 1.25 * x, "approx", raw)
        return ParsedNumber(x, x, x, "exact", raw)
    return ParsedNumber(None, None, None, "qualitative", raw)


def parse_time_hours(value: Any) -> tuple[float | None, float | None, float | None, str]:
    raw = _clean_text(value)
    if raw == "":
        return None, None, None, "none"
    s = raw.lower()
    if "last 4 months" in s or "last four months" in s:
        mid = 4.0 * 30.44 * 24.0
        return 0.8 * mid, mid, 1.25 * mid, "approx"
    parsed = parse_numeric_interval(value)
    return parsed.low, parsed.mid, parsed.high, parsed.relation


def parse_temperature_c(value: Any, default: float | None = None) -> float | None:
    raw = _clean_text(value)
    if raw == "":
        return default
    if isinstance(value, (int, float)) and not (isinstance(value, float) and math.isnan(value)):
        return float(value)
    s = raw.lower().replace(",", "")
    nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)]
    if not nums:
        return default
    # For strings such as '550 and 650 static', use the average as the nominal condition.
    if "and" in s and len(nums) >= 2:
        return float(np.mean(nums[:2]))
    return nums[0]


def parse_delta_t_c(value: Any) -> float:
    raw = _clean_text(value)
    if raw == "":
        return 0.0
    # MSRE sometimes reports 40 F core/HX gradient; convert if explicit F is present.
    if "40 f" in raw.lower():
        return 40.0 * 5.0 / 9.0
    x = parse_numeric_interval(raw).mid
    return float(x) if x is not None else 0.0


def flow_class_and_factor(row: pd.Series) -> tuple[str, float]:
    text = " ".join(
        _clean_text(row.get(c))
        for c in [
            "circulation_or_flow",
            "velocity_or_flow_rate",
            "experiment_class",
            "location_or_condition",
            "facility_or_test_id",
        ]
    ).lower()
    if "no salt" in text or "gas exposure" in text or "cell gas" in text:
        return "gas_control", 0.02
    if "off-gas" in text or "helium purge" in text:
        return "offgas", 0.20
    if "turbulent" in text:
        return "turbulent_forced", 2.80
    if "laminar" in text:
        return "laminar_forced", 1.25
    if "1200 gpm" in text or "operating msre" in text or "msre fuel loop" in text or "circulating primary" in text:
        return "reactor_circulation", 2.00
    if "forced" in text or "circulation loop" in text or "molten salt circulation" in text:
        return "forced_circulation", 1.80
    if "natural circulation" in text or "natural convection" in text:
        return "natural_convection", 1.10
    if "thermal convection" in text or "tcl" in text or "density-driven" in text:
        return "thermal_convection", 1.00
    if "static" in text or "capsule" in text or "pot" in text or "crucible" in text:
        return "static", 0.35
    if "electrochemical" in text:
        return "electrochemical", 0.50
    return "unspecified_flow", 0.75


def redox_class(row: pd.Series) -> str:
    text = " ".join(
        _clean_text(row.get(c))
        for c in [
            "redox_impurity_chemistry",
            "salt_purification_initial_condition",
            "location_or_condition",
            "observable",
            "facility_or_test_id",
            "validation_notes",
        ]
    ).lower()
    salt = canonical_salt(row.get("salt_system"))
    case_id = _clean_text(row.get("case_id"))
    observable = _clean_text(row.get("observable")).lower()

    if "no salt" in text or "oxidizing cell atmosphere" in text or "n2 +" in text:
        return "gas_control"

    # Reactor surveillance rows represent integrated MSRE fuel-salt exposure.
    # Mentions of later Be adjustments in chemistry notes should not turn the
    # whole in-reactor corrosion history into a Be-reduced bench test.
    if case_id.startswith("MSRE-SURV"):
        return "msre_or_fuel_baseline"

    # Explicit pre/post redox perturbation handling for NCL-16.
    if "pre-fef2" in text or "before fef2" in text:
        return "purified_baseline"
    if "after fef2" in text or "fef2 addition" in text or "oxidizing additions" in text or " fef2 " in f" {text} ":
        return "oxidizing_fef2"

    # NCL-31 as-received FLiBe-like salt was specifically described as
    # relatively oxidizing/impure, distinct from moisture-bearing chlorides but
    # sharing the same effective oxidizing overpotential class in this sparse fit.
    if "as-received" in text or "relatively oxidizing" in text or "fe and cr impurities" in text or "initial equivalent corrosion" in observable:
        return "impure_moisture"

    if ("be addition" in text or "be metal addition" in text or "be redox" in text or "be additions" in text or "reduced dissolution" in text) and "no be addition" not in text:
        return "reducing_be"
    if "impure" in text or "moisture" in text or "water" in text or "mgcl2.6h2o" in text:
        return "impure_moisture"
    if "tellur" in text or "ni3te2" in text or "li2te" in text:
        return "tellurium"

    # Use stress only when the specific measurement is loaded/stressed, not for
    # companion unloaded/literature rows from the same case.
    if ("loaded" in observable or "stressed" in observable or "stress" in observable) and "unloaded" not in observable:
        return "stressed"
    if "unloaded" in observable:
        return "purified_baseline"
    if "stressed" in text and "unloaded" not in text:
        return "stressed"

    if "multi-alloy" in text or "plating from other alloys" in text or "galvanic" in text:
        return "multi_alloy"
    if "fission" in text or "noble metal" in text or "uf3" in text or "nb behavior" in text:
        return "fission_product"
    if "purified" in text or "high-purity" in text or "baseline" in text:
        return "purified_baseline"
    if salt == "chloride":
        return "chloride_unspecified"
    if "msre" in text or salt == "fluoride_fuel":
        return "msre_or_fuel_baseline"
    return "purified_baseline"


def position_class(row: pd.Series) -> str:
    text = " ".join(
        _clean_text(row.get(c))
        for c in ["location_or_condition", "observable", "facility_or_test_id", "primary_model_use"]
    ).lower()
    if "cold" in text or "radiator" in text:
        return "cold_leg"
    if "hot" in text or "hottest" in text or "under heaters" in text:
        return "hot_leg"
    if "core" in text or "msre" in text:
        return "core_or_reactor"
    if "off-gas" in text or "gas" in text:
        return "gas_or_offgas"
    return "nominal"


def surface_class(row: pd.Series) -> str:
    text = " ".join(
        _clean_text(row.get(c))
        for c in ["material_or_species", "materials_exposed", "location_or_condition", "observable", "value"]
    ).lower()
    if "turbulent" in text and "metal" in text:
        return "turbulent_metal"
    if "laminar" in text and "metal" in text:
        return "laminar_metal"
    if "graphite" in text and "hastelloy" not in text:
        return "graphite"
    if "metal" in text or "hastelloy" in text or "steel" in text or "alloy" in text:
        return "metal"
    return "generic_surface"


def quality_weight(row: pd.Series) -> float:
    q = _clean_text(row.get("data_quality")).upper()
    source = _clean_text(row.get("source_id"))
    status = _clean_text(row.get("extraction_status")).lower()
    weight = {"A": 1.0, "B": 0.8, "C": 0.55, "D": 0.35}.get(q, 0.65)
    if "abstract" in status or "indexed" in status:
        weight *= 0.65
    if "ocr" in status:
        weight *= 0.65
    if source.startswith("DB-"):
        weight *= 0.35
    return weight


def classify_response(row: pd.Series) -> dict[str, Any]:
    obs = _clean_text(row.get("observable")).lower()
    mapping = _clean_text(row.get("model_mapping")).lower()
    units = _clean_text(row.get("units")).lower()
    value = row.get("value")
    material_class = row.get("material_class", canonical_material(row.get("material_or_species")))
    parsed = parse_numeric_interval(value)

    result: dict[str, Any] = {
        "response_kind": "unclassified",
        "fit_role": "input_only",
        "target_low": parsed.low,
        "target_mid": parsed.mid,
        "target_high": parsed.high,
        "target_relation": parsed.relation,
        "target_units_model": None,
        "usage_reason": "Not a direct corrosion/plating response target.",
        "default_sigma_ln": 0.75,
    }

    # Source-term inputs and experimental method descriptors.
    if "effective full-power hours" in obs:
        result.update(response_kind="source_duration_h", fit_role="input_only", target_units_model="h", usage_reason="Used as a fission-product source-term input, not fitted as corrosion response.")
        return result
    if "fef2 addition" in obs:
        result.update(response_kind="redox_perturbation_input", fit_role="input_only", target_units_model="ppm", usage_reason="Used to classify redox perturbation; not a response variable.")
        return result
    if "tracer radionuclides" in obs or "mev protons" in units:
        result.update(response_kind="tracer_method_input", fit_role="input_only", target_units_model="MeV", usage_reason="Describes tracer generation method; no corrosion response magnitude.")
        return result
    if "carbon-rich layer" in obs:
        # Keep in audit, but do not fit the molten-salt BV corrosion/plating parameters.
        low, mid, high = parsed.low, parsed.mid, parsed.high
        if units == "mil":
            low = mil_to_um(low) if low is not None else None
            mid = mil_to_um(mid) if mid is not None else None
            high = mil_to_um(high) if high is not None else None
        result.update(response_kind="carbon_layer_depth_um", fit_role="excluded_auxiliary", target_low=low, target_mid=mid, target_high=high, target_units_model="um", usage_reason="Graphite-contact carbon transport is auxiliary to corrosion/plating BV model.")
        return result

    # Diffusion and redox thresholds.
    if "diffusion coefficient" in obs:
        result.update(response_kind="cr_diffusion_cm2_s", fit_role="direct", target_units_model="cm2/s", usage_reason="Directly calibrates chromium diffusion submodel.", default_sigma_ln=0.35)
        return result
    if "critical redox" in obs or "u(iv)/u(iii)" in _clean_text(row.get("material_or_species")).lower():
        result.update(response_kind="te_redox_threshold_ratio", fit_role="direct", target_units_model="ratio", usage_reason="Calibrates tellurium redox threshold submodel.", default_sigma_ln=0.45)
        return result

    # Salt chemistry inventories.
    if "salt chromium" in obs:
        result.update(response_kind="salt_cr_ppm", fit_role="range" if parsed.relation == "range" else "direct", target_units_model="ppm", usage_reason="Calibrates salt-inventory geometry/source-term closure.", default_sigma_ln=0.70)
        return result
    if "salt iron" in obs:
        # Use absolute magnitude because model predicts depletion sign separately.
        low = abs(parsed.low) if parsed.low is not None else None
        mid = abs(parsed.mid) if parsed.mid is not None else None
        high = abs(parsed.high) if parsed.high is not None else None
        result.update(response_kind="salt_fe_decrease_ppm", fit_role="direct", target_low=low, target_mid=mid, target_high=high, target_units_model="ppm", usage_reason="Species-balance check tied to chromium inventory submodel.", default_sigma_ln=0.85)
        return result
    if "tellurium concentration" in obs:
        result.update(response_kind="te_soluble_ppm", fit_role="upper", target_units_model="ppm", usage_reason="Upper-bound constraint on soluble tellurium activity.", default_sigma_ln=0.80)
        return result

    # Qualitative fission-product/plating constraints.
    if "relative deposition" in obs:
        result.update(response_kind="noble_metal_deposition_ranking", fit_role="ranking", target_units_model="ordinal", usage_reason="Qualitative constraint: metal > graphite and turbulent metal > laminar metal.", default_sigma_ln=1.0)
        return result
    if "off-gas" in obs and "fraction" in obs:
        result.update(response_kind="offgas_fraction_percent", fit_role="upper", target_units_model="percent", usage_reason="Upper-bound constraint on off-gas noble-metal escape fraction.", default_sigma_ln=0.85)
        return result

    # Qualitative no/minimal/attack constraints.
    text_value = _clean_text(value).lower()
    if "no visible corrosion" in obs or "no corrosion observed" in text_value:
        result.update(response_kind="corrosion_depth_um", fit_role="upper", target_low=None, target_mid=12.5, target_high=25.0, target_relation="upper", target_units_model="um", usage_reason="No visible corrosion represented as an upper-bound depth constraint.", default_sigma_ln=0.85)
        return result
    if "minimal" in text_value:
        result.update(response_kind="corrosion_depth_um", fit_role="upper", target_low=None, target_mid=2.5, target_high=5.0, target_relation="upper", target_units_model="um", usage_reason="Minimal attack represented as an upper-bound depth constraint.", default_sigma_ln=0.85)
        return result
    if "below detection" in text_value and "weight" in obs:
        result.update(response_kind="mass_loss_mg_cm2", fit_role="upper", target_low=None, target_mid=0.05, target_high=0.10, target_relation="upper", target_units_model="mg/cm2", usage_reason="Below-detection weight change represented as an upper-bound mass-change constraint.", default_sigma_ln=0.90)
        return result
    if "no detectable weight" in text_value:
        result.update(response_kind="mass_loss_mg_cm2", fit_role="upper", target_low=None, target_mid=0.05, target_high=0.10, target_relation="upper", target_units_model="mg/cm2", usage_reason="No detectable weight change represented as an upper bound.", default_sigma_ln=0.90)
        return result
    if "extensive" in text_value and ("igc" in obs or "cracking" in obs):
        result.update(response_kind="igc_depth_um", fit_role="lower", target_low=50.0, target_mid=100.0, target_high=None, target_relation="lower", target_units_model="um", usage_reason="Extensive IGC represented as a lower-bound damage-depth constraint.", default_sigma_ln=0.95)
        return result

    # Numeric rates/depths/mass changes.
    low, mid, high = parsed.low, parsed.mid, parsed.high
    relation = parsed.relation
    fit_role = "range" if relation == "range" else ("upper" if relation == "upper" else ("lower" if relation == "lower" else "direct"))

    if units == "mil/year":
        low = mil_to_um(low) if low is not None else None
        mid = mil_to_um(mid) if mid is not None else None
        high = mil_to_um(high) if high is not None else None
        result.update(response_kind="corrosion_rate_um_y", fit_role=fit_role, target_low=low, target_mid=mid, target_high=high, target_relation=relation, target_units_model="um/year", usage_reason="Direct corrosion-rate validation target.", default_sigma_ln=0.55)
        return result
    if units == "mg/cm2 h":
        low = mg_cm2_h_to_um_y(low, material_class) if low is not None else None
        mid = mg_cm2_h_to_um_y(mid, material_class) if mid is not None else None
        high = mg_cm2_h_to_um_y(high, material_class) if high is not None else None
        result.update(response_kind="corrosion_rate_um_y", fit_role=fit_role, target_low=low, target_mid=mid, target_high=high, target_relation=relation, target_units_model="um/year", usage_reason="Areal mass-loss rate converted to equivalent penetration rate.", default_sigma_ln=0.65)
        return result
    if units == "mg/cm2":
        kind = "mass_gain_mg_cm2" if "gain" in obs or "deposit" in mapping or "cold" in _clean_text(row.get("location_or_condition")).lower() else "mass_loss_mg_cm2"
        result.update(response_kind=kind, fit_role=fit_role, target_low=low, target_mid=mid, target_high=high, target_relation=relation, target_units_model="mg/cm2", usage_reason="Integrated areal mass-change target.", default_sigma_ln=0.55)
        return result
    if units == "mil":
        low = mil_to_um(low) if low is not None else None
        mid = mil_to_um(mid) if mid is not None else None
        high = mil_to_um(high) if high is not None else None
        kind = "igc_depth_um" if "crack" in obs or "void" in obs or "defect" in obs else "corrosion_depth_um"
        result.update(response_kind=kind, fit_role=fit_role, target_low=low, target_mid=mid, target_high=high, target_relation=relation, target_units_model="um", usage_reason="Attack/depletion depth target converted from mil.", default_sigma_ln=0.65)
        return result
    if units == "um":
        kind = "igc_depth_um" if "igc" in obs or "crack" in obs or "void" in obs or "damage" in obs else "corrosion_depth_um"
        result.update(response_kind=kind, fit_role=fit_role, target_low=low, target_mid=mid, target_high=high, target_relation=relation, target_units_model="um", usage_reason="Attack/depletion depth target.", default_sigma_ln=0.70)
        return result
    if units == "times baseline":
        result.update(response_kind="redox_acceleration_ratio", fit_role=fit_role, target_low=low, target_mid=mid, target_high=high, target_relation=relation, target_units_model="ratio", usage_reason="Calibrates oxidizing-redox acceleration relative to pre-FeF2 baseline.", default_sigma_ln=0.45)
        return result
    if units == "qualitative ratio" and "450 h change" in text_value:
        result.update(response_kind="redox_acceleration_qualitative", fit_role="lower", target_low=6.0, target_mid=10.0, target_high=None, target_relation="lower", target_units_model="ratio", usage_reason="Narrative says 450 h after FeF2 matched roughly 10000 h baseline; represented as a lower-bound acceleration.", default_sigma_ln=0.90)
        return result

    # Ranges in the stress-free literature are useful even with limited conditions.
    if units in ["qualitative", "percent", "ratio"]:
        result.update(response_kind="qualitative_constraint", fit_role="qualitative", target_units_model=units, usage_reason="Qualitative validation-only constraint.")
        return result

    return result


def build_model_tables(workbook_path: str | Path) -> dict[str, pd.DataFrame]:
    """Load and normalize the local calibration workbook.

    Returns the original case, measurement, and source sheets together with the
    normalized targets, all-case feature matrix, and data-use summary.  A clear
    error is raised if the workbook is missing a required sheet or a measurement
    references a case that is absent from ``Validation_Cases``.
    """
    workbook_path = Path(workbook_path)
    if not workbook_path.is_file():
        raise FileNotFoundError(f"Calibration workbook does not exist: {workbook_path}")

    workbook = pd.ExcelFile(workbook_path)
    missing_sheets = sorted(REQUIRED_WORKBOOK_SHEETS - set(workbook.sheet_names))
    if missing_sheets:
        raise ValueError(
            f"Calibration workbook {workbook_path} is missing required sheets: "
            + ", ".join(missing_sheets)
        )

    validation_cases = pd.read_excel(workbook, sheet_name="Validation_Cases")
    measurements = pd.read_excel(workbook, sheet_name="Detailed_Measurements")
    source_index = pd.read_excel(workbook, sheet_name="Source_Index")

    # Merge case context into measurement rows.
    merged = measurements.merge(validation_cases, on="case_id", how="left", suffixes=("", "_case"))
    known_cases = set(validation_cases["case_id"].dropna())
    orphan_measurements = sorted(
        set(measurements.loc[~measurements["case_id"].isin(known_cases), "case_id"].dropna())
    )
    if orphan_measurements:
        raise ValueError(
            "Detailed_Measurements contains case IDs absent from Validation_Cases: "
            + ", ".join(str(case_id) for case_id in orphan_measurements)
        )

    feature_rows: list[dict[str, Any]] = []
    for _, row in merged.iterrows():
        row = row.copy()
        mat_source = row.get("material_or_species") if _clean_text(row.get("material_or_species")) else row.get("materials_exposed")
        material_class = canonical_material(mat_source)
        salt_class = canonical_salt(row.get("salt_system"))
        t_meas_low, t_meas_mid, t_meas_high, t_rel = parse_time_hours(row.get("time_h"))
        t_case_low, t_case_mid, t_case_high, _ = parse_time_hours(row.get("exposure_h"))
        time_low = t_meas_low if t_meas_low is not None else t_case_low
        time_h = t_meas_mid if t_meas_mid is not None else t_case_mid
        time_high = t_meas_high if t_meas_high is not None else t_case_high
        T = row.get("temperature_C")
        if T is None or (isinstance(T, float) and math.isnan(T)):
            T = parse_temperature_c(row.get("temperature_C_hot_or_nominal"), default=None)
        else:
            T = float(T)
        if T is None:
            # Sparse modern cases with '650 to 750 context' are parsed from the case field.
            T = parse_temperature_c(row.get("temperature_C_hot_or_nominal"), default=650.0)
        T_cold = parse_temperature_c(row.get("temperature_C_cold"), default=None)
        delta_T = parse_delta_t_c(row.get("delta_T_C"))
        if delta_T == 0.0 and T is not None and T_cold is not None:
            delta_T = max(float(T) - float(T_cold), 0.0)
        f_class, f_factor = flow_class_and_factor(row)
        p_class = position_class(row)
        r_class = redox_class(row)
        surf_class = surface_class(row)
        row["material_class"] = material_class
        response = classify_response(row)
        weight = quality_weight(row)
        if row.get("validation_priority") == "Critical":
            weight *= 1.15
        elif row.get("validation_priority") == "Medium":
            weight *= 0.85
        elif row.get("validation_priority") == "Low":
            weight *= 0.65

        # Wider uncertainty for OCR/abstract extracted values.
        sigma = response.get("default_sigma_ln", 0.75)
        status = _clean_text(row.get("extraction_status")).lower()
        note = _clean_text(row.get("uncertainty_or_note")).lower()
        if "ocr" in status or "verify" in note:
            sigma += 0.35
        if "abstract" in status or "indexed" in status or "search-result" in _clean_text(row.get("key_numeric_outputs_extracted")).lower():
            sigma += 0.25
        if response.get("target_relation") == "range" and response.get("target_low") and response.get("target_high"):
            low = max(float(response["target_low"]), 1.0e-30)
            high = max(float(response["target_high"]), low * 1.01)
            sigma = max(sigma, 0.5 * math.log(high / low) + 0.20)

        feature_rows.append(
            {
                "measurement_id": row.get("measurement_id"),
                "case_id": row.get("case_id"),
                "source_id": row.get("source_id"),
                "validation_priority": row.get("validation_priority"),
                "experiment_family": row.get("experiment_family"),
                "experiment_class": row.get("experiment_class"),
                "observable": row.get("observable"),
                "material_or_species": row.get("material_or_species"),
                "location_or_condition": row.get("location_or_condition"),
                "raw_value": row.get("value"),
                "raw_units": row.get("units"),
                "material_class": material_class,
                "salt_class": salt_class,
                "flow_class": f_class,
                "flow_factor": f_factor,
                "position_class": p_class,
                "redox_class": r_class,
                "surface_class": surf_class,
                "temperature_C": T,
                "temperature_K": float(T) + 273.15 if T is not None else np.nan,
                "time_h_low": time_low,
                "time_h": time_h,
                "time_h_high": time_high,
                "time_years": float(time_h) / HOURS_PER_YEAR if time_h is not None else np.nan,
                "delta_T_C": delta_T,
                "temperature_cold_C": T_cold,
                "source_url": row.get("source_url"),
                "data_quality": row.get("data_quality"),
                "extraction_status": row.get("extraction_status"),
                "quality_weight": weight,
                "model_mapping": row.get("model_mapping"),
                **response,
            }
        )

    targets = pd.DataFrame(feature_rows)

    # Apply final role corrections where model can evaluate only after features are available.
    missing_time = targets["time_h"].isna()
    time_needed = targets["response_kind"].isin(
        ["corrosion_depth_um", "igc_depth_um", "mass_loss_mg_cm2", "mass_gain_mg_cm2", "salt_cr_ppm", "salt_fe_decrease_ppm"]
    )
    targets.loc[missing_time & time_needed & targets["fit_role"].isin(["direct", "range", "upper", "lower"]), "fit_role"] = "prediction_only"
    targets.loc[missing_time & time_needed, "usage_reason"] = targets.loc[missing_time & time_needed, "usage_reason"].astype(str) + " Time is not sufficiently specified, so row is prediction-only until digitized."

    # Add model-readable case features for all validation cases, including those without detailed target rows.
    case_rows: list[dict[str, Any]] = []
    for _, row in validation_cases.iterrows():
        T = parse_temperature_c(row.get("temperature_C_hot_or_nominal"), default=650.0)
        T_cold = parse_temperature_c(row.get("temperature_C_cold"), default=None)
        delta_T = parse_delta_t_c(row.get("delta_T_C"))
        if delta_T == 0 and T is not None and T_cold is not None:
            delta_T = max(T - T_cold, 0.0)
        t_low, t_mid, t_high, _ = parse_time_hours(row.get("exposure_h"))
        f_class, f_factor = flow_class_and_factor(row)
        r_class = redox_class(row)
        p_class = position_class(row)
        mat = canonical_material(row.get("materials_exposed"))
        salt = canonical_salt(row.get("salt_system"))
        case_rows.append(
            {
                "case_id": row.get("case_id"),
                "source_id": row.get("source_id"),
                "validation_priority": row.get("validation_priority"),
                "experiment_family": row.get("experiment_family"),
                "experiment_class": row.get("experiment_class"),
                "primary_model_use": row.get("primary_model_use"),
                "material_class": mat,
                "salt_class": salt,
                "flow_class": f_class,
                "flow_factor": f_factor,
                "position_class": p_class,
                "redox_class": r_class,
                "surface_class": surface_class(row),
                "temperature_C": T,
                "temperature_K": T + 273.15 if T is not None else np.nan,
                "temperature_cold_C": T_cold,
                "delta_T_C": delta_T,
                "time_h": t_mid,
                "time_years": t_mid / HOURS_PER_YEAR if t_mid is not None else np.nan,
                "source_url": row.get("source_url"),
                "extraction_status": row.get("extraction_status"),
                "data_quality": row.get("data_quality"),
            }
        )
    case_features = pd.DataFrame(case_rows)

    usage = (
        targets.groupby(["fit_role", "response_kind"], dropna=False)
        .size()
        .reset_index(name="n_measurements")
        .sort_values(["fit_role", "response_kind"])
    )

    return {
        "validation_cases": validation_cases,
        "measurements": measurements,
        "source_index": source_index,
        "targets": targets,
        "case_features": case_features,
        "usage_summary": usage,
    }
