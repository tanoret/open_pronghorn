"""Explicit, measurement-keyed geometry/context inputs for advanced models."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ADVANCED_CONTEXT_COLUMNS: tuple[str, ...] = (
    "context_revision",
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


def merge_advanced_case_context(
    targets: pd.DataFrame,
    context: pd.DataFrame | str | Path,
) -> pd.DataFrame:
    """Merge explicit advanced-model context into targets with strict checks."""
    frame = pd.read_csv(context) if isinstance(context, (str, Path)) else context.copy()
    required = {"measurement_id", *ADVANCED_CONTEXT_COLUMNS}
    missing_columns = sorted(required.difference(frame.columns))
    if missing_columns:
        raise ValueError(f"Advanced case context is missing columns: {missing_columns}")
    if frame["measurement_id"].isna().any() or frame["measurement_id"].duplicated().any():
        raise ValueError("Advanced case context measurement_id values must be non-null and unique")
    overlapping = sorted(set(ADVANCED_CONTEXT_COLUMNS).intersection(targets.columns))
    if overlapping:
        raise ValueError(f"Targets already define advanced context columns: {overlapping}")
    merged = targets.merge(
        frame[["measurement_id", *ADVANCED_CONTEXT_COLUMNS]],
        on="measurement_id",
        how="left",
        validate="many_to_one",
    )
    missing_rows = merged.loc[
        merged[list(ADVANCED_CONTEXT_COLUMNS)].isna().any(axis=1), "measurement_id"
    ].astype(str).tolist()
    if missing_rows:
        raise ValueError(f"Advanced case context has no complete row for: {missing_rows}")
    return merged
