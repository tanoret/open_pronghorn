"""Publication-quality validation and dynamic-model figures.

All figures are generated with Matplotlib and exported as 600-dpi PNG plus
editable SVG and publication-ready vector PDF.  The palette is color-vision-deficiency friendly and the plotting
functions consume the checked validation tables so every point is auditable.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

from .dynamic_network import DynamicRedoxInventoryDepletionModel, ELEMENTS

MODEL_ORDER = ("reduced_mechanistic", "thermochemical", "dynamic_network")
MODEL_LABELS = {
    "effective": "Effective BV reference",
    "reduced_mechanistic": "Reduced mechanistic",
    "thermochemical": "MSTDB-TC thermochemical",
    "dynamic_network": "DRIDN dynamic simulation",
}
MODEL_COLORS = {
    "effective": "#6B6B6B",
    "reduced_mechanistic": "#0072B2",
    "thermochemical": "#D55E00",
    "dynamic_network": "#009E73",
}
MODEL_MARKERS = {
    "effective": "D",
    "reduced_mechanistic": "o",
    "thermochemical": "^",
    "dynamic_network": "s",
}
SPECIES_COLORS = {"Cr": "#7B3294", "Fe": "#E69F00", "Ni": "#56B4E9"}


def apply_publication_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.0,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.0,
            "axes.linewidth": 0.8,
            "legend.fontsize": 7.2,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.minor.width": 0.55,
            "ytick.minor.width": 0.55,
            "lines.linewidth": 1.35,
            "lines.markersize": 4.5,
            "legend.frameon": False,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> list[Path]:
    """Write one figure in raster and two vector publication formats."""
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = [
        figure_dir / f"{stem}.png",
        figure_dir / f"{stem}.svg",
        figure_dir / f"{stem}.pdf",
    ]
    fig.savefig(paths[0], dpi=600, facecolor="white")
    fig.savefig(paths[1], facecolor="white")
    fig.savefig(paths[2], facecolor="white")
    plt.close(fig)
    return paths


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.5,
        fontweight="bold",
    )


def _response_group(kind: str) -> str:
    mapping = {
        "corrosion_rate_um_y": "Corrosion rate",
        "corrosion_depth_um": "Selective depletion",
        "igc_depth_um": "Localized damage",
        "mass_loss_mg_cm2": "Areal mass change",
        "mass_gain_mg_cm2": "Areal mass change",
        "salt_cr_ppm": "Salt inventory",
        "salt_fe_decrease_ppm": "Salt inventory",
        "redox_acceleration_ratio": "Redox / diffusion",
        "redox_acceleration_qualitative": "Redox / diffusion",
        "cr_diffusion_cm2_s": "Redox / diffusion",
    }
    return mapping.get(str(kind), "Other")


def _direct_range(comparison: pd.DataFrame) -> pd.DataFrame:
    frame = comparison[comparison["fit_role"].isin(["direct", "range"])].copy()
    frame = frame[pd.to_numeric(frame["target_mid"], errors="coerce").gt(0)]
    frame["response_group"] = frame["response_kind"].map(_response_group)
    order = {
        "Corrosion rate": 0,
        "Selective depletion": 1,
        "Localized damage": 2,
        "Areal mass change": 3,
        "Salt inventory": 4,
        "Redox / diffusion": 5,
        "Other": 6,
    }
    frame["_group_order"] = frame["response_group"].map(order).fillna(99)
    frame["_id_order"] = frame["measurement_id"].astype(str).str.extract(r"(\d+)")[0].astype(int)
    return frame.sort_values(["_group_order", "_id_order"]).reset_index(drop=True)


def _group_boundaries(frame: pd.DataFrame) -> tuple[list[float], list[tuple[float, str]]]:
    boundaries: list[float] = []
    labels: list[tuple[float, str]] = []
    start = 0
    for group, subset in frame.groupby("response_group", sort=False):
        end = start + len(subset)
        labels.append(((start + end - 1) / 2.0, str(group)))
        if end < len(frame):
            boundaries.append(end - 0.5)
        start = end
    return boundaries, labels


def figure_validation_ratios(comparison: pd.DataFrame, figure_dir: Path) -> list[Path]:
    frame = _direct_range(comparison)
    x = np.arange(len(frame), dtype=float)
    target = frame["target_mid"].astype(float).to_numpy()
    fig, ax = plt.subplots(figsize=(7.6, 4.5), constrained_layout=True)
    ax.axhspan(0.5, 2.0, color="#D9D9D9", alpha=0.32, zorder=0, label="factor-of-two band")
    ax.axhline(1.0, color="black", linewidth=1.0, zorder=1)

    low = pd.to_numeric(frame["target_low"], errors="coerce").to_numpy(dtype=float)
    high = pd.to_numeric(frame["target_high"], errors="coerce").to_numpy(dtype=float)
    for index in range(len(frame)):
        if frame.loc[index, "fit_role"] == "range" and np.isfinite(low[index]) and np.isfinite(high[index]):
            ax.vlines(
                x[index],
                low[index] / target[index],
                high[index] / target[index],
                color="black",
                linewidth=2.0,
                alpha=0.75,
                zorder=2,
            )
            ax.plot(x[index], 1.0, marker="_", color="black", markersize=6, zorder=3)

    offsets = (-0.20, 0.0, 0.20)
    for model, offset in zip(MODEL_ORDER, offsets):
        ratios = frame[f"{model}_prediction"].astype(float).to_numpy() / target
        ax.scatter(
            x + offset,
            ratios,
            s=24,
            marker=MODEL_MARKERS[model],
            facecolor=MODEL_COLORS[model] if model != "reduced_mechanistic" else "white",
            edgecolor=MODEL_COLORS[model],
            linewidth=0.9,
            label=MODEL_LABELS[model],
            zorder=4,
        )

    boundaries, group_labels = _group_boundaries(frame)
    for boundary in boundaries:
        ax.axvline(boundary, color="#BDBDBD", linewidth=0.65, zorder=0)
    for midpoint, label in group_labels:
        ax.text(
            midpoint,
            1.02,
            label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=6.6,
        )

    ax.set_yscale("log")
    finite = []
    for model in MODEL_ORDER:
        finite.extend((frame[f"{model}_prediction"].astype(float) / frame["target_mid"].astype(float)).tolist())
    ymin = max(0.12, min(finite) / 1.35)
    ymax = min(8.0, max(finite) * 1.35)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(-0.65, len(frame) - 0.35)
    ax.set_ylabel("Prediction / experimental target midpoint")
    ax.set_xlabel("Validation measurement")
    ax.set_xticks(x)
    ax.set_xticklabels(frame["measurement_id"].astype(str), rotation=58, ha="right")
    ax.legend(loc="lower left", ncol=2, columnspacing=1.0, handletextpad=0.45)
    ax.set_title("Like-for-like validation across the three mechanistic approaches", pad=18)
    _panel_label(ax, "a")
    return _save_figure(fig, figure_dir, "fig10_three_model_validation_ratios")


def _add_parity_guides(ax: plt.Axes, values: np.ndarray) -> None:
    positive = values[np.isfinite(values) & (values > 0)]
    if not len(positive):
        return
    lo = 10 ** math.floor(math.log10(positive.min() / 1.35))
    hi = 10 ** math.ceil(math.log10(positive.max() * 1.35))
    xx = np.geomspace(lo, hi, 200)
    ax.fill_between(xx, 0.5 * xx, 2.0 * xx, color="#D9D9D9", alpha=0.28, zorder=0)
    ax.plot(xx, xx, color="black", linewidth=0.9, zorder=1)
    ax.plot(xx, 0.5 * xx, color="#8C8C8C", linewidth=0.65, linestyle="--", zorder=1)
    ax.plot(xx, 2.0 * xx, color="#8C8C8C", linewidth=0.65, linestyle="--", zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)


def figure_parity_by_response(comparison: pd.DataFrame, figure_dir: Path) -> list[Path]:
    frame = _direct_range(comparison)
    groups = [
        ("Corrosion rate", ["corrosion_rate_um_y"], r"Experimental rate ($\mu$m y$^{-1}$)", r"Predicted rate ($\mu$m y$^{-1}$)"),
        ("Selective depletion", ["corrosion_depth_um"], r"Experimental depth ($\mu$m)", r"Predicted depth ($\mu$m)"),
        ("Localized damage", ["igc_depth_um"], r"Experimental depth ($\mu$m)", r"Predicted depth ($\mu$m)"),
        ("Areal mass change", ["mass_loss_mg_cm2", "mass_gain_mg_cm2"], r"Experimental mass change (mg cm$^{-2}$)", r"Predicted mass change (mg cm$^{-2}$)"),
        ("Salt inventory", ["salt_cr_ppm", "salt_fe_decrease_ppm"], "Experimental inventory (ppm)", "Predicted inventory (ppm)"),
        ("Redox and diffusion", ["redox_acceleration_ratio", "redox_acceleration_qualitative", "cr_diffusion_cm2_s"], "Normalized experimental target", "Normalized prediction"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(7.6, 5.8), constrained_layout=False)
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.10, top=0.87, wspace=0.30, hspace=0.36)
    letters = "abcdef"
    for ax, (title, kinds, xlabel, ylabel), letter in zip(axes.flat, groups, letters):
        subset = frame[frame["response_kind"].isin(kinds)].copy()
        if subset.empty:
            ax.set_visible(False)
            continue
        target = subset["target_mid"].astype(float).to_numpy()
        normalized = title == "Redox and diffusion"
        all_values: list[float] = []
        for model in MODEL_ORDER:
            prediction = subset[f"{model}_prediction"].astype(float).to_numpy()
            xx = np.ones_like(target) if normalized else target
            yy = prediction / target if normalized else prediction
            all_values.extend(xx.tolist())
            all_values.extend(yy.tolist())
            ax.scatter(
                xx,
                yy,
                s=25,
                marker=MODEL_MARKERS[model],
                facecolor=MODEL_COLORS[model] if model != "reduced_mechanistic" else "white",
                edgecolor=MODEL_COLORS[model],
                linewidth=0.9,
                label=MODEL_LABELS[model],
                zorder=3,
            )
        _add_parity_guides(ax, np.asarray(all_values, dtype=float))
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        _panel_label(ax, letter)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.985))
    return _save_figure(fig, figure_dir, "fig11_three_model_parity_by_response")


def figure_residual_heatmap(comparison: pd.DataFrame, figure_dir: Path) -> list[Path]:
    frame = _direct_range(comparison)
    matrix = np.asarray(
        [
            frame[f"{model}_range_aware_log_error"].astype(float).to_numpy() / math.log(2.0)
            for model in MODEL_ORDER
        ]
    )
    vmax = max(1.0, float(np.nanmax(np.abs(matrix))))
    vmax = min(2.5, math.ceil(vmax * 2.0) / 2.0)
    fig, ax = plt.subplots(figsize=(7.8, 2.4), constrained_layout=True)
    image = ax.imshow(
        matrix,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
        interpolation="nearest",
    )
    ax.set_yticks(np.arange(len(MODEL_ORDER)))
    ax.set_yticklabels([MODEL_LABELS[model] for model in MODEL_ORDER])
    ax.set_xticks(np.arange(len(frame)))
    ax.set_xticklabels(frame["measurement_id"].astype(str), rotation=58, ha="right")
    boundaries, _ = _group_boundaries(frame)
    for boundary in boundaries:
        ax.axvline(boundary, color="black", linewidth=0.65)
    for row_index, model in enumerate(MODEL_ORDER):
        passed = frame[f"{model}_constraint_pass"].astype(bool).to_numpy()
        for column_index, is_pass in enumerate(passed):
            if not is_pass:
                ax.plot(column_index, row_index, marker="x", color="black", markersize=5.5, mew=1.0)
    colorbar = fig.colorbar(image, ax=ax, pad=0.015, fraction=0.035)
    colorbar.set_label(r"Interval-aware error, $\log_2$ scale")
    ax.set_title("Case-resolved interval-aware validation residuals")
    _panel_label(ax, "b")
    return _save_figure(fig, figure_dir, "fig12_three_model_residual_heatmap")


def figure_performance_summary(metrics: Mapping[str, Any], figure_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0), constrained_layout=True)
    x = np.arange(len(MODEL_ORDER))
    panels = [
        ("log_rmse_direct", "Midpoint log RMSE", "lower is better", "a"),
        ("range_aware_log_rmse_direct", "Interval-aware log RMSE", "lower is better", "b"),
        ("range_aware_median_factor_error_direct", "Median factor error", "lower is better", "c"),
        ("constraint_pass_fraction", "Active-constraint pass fraction", "higher is better", "d"),
    ]
    for ax, (suffix, title, subtitle, letter) in zip(axes.flat, panels):
        values = [float(metrics[f"{model}_{suffix}"]) for model in MODEL_ORDER]
        bars = ax.bar(
            x,
            values,
            width=0.68,
            color=[MODEL_COLORS[model] for model in MODEL_ORDER],
            edgecolor="black",
            linewidth=0.45,
        )
        reference = float(metrics[f"effective_{suffix}"])
        ax.axhline(reference, color=MODEL_COLORS["effective"], linestyle="--", linewidth=1.0)
        ax.text(
            0.99,
            reference,
            " effective BV",
            color=MODEL_COLORS["effective"],
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="bottom",
            fontsize=6.8,
        )
        for bar, value in zip(bars, values):
            text = f"{100.0 * value:.1f}%" if suffix == "constraint_pass_fraction" else f"{value:.3f}"
            ax.text(bar.get_x() + bar.get_width() / 2.0, value, text, ha="center", va="bottom", fontsize=7.0)
        ax.set_xticks(x)
        ax.set_xticklabels(["Reduced\nmechanistic", "MSTDB-TC", "DRIDN"])
        ax.set_title(f"{title}\n{subtitle}")
        if suffix == "constraint_pass_fraction":
            ax.set_ylim(0.82, 1.02)
            ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
        else:
            ax.set_ylim(0.0, max(max(values), reference) * 1.25)
        _panel_label(ax, letter)
    return _save_figure(fig, figure_dir, "fig13_model_performance_summary")


def _find_target(targets: pd.DataFrame, measurement_id: str) -> pd.Series:
    selected = targets[targets["measurement_id"].astype(str) == measurement_id]
    if selected.empty:
        raise KeyError(f"Validation measurement {measurement_id} not found")
    return selected.iloc[0]


def figure_dynamic_trajectories(
    dynamic_model: DynamicRedoxInventoryDepletionModel,
    targets: pd.DataFrame,
    figure_dir: Path,
    results_dir: Path,
) -> list[Path]:
    case_ids = ("M-003", "M-014", "M-018", "M-030")
    simulations = {
        measurement_id: dynamic_model.simulate(
            _find_target(targets, measurement_id),
            return_trajectory=True,
        )
        for measurement_id in case_ids
    }
    trajectory_frames: list[pd.DataFrame] = []
    for measurement_id, result in simulations.items():
        frame = result.trajectory.copy() if result.trajectory is not None else pd.DataFrame()
        frame.insert(0, "measurement_id", measurement_id)
        trajectory_frames.append(frame)
    combined = pd.concat(trajectory_frames, ignore_index=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    combined.to_csv(results_dir / "dynamic_network_representative_trajectories.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(7.6, 5.7), constrained_layout=True)

    # M-003: selective depletion and surface composition.
    ax = axes[0, 0]
    trajectory = simulations["M-003"].trajectory
    hours = trajectory["time_years"] * 365.25 * 24.0
    ax.plot(hours, trajectory["front_depth_um"], color=MODEL_COLORS["dynamic_network"], label="reaction-front depth")
    target = _find_target(targets, "M-003")
    ax.axhline(float(target["target_mid"]), color="black", linestyle="--", linewidth=0.8, label="measured depth")
    ax.set_xlabel("Exposure time (h)")
    ax.set_ylabel(r"Reaction-front depth ($\mu$m)")
    twin = ax.twinx()
    twin.plot(hours, trajectory["surface_Cr_availability"], color=SPECIES_COLORS["Cr"], label="surface Cr availability")
    twin.set_ylabel("Normalized surface Cr availability", color=SPECIES_COLORS["Cr"])
    twin.tick_params(axis="y", colors=SPECIES_COLORS["Cr"])
    ax.set_title("M-003: MSRE selective Cr depletion")
    _panel_label(ax, "a")

    # M-014: salt inventories.
    ax = axes[0, 1]
    trajectory = simulations["M-014"].trajectory
    hours = trajectory["time_years"] * 365.25 * 24.0
    for element in ELEMENTS:
        ax.plot(hours, trajectory[f"dissolved_{element}_ppm"], color=SPECIES_COLORS[element], label=element)
    target = _find_target(targets, "M-014")
    ax.axhline(float(target["target_mid"]), color="black", linestyle="--", linewidth=0.8, label="measured Cr")
    ax.set_xlabel("Exposure time (h)")
    ax.set_ylabel("Dissolved inventory (ppm)")
    ax.set_title("M-014: NCL-16 dissolved-metal inventory")
    ax.legend(ncol=2, loc="best")
    _panel_label(ax, "b")

    # M-018: redox perturbation and instantaneous rate.
    ax = axes[1, 0]
    trajectory = simulations["M-018"].trajectory
    hours = trajectory["time_years"] * 365.25 * 24.0
    ax.plot(hours, trajectory["instantaneous_corrosion_rate_um_y"], color=MODEL_COLORS["dynamic_network"])
    baseline_result = dynamic_model.simulate(_find_target(targets, "M-018"), redox_override="purified_baseline")
    ax.axhline(baseline_result.corrosion_rate_um_y, color="black", linestyle="--", linewidth=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("Exposure time (h)")
    ax.set_ylabel(r"Instantaneous rate ($\mu$m y$^{-1}$)")
    twin = ax.twinx()
    twin.plot(hours, trajectory["redox_log_shift"], color="#CC79A7")
    twin.set_ylabel("Redox log-driving-force shift", color="#CC79A7")
    twin.tick_params(axis="y", colors="#CC79A7")
    ax.set_title("M-018: FeF$_2$ redox perturbation")
    _panel_label(ax, "c")

    # M-030: species-resolved cold-leg deposits.
    ax = axes[1, 1]
    trajectory = simulations["M-030"].trajectory
    hours = trajectory["time_years"] * 365.25 * 24.0
    for element in ELEMENTS:
        ax.plot(hours, trajectory[f"deposit_{element}_mg_cm2"], color=SPECIES_COLORS[element], label=element)
    total = sum(trajectory[f"deposit_{element}_mg_cm2"] for element in ELEMENTS)
    ax.plot(hours, total, color="black", linewidth=1.1, label="total")
    target = _find_target(targets, "M-030")
    ax.axhline(float(target["target_mid"]), color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Exposure time (h)")
    ax.set_ylabel(r"Cumulative deposit (mg cm$^{-2}$)")
    ax.set_title("M-030: 316H FLiNaK cold-leg deposition")
    ax.legend(ncol=2, loc="best")
    _panel_label(ax, "d")

    return _save_figure(fig, figure_dir, "fig14_dridn_representative_trajectories")


def figure_species_partitioning(
    dynamic_model: DynamicRedoxInventoryDepletionModel,
    targets: pd.DataFrame,
    figure_dir: Path,
    results_dir: Path,
) -> list[Path]:
    thermochemical = dynamic_model.thermochemical_model
    dissolution_ids = ("M-003", "M-029", "M-038")
    deposition_ids = ("M-012", "M-030")
    records: list[dict[str, Any]] = []

    def append_fraction(measurement_id: str, process: str, model: str, fractions: Mapping[str, float]) -> None:
        for element in ELEMENTS:
            records.append(
                {
                    "measurement_id": measurement_id,
                    "process": process,
                    "model": model,
                    "element": element,
                    "fraction": float(fractions[element]),
                }
            )

    for measurement_id in dissolution_ids:
        row = _find_target(targets, measurement_id)
        append_fraction(measurement_id, "dissolution source", "MSTDB-TC", thermochemical.species_flux_fractions(row))
        append_fraction(
            measurement_id,
            "dissolution source",
            "DRIDN",
            dynamic_model.simulate(row).cumulative_source_fractions,
        )
    for measurement_id in deposition_ids:
        row = _find_target(targets, measurement_id)
        append_fraction(measurement_id, "cold-leg deposit", "MSTDB-TC", thermochemical.deposition_species_fractions(row))
        append_fraction(measurement_id, "cold-leg deposit", "DRIDN", dynamic_model.simulate(row).deposit_fractions)

    table = pd.DataFrame(records)
    results_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(results_dir / "thermochemical_dynamic_species_partitions.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6), constrained_layout=True)
    for ax, process, measurement_ids, title, letter in (
        (axes[0], "dissolution source", dissolution_ids, "Dissolved-source composition", "a"),
        (axes[1], "cold-leg deposit", deposition_ids, "Cold-leg deposit composition", "b"),
    ):
        labels: list[str] = []
        values: list[list[float]] = []
        for measurement_id in measurement_ids:
            for model in ("MSTDB-TC", "DRIDN"):
                subset = table[
                    (table["measurement_id"] == measurement_id)
                    & (table["process"] == process)
                    & (table["model"] == model)
                ]
                labels.append(f"{measurement_id}\n{model}")
                values.append(
                    [
                        float(subset.loc[subset["element"] == element, "fraction"].iloc[0])
                        for element in ELEMENTS
                    ]
                )
        array = np.asarray(values, dtype=float)
        bottom = np.zeros(len(labels))
        for index, element in enumerate(ELEMENTS):
            ax.bar(
                np.arange(len(labels)),
                array[:, index],
                bottom=bottom,
                color=SPECIES_COLORS[element],
                edgecolor="white",
                linewidth=0.4,
                label=element,
            )
            bottom += array[:, index]
        ax.set_ylim(0.0, 1.0)
        ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylabel("Element mass fraction")
        ax.set_title(title)
        _panel_label(ax, letter)
    axes[1].legend(loc="upper right", ncol=3)
    return _save_figure(fig, figure_dir, "fig15_species_partitioning_static_vs_dynamic")


def _constraint_margin_log2(row: pd.Series, prediction: float) -> float:
    """Return a signed validation margin; non-negative values satisfy the criterion."""
    prediction = max(float(prediction), 1.0e-300)
    role = str(row.get("fit_role", ""))
    low = pd.to_numeric(pd.Series([row.get("target_low")]), errors="coerce").iloc[0]
    mid = pd.to_numeric(pd.Series([row.get("target_mid")]), errors="coerce").iloc[0]
    high = pd.to_numeric(pd.Series([row.get("target_high")]), errors="coerce").iloc[0]
    log_two = math.log(2.0)
    if role == "direct" and np.isfinite(mid) and mid > 0:
        return 1.0 - abs(math.log(prediction / float(mid))) / log_two
    if role == "range" and np.isfinite(low) and low > 0 and np.isfinite(high) and high > 0:
        return min(math.log(prediction / float(low)), math.log(float(high) / prediction)) / log_two
    if role == "upper":
        limit = high if np.isfinite(high) and high > 0 else mid
        return math.log(float(limit) / prediction) / log_two if np.isfinite(limit) and limit > 0 else np.nan
    if role == "lower":
        limit = low if np.isfinite(low) and low > 0 else mid
        return math.log(prediction / float(limit)) / log_two if np.isfinite(limit) and limit > 0 else np.nan
    return np.nan


def figure_all_constraint_margins(comparison: pd.DataFrame, figure_dir: Path) -> list[Path]:
    """Show every supported active numerical/inequality validation constraint."""
    frame = comparison.copy()
    frame["response_group"] = frame["response_kind"].map(_response_group)
    order = {
        "Corrosion rate": 0,
        "Selective depletion": 1,
        "Localized damage": 2,
        "Areal mass change": 3,
        "Salt inventory": 4,
        "Redox / diffusion": 5,
        "Other": 6,
    }
    frame["_group_order"] = frame["response_group"].map(order).fillna(99)
    frame["_id_order"] = frame["measurement_id"].astype(str).str.extract(r"(\d+)")[0].astype(int)
    frame = frame.sort_values(["_group_order", "_id_order"]).reset_index(drop=True)
    matrix = np.asarray(
        [
            [
                _constraint_margin_log2(row, float(row[f"{model}_prediction"]))
                for _, row in frame.iterrows()
            ]
            for model in MODEL_ORDER
        ],
        dtype=float,
    )
    finite = matrix[np.isfinite(matrix)]
    scale = min(3.0, max(1.0, float(np.nanmax(np.abs(finite))) if finite.size else 1.0))
    fig, ax = plt.subplots(figsize=(8.2, 2.75), constrained_layout=True)
    image = ax.imshow(
        np.clip(matrix, -scale, scale),
        aspect="auto",
        cmap="RdYlGn",
        norm=TwoSlopeNorm(vmin=-scale, vcenter=0.0, vmax=scale),
        interpolation="nearest",
    )
    ax.set_yticks(np.arange(len(MODEL_ORDER)))
    ax.set_yticklabels([MODEL_LABELS[model] for model in MODEL_ORDER])
    ax.set_xticks(np.arange(len(frame)))
    ax.set_xticklabels(frame["measurement_id"].astype(str), rotation=58, ha="right")
    boundaries, group_labels = _group_boundaries(frame)
    for boundary in boundaries:
        ax.axvline(boundary, color="black", linewidth=0.65)
    group_display = {
        "Salt inventory": "Salt\ninventory",
        "Redox / diffusion": "Redox /\ndiffusion",
    }
    for midpoint, label in group_labels:
        ax.text(
            midpoint,
            1.03,
            group_display.get(label, label),
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=6.3,
            linespacing=0.90,
        )
    for row_index, model in enumerate(MODEL_ORDER):
        passed = frame[f"{model}_constraint_pass"].astype(bool).to_numpy()
        for column_index, is_pass in enumerate(passed):
            if not is_pass:
                ax.plot(column_index, row_index, marker="x", color="black", markersize=5.5, mew=1.0)
    colorbar = fig.colorbar(image, ax=ax, pad=0.015, fraction=0.035)
    colorbar.set_label(r"Signed constraint margin ($\log_2$ factor units)")
    ax.set_title("All 32 supported active validation constraints (green: pass; red: violation)", pad=19)
    _panel_label(ax, "c")
    return _save_figure(fig, figure_dir, "fig17_all_active_constraint_margins")


def figure_error_distributions(comparison: pd.DataFrame, figure_dir: Path) -> list[Path]:
    frame = _direct_range(comparison)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2), constrained_layout=True)
    for ax, error_suffix, title, letter in (
        (axes[0], "midpoint_log_error", "Midpoint-referenced factor error", "a"),
        (axes[1], "range_aware_log_error", "Interval-aware factor error", "b"),
    ):
        for model in MODEL_ORDER:
            errors = np.abs(frame[f"{model}_{error_suffix}"].astype(float).to_numpy())
            factors = np.sort(np.exp(errors))
            cumulative = np.arange(1, len(factors) + 1) / len(factors)
            ax.step(
                factors,
                cumulative,
                where="post",
                color=MODEL_COLORS[model],
                label=MODEL_LABELS[model],
            )
        ax.axvline(2.0, color="black", linestyle="--", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_xlim(1.0, max(3.2, ax.get_xlim()[1]))
        ax.set_ylim(0.0, 1.02)
        ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
        ax.set_xlabel("Factor error")
        ax.set_ylabel("Cumulative fraction of targets")
        ax.set_title(title)
        _panel_label(ax, letter)
    axes[1].legend(loc="lower right")
    return _save_figure(fig, figure_dir, "fig16_factor_error_distributions")


def make_all_publication_figures(
    comparison: pd.DataFrame,
    metrics: Mapping[str, Any],
    dynamic_model: DynamicRedoxInventoryDepletionModel,
    targets: pd.DataFrame,
    figure_dir: str | Path,
    results_dir: str | Path,
) -> dict[str, list[str]]:
    apply_publication_style()
    figure_dir = Path(figure_dir)
    results_dir = Path(results_dir)
    generated = {
        "validation_ratios": figure_validation_ratios(comparison, figure_dir),
        "parity_by_response": figure_parity_by_response(comparison, figure_dir),
        "residual_heatmap": figure_residual_heatmap(comparison, figure_dir),
        "performance_summary": figure_performance_summary(metrics, figure_dir),
        "dynamic_trajectories": figure_dynamic_trajectories(dynamic_model, targets, figure_dir, results_dir),
        "species_partitioning": figure_species_partitioning(dynamic_model, targets, figure_dir, results_dir),
        "factor_error_distributions": figure_error_distributions(comparison, figure_dir),
        "all_constraint_margins": figure_all_constraint_margins(comparison, figure_dir),
    }
    repository_root = results_dir.parent
    manifest: dict[str, list[str]] = {}
    for key, paths in generated.items():
        values: list[str] = []
        for path in paths:
            try:
                values.append(str(path.relative_to(repository_root)))
            except ValueError:
                values.append(str(path))
        manifest[key] = values
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "publication_figure_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    _direct_range(comparison).drop(columns=["_group_order", "_id_order"]).to_csv(
        results_dir / "publication_validation_case_order.csv",
        index=False,
    )
    return manifest
