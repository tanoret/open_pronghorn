#!/usr/bin/env python3
"""
Plot chromium tracer diffusivity correlations from chromium_diffusivity_library_v5.json.

Produces FOUR plots:
    1. chromium_bulk_lattice_comparison
    2. chromium_grain_boundary_comparison
    3. stainless_316_bulk_vs_grain_boundary
    4. incoloy_800_bulk_vs_grain_boundary

Plot conventions:
    - Solid line  = within the experimentally fitted temperature range.
    - Dashed line = extrapolation outside the fitted range.
    - A given dataset ALWAYS keeps the same color between solid and dashed parts.
    - Material colors are kept consistent across all four figures.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


R = 8.31446261815324  # J/(mol K)


# Use Matplotlib's normal default palette, but assign it ONCE globally so
# materials do not change color from one figure to another.
_DEFAULT_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]

COLOR_HASTELLOY_N = _DEFAULT_COLORS[0]
COLOR_304 = _DEFAULT_COLORS[1]
COLOR_316 = _DEFAULT_COLORS[2]
COLOR_INCOLOY_800 = _DEFAULT_COLORS[3]
COLOR_316L = _DEFAULT_COLORS[4]


def diffusivity_arrhenius(temperature_k, D0_cm2_s, Q_j_mol):
    temperature_k = np.asarray(temperature_k, dtype=float)
    if np.any(temperature_k <= 0.0):
        raise ValueError("Temperature must be greater than 0 K.")
    return D0_cm2_s * np.exp(-Q_j_mol / (R * temperature_k))


def arrhenius_to_celsius(value):
    value = np.asarray(value, dtype=float)
    return 1000.0 / value - 273.15


def celsius_to_arrhenius(value):
    value = np.asarray(value, dtype=float)
    return 1000.0 / (value + 273.15)


def load_database(path):
    with Path(path).open("r", encoding="utf-8") as f:
        return json.load(f)


def diffusion_records(database):
    return database["diffusion_coefficients"]["records"]


def select_records(database, **criteria):
    matches = []
    for record in diffusion_records(database):
        if all(record.get(key) == value for key, value in criteria.items()):
            matches.append(record)
    return matches


def get_single_record(database, **criteria):
    matches = select_records(database, **criteria)

    if len(matches) == 0:
        text = ", ".join(f"{k}={v!r}" for k, v in criteria.items())
        raise KeyError(f"No diffusivity record found for {text}")

    if len(matches) > 1:
        text = ", ".join(f"{k}={v!r}" for k, v in criteria.items())
        raise ValueError(
            f"Expected one model-ready record for {text}, found {len(matches)}."
        )

    return matches[0]


def add_temperature_axis(ax):
    top_axis = ax.secondary_xaxis(
        "top",
        functions=(arrhenius_to_celsius, celsius_to_arrhenius),
    )
    top_axis.set_xlabel("Temperature (degC)", fontsize=12, labelpad=9)
    top_axis.tick_params(labelsize=10, direction="in")


def format_axis(ax, title, ylabel):
    ax.set_xlabel(r"$1000/T$ (K$^{-1}$)", fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(title, fontsize=15, pad=14)
    ax.grid(True, which="major", alpha=0.28)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=11, direction="in")
    ax.tick_params(axis="both", which="minor", direction="in")
    add_temperature_axis(ax)


def plot_record_with_extrapolation(
    ax,
    record,
    temperature_plot_k,
    label,
    color,
    marker=None,
    linewidth=2.2,
):
    """
    Plot one Arrhenius correlation.

    The solid measured segment and dashed extrapolated segment ALWAYS use
    exactly the same color. Optional markers distinguish conditions without
    creating a rainbow of extra colors.
    """
    T_min = float(record["T_min_K"])
    T_max = float(record["T_max_K"])

    values = diffusivity_arrhenius(
        temperature_plot_k,
        float(record["D0_cm2_s"]),
        float(record["Q_J_mol"]),
    )

    x = 1000.0 / temperature_plot_k
    y = np.log10(values)

    measured = (temperature_plot_k >= T_min) & (temperature_plot_k <= T_max)
    extrapolated = ~measured

    common = dict(
        color=color,
        linewidth=linewidth,
    )

    if marker is not None:
        common.update(
            marker=marker,
            markevery=65,
            markersize=4.5,
            markerfacecolor="none",
        )

    ax.plot(
        x,
        np.where(measured, y, np.nan),
        linestyle="-",
        label=label,
        **common,
    )

    ax.plot(
        x,
        np.where(extrapolated, y, np.nan),
        linestyle="--",
        alpha=0.85,
        **common,
    )


def save_figure(fig, output_directory, stem):
    output_directory.mkdir(parents=True, exist_ok=True)

    png_path = output_directory / f"{stem}.png"
    pdf_path = output_directory / f"{stem}.pdf"

    fig.tight_layout()
    fig.savefig(png_path, dpi=400, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


def make_bulk_lattice_comparison(database, output_directory):
    records = [
        (
            get_single_record(
                database,
                normalized_name="hastelloy_n",
                diffusion_type="bulk",
            ),
            "INOR-8 / Hastelloy N",
            COLOR_HASTELLOY_N,
        ),
        (
            get_single_record(
                database,
                normalized_name="stainless_304",
                diffusion_type="volume",
            ),
            "AISI 304",
            COLOR_304,
        ),
        (
            get_single_record(
                database,
                normalized_name="stainless_316",
                diffusion_type="volume",
                condition="solution_treated_grain_grown",
            ),
            "SUS316",
            COLOR_316,
        ),
        (
            get_single_record(
                database,
                normalized_name="incoloy_800",
                diffusion_type="volume",
            ),
            "Incoloy-800",
            COLOR_INCOLOY_800,
        ),
    ]

    temperature_c = np.linspace(500.0, 1100.0, 900)
    temperature_k = temperature_c + 273.15

    fig, ax = plt.subplots(figsize=(9.0, 6.6))

    for record, label, color in records:
        plot_record_with_extrapolation(
            ax,
            record,
            temperature_k,
            label,
            color=color,
        )

    format_axis(
        ax,
        "Chromium Bulk / Lattice Diffusivity",
        r"$\log_{10}\!\left[D\;(\mathrm{cm^2\,s^{-1}})\right]$",
    )

    ax.legend(
        frameon=True,
        fontsize=10,
        loc="best",
        title="Solid: measured range\nDashed: extrapolation",
        title_fontsize=9,
    )

    save_figure(fig, output_directory, "chromium_bulk_lattice_comparison")


def make_grain_boundary_comparison(database, output_directory):
    """
    Compare direct D_gb values only.

    All SUS316 conditions use the SAME 316 material color; markers distinguish
    processing condition. SUS316L and Incoloy-800 retain their own material
    colors. This keeps the figure readable without introducing many colors.
    """
    records = [
        (
            get_single_record(
                database,
                normalized_name="stainless_316",
                diffusion_type="grain_boundary",
                condition="solution_treated",
            ),
            "SUS316, solution treated",
            COLOR_316,
            "o",
        ),
        (
            get_single_record(
                database,
                normalized_name="stainless_316",
                diffusion_type="grain_boundary",
                condition="sensitized_923K_172.8ks_air",
            ),
            "SUS316, sensitized",
            COLOR_316,
            "s",
        ),
        (
            get_single_record(
                database,
                normalized_name="stainless_316",
                diffusion_type="grain_boundary",
                condition="cold_worked_20pct_thickness_reduction",
            ),
            "SUS316, 20% cold worked",
            COLOR_316,
            "^",
        ),
        (
            get_single_record(
                database,
                normalized_name="stainless_316l",
                diffusion_type="grain_boundary",
                condition="solution_treated",
            ),
            "SUS316L, solution treated",
            COLOR_316L,
            "D",
        ),
        (
            get_single_record(
                database,
                normalized_name="incoloy_800",
                diffusion_type="grain_boundary",
            ),
            "Incoloy-800",
            COLOR_INCOLOY_800,
            None,
        ),
    ]

    temperature_c = np.linspace(200.0, 900.0, 1000)
    temperature_k = temperature_c + 273.15

    fig, ax = plt.subplots(figsize=(9.0, 6.6))

    for record, label, color, marker in records:
        plot_record_with_extrapolation(
            ax,
            record,
            temperature_k,
            label,
            color=color,
            marker=marker,
        )

    format_axis(
        ax,
        "Chromium Grain-Boundary Diffusivity",
        r"$\log_{10}\!\left[D_{gb}\;(\mathrm{cm^2\,s^{-1}})\right]$",
    )

    ax.legend(
        frameon=True,
        fontsize=9.5,
        loc="best",
        title="Solid: measured range\nDashed: extrapolation",
        title_fontsize=9,
    )

    save_figure(fig, output_directory, "chromium_grain_boundary_comparison")


def make_stainless_316_bulk_vs_grain_boundary(database, output_directory):
    bulk_record = get_single_record(
        database,
        normalized_name="stainless_316",
        diffusion_type="volume",
        condition="solution_treated_grain_grown",
    )

    gb_record = get_single_record(
        database,
        normalized_name="stainless_316",
        diffusion_type="grain_boundary",
        condition="solution_treated",
    )

    temperature_c = np.linspace(200.0, 950.0, 1100)
    temperature_k = temperature_c + 273.15

    fig, ax = plt.subplots(figsize=(8.8, 6.6))

    # Same material = same material color. Marker distinguishes D_gb.
    plot_record_with_extrapolation(
        ax,
        bulk_record,
        temperature_k,
        "SUS316 bulk / lattice",
        color=COLOR_316,
    )

    plot_record_with_extrapolation(
        ax,
        gb_record,
        temperature_k,
        r"SUS316 grain boundary $D_{gb}$",
        color=COLOR_316,
        marker="o",
    )

    format_axis(
        ax,
        "Chromium Diffusivity in SUS316: Bulk vs Grain Boundary",
        r"$\log_{10}\!\left[D\;(\mathrm{cm^2\,s^{-1}})\right]$",
    )

    ax.legend(
        frameon=True,
        fontsize=10,
        loc="best",
        title="Solid: measured range\nDashed: extrapolation",
        title_fontsize=9,
    )

    ax.text(
        0.02,
        0.03,
        "Measured ranges do not overlap:\n"
        "bulk 888-1173 K; GB 518-698 K",
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        ha="left",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

    save_figure(fig, output_directory, "stainless_316_bulk_vs_grain_boundary")


def make_incoloy_800_bulk_vs_grain_boundary(database, output_directory):
    bulk_record = get_single_record(
        database,
        normalized_name="incoloy_800",
        diffusion_type="volume",
    )

    gb_record = get_single_record(
        database,
        normalized_name="incoloy_800",
        diffusion_type="grain_boundary",
    )

    temperature_c = np.linspace(450.0, 1250.0, 1100)
    temperature_k = temperature_c + 273.15

    fig, ax = plt.subplots(figsize=(8.8, 6.6))

    # Same material = same material color. Marker distinguishes D_gb.
    plot_record_with_extrapolation(
        ax,
        bulk_record,
        temperature_k,
        "Incoloy-800 bulk / lattice",
        color=COLOR_INCOLOY_800,
    )

    plot_record_with_extrapolation(
        ax,
        gb_record,
        temperature_k,
        r"Incoloy-800 grain boundary $D_{gb}$",
        color=COLOR_INCOLOY_800,
        marker="o",
    )

    format_axis(
        ax,
        "Chromium Diffusivity in Incoloy-800: Bulk vs Grain Boundary",
        r"$\log_{10}\!\left[D\;(\mathrm{cm^2\,s^{-1}})\right]$",
    )

    ax.legend(
        frameon=True,
        fontsize=10,
        loc="best",
        title="Solid: measured range\nDashed: extrapolation",
        title_fontsize=9,
    )

    ax.text(
        0.02,
        0.03,
        "Both fits are measured over 1060-1170 K",
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        ha="left",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

    save_figure(fig, output_directory, "incoloy_800_bulk_vs_grain_boundary")


def parse_args():
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="Plot chromium tracer diffusivity correlations."
    )
    parser.add_argument(
        "--database",
        type=Path,
        default=script_dir / "chromium_diffusivity_library_v5.json",
        help="Path to the chromium diffusivity JSON database.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir / "diffusivity_plots",
        help="Directory for PNG/PDF plot files.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    database = load_database(args.database)

    make_bulk_lattice_comparison(database, args.output_dir)
    make_grain_boundary_comparison(database, args.output_dir)
    make_stainless_316_bulk_vs_grain_boundary(database, args.output_dir)
    make_incoloy_800_bulk_vs_grain_boundary(database, args.output_dir)

    print()
    print("Generated 4 plot sets:")
    print("  1. chromium_bulk_lattice_comparison")
    print("  2. chromium_grain_boundary_comparison")
    print("  3. stainless_316_bulk_vs_grain_boundary")
    print("  4. incoloy_800_bulk_vs_grain_boundary")


if __name__ == "__main__":
    main()
