"""Publication-quality plotting for the MSR corrosion/plating validation package."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .model import MoltenSaltBVModel


def set_publication_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 450,
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.8,
            "lines.markersize": 5.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def _positive_df(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    mask = np.ones(len(out), dtype=bool)
    for col in cols:
        vals = pd.to_numeric(out[col], errors="coerce")
        mask &= np.isfinite(vals) & (vals > 0)
    return out.loc[mask].copy()


def plot_parity(predictions: pd.DataFrame, outdir: str | Path) -> None:
    set_publication_style()
    outdir = Path(outdir)
    df = predictions[predictions["fit_role"].isin(["direct", "range"])].copy()
    # Keep parity to corrosion, plating, salt-inventory, and depth-like responses.
    # Diffusion coefficients and redox ratios differ by many orders of magnitude
    # and are better reported in parameter/residual tables.
    df = df[~df["response_kind"].isin(["cr_diffusion_cm2_s", "te_redox_threshold_ratio"])]
    df = _positive_df(df, ["target_mid", "prediction"])
    if df.empty:
        return
    kinds = sorted(df["response_kind"].dropna().unique())
    markers = ["o", "s", "^", "D", "v", "P", "X", "<", ">"]
    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    for i, kind in enumerate(kinds):
        sub = df[df["response_kind"] == kind]
        ax.scatter(
            sub["target_mid"].astype(float),
            sub["prediction"].astype(float),
            marker=markers[i % len(markers)],
            edgecolor="black",
            linewidth=0.35,
            alpha=0.85,
            label=kind.replace("_", " "),
        )
    all_values = np.concatenate([df["target_mid"].astype(float).to_numpy(), df["prediction"].astype(float).to_numpy()])
    lo = 10 ** math.floor(math.log10(max(np.nanmin(all_values), 1e-16)))
    hi = 10 ** math.ceil(math.log10(np.nanmax(all_values)))
    ax.plot([lo, hi], [lo, hi], linestyle="--", color="black", linewidth=1.0, label="1:1")
    ax.plot([lo, hi], [2 * lo, 2 * hi], linestyle=":", color="0.35", linewidth=0.9)
    ax.plot([lo, hi], [0.5 * lo, 0.5 * hi], linestyle=":", color="0.35", linewidth=0.9, label="factor 2")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Observed target in model-native units")
    ax.set_ylabel("Predicted value in same units")
    ax.set_title("Calibration parity across extracted numerical targets")
    ax.legend(frameon=False, loc="best")
    ax.grid(True, which="major", linewidth=0.35, alpha=0.35)
    _save(fig, outdir / "fig01_calibration_parity")


def plot_residuals(predictions: pd.DataFrame, outdir: str | Path) -> None:
    set_publication_style()
    outdir = Path(outdir)
    df = predictions[predictions["fit_role"].isin(["direct", "range"])].copy()
    df = _positive_df(df, ["target_mid", "prediction"])
    if df.empty:
        return
    df["ln_error"] = np.log(df["prediction"].astype(float) / df["target_mid"].astype(float))
    df = df.sort_values("ln_error")
    labels = df["measurement_id"].astype(str).to_list()
    y = np.arange(len(df))
    fig_h = max(4.8, 0.26 * len(df) + 1.4)
    fig, ax = plt.subplots(figsize=(7.2, fig_h))
    ax.axvline(0.0, color="black", linewidth=1.0)
    ax.axvline(math.log(2.0), color="0.35", linestyle=":", linewidth=0.9)
    ax.axvline(-math.log(2.0), color="0.35", linestyle=":", linewidth=0.9)
    ax.scatter(df["ln_error"], y, edgecolor="black", linewidth=0.3, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("log(predicted / observed)")
    ax.set_ylabel("Measurement")
    ax.set_title("Signed calibration residuals for direct/range targets")
    ax.grid(True, axis="x", linewidth=0.35, alpha=0.35)
    _save(fig, outdir / "fig02_signed_residuals")


def plot_temperature_flow_sensitivity(model: MoltenSaltBVModel, outdir: str | Path) -> None:
    set_publication_style()
    outdir = Path(outdir)
    temps = np.linspace(520, 760, 121)
    flow_scenarios = [
        ("static", 0.35),
        ("thermal convection", 1.0),
        ("reactor circulation", 2.0),
        ("turbulent metal", 2.8),
    ]
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    for label, flow in flow_scenarios:
        rates = []
        for T in temps:
            row = {
                "temperature_C": T,
                "temperature_K": T + 273.15,
                "material_class": "hastelloy_n",
                "salt_class": "fluoride_fuel",
                "redox_class": "purified_baseline",
                "flow_factor": flow,
                "delta_T_C": 100.0 if flow >= 1.0 else 0.0,
                "position_class": "hot_leg",
                "surface_class": "metal",
            }
            rates.append(model.corrosion_rate_um_y(row))
        ax.plot(temps, rates, label=label)
    ax.set_yscale("log")
    ax.set_xlabel("Temperature (C)")
    ax.set_ylabel("Predicted corrosion rate (um/year)")
    ax.set_title("Temperature and circulation sensitivity: Hastelloy N / fuel fluoride")
    ax.legend(frameon=False)
    ax.grid(True, which="both", linewidth=0.35, alpha=0.35)
    _save(fig, outdir / "fig03_temperature_flow_sensitivity")


def plot_redox_effects(redox_table: pd.DataFrame, outdir: str | Path) -> None:
    set_publication_style()
    outdir = Path(outdir)
    df = redox_table.copy()
    df = df.sort_values("corrosion_multiplier_vs_baseline")
    fig_h = max(4.8, 0.33 * len(df) + 1.4)
    fig, ax = plt.subplots(figsize=(7.0, fig_h))
    y = np.arange(len(df))
    ax.barh(y, df["corrosion_multiplier_vs_baseline"].astype(float))
    ax.axvline(1.0, color="black", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(df["redox_class"].str.replace("_", " "))
    ax.set_xlabel("Corrosion multiplier relative to purified baseline")
    ax.set_title("Fitted effective Butler-Volmer redox terms")
    ax.grid(True, axis="x", which="both", linewidth=0.35, alpha=0.35)
    _save(fig, outdir / "fig04_redox_multipliers")


def plot_plating_ranking(model: MoltenSaltBVModel, outdir: str | Path) -> None:
    set_publication_style()
    outdir = Path(outdir)
    base = {
        "temperature_K": 650.0 + 273.15,
        "material_class": "hastelloy_n",
        "salt_class": "fluoride_fuel",
        "redox_class": "fission_product",
        "flow_factor": 2.0,
        "delta_T_C": 22.0,
        "position_class": "core_or_reactor",
        "surface_class": "metal",
        "time_years": 1.0,
    }
    ranking = model.deposition_ranking(base)
    df = pd.DataFrame({"surface": list(ranking.keys()), "deposition_rate_um_y": list(ranking.values())})
    df = df.sort_values("deposition_rate_um_y")
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.bar(df["surface"].str.replace("_", " "), df["deposition_rate_um_y"])
    ax.set_yscale("log")
    ax.set_ylabel("Relative deposition rate (um/year equivalent)")
    ax.set_xlabel("Surface / flow regime")
    ax.set_title("Noble-metal plating branch reproduces MSRE ranking")
    ax.tick_params(axis="x", rotation=25)
    ax.grid(True, axis="y", which="both", linewidth=0.35, alpha=0.35)
    _save(fig, outdir / "fig05_plating_surface_ranking")


def plot_usage_audit(predictions: pd.DataFrame, outdir: str | Path) -> None:
    set_publication_style()
    outdir = Path(outdir)
    counts = predictions.groupby("fit_role", dropna=False).size().sort_values()
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.barh(counts.index.astype(str), counts.values)
    ax.set_xlabel("Number of measurements")
    ax.set_ylabel("Data-use category")
    ax.set_title("All detailed measurements are retained in the validation audit")
    for y, v in enumerate(counts.values):
        ax.text(v + 0.2, y, str(int(v)), va="center", fontsize=8)
    ax.grid(True, axis="x", linewidth=0.35, alpha=0.35)
    _save(fig, outdir / "fig06_data_usage_audit")


def plot_loop_simulation(model: MoltenSaltBVModel, outdir: str | Path) -> None:
    set_publication_style()
    outdir = Path(outdir)
    segments = pd.DataFrame(
        [
            {
                "segment": "NCL-16 hot leg",
                "temperature_C": 704.0,
                "temperature_K": 704.0 + 273.15,
                "material_class": "hastelloy_n",
                "salt_class": "fluoride_fuel",
                "redox_class": "purified_baseline",
                "flow_factor": 1.0,
                "delta_T_C": 166.0,
                "position_class": "hot_leg",
                "surface_class": "metal",
                "surface_area_cm2": 100.0,
            },
            {
                "segment": "NCL-16 cold leg",
                "temperature_C": 538.0,
                "temperature_K": 538.0 + 273.15,
                "material_class": "hastelloy_n",
                "salt_class": "fluoride_fuel",
                "redox_class": "purified_baseline",
                "flow_factor": 1.0,
                "delta_T_C": 166.0,
                "position_class": "cold_leg",
                "surface_class": "metal",
                "surface_area_cm2": 100.0,
            },
        ]
    )
    sim = model.simulate_loop(segments, duration_h=29500.0, dt_h=250.0, salt_volume_cm3=250.0, initial_cr_ppm=0.0)
    sim.to_csv(Path(outdir).parent / "results" / "ncl16_simplified_loop_simulation.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(sim["time_h"] / 1000.0, sim["salt_cr_ppm"])
    ax.set_xlabel("Time (10^3 h)")
    ax.set_ylabel("Predicted salt Cr concentration (ppm, normalized geometry)")
    ax.set_title("Example dynamic loop simulation: NCL-16 baseline")
    ax.grid(True, linewidth=0.35, alpha=0.35)
    _save(fig, outdir / "fig07_loop_simulation_ncl16")


def make_all_plots(outputs: dict[str, Any], outdir: str | Path) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    model: MoltenSaltBVModel = outputs["model"]
    predictions: pd.DataFrame = outputs["target_predictions"]
    plot_parity(predictions, outdir)
    plot_residuals(predictions, outdir)
    plot_temperature_flow_sensitivity(model, outdir)
    plot_redox_effects(outputs["redox_table"], outdir)
    plot_plating_ranking(model, outdir)
    plot_usage_audit(predictions, outdir)
    plot_loop_simulation(model, outdir)
