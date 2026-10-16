"""Shared helpers for the chemistry training plotting scripts.

Every case directory has a plot.py that imports this module:

    import sys, pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
    import plotting

The helpers read the MOOSE CSV outputs with the standard library only (no pandas), apply a
common style and save the figure next to the input and into the shared figures/ directory used by
the slides.
"""

import csv
import glob
import os
import pathlib
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TRAINING_DIR = pathlib.Path(__file__).resolve().parents[1]
FIGURES_DIR = TRAINING_DIR / "figures"

# Categorical series colours in fixed order (never cycled past the list); analytic/reference curves
# are drawn in INK_MUTED with a dashed style so that the model result keeps the series colour.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8985"
SURFACE = "#fcfcfb"

plt.rcParams.update(
    {
        "figure.figsize": (6.4, 4.2),
        "figure.dpi": 120,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.prop_cycle": matplotlib.cycler(color=SERIES),
        "axes.edgecolor": INK_MUTED,
        "axes.labelcolor": INK_SECONDARY,
        "axes.titlecolor": INK,
        "xtick.color": INK_SECONDARY,
        "ytick.color": INK_SECONDARY,
        "text.color": INK,
        "axes.grid": True,
        "grid.color": "#e4e3df",
        "grid.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 11,
        "legend.frameon": False,
        "lines.linewidth": 2.0,
        "lines.markersize": 7,
    }
)


def reference(ax, x, y, label="analytic"):
    """Draw an analytic / reference curve (dashed, muted ink) under the model curves."""
    ax.plot(x, y, ls="--", lw=1.6, color=INK_MUTED, label=label, zorder=1)


def read_csv(path):
    """Read a MOOSE CSV file into a dict of float lists keyed by column name."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"{path} is empty")
    return {key: [float(r[key]) for r in rows] for key in rows[0]}


def latest_vpp(case_dir, prefix):
    """Return the last time-step file '<prefix>_NNNN.csv' of a VectorPostprocessor (or the
    '<prefix>_FINAL.csv' file when present)."""
    final = pathlib.Path(case_dir) / f"{prefix}_FINAL.csv"
    if final.exists():
        return final
    files = sorted(glob.glob(str(pathlib.Path(case_dir) / f"{prefix}_[0-9]*.csv")))
    if not files:
        raise FileNotFoundError(f"no VectorPostprocessor output {prefix}_*.csv in {case_dir}")
    return pathlib.Path(files[-1])


def ensure_output(case_dir, input_file, csv_file, extra_args=(), ranks=1):
    """Run the input with open_pronghorn-opt if its CSV output is missing.

    The executable is taken from $OPEN_PRONGHORN_EXE, else from the repository root. With ranks > 1
    the input is run with 'mpiexec -n <ranks>'."""
    case_dir = pathlib.Path(case_dir)
    if (case_dir / csv_file).exists():
        return
    exe = os.environ.get("OPEN_PRONGHORN_EXE", str(TRAINING_DIR.parents[2] / "open_pronghorn-opt"))
    print(f"running {input_file} on {ranks} rank(s) ...")
    launcher = ["mpiexec", "-n", str(ranks)] if ranks > 1 else []
    subprocess.run([*launcher, exe, "-i", input_file, *extra_args], cwd=case_dir, check=True,
                   stdout=subprocess.DEVNULL)


def save(fig, case_dir, name):
    """Save a figure as figures/<name>.png (shared by the READMEs and the slides)."""
    fig.tight_layout()
    FIGURES_DIR.mkdir(exist_ok=True)
    fig.savefig(FIGURES_DIR / f"{name}.png")
    print(f"wrote figures/{name}.png")
    plt.close(fig)


def case_setup(case_file, input_file, csv_file, extra_args=(), ranks=1):
    """Common start of a plot.py: return (case_dir, data) after running the input if needed."""
    case_dir = pathlib.Path(case_file).resolve().parent
    ensure_output(case_dir, input_file, csv_file, extra_args, ranks)
    return case_dir, read_csv(case_dir / csv_file)
