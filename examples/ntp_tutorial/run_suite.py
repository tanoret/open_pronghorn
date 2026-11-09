#!/usr/bin/env python3
"""Demonstration suite of the NTP fuel-element tutorial.

Runs the documented commands of every stage (the stage README.md files), from each stage
directory and with the output file bases the READMEs use, then extracts the figures of merit
from the CSV (and, for the 3-D transverse spread, Exodus) outputs and compares each one with the
value the READMEs document. The comparison is written to suite_results.md and suite_results.csv
next to this script; every run's log goes to suite_logs/<job>.log.

Groups (--only takes one or more, comma separated):
  01    closures: check_hydrogen.py, validate.py (46 Texas A&M points through both solvers)
  02    element: mesh and mesh.png, Newton and SIMPLE cold and 5 % power, cold_flow.png
  03    anchor: 1 GW anchor (continuation + implicit orifice), uniform-frit flow-excursion folds
  04    segregated Navier-Stokes (linear FV) at 1 GW in RZ: SIMPLE and PIMPLE with the orifice
        designed on the slice flows; extract_orifice.py. The findings runs (a) Newton orifice and
        (b) lagged temperature design run only with --placeholders
  05    3-D: mesh, uniform Newton anchor, segregated SIMPLE cold flow
  long  3-D: Newton power tilts 0.01, 0.05 and 0.15 and transverse.png; segregated SIMPLE and
        PIMPLE anchors with the column orifice, and the ring orifice that loses the flow around
        the element (minutes on 8 ranks each); run with --long or --only long
The default runs 01-05 without 'long'. Each run records its jobs' status in suite_logs/state.json;
the report covers every job recorded there, so a partial run (e.g. --only 04 --placeholders)
updates its rows of the full table and leaves the others as last run.

Tolerances: a documented number such as "2989.4 K" is matched when the fresh value rounds to it
(|fresh - doc| <= half a unit of its last digit). Round-off quantities ("1e-13") are checked as
upper bounds, run times within a factor of three.

Usage (inside the MOOSE environment, e.g. moose-dev-shell, from this directory):
    python3 run_suite.py --list
    python3 run_suite.py --dry-run --long
    python3 run_suite.py                                   # 01-05, ~8 minutes
    python3 run_suite.py --only 03,05 --jobs 3
    python3 run_suite.py --long --max-ranks 24
    python3 run_suite.py --report-only --long              # re-extract from existing outputs
Options:
    --exe PATH         OpenPronghorn executable (default ../../open_pronghorn-opt)
    --mpirun CMD       MPI launcher (default 'mpirun --bind-to none': Open MPI binds every job
                       to cores 0..n-1 otherwise, so concurrent jobs would share cores)
    --python PATH      interpreter of the tutorial scripts (default: this one); the plotting
                       scripts need numpy, matplotlib and netCDF4
    --jobs N           concurrent jobs (default 3); --max-ranks R total MPI ranks (default 24)
    --ranks-scale F    scale the documented rank counts (e.g. 0.5 on a small machine)
The process exit code is 0 when every selected job ran and every check matched.
"""
import argparse
import csv
import json
import math
import os
import re
import shlex
import subprocess
import sys
import time

sys.dont_write_bytecode = True  # no __pycache__ in the tutorial directories
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "common"))
import hydrogen as H  # noqa: E402

LOGS = os.path.join(HERE, "suite_logs")
GROUPS = ["01", "02", "03", "04", "05", "long"]
STAGE_DIRS = {"01": "01_closures", "02": "02_element", "03": "03_ludewig_anchor",
              "04": "04_navier_stokes", "05": "05_element_3d", "long": "05_element_3d"}
# stop a power continuation at its turning point: without it the time step creeps at the default
# dtmin of 1e-12 for a long time (the run ends with a non-zero code: its last solve failed)
FOLD_STOP = ["Executioner/dtmin=1e-6", "Executioner/error_on_dtmin=false"]
FOLD_MARK = "Aborting as timestep already at or below dtmin"
# the stage-4 findings (a) and (b) end in the runaway the README documents: the enthalpy of a
# starved cell leaves the range of the hydrogen property inversion
RUNAWAY_MARK = "T_from_h did not converge"


# ----------------------------------------------------------------------------------------------
# checks: documented value, tolerance, source
# ----------------------------------------------------------------------------------------------
def doc(text, unit="", scale=1.0, src=""):
    """A number as documented: matched when the fresh value rounds to it."""
    m, _, e = text.lower().partition("e")
    decimals = len(m.split(".")[1]) if "." in m else 0
    half = 0.5 * 10.0 ** (-decimals) * 10.0 ** (int(e) if e else 0)
    return dict(kind="abs", doc=float(text) * scale, tol=half * abs(scale), text=f"{text} {unit}".strip(),
                src=src, scale=scale, unit=unit)


def below(text, bound, unit="", src=""):
    """A round-off quantity: matched when |fresh| <= bound."""
    return dict(kind="below", doc=float(text), tol=bound, text=f"{text} {unit}".strip(), src=src, scale=1.0,
                unit=unit)


def above(text, bound, unit="", scale=1.0, src=""):
    """A lower bound (e.g. a runaway): matched when fresh >= bound."""
    return dict(kind="above", doc=float(text) * scale, tol=bound * scale, text=f"> {text} {unit}".strip(), src=src,
                scale=scale, unit=unit)


def factor(text, f=3.0, unit="s", src=""):
    """A run time: matched within a factor f."""
    return dict(kind="factor", doc=float(text), tol=f, text=f"{text} {unit}".strip(), src=src, scale=1.0,
                unit=unit)


def evaluate(check, value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "missing"
    if check["kind"] == "abs":
        ok = abs(value - check["doc"]) <= check["tol"] * (1 + 1e-9)
    elif check["kind"] == "below":
        ok = abs(value) <= check["tol"]
    elif check["kind"] == "above":
        ok = value >= check["tol"]
    else:
        ok = check["doc"] / check["tol"] <= value <= check["doc"] * check["tol"]
    return "ok" if ok else "MISMATCH"


def tol_text(check):
    if check["kind"] == "abs":
        return f"±{check['tol'] / abs(check['scale']):.3g} {check['unit']}".strip()
    if check["kind"] == "below":
        return f"|x| <= {check['tol']:.0e}"
    if check["kind"] == "above":
        return f">= {check['tol'] / check['scale']:g} {check['unit']}".strip()
    return f"×/÷{check['tol']:g}"


# ----------------------------------------------------------------------------------------------
# figures of merit from the outputs
# ----------------------------------------------------------------------------------------------
def read_csv(path):
    with open(path) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def newton_metrics(path):
    """Darcy-Forchheimer model (anchor_body.i postprocessors): last step."""
    rows = read_csv(path)
    r = rows[-1]
    out = dict(T_mixed=float(H.T_from_h(r["h_out_mixed"])) if r["mdot_out"] > 0 else float("nan"),
               T_fuel_max=r["T_fuel_max"], p_in=r["p_in"], dp_in=r["p_in"] - 7e6,
               energy_error=r["energy_error"], mass_error=r["mass_error"], power=r["power"],
               s_power=r.get("s_power", float("nan")), s_max=max(x.get("s_power", 0) for x in rows),
               steps=len(rows), faceting_pc=100 * (1 - 1 / r["volume_ratio"]) if r.get("volume_ratio") else 0.0)
    lam = [r[k] for k in sorted((k for k in r if re.fullmatch(r"log_chi_\d+", k)), key=lambda k: int(k[8:]))]
    if lam:
        out.update(chi_min=math.exp(min(lam)), chi_max=math.exp(max(lam)), chi_inlet=math.exp(lam[0]),
                   chi_exit=math.exp(lam[-1]))
    return out


def segregated_metrics(path):
    """Segregated model (element_segregated_postprocessors.i): last step or solve."""
    rows = read_csv(path)
    r = rows[-1]
    out = dict(T_mixed=r["T_mixed"], T_fuel_max=r["T_fuel_max"], p_in=r["p_in"], dp_in=r["p_in"] - 7e6,
               energy_error=r["energy_error"], mass_error=r["mass_error"], power=r["power"], mdot=r["mdot"],
               steps=len(rows))
    for k in ("chi_min", "chi_max", "s_power", "T_slice_spread", "T_slice_min", "T_slice_max"):
        if k in r:
            out[k] = r[k]
    if "s_power" in r:
        out["s_max"] = max(x.get("s_power", 0) for x in rows)
    out["T_fuel_peak"] = max(x["T_fuel_max"] for x in rows)
    return out


def asymmetry_metrics(exodus, python=None):
    """Azimuthal asymmetry of the last step of a 3-D segregated run: the largest spread of the
    fuel temperature around the element (over the axial rows) relative to the mean rise."""
    code = f"""
import json, sys
import numpy as np, netCDF4
d = netCDF4.Dataset({exodus!r})
names = [b"".join(n).decode().strip(chr(0)) for n in d["name_elem_var"][:]]
blocks = [b"".join(n).decode().strip(chr(0)) for n in d["eb_names"][:]]
X = np.vstack([np.asarray(d["coord" + c][:]) for c in "xyz"]).T
b = blocks.index("fuel_bed")
c = np.asarray(d["connect%d" % (b + 1)][:]) - 1
y = np.round(X[c].mean(1)[:, 1], 4)
T = np.asarray(d["vals_elem_var%deb%d" % (names.index("T_solid") + 1, b + 1)][-1])
print(json.dumps(dict(asymmetry=max(np.ptp(T[y == r]) for r in np.unique(y)) / max(T.mean() - 30, 1))))
"""
    try:
        import netCDF4  # noqa: F401
        cmd = [sys.executable, "-c", code]
    except ImportError:
        cmd = shlex.split(python or sys.executable) + ["-c", code]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("asymmetry reduction failed: " + out.stderr.strip().splitlines()[-1])
    return {k: float(v) for k, v in json.loads(out.stdout.strip().splitlines()[-1]).items()}


def with_asymmetry(csv_path, exodus, python):
    """segregated_metrics plus the azimuthal asymmetry (left out if it cannot be read)."""
    out = segregated_metrics(csv_path)
    try:
        out.update(asymmetry_metrics(exodus, python))
    except Exception as e:
        print(f"  asymmetry of {os.path.basename(exodus)} not extracted: {e}", flush=True)
    return out


def transverse_metrics(exodus, python=None):
    """Azimuthal spread of the last step of a 3-D run (plot_3d.py's reduction).

    Needs netCDF4; without it in this interpreter, the reduction runs in 'python' (--python)."""
    code = f"""
import json, sys
sys.dont_write_bytecode = True
sys.path.insert(0, {os.path.join(HERE, '05_element_3d')!r})
import numpy as np
import plot_3d
_, T_gas = plot_3d.azimuthal({exodus!r}, "hot_frit", "T_fluid", np.mean)
_, T_fuel = plot_3d.azimuthal({exodus!r}, "fuel_bed", "T_solid", np.max)
rise = T_gas.mean() - plot_3d.T_IN
print(json.dumps(dict(gas_rise_max=(T_gas.max() - plot_3d.T_IN) / rise - 1,
                      gas_rise_min=(T_gas.min() - plot_3d.T_IN) / rise - 1,
                      fuel_rise_ratio=(T_fuel.max() - plot_3d.T_IN) / rise)))
"""
    try:
        import netCDF4  # noqa: F401
        cmd = [sys.executable, "-c", code]
    except ImportError:
        cmd = shlex.split(python or sys.executable) + ["-c", code]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("transverse reduction failed: " + out.stderr.strip().splitlines()[-1])
    return {k: float(v) for k, v in json.loads(out.stdout.strip().splitlines()[-1]).items()}


def with_transverse(csv_path, exodus, python):
    """newton_metrics plus the transverse spread; the spread is left out if it cannot be read."""
    out = newton_metrics(csv_path)
    try:
        out.update(transverse_metrics(exodus, python))
    except Exception as e:  # e.g. no netCDF4: the checks of the spread are reported as missing
        print(f"  transverse spread of {os.path.basename(exodus)} not extracted: {e}", flush=True)
    return out


def stability_metrics(cwd):
    """stability.i: largest rise of a slice outlet temperature above the anchor after the kick and
    its e-folding time (from the slice vectors), the peak fuel temperature, the frozen orifice."""
    import glob
    rows = read_csv(os.path.join(cwd, "stability_out.csv"))
    files = sorted(glob.glob(os.path.join(cwd, "stability_out_slice_outlet_*.csv")))
    t_freeze = 1.05
    series = []
    for r, f in zip(rows, files):
        with open(f) as fh:
            series.append((r["time"], [float(x["T_fluid"]) for x in csv.DictReader(fh)]))
    base = min(series, key=lambda x: abs(x[0] - t_freeze))[1]
    dev = [(t, max(abs(a - b) for a, b in zip(T, base))) for t, T in series if t >= t_freeze]
    peak_t, peak = max(dev, key=lambda x: x[1])
    late = [x for x in dev if x[0] > peak_t + 0.1]
    efold = (late[-1][0] - late[0][0]) / math.log(late[0][1] / late[-1][1]) if len(late) > 1 else float("nan")
    return dict(peak_rise=peak, efold=efold, T_fuel_peak=max(r["T_fuel_max"] for r in rows if r["time"] >= t_freeze),
                chi_exit=rows[-1]["chi_exit"])


def log_floats(path, pattern):
    with open(path) as f:
        return [float(x) for x in re.findall(pattern, f.read())]


def at_power(path, fractions, keys):
    """Values of 'keys' at the rows whose s_power is closest to each of 'fractions'."""
    rows = [r for r in read_csv(path) if "s_power" in r]
    out = {}
    for s in fractions:
        r = min(rows, key=lambda x: abs(x["s_power"] - s))
        if abs(r["s_power"] - s) > 0.01:
            continue
        for k in keys:
            if k in r:
                out[f"{k}@{round(100 * s)}%"] = r[k]
    return out


def hydrogen_metrics(log):
    with open(log) as f:
        text = f.read()
    diffs = [float(x) for x in re.findall(r"max relative difference ([0-9.eE+-]+)", text)]
    out = dict(props_max_rel_diff=max(diffs) if diffs else float("nan"))
    m = re.search(r"h\(3000 K\) - h\(30 K\) = ([0-9.eE+-]+)", text)
    if m:
        out["dh_30_3000"] = float(m.group(1))
    m = re.search(r"bundle flow Q / dh = ([0-9.]+) kg/s, element flow = ([0-9.]+)", text)
    if m:
        out.update(mdot_bundle=float(m.group(1)), mdot_element=float(m.group(2)))
    return out


def validation_metrics(summary):
    with open(summary) as f:
        text = f.read()
    out = {}
    m = re.search(r"mean \|error\| ([0-9.]+) %, max \|error\| ([0-9.]+) %", text)
    if m:
        out.update(data_mean_err_pc=float(m.group(1)), data_max_err_pc=float(m.group(2)))
    m = re.search(r"\(Newton\) vs closed form: max \|relative difference\| ([0-9.eE+-]+)", text)
    if m:
        out["newton_vs_closed_form"] = float(m.group(1))
    m = re.search(r"\(SIMPLE\) vs closed form: max \|relative difference\| ([0-9.eE+-]+)", text)
    if m:
        out["simple_vs_closed_form"] = float(m.group(1))
    m = re.search(r"Texas A&M air data, (\d+) points", text)
    if m:
        out["points"] = float(m.group(1))
    return out


def compare_3d_2d(csv_3d, csv_2d):
    m3, m2 = newton_metrics(csv_3d), newton_metrics(csv_2d)
    return dict(fuel_diff_pc=100 * (m3["T_fuel_max"] / m2["T_fuel_max"] - 1),
                p_in_diff_pc=100 * abs(m3["p_in"] / m2["p_in"] - 1))


def extract_orifice_metrics(log):
    with open(log) as f:
        m = re.search(r"chi_cf from ([0-9.]+) to ([0-9.]+) over (\d+) slices", f.read())
    return dict(chi_min=float(m.group(1)), chi_max=float(m.group(2))) if m else {}


# ----------------------------------------------------------------------------------------------
# the jobs: the documented commands of the stage READMEs
# ----------------------------------------------------------------------------------------------
def moose(args, ranks, *inputs_and_options):
    cmd = [args.exe] + list(inputs_and_options) + ["--color", "off"]
    if ranks > 1:
        cmd = shlex.split(args.mpirun) + ["-n", str(ranks)] + cmd
    return cmd


def job_list(args):
    py, J = shlex.split(args.python), []

    def add(name, group, cmd, ranks=1, needs=(), outputs=(), metrics=None, checks=None, fold=False,
            placeholder=False, documented=None, note="", exodus_reader=False, runaway=False):
        J.append(dict(name=name, group=group, cwd=os.path.join(HERE, STAGE_DIRS[group]), cmd=cmd,
                      ranks=ranks, needs=list(needs), outputs=list(outputs), metrics=metrics,
                      checks=checks or {}, fold=fold, placeholder=placeholder,
                      documented=documented, note=note, exodus_reader=exodus_reader, runaway=runaway))

    a3 = ["-i", "../03_ludewig_anchor/anchor.i", "../03_ludewig_anchor/continuation.i"]
    o3 = ["../03_ludewig_anchor/orifice_implicit.i"]
    S1, S2, S3, S4, S5, R = ("01_closures/README.md", "02_element/README.md", "03_ludewig_anchor/README.md",
                             "04_navier_stokes/README.md", "05_element_3d/README.md", "README.md")

    # ---- 01 closures
    add("01_hydrogen", "01", [*py, "check_hydrogen.py", "--exe", args.exe],
        metrics=lambda j: hydrogen_metrics(j["log"]),
        documented="python3 check_hydrogen.py",
        checks=dict(props_max_rel_diff=[below("7e-14", 1e-12, src=S1), below("7e-14", 1e-12, src=R)],
                    dh_30_3000=[doc("4.749625e7", "J/kg", src=S1)],
                    mdot_bundle=[doc("21.054295", "kg/s", src=S1)],
                    mdot_element=[doc("0.569035", "kg/s", src=S1)]))
    add("01_validate", "01", [*py, "validate.py", "--exe", args.exe],
        outputs=["validation_points.csv", "validation_summary.txt", "validation.png"],
        metrics=lambda j: validation_metrics(os.path.join(j["cwd"], "validation_summary.txt")),
        documented="python3 validate.py",
        checks=dict(points=[doc("46", src=S1)],
                    data_mean_err_pc=[doc("5.32", "%", src=S1), doc("5.3", "%", src=R)],
                    data_max_err_pc=[doc("13.73", "%", src=S1), doc("13.7", "%", src=R)],
                    newton_vs_closed_form=[below("8e-12", 1e-10, src=S1)],
                    simple_vs_closed_form=[below("2e-13", 1e-11, src=S1)]))

    # ---- 02 element
    add("02_mesh", "02", moose(args, 1, "-i", "mesh.i", "--mesh-only", "--allow-unused"), outputs=["mesh_in.e"],
        documented="open_pronghorn-opt -i mesh.i --mesh-only --allow-unused")
    add("02_plot_mesh", "02", [*py, "plot_mesh.py"], needs=["02_mesh"], outputs=["mesh.png"], exodus_reader=True,
        documented="python3 plot_mesh.py")
    add("02_newton_cold", "02", moose(args, 1, *a3, "Executioner/end_time=0", "Outputs/file_base=newton_cold"),
        outputs=["newton_cold.csv", "newton_cold.e"], metrics=lambda j: newton_metrics(j["csv"]),
        documented="open_pronghorn-opt -i ../03_ludewig_anchor/anchor.i ../03_ludewig_anchor/continuation.i "
                   "Executioner/end_time=0 Outputs/file_base=newton_cold",
        checks=dict(dp_in=[doc("16.2", "kPa", 1e3, S2)], mass_error=[doc("-9e-7", src=S2)]))
    add("02_newton_5pc", "02", moose(args, 1, *a3, "Executioner/end_time=0.1", "Outputs/file_base=newton_5pc"),
        outputs=["newton_5pc.csv", "newton_5pc.e"], metrics=lambda j: newton_metrics(j["csv"]),
        documented="open_pronghorn-opt -i ../03_ludewig_anchor/anchor.i ../03_ludewig_anchor/continuation.i "
                   "Executioner/end_time=0.1 Outputs/file_base=newton_5pc",
        checks=dict(T_mixed=[doc("218.9", "K", src=S2)], energy_error=[doc("8e-8", src=S2)],
                    T_fuel_max=[doc("325", "K", src=S2)]))
    add("02_simple_cold", "02", moose(args, 1, "-i", "element_simple_rz.i", "Outputs/file_base=simple_cold"),
        outputs=["simple_cold.csv", "simple_cold.e"], metrics=lambda j: segregated_metrics(j["csv"]),
        documented="open_pronghorn-opt -i element_simple_rz.i Outputs/file_base=simple_cold",
        checks=dict(dp_in=[doc("11.6", "kPa", 1e3, S2)], mass_error=[below("5e-13", 1e-11, src=S2)]))
    add("02_simple_5pc", "02", moose(args, 1, "-i", "element_simple_rz.i", "s_power=0.05", "simple_fixed_point=40"),
        outputs=["element_simple_rz_out.csv", "element_simple_rz_out.e"],
        metrics=lambda j: segregated_metrics(j["csv"]),
        documented="open_pronghorn-opt -i element_simple_rz.i s_power=0.05 simple_fixed_point=40",
        checks=dict(T_mixed=[doc("218.9", "K", src=S2)], energy_error=[doc("-5e-8", src=S2)],
                    T_fuel_max=[doc("1085", "K", src=S2)]))
    add("02_compare_cold_flow", "02", [*py, "compare_cold_flow.py"], needs=["02_newton_cold", "02_simple_cold"],
        outputs=["cold_flow.png"], exodus_reader=True, documented="python3 compare_cold_flow.py")

    # ---- 03 anchor
    a = ["-i", "anchor.i", "continuation.i"]
    add("03_anchor", "03", moose(args, 4, *a, "orifice_implicit.i", "Outputs/file_base=anchor_out"), ranks=4,
        outputs=["anchor_out.csv", "anchor_out.e"], metrics=lambda j: newton_metrics(j["csv"]),
        documented="mpirun -n 4 open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i "
                   "Outputs/file_base=anchor_out",
        checks=dict(T_mixed=[doc("2989.4", "K", src=S3), doc("2989.4", "K", src=R), doc("2989.4", "K", src=S5)],
                    T_fuel_max=[doc("3193", "K", src=S3), doc("3193", "K", src=R), doc("3193.2", "K", src=S5)],
                    p_in=[doc("7.84", "MPa", 1e6, S3), doc("7.84", "MPa", 1e6, R), doc("7.841", "MPa", 1e6, S5)],
                    chi_inlet=[doc("0.15", src=S3)], chi_exit=[doc("20.6", src=S3), doc("20.63", src=S5)],
                    energy_error=[doc("3e-7", src=S3), doc("3e-7", src=R), doc("3e-7", src=S5)],
                    mass_error=[below("1e-13", 1e-12, src=S3), below("1e-13", 1e-12, src=R)],
                    power=[doc("26.92", "MW", 1e6, S3)], s_power=[doc("1", src=S3)],
                    run_time=[factor("10", src=S3 + " ('about 10 s, 4 ranks')")]))
    add("03_stability", "03", moose(args, 4, *a, "orifice_implicit.i", "stability.i"), ranks=4,
        outputs=["stability_out.csv"], metrics=lambda j: stability_metrics(j["cwd"]),
        documented="open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i stability.i",
        note="frozen orifice in physical time, one slice kicked",
        checks=dict(peak_rise=[doc("142", "K", src=S3)], efold=[doc("2.0", "s", src=S3)],
                    T_fuel_peak=[doc("3345", "K", src=S3)], chi_exit=[doc("20.63", src=S3)]))
    add("03_fold_uniform", "03", moose(args, 1, *a, "s_design=0.3", "t_design=0.6", *FOLD_STOP), fold=True,
        outputs=["continuation_out.csv", "continuation_out.e"], metrics=lambda j: newton_metrics(j["csv"]),
        documented="open_pronghorn-opt -i anchor.i continuation.i s_design=0.3 t_design=0.6",
        note="stopped at the turning point by " + " ".join(FOLD_STOP),
        checks=dict(s_max=[doc("0.085", src=S3 + " (8.5 %)"), doc("0.085", src=R + " (8.5 %)")]))
    add("03_fold_chi30", "03", moose(args, 1, *a, "s_design=0.3", "t_design=0.6", "cold_frit_multiplier=30",
                                     "Outputs/file_base=fold_chi30_out", *FOLD_STOP), fold=True,
        outputs=["fold_chi30_out.csv"], metrics=lambda j: newton_metrics(j["csv"]),
        note="uniform throttling chi = 30 (README: fold at 3.3 %)",
        checks=dict(s_max=[doc("0.033", src=S3 + " (3.3 % at chi = 30)")]))
    add("03_fold_chi100", "03", moose(args, 1, *a, "s_design=0.3", "t_design=0.6", "cold_frit_multiplier=100",
                                      "Outputs/file_base=fold_chi100_out", *FOLD_STOP), fold=True,
        outputs=["fold_chi100_out.csv"], metrics=lambda j: newton_metrics(j["csv"]),
        note="uniform throttling chi = 100 (README: fold at 2.0 %)",
        checks=dict(s_max=[doc("0.020", src=S3 + " (2.0 % at chi = 100)")]))
    add("03_anchor_nz24", "03", moose(args, 4, *a, "orifice_implicit.i", "nz=24", "log_chi_order=TWENTYFOURTH",
                                      "Outputs/file_base=anchor_nz24_out"), ranks=4,
        outputs=["anchor_nz24_out.csv"], metrics=lambda j: newton_metrics(j["csv"]),
        documented="open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i nz=24 "
                   "log_chi_order=TWENTYFOURTH", note="mesh variant of the README; no documented values",
        checks=dict(s_power=[doc("1", src=S3 + " (reaches the design power)")]))

    # ---- 04 segregated Navier-Stokes (linear FV) at the design power, RZ
    seg = lambda j: segregated_metrics(j["csv"])  # noqa: E731
    add("04_simple_flow", "04", moose(args, 1, "-i", "anchor_simple.i", "orifice_flow.i", "orifice_verbose=false",
                                      "Outputs/file_base=anchor_simple_out"),
        outputs=["anchor_simple_out.csv", "anchor_simple_out.e"], metrics=seg,
        documented="open_pronghorn-opt -i anchor_simple.i orifice_flow.i",
        note="(c) SIMPLE, orifice designed on the slice flows",
        checks=dict(T_mixed=[doc("2989.81", "K", src=S4), doc("2989.8", "K", src=R)],
                    T_fuel_max=[doc("3167", "K", src=S4), doc("3167", "K", src=R)],
                    p_in=[doc("8.873", "MPa", 1e6, S4)], chi_min=[doc("0.41", src=S4)],
                    chi_max=[doc("46.3", src=S4)], T_slice_spread=[doc("1.1e-3", src=S4)],
                    energy_error=[below("-5e-9", 1e-8, src=S4)], mass_error=[below("3e-14", 1e-12, src=S4)],
                    s_power=[doc("1", src=S4)], run_time=[factor("112", src=S4 + " (1 rank)")]))
    add("04_pimple_flow", "04", moose(args, 1, "-i", "anchor_pimple.i", "orifice_flow.i", "orifice_verbose=false",
                                      "Outputs/file_base=anchor_pimple_out"),
        outputs=["anchor_pimple_out.csv", "anchor_pimple_out.e"], metrics=seg,
        documented="open_pronghorn-opt -i anchor_pimple.i orifice_flow.i",
        note="(c) PIMPLE, orifice designed on the slice flows",
        checks=dict(T_mixed=[doc("2989.81", "K", src=S4)], T_fuel_max=[doc("3167", "K", src=S4)],
                    p_in=[doc("8.872", "MPa", 1e6, S4)], chi_min=[doc("0.41", src=S4)],
                    chi_max=[doc("46.1", src=S4)], T_slice_spread=[doc("1.1e-3", src=S4)],
                    energy_error=[below("-7e-11", 1e-8, src=S4)], mass_error=[below("8e-14", 1e-12, src=S4)],
                    s_power=[doc("1", src=S4)], run_time=[factor("165", src=S4 + " (1 rank)")]))
    add("04_extract_orifice", "04", [*py, "extract_orifice.py", "../03_ludewig_anchor/anchor_out.csv"],
        needs=["03_anchor"], outputs=["orifice_profile.i"], metrics=lambda j: extract_orifice_metrics(j["log"]),
        documented="python3 extract_orifice.py ../03_ludewig_anchor/anchor_out.csv",
        checks=dict(chi_min=[doc("0.15", src=S3)], chi_max=[doc("20.6", src=S3)]))
    add("04_pimple_newton_orifice", "04", moose(args, 1, "-i", "anchor_pimple.i",
                                                "Outputs/file_base=anchor_pimple_newton_out"), placeholder=True,
        needs=["04_extract_orifice"], outputs=["anchor_pimple_newton_out.csv"], metrics=seg, runaway=True,
        documented="open_pronghorn-opt -i anchor_pimple.i",
        note="(a) finding: the Newton orifice starves the middle slices (README)",
        checks=dict(T_mixed=[doc("2990", "K", src=S4)], T_fuel_peak=[above("1e4", 1e4, "K", src=S4)]))
    lag = lambda j: dict(segregated_metrics(j["csv"]),  # noqa: E731
                         **at_power(j["csv"], (0.10, 0.25, 0.50), ("T_mixed", "T_fuel_max", "p_in", "chi_min",
                                                                   "chi_max")))
    add("04_pimple_lagged", "04", moose(args, 1, "-i", "anchor_pimple.i", "orifice_lagged.i", "s_final=0.5",
                                        "t_ramp=30", "t_end=60", "pimple_steady_state_detection=false",
                                        "Outputs/file_base=orifice_lagged_out"),
        placeholder=True, outputs=["orifice_lagged_out.csv"], metrics=lag, runaway=True,
        documented="open_pronghorn-opt -i anchor_pimple.i orifice_lagged.i s_final=0.5 t_ramp=30 t_end=60 "
                   "pimple_steady_state_detection=false",
        note="(b) finding: lagged temperature design along the ramp to 50 %",
        checks={"T_mixed@10%": [doc("384", "K", src=S4)], "T_mixed@25%": [doc("872", "K", src=S4)],
                "T_mixed@50%": [doc("1634", "K", src=S4)], "T_fuel_max@10%": [doc("414", "K", src=S4)],
                "T_fuel_max@25%": [doc("935", "K", src=S4)], "T_fuel_max@50%": [doc("1739", "K", src=S4)],
                "p_in@10%": [doc("7.26", "MPa", 1e6, S4)], "p_in@25%": [doc("7.66", "MPa", 1e6, S4)],
                "p_in@50%": [doc("8.34", "MPa", 1e6, S4)], "chi_min@10%": [doc("17.9", src=S4)],
                "chi_max@10%": [doc("94", src=S4)], "chi_min@50%": [doc("24.3", src=S4)],
                "chi_max@50%": [doc("77", src=S4)], "T_fuel_peak": [above("1e4", 1e4, "K", src=S4)]})

    # ---- 05 3-D
    t3 = ["-i", "anchor_3d.i", "../03_ludewig_anchor/continuation.i", "../03_ludewig_anchor/orifice_implicit.i",
          "tolerance_3d.i"]
    add("05_mesh", "05", moose(args, 1, "-i", "mesh_3d.i", "--mesh-only", "--allow-unused"),
        outputs=["mesh_3d_in.e"], documented="open_pronghorn-opt -i mesh_3d.i --mesh-only --allow-unused")
    add("05_anchor_3d", "05", moose(args, 12, *t3, "Outputs/file_base=anchor_3d_out"), ranks=12,
        outputs=["anchor_3d_out.csv", "anchor_3d_out.e"],
        metrics=lambda j: with_transverse(j["csv"], j["exodus"], args.python),
        documented="mpirun -n 12 open_pronghorn-opt -i anchor_3d.i ../03_ludewig_anchor/continuation.i "
                   "../03_ludewig_anchor/orifice_implicit.i tolerance_3d.i Outputs/file_base=anchor_3d_out",
        checks=dict(T_mixed=[doc("2989.4", "K", src=S5)], T_fuel_max=[doc("3197.1", "K", src=S5)],
                    p_in=[doc("7.836", "MPa", 1e6, S5)], chi_exit=[doc("20.63", src=S5)],
                    energy_error=[doc("3e-7", src=S5)], s_power=[doc("1", src=S5 + " (tilt 0: 1.00)")],
                    fuel_rise_ratio=[doc("1.07", src=S5 + " (no tilt)")],
                    faceting_pc=[doc("2.5", "%", src=S5 + " (faceted bed 2.5 % smaller, n_theta = 8)")],
                    run_time=[factor("40", src=S5 + " ('~40 s', 12 ranks)")]))
    add("05_vs_03", "05", None, needs=["03_anchor", "05_anchor_3d"], note="3-D vs axisymmetric anchor",
        metrics=lambda j: compare_3d_2d(os.path.join(HERE, "05_element_3d", "anchor_3d_out.csv"),
                                        os.path.join(HERE, "03_ludewig_anchor", "anchor_out.csv")),
        checks=dict(fuel_diff_pc=[doc("0.12", "%", src=S5), doc("0.12", "%", src=R)],
                    p_in_diff_pc=[doc("0.06", "%", src=S5)]))
    add("05_tilt_005", "long", moose(args, 8, *t3, "tilt=0.05", "Outputs/file_base=anchor_3d_tilt_out",
                                     *FOLD_STOP), ranks=8, fold=True,
        outputs=["anchor_3d_tilt_out.csv", "anchor_3d_tilt_out.e"],
        metrics=lambda j: with_transverse(j["csv"], j["exodus"], args.python),
        documented="mpirun -n 8 open_pronghorn-opt -i anchor_3d.i ../03_ludewig_anchor/continuation.i "
                   "../03_ludewig_anchor/orifice_implicit.i tolerance_3d.i tilt=0.05 "
                   "Outputs/file_base=anchor_3d_tilt_out",
        note="stopped at the turning point by " + " ".join(FOLD_STOP),
        checks=dict(s_max=[doc("0.178", src=S5), doc("0.18", src=R + " (18 %)")],
                    gas_rise_max=[doc("0.47", src=S5 + " (+47 %)"), doc("0.47", src=R + " (+47 %)")],
                    gas_rise_min=[doc("-0.24", src=S5 + " (-24 %)"), doc("-0.24", src=R + " (-24 %)")],
                    fuel_rise_ratio=[doc("2.2", src=S5)]))
    add("05_tilt_015", "long", moose(args, 8, *t3, "tilt=0.15", "Outputs/file_base=anchor_3d_tilt015_out",
                                     *FOLD_STOP), ranks=8, fold=True,
        outputs=["anchor_3d_tilt015_out.csv", "anchor_3d_tilt015_out.e"],
        metrics=lambda j: newton_metrics(j["csv"]),
        documented="(table row only; same command with tilt=0.15)",
        note="stopped at the turning point by " + " ".join(FOLD_STOP),
        checks=dict(s_max=[doc("0.050", src=S5)]))
    add("05_tilt_001", "long", moose(args, 8, *t3, "tilt=0.01", "Outputs/file_base=anchor_3d_tilt001_out",
                                     *FOLD_STOP), ranks=8, fold=True,
        outputs=["anchor_3d_tilt001_out.csv"], metrics=lambda j: newton_metrics(j["csv"]),
        documented="(table row only; same command with tilt=0.01)",
        note="stopped at the turning point by " + " ".join(FOLD_STOP),
        checks=dict(s_max=[doc("0.513", src=S5)]))
    # segregated (linear FV) model in 3-D
    seg3 = lambda j: with_asymmetry(j["csv"], j["exodus"], args.python)  # noqa: E731
    add("05_simple_cold_3d", "05", moose(args, 8, "-i", "element_simple_3d.i",
                                         "Outputs/file_base=element_simple_3d_out"), ranks=8,
        outputs=["element_simple_3d_out.csv", "element_simple_3d_out.e"], metrics=seg,
        documented="mpirun -n 8 open_pronghorn-opt -i element_simple_3d.i",
        note="cold flow from the plug-flow start",
        checks=dict(dp_in=[doc("11.97", "kPa", 1e3, S5)], mass_error=[below("1e-10", 1e-9, src=S5)]))
    add("05_simple_3d", "long", moose(args, 8, "-i", "anchor_simple_3d.i", "../04_navier_stokes/orifice_flow.i",
                                      "orifice_sectors=8", "orifice_verbose=false",
                                      "Outputs/file_base=anchor_simple_3d_out"), ranks=8,
        outputs=["anchor_simple_3d_out.csv", "anchor_simple_3d_out.e"], metrics=seg3,
        documented="mpirun -n 8 open_pronghorn-opt -i anchor_simple_3d.i ../04_navier_stokes/orifice_flow.i "
                   "orifice_sectors=8",
        note="SIMPLE, column orifice (12 slices x 8 sectors)",
        checks=dict(T_mixed=[doc("2989.42", "K", src=S5)], T_fuel_max=[doc("3170", "K", src=S5)],
                    p_in=[doc("8.900", "MPa", 1e6, S5)], T_slice_spread=[doc("1.1e-3", src=S5)],
                    energy_error=[below("1e-8", 1e-7, src=S5)], mass_error=[below("6e-12", 1e-10, src=S5)],
                    asymmetry=[below("2e-7", 1e-6, src=S5)], s_power=[doc("1", src=S5)],
                    run_time=[factor("300", src=S5 + " (8 ranks)")]))
    add("05_pimple_3d", "long", moose(args, 8, "-i", "anchor_pimple_3d.i", "../04_navier_stokes/orifice_flow.i",
                                      "orifice_sectors=8", "orifice_verbose=false",
                                      "Outputs/file_base=anchor_pimple_3d_out"), ranks=8,
        outputs=["anchor_pimple_3d_out.csv", "anchor_pimple_3d_out.e"], metrics=seg3,
        documented="mpirun -n 8 open_pronghorn-opt -i anchor_pimple_3d.i ../04_navier_stokes/orifice_flow.i "
                   "orifice_sectors=8",
        note="PIMPLE from rest, column orifice",
        checks=dict(T_mixed=[doc("2989.42", "K", src=S5)], T_fuel_max=[doc("3170", "K", src=S5)],
                    p_in=[doc("8.899", "MPa", 1e6, S5)], T_slice_spread=[doc("1.0e-3", src=S5)],
                    energy_error=[below("1e-8", 1e-7, src=S5)], asymmetry=[below("2e-7", 1e-6, src=S5)],
                    s_power=[doc("1", src=S5)], run_time=[factor("470", src=S5 + " (8 ranks)")]))
    add("05_pimple_3d_rings", "long", moose(args, 8, "-i", "anchor_pimple_3d.i",
                                            "../04_navier_stokes/orifice_flow.i", "t_end=10",
                                            "pimple_steady_state_detection=false", "orifice_verbose=false",
                                            "Outputs/file_base=anchor_pimple_3d_rings_out"), ranks=8,
        outputs=["anchor_pimple_3d_rings_out.csv", "anchor_pimple_3d_rings_out.e"], metrics=seg3,
        documented="mpirun -n 8 open_pronghorn-opt -i anchor_pimple_3d.i ../04_navier_stokes/orifice_flow.i "
                   "t_end=10 pimple_steady_state_detection=false",
        note="ring orifice: the azimuthal flow excursion",
        checks=dict(asymmetry=[above("1", 1.0, src=S5)], T_fuel_peak=[above("1e4", 1e4, "K", src=S5)]))
    add("05_plot_3d", "long", [*py, "plot_3d.py", "anchor_3d_out.e:no tilt", "anchor_3d_tilt_out.e:tilt 0.05"],
        needs=["05_anchor_3d", "05_tilt_005"], outputs=["transverse.png"], exodus_reader=True,
        documented='python3 plot_3d.py anchor_3d_out.e:"no tilt" anchor_3d_tilt_out.e:"tilt 0.05"')

    for j in J:
        j["log"] = os.path.join(LOGS, j["name"] + ".log")
        fb = next((c.split("=", 1)[1] for c in (j["cmd"] or []) if c.startswith("Outputs/file_base=")), None)
        inputs = [c for c in (j["cmd"] or []) if c.endswith(".i")]
        if fb is None and inputs and "--mesh-only" not in j["cmd"]:
            fb = os.path.splitext(os.path.basename(inputs[-1]))[0] + "_out"
        j["csv"] = os.path.join(j["cwd"], fb + ".csv") if fb else None
        j["exodus"] = os.path.join(j["cwd"], fb + ".e") if fb else None
        j["ranks"] = max(1, round(j["ranks"] * args.ranks_scale)) if j["ranks"] > 1 else 1
        if j["ranks"] > 1 and j["cmd"]:
            k = j["cmd"].index("-n")
            j["cmd"][k + 1] = str(j["ranks"])
    return J


# ----------------------------------------------------------------------------------------------
# running
# ----------------------------------------------------------------------------------------------
def shell_line(j, args=None):
    """The job's command from its stage directory; with 'args', the executable and the interpreter
    are shown by name (open_pronghorn-opt, python3) so that the report holds no local paths."""
    cmd = list(j["cmd"])
    if args is not None:
        names = {args.exe: "open_pronghorn-opt"}
        names.update({p: "python3" for p in shlex.split(args.python)})
        cmd = [names.get(c, c) for c in cmd]
    return f"cd {os.path.relpath(j['cwd'], HERE)} && " + " ".join(shlex.quote(c) for c in cmd)


def has_netcdf4(python):
    cmd = shlex.split(python) + ["-c", "import netCDF4, matplotlib, numpy"]
    return subprocess.run(cmd, capture_output=True).returncode == 0


def run_jobs(jobs, args, state):
    pending = [j for j in jobs if j["cmd"]]
    if not has_netcdf4(args.python):
        # the plotting scripts read Exodus files with netCDF4
        for j in [j for j in pending if j["exodus_reader"]]:
            state[j["name"]] = dict(status="skipped", wall=0.0, rc=None,
                                    why=f"{args.python} has no netCDF4 (pass --python)")
            print(f"[skip] {j['name']}: {args.python} cannot import netCDF4; pass --python", flush=True)
            pending.remove(j)
    running, done = {}, set(j["name"] for j in jobs if not j["cmd"] or j["name"] in state)
    names = set(j["name"] for j in jobs)
    os.makedirs(LOGS, exist_ok=True)
    while pending or running:
        for j in list(pending):
            unmet = [n for n in j["needs"] if n in names and n not in done]
            failed = [n for n in j["needs"] if state.get(n, {}).get("status") in ("failed", "skipped")]
            if failed:
                state[j["name"]] = dict(status="skipped", wall=0.0, rc=None, why="needs " + " ".join(failed))
                print(f"[skip] {j['name']}: needs {' '.join(failed)}", flush=True)
                pending.remove(j)
                done.add(j["name"])
                continue
            ranks_used = sum(r["ranks"] for r in running.values())
            if unmet or len(running) >= args.jobs or (running and ranks_used + j["ranks"] > args.max_ranks):
                continue
            log = open(j["log"], "w")
            log.write(f"# {shell_line(j)}\n")
            log.flush()
            p = subprocess.Popen(j["cmd"], cwd=j["cwd"], stdout=log, stderr=subprocess.STDOUT)
            running[j["name"]] = dict(job=j, proc=p, start=time.time(), logf=log, ranks=j["ranks"])
            pending.remove(j)
            print(f"[start] {j['name']} ({j['ranks']} rank{'s' if j['ranks'] > 1 else ''})", flush=True)
        for name, r in list(running.items()):
            j, p = r["job"], r["proc"]
            wall = time.time() - r["start"]
            if p.poll() is None and wall > args.timeout:
                p.kill()
            if p.poll() is None:
                continue
            r["logf"].close()
            with open(j["log"]) as f:
                text = f.read()
            status = "ok" if p.returncode == 0 else "failed"
            if j["fold"] and FOLD_MARK in text and j["csv"] and os.path.exists(j["csv"]):
                status = "ok (turning point)"
            if j.get("runaway") and RUNAWAY_MARK in text and j["csv"] and os.path.exists(j["csv"]):
                status = "ok (runaway, as documented)"
            if wall > args.timeout:
                status = "failed"
            state[name] = dict(status=status, wall=wall, rc=p.returncode)
            print(f"[{status}] {name} in {wall:.1f} s (rc {p.returncode}) -> {os.path.relpath(j['log'], HERE)}",
                  flush=True)
            del running[name]
            done.add(name)
        time.sleep(0.5)
    for j in jobs:
        if not j["cmd"]:
            state[j["name"]] = dict(status="ok", wall=0.0, rc=0)


# ----------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------
def fmt(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    if x == 0 or 1e-3 <= abs(x) < 1e5:
        return f"{x:.6g}"
    return f"{x:.3e}"


def report(jobs, state, args):
    rows, table = [], []
    for j in jobs:
        st = state.get(j["name"], dict(status="not run", wall=float("nan")))
        metrics = {}
        if j["metrics"] and st["status"].startswith("ok"):
            try:
                metrics = j["metrics"](j)
            except Exception as e:  # a missing output is reported, not fatal
                metrics = {}
                st = dict(st, status=f"no metrics ({type(e).__name__}: {e})")
        if j["cmd"]:
            metrics["run_time"] = st.get("wall", float("nan"))
        j["status"] = st["status"]
        for q in sorted(set(metrics) | set(j["checks"])):
            value = metrics.get(q)
            for c in j["checks"].get(q, [None]):
                if c is None:
                    rows.append(dict(job=j["name"], quantity=q, value=value, documented="", tolerance="",
                                     status="", source=""))
                    continue
                status = evaluate(c, value) if st["status"].startswith("ok") else "not run"
                rows.append(dict(job=j["name"], quantity=q, value=value, documented=c["text"], tolerance=tol_text(c),
                                 status=status, source=c["src"], scale=c["scale"]))
                table.append(rows[-1])
    with open(os.path.join(HERE, "suite_results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["job", "quantity", "value", "documented", "tolerance", "status", "source"],
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    n_bad = sum(r["status"] in ("MISMATCH", "missing") for r in table)
    n_failed = sum(not j["status"].startswith(("ok", "not run", "skipped")) for j in jobs)
    n_skipped = sum(j["status"].startswith("skipped") for j in jobs)
    lines = ["# NTP tutorial: demonstration suite results", "",
             f"Written by `run_suite.py` on {time.strftime('%Y-%m-%d %H:%M')}, executable "
             f"`{os.path.relpath(args.exe, HERE)}`. "
             f"{len(table)} checks against the READMEs: "
             f"{sum(r['status'] == 'ok' for r in table)} match, {n_bad} mismatch or missing, "
             f"{sum(r['status'] == 'not run' for r in table)} not run; {n_failed} job(s) failed, "
             f"{n_skipped} skipped.", "",
             "Tolerance: a documented number matches when the fresh value rounds to it; round-off "
             "quantities are upper bounds; run times within a factor of three. The CSV "
             "`suite_results.csv` holds every extracted figure (SI units) including the unchecked ones.",
             "", "## Jobs", "", "| job | group | status | ranks | wall [s] | command (from the stage directory) |",
             "|---|---|---|---|---|---|"]
    for j in jobs:
        st = state.get(j["name"], {})
        cmd = shell_line(j, args) if j["cmd"] else "(comparison of outputs)"
        note = f" — {j['note']}" if j["note"] else ""
        lines.append(f"| {j['name']} | {j['group']} | {j['status']} | {j['ranks']} | {fmt(st.get('wall'))} | "
                     f"`{cmd.replace('|', '/')}`{note} |")
    lines += ["", "## Figures of merit vs the READMEs", "",
              "| job | quantity | this run | documented | tolerance | status | source |", "|---|---|---|---|---|---|---|"]
    for r in table:
        v = r["value"] / r["scale"] if isinstance(r["value"], float) and r["scale"] else r["value"]
        flag = f"**{r['status']}**" if r["status"] in ("MISMATCH", "missing") else r["status"]
        lines.append(f"| {r['job']} | {r['quantity']} | {fmt(v)} | {r['documented']} | {r['tolerance']} | "
                     f"{flag} | {r['source']} |")
    with open(os.path.join(HERE, "suite_results.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n{len(table)} checks: {sum(r['status'] == 'ok' for r in table)} ok, {n_bad} mismatch/missing; "
          f"{n_failed} failed, {n_skipped} skipped job(s) -> suite_results.md, suite_results.csv")
    for r in table:
        if r["status"] in ("MISMATCH", "missing"):
            v = r["value"] / r["scale"] if isinstance(r["value"], float) else r["value"]
            print(f"  {r['status']}: {r['job']} {r['quantity']} = {fmt(v)}, documented {r['documented']} "
                  f"({r['source']})")
    return n_bad == 0 and n_failed == 0 and n_skipped == 0


def select(jobs, args):
    groups = set(GROUPS[:5])
    if args.only:
        groups = set(g.strip() for g in args.only.split(","))
        unknown = groups - set(GROUPS)
        if unknown:
            sys.exit(f"unknown group(s) {' '.join(sorted(unknown))}; groups: {' '.join(GROUPS)}")
    if args.long:
        groups.add("long")
    chosen = [j for j in jobs if j["group"] in groups and (args.placeholders or not j["placeholder"])]
    if args.jobs_re:
        chosen = [j for j in chosen if re.search(args.jobs_re, j["name"])]
    return chosen


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--exe", default=os.path.join(HERE, "..", "..", "open_pronghorn-opt"))
    parser.add_argument("--mpirun", default="mpirun --bind-to none")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--only", help="comma-separated groups: " + ",".join(GROUPS))
    parser.add_argument("--long", action="store_true", help="add the 'long' group (3-D power tilts)")
    parser.add_argument("--placeholders", action="store_true", help="also run the open-work stage-4 PIMPLE jobs")
    parser.add_argument("--select", dest="jobs_re", help="regular expression on the job names")
    parser.add_argument("--jobs", type=int, default=3, help="concurrent jobs")
    parser.add_argument("--max-ranks", type=int, default=24, help="largest total number of MPI ranks at once")
    parser.add_argument("--ranks-scale", type=float, default=1.0, help="scale the documented rank counts")
    parser.add_argument("--timeout", type=float, default=3600, help="seconds per job")
    parser.add_argument("--list", action="store_true", help="list the jobs and exit")
    parser.add_argument("--dry-run", action="store_true", help="print the commands and exit")
    parser.add_argument("--report-only", action="store_true", help="only extract from existing outputs")
    args = parser.parse_args()
    args.exe = os.path.abspath(args.exe)
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"  # no __pycache__ in the stage directories

    jobs = select(job_list(args), args)
    if args.list:
        for j in job_list(args):
            chosen = j["name"] in {c["name"] for c in jobs}
            tags = [t for t, on in (("placeholder", j["placeholder"]), ("selected", chosen)) if on]
            print(f"{j['name']:26s} {j['group']:5s} {j['ranks']:3d} rank(s)  {', '.join(tags)}")
        return 0
    if args.dry_run:
        for j in jobs:
            print(f"# {j['name']} ({j['group']}, {j['ranks']} ranks"
                  f"{', needs ' + ' '.join(j['needs']) if j['needs'] else ''})")
            print(shell_line(j) if j["cmd"] else "(comparison of outputs)")
        return 0
    if not os.path.exists(args.exe):
        sys.exit(f"executable not found: {args.exe} (use --exe)")

    state_path = os.path.join(LOGS, "state.json")
    state = {}
    if os.path.exists(state_path):
        with open(state_path) as f:
            state = json.load(f)
    if not args.report_only:
        t0 = time.time()
        for j in jobs:  # the selected jobs run again
            state.pop(j["name"], None)
        run_jobs(jobs, args, state)
        print(f"suite wall time {time.time() - t0:.0f} s")
        os.makedirs(LOGS, exist_ok=True)
        with open(state_path, "w") as f:
            json.dump(state, f, indent=1)
    else:
        for j in jobs:
            if j["name"] not in state:
                have = not j["cmd"] or all(os.path.exists(os.path.join(j["cwd"], o)) for o in j["outputs"])
                state[j["name"]] = dict(status="ok" if have else "not run", wall=float("nan"))
    # the report covers every job with a recorded state, so that a partial run (e.g. --only 04)
    # updates its rows of the full table
    chosen = {j["name"] for j in jobs}
    reported = [j for j in job_list(args) if j["name"] in chosen or j["name"] in state]
    return 0 if report(reported, state, args) else 1


if __name__ == "__main__":
    sys.exit(main())
