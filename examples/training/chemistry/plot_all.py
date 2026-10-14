#!/usr/bin/env python3
"""Regenerate every training figure (runs each plot.py; inputs are run first if outputs are missing).

Usage:  python plot_all.py   (set OPEN_PRONGHORN_EXE if the executable is not at the repo root)
"""
import pathlib
import runpy

for script in sorted(pathlib.Path(__file__).resolve().parent.glob("[0-9][0-9]_*/plot.py")):
    print(f"--- {script.parent.name}")
    runpy.run_path(str(script), run_name="__main__")
