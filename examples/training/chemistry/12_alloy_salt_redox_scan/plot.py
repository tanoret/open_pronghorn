#!/usr/bin/env python3
"""Case 12: the parity plot is produced by run_scan.py (it reruns the input seven times)."""
import runpy
import pathlib

runpy.run_path(str(pathlib.Path(__file__).resolve().parent / "run_scan.py"), run_name="__main__")
