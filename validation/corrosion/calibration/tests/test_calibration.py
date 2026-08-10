"""Reproducibility tests for the in-repository corrosion/plating calibration."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

import numpy as np


CALIBRATION_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = CALIBRATION_ROOT.parents[2]
SRC = CALIBRATION_ROOT / "src"
WORKBOOK = CALIBRATION_ROOT / "data" / "msr_corrosion_plating_validation_data.xlsx"
REFERENCE = CALIBRATION_ROOT / "reference"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from msr_corrosion_bv.calibrate import (  # noqa: E402
    ACTIVE_FIT_ROLES,
    fit_model,
    objective_vector,
)
from msr_corrosion_bv.ingest import build_model_tables  # noqa: E402
from msr_corrosion_bv.model import (  # noqa: E402
    PARAMETER_SPECS,
    initial_parameter_vector,
)


def read_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


class CalibrationDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tables = build_model_tables(WORKBOOK)

    def test_workbook_inventory_and_active_sources(self) -> None:
        targets = self.tables["targets"]
        active = targets[targets["fit_role"].isin(ACTIVE_FIT_ROLES)]

        self.assertEqual(len(self.tables["case_features"]), 76)
        self.assertEqual(len(targets), 43)
        self.assertEqual(len(active), 37)
        self.assertEqual(active["source_id"].nunique(), 14)

    def test_objective_contains_data_and_prior_residuals(self) -> None:
        targets = self.tables["targets"]
        objective = objective_vector(initial_parameter_vector(), targets)

        # Thirty-seven active rows produce 38 data residuals because the one
        # deposition-ranking row encodes two inequalities. Every one of the
        # 59 parameters contributes a regularizing prior residual.
        self.assertEqual(len(PARAMETER_SPECS), 59)
        self.assertEqual(len(objective), 38 + 59)

    def test_frozen_parameters_match_existing_vendored_copy(self) -> None:
        frozen = read_json(REFERENCE / "parameters.json")
        vendored = read_json(
            REPOSITORY_ROOT / "validation" / "corrosion" / "data" / "parameters.json"
        )
        self.assertEqual(frozen, vendored)

    def test_production_database_has_only_documented_extension(self) -> None:
        frozen = read_json(REFERENCE / "parameters.json")
        production = read_json(REPOSITORY_ROOT / "data" / "corrosion_database.json")[
            "calibrated_parameters"
        ]

        self.assertEqual(set(production) - set(frozen), {"log_ncl16_cr_inventory_bonus"})
        self.assertEqual(set(frozen) - set(production), set())
        for name, value in frozen.items():
            self.assertEqual(value, production[name], msg=name)


class CalibrationFitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        # Run the optimizer once for this class; individual tests inspect the
        # shared in-memory result without writing repository files.
        cls.outputs = fit_model(WORKBOOK, output_dir=None)
        cls.expected_parameters = read_json(REFERENCE / "parameters.json")
        cls.expected_metrics = read_json(REFERENCE / "metrics.json")

    def test_optimizer_converges(self) -> None:
        metrics = self.outputs["metrics"]
        self.assertTrue(metrics["optimizer_success"])
        self.assertEqual(metrics["optimizer_nfev"], 21)

    def test_fitted_parameters_reproduce_frozen_vector(self) -> None:
        observed = self.outputs["model"].params
        differences = {
            name: abs(float(observed[name]) - float(expected))
            for name, expected in self.expected_parameters.items()
        }
        worst = max(differences, key=differences.get)
        self.assertLessEqual(differences[worst], 1.0e-3, msg=f"largest drift: {worst}")

    def test_core_metrics_reproduce_frozen_run(self) -> None:
        observed = self.outputs["metrics"]
        expected = self.expected_metrics

        for name in (
            "n_measurements_total",
            "n_active_constraints",
            "n_direct_or_range_targets",
            "n_upper_lower_ranking_constraints",
        ):
            self.assertEqual(observed[name], expected[name], msg=name)

        for name in (
            "median_factor_error_direct",
            "within_factor_2_direct",
            "within_factor_5_direct",
            "constraint_pass_fraction",
        ):
            self.assertTrue(
                np.isclose(observed[name], expected[name], atol=1.0e-6, rtol=0.0),
                msg=name,
            )


if __name__ == "__main__":
    unittest.main()
