from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path

import pandas as pd


CALIBRATION_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = CALIBRATION_ROOT.parents[2]
ADVANCED_RESULTS = CALIBRATION_ROOT / "results" / "advanced"


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class AdvancedArtifactConsistencyTests(unittest.TestCase):
    def test_combined_database_matches_scalar_parameter_artifacts(self) -> None:
        combined = _json(REPOSITORY_ROOT / "data" / "advanced_corrosion_models.json")
        thermochemical = _json(ADVANCED_RESULTS / "thermochemical_parameters.json")
        dynamic = _json(ADVANCED_RESULTS / "dynamic_network_parameters.json")
        self.assertEqual(combined["schema_version"], "1.0")
        self.assertEqual(
            combined["model_revision"],
            "mstdb-nst-v2-fe-identity_dridn-v2-explicit-geometry",
        )
        self.assertEqual(
            combined["calibration_id"],
            "pr-corrosion-advanced-explicit-context-2026-08-28-v5",
        )
        self.assertEqual(
            combined["calibration_definition"]["case_context_revision"],
            "explicit-geometry-v2",
        )
        self.assertEqual(
            combined["calibration_definition"]["case_context_sha256"],
            _sha256(ADVANCED_RESULTS / "advanced_case_context.csv"),
        )
        self.assertEqual(combined["thermochemical_parameters"], thermochemical)
        self.assertEqual(combined["dynamic_parameters"], dynamic)
        self.assertEqual(
            combined["base_model_provenance"]["source_sha256"],
            "bfad8e09ada73900b800c4df262ae70192296e6c7b95cf5fc20943272637d0e6",
        )
        self.assertIn(
            "base-sha256=bfad8e09ada73900b800c4df262ae70192296e6c7b95cf5fc20943272637d0e6",
            combined["calibration_data_revision"],
        )
        self.assertEqual(
            set(combined["base_model_provenance"]["required_element_properties"]),
            {"Cr", "Fe", "Ni"},
        )

        thermo_csv = pd.read_csv(ADVANCED_RESULTS / "thermochemical_parameters.csv")
        dynamic_csv = pd.read_csv(ADVANCED_RESULTS / "dynamic_network_parameters.csv")
        thermo_from_csv = dict(zip(thermo_csv["parameter"], thermo_csv["value"].astype(float)))
        dynamic_from_csv = dict(zip(dynamic_csv["parameter"], dynamic_csv["value"].astype(float)))
        self.assertEqual(set(thermo_from_csv), set(thermochemical))
        self.assertEqual(set(dynamic_from_csv), set(dynamic))
        for name, value in thermochemical.items():
            self.assertAlmostEqual(thermo_from_csv[name], value, places=14, msg=name)
        for name, value in dynamic.items():
            self.assertAlmostEqual(dynamic_from_csv[name], value, places=14, msg=name)

    def test_combined_metrics_match_checked_reports(self) -> None:
        combined = _json(REPOSITORY_ROOT / "data" / "advanced_corrosion_models.json")
        thermo = _json(ADVANCED_RESULTS / "thermochemical_comparison_metrics.json")
        dynamic = _json(ADVANCED_RESULTS / "three_model_comparison_metrics.json")
        checked = combined["checked_metrics"]
        self.assertEqual(
            checked["thermochemical_midpoint_log_rmse"],
            thermo["thermochemical_log_rmse_direct"],
        )
        self.assertEqual(
            checked["thermochemical_relation_aware_log_rmse"],
            thermo["thermochemical_relation_aware_log_rmse_direct"],
        )
        self.assertEqual(
            checked["dridn_midpoint_log_rmse"],
            dynamic["dynamic_network_log_rmse_direct"],
        )
        self.assertEqual(
            checked["dridn_relation_aware_log_rmse"],
            dynamic["dynamic_network_range_aware_log_rmse_direct"],
        )
        self.assertEqual(
            checked["dridn_maximum_mass_balance_relative_error"],
            dynamic["endpoint_mass_balance"]["maximum_relative_error"],
        )

    def test_all_provenance_hashes_match_bytes(self) -> None:
        provenance = _json(CALIBRATION_ROOT / "advanced_provenance.json")
        for relative, expected in provenance["inputs"].items():
            path = REPOSITORY_ROOT / relative
            if relative == "data/corrosion_database.json" and not path.is_file():
                # This overlay omits the unchanged base-branch file. A merged-tree
                # run verifies its bytes; an overlay-only run verifies the exact
                # semantic binding embedded in the advanced database.
                combined = _json(REPOSITORY_ROOT / "data" / "advanced_corrosion_models.json")
                self.assertEqual(
                    combined["base_model_provenance"]["source_sha256"], expected
                )
                continue
            self.assertEqual(_sha256(path), expected, relative)
        for relative, expected in provenance["curated_outputs"].items():
            path = REPOSITORY_ROOT / relative
            self.assertEqual(_sha256(path), expected, relative)

    def test_curated_manifest_covers_every_checked_advanced_result(self) -> None:
        provenance = _json(CALIBRATION_ROOT / "advanced_provenance.json")
        prefix = "validation/corrosion/calibration/results/advanced/"
        checked_results = {
            f"{prefix}{path.name}"
            for path in ADVANCED_RESULTS.iterdir()
            if path.is_file() and path.suffix in {".csv", ".json"}
        }
        manifested_results = {
            relative
            for relative in provenance["curated_outputs"]
            if relative.startswith(prefix)
        }
        self.assertEqual(len(checked_results), 18)
        self.assertEqual(manifested_results, checked_results)

    def test_cpp_parity_case_set_is_frozen(self) -> None:
        parity = pd.read_csv(ADVANCED_RESULTS / "cpp_parity_cases.csv")
        context = pd.read_csv(ADVANCED_RESULTS / "advanced_case_context.csv")
        explicit_fields = [
            "initial_dissolved_Cr_ppm",
            "initial_dissolved_Fe_ppm",
            "initial_dissolved_Ni_ppm",
            "explicit_inventory_scale",
            "transient_redox",
            "stress_interfacial_activation",
            "fluoride_impurity_interfacial_activation",
            "chloride_salt",
        ]
        self.assertTrue(set(explicit_fields).issubset(parity.columns))
        self.assertEqual(
            parity["measurement_id"].tolist(),
            [
                "M-003", "M-005", "M-014", "M-018", "M-027",
                "M-029", "M-030", "M-038", "M-041", "boundary_zero_elements",
            ],
        )
        self.assertTrue(parity["mass_balance_relative_error"].lt(1.0e-10).all())
        measured = parity[parity["measurement_id"] != "boundary_zero_elements"].copy()
        expected = context.set_index("measurement_id").loc[measured["measurement_id"]]
        expected.index = measured.index
        for field in explicit_fields:
            if field in {
                "transient_redox",
                "stress_interfacial_activation",
                "fluoride_impurity_interfacial_activation",
                "chloride_salt",
            }:
                actual_bool = measured[field].astype(str).str.strip().str.lower()
                expected_bool = expected[field].astype(str).str.strip().str.lower()
                self.assertTrue(actual_bool.eq(expected_bool).all(), field)
            else:
                difference = (measured[field].astype(float) - expected[field].astype(float)).abs()
                self.assertTrue(difference.le(1.0e-15).all(), field)
        boundary = parity.iloc[-1]
        self.assertEqual(
            [float(boundary[f"mass_fraction_{element}"]) for element in ("Cr", "Fe", "Ni")],
            [1.0, 0.0, 0.0],
        )
        self.assertEqual(float(boundary["cumulative_source_Fe_ppm"]), 0.0)
        self.assertEqual(float(boundary["cumulative_source_Ni_ppm"]), 0.0)


if __name__ == "__main__":
    unittest.main()
