from __future__ import annotations

import json
import math
import os
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from msr_corrosion_bv.calibrate import load_model_from_json
from msr_corrosion_bv.advanced_context import merge_advanced_case_context
from msr_corrosion_bv.dynamic_network import DynamicRedoxInventoryDepletionModel
from msr_corrosion_bv.dynamic_network_fit import dynamic_physics_residuals
from msr_corrosion_bv.mechanistic import mechanistic_is_supported
from msr_corrosion_bv.mstdb import MSTDBPair
from msr_corrosion_bv.thermochemical_model import MSTDBThermochemicalCorrosionModel


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(os.environ.get("MSTDB_TC_DIR"), "set MSTDB_TC_DIR for DRIDN integration tests")
class DynamicNetworkIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        effective = load_model_from_json(ROOT / "data" / "parameters.json")
        mstdb = MSTDBPair.from_calibration_directory(os.environ["MSTDB_TC_DIR"])
        thermo_params = json.loads(
            (ROOT / "results" / "advanced" / "thermochemical_parameters.json").read_text(encoding="utf-8")
        )
        dynamic_params = json.loads(
            (ROOT / "results" / "advanced" / "dynamic_network_parameters.json").read_text(encoding="utf-8")
        )
        cls.thermochemical = MSTDBThermochemicalCorrosionModel(effective, mstdb, thermo_params)
        cls.model = DynamicRedoxInventoryDepletionModel(
            cls.thermochemical,
            dynamic_params,
            integration_steps=120,
        )
        cls.targets = merge_advanced_case_context(
            pd.read_csv(ROOT / "data" / "targets.csv"),
            ROOT / "results" / "advanced" / "advanced_case_context.csv",
        )
        cls.supported = cls.targets[cls.targets.apply(mechanistic_is_supported, axis=1)].copy()

    def test_all_supported_predictions_are_positive_and_element_balanced(self) -> None:
        self.assertEqual(len(self.supported), 32)
        maximum_balance_error = 0.0
        for _, row in self.supported.iterrows():
            prediction = float(self.model.predict_response(row))
            self.assertTrue(math.isfinite(prediction), row.get("measurement_id"))
            self.assertGreater(prediction, 0.0, row.get("measurement_id"))
            result = self.model.simulate(row)
            maximum_balance_error = max(maximum_balance_error, result.mass_balance_relative_error)
            for collection in (
                result.dissolved_ppm,
                result.cumulative_source_ppm,
                result.deposit_mg_cm2,
                result.bulk_captured_ppm,
                result.surface_availability,
            ):
                self.assertTrue(all(math.isfinite(float(value)) for value in collection.values()))
                self.assertTrue(all(float(value) >= 0.0 for value in collection.values()))
        self.assertLess(maximum_balance_error, 1.0e-10)

    def test_species_mechanism_constraints_pass(self) -> None:
        _, details = dynamic_physics_residuals(self.model, self.targets)
        self.assertEqual(len(details), 3)
        self.assertTrue(all(bool(item["constraint_pass"]) for item in details))

    def test_integration_resolution_is_converged(self) -> None:
        row = self.targets[self.targets["measurement_id"] == "M-014"].iloc[0]
        params = self.model.params
        coarse = DynamicRedoxInventoryDepletionModel(
            self.thermochemical,
            params,
            integration_steps=60,
        ).predict_response(row)
        fine = DynamicRedoxInventoryDepletionModel(
            self.thermochemical,
            params,
            integration_steps=240,
        ).predict_response(row)
        relative = abs(float(coarse) - float(fine)) / max(abs(float(fine)), 1.0e-14)
        self.assertLess(relative, 1.0e-4)

    def test_checked_metrics_preserve_comparative_result(self) -> None:
        metrics = json.loads(
            (ROOT / "results" / "advanced" / "three_model_comparison_metrics.json").read_text(encoding="utf-8")
        )
        self.assertEqual(metrics["n_supported_constraints"], 32)
        self.assertEqual(metrics["n_direct_or_range_targets"], 26)
        self.assertTrue(metrics["optimizer_success"])
        self.assertEqual(metrics["dynamic_network_within_factor_2_direct"], 1.0)
        self.assertLess(
            metrics["dynamic_network_range_aware_log_rmse_direct"],
            metrics["reduced_mechanistic_range_aware_log_rmse_direct"],
        )
        self.assertLess(
            metrics["dynamic_network_range_aware_log_rmse_direct"],
            metrics["thermochemical_range_aware_log_rmse_direct"],
        )
        self.assertLess(metrics["endpoint_mass_balance"]["maximum_relative_error"], 1.0e-10)
        self.assertLess(
            metrics["numerical_convergence"]["comparisons"]["120_vs_240"][
                "maximum_relative_difference"
            ],
            1.0e-4,
        )

    def test_corrected_branch_targets_are_in_force(self) -> None:
        modern = self.targets.set_index("measurement_id")
        for measurement_id in ("M-035", "M-036", "M-037"):
            self.assertEqual(float(modern.loc[measurement_id, "temperature_K"]), 973.15)
        self.assertEqual(str(modern.loc["M-041", "fit_role"]), "direct")
        self.assertEqual(float(modern.loc["M-041", "target_mid"]), 19.7)
        self.assertAlmostEqual(float(modern.loc["M-041", "time_h"]), 1079.0)

    def test_physics_context_is_invariant_to_metadata_renaming(self) -> None:
        original = self.targets[self.targets["measurement_id"] == "M-014"].iloc[0].copy()
        renamed = original.copy()
        renamed["response_kind"] = "renamed_observable"
        renamed["source_id"] = "renamed-source"
        renamed["experiment_family"] = "renamed family"
        first = self.model.build_context(original)
        second = self.model.build_context(renamed)
        for name in (
            "material",
            "inventory_scale",
            "deposition_closure",
            "area_to_salt_mass_cm2_g",
            "inventory_coupling_factor",
            "deposit_area_factor",
            "explicit_inventory_scale",
            "transient_redox",
            "stress_interfacial_activation",
            "fluoride_impurity_interfacial_activation",
            "chloride_salt",
        ):
            self.assertEqual(getattr(first, name), getattr(second, name))
        np.testing.assert_array_equal(first.initial_dissolved_ppm, second.initial_dissolved_ppm)
        first_result = self.model.simulate_context(first)
        second_result = self.model.simulate_context(second)
        relabeled_result = self.model.simulate_context(
            replace(first, redox_class="renamed-redox", salt_class="renamed-salt")
        )
        for name in (
            "front_depth_um",
            "mass_recession_um",
            "mass_loss_mg_cm2",
            "mass_gain_mg_cm2",
            "igc_depth_um",
            "corrosion_rate_um_y",
            "mass_balance_relative_error",
        ):
            self.assertEqual(getattr(first_result, name), getattr(second_result, name))
            self.assertEqual(getattr(first_result, name), getattr(relabeled_result, name))

    def test_dynamic_context_rejects_each_missing_explicit_parity_field(self) -> None:
        row = self.targets[self.targets["measurement_id"] == "M-014"].iloc[0]
        for field in (
            "initial_dissolved_Cr_ppm",
            "initial_dissolved_Fe_ppm",
            "initial_dissolved_Ni_ppm",
            "explicit_inventory_scale",
            "transient_redox",
            "stress_interfacial_activation",
            "fluoride_impurity_interfacial_activation",
            "chloride_salt",
        ):
            incomplete = row.drop(labels=[field])
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, field):
                self.model.build_context(incomplete)

    def test_dynamic_context_rejects_missing_normalized_physics_inputs(self) -> None:
        row = self.targets[self.targets["measurement_id"] == "M-014"].iloc[0]
        for field in (
            "measurement_id",
            "salt_class",
            "redox_class",
            "position_class",
            "temperature_K",
            "delta_T_C",
            "flow_factor",
            "time_years",
        ):
            incomplete = row.drop(labels=[field])
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, field):
                self.model.build_context(incomplete)

    def test_explicit_inventory_scale_mode_matches_the_frozen_loop_scale(self) -> None:
        row = self.targets[self.targets["measurement_id"] == "M-014"].iloc[0]
        loop_context = self.model.build_context(row)
        explicit_context = replace(
            loop_context,
            inventory_scale="explicit",
            explicit_inventory_scale=math.exp(self.model.params["log_inventory_scale_loop"]),
        )
        loop_result = self.model.simulate_context(loop_context)
        explicit_result = self.model.simulate_context(explicit_context)
        self.assertEqual(loop_result.front_depth_um, explicit_result.front_depth_um)
        self.assertEqual(loop_result.dissolved_ppm, explicit_result.dissolved_ppm)

    def test_zero_inventory_coupling_has_no_salt_source_and_conserves_inventory(self) -> None:
        row = self.targets[self.targets["measurement_id"] == "M-014"].iloc[0].copy()
        row["inventory_coupling_factor"] = 0.0
        result = self.model.simulate(row)
        self.assertTrue(all(value == 0.0 for value in result.cumulative_source_ppm.values()))
        self.assertLess(result.mass_balance_relative_error, 1.0e-10)

    def test_zero_exposure_returns_the_initial_state_without_integration(self) -> None:
        row = self.targets[self.targets["measurement_id"] == "M-014"].iloc[0].copy()
        row["time_years"] = 0.0
        context = self.model.build_context(row)
        result = self.model.simulate_context(context, return_trajectory=True)
        self.assertEqual(result.front_depth_um, 0.0)
        self.assertEqual(result.mass_recession_um, 0.0)
        self.assertEqual(result.mass_loss_mg_cm2, 0.0)
        self.assertEqual(result.mass_gain_mg_cm2, 0.0)
        self.assertEqual(result.igc_depth_um, 0.0)
        self.assertEqual(result.corrosion_rate_um_y, 0.0)
        self.assertEqual(result.mass_balance_relative_error, 0.0)
        self.assertEqual(result.trajectory.shape[0], 1)
        np.testing.assert_array_equal(
            np.asarray(list(result.dissolved_ppm.values())),
            context.initial_dissolved_ppm,
        )
        self.assertTrue(all(value == 0.0 for value in result.cumulative_source_ppm.values()))

    def test_dynamic_context_rejects_mismatched_geometry_labels_and_boundaries(self) -> None:
        original = self.targets[self.targets["measurement_id"] == "M-014"].iloc[0]
        mismatched = original.copy()
        mismatched["delta_T_C"] = 0.0
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            self.model.build_context(mismatched)
        typo = original.copy()
        typo["position_class"] = "cold-leg-typo"
        with self.assertRaisesRegex(ValueError, "position_class"):
            self.model.build_context(typo)
        for field, invalid in (
            ("inventory_source_material", "graphite"),
            ("salt_class", "no_salt"),
            ("redox_class", "unknown-redox"),
        ):
            row = original.copy()
            row[field] = invalid
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, field):
                self.model.build_context(row)
        for field, invalid in (
            ("temperature_K", 0.0),
            ("flow_factor", 0.0),
            ("time_years", -1.0),
        ):
            row = original.copy()
            row[field] = invalid
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, field):
                self.model.build_context(row)

    def test_small_positive_flow_builds_a_finite_context(self) -> None:
        row = self.targets[self.targets["measurement_id"] == "M-014"].iloc[0].copy()
        row["flow_factor"] = 1.0e-4
        context = self.model.build_context(row)
        self.assertEqual(context.flow_factor, 1.0e-4)
        self.assertTrue(math.isfinite(context.mass_transfer_rate_um_y))
        self.assertGreater(context.mass_transfer_rate_um_y, 0.0)

    def test_absent_alloy_elements_have_exactly_zero_source_fraction(self) -> None:
        row = self.targets[self.targets["measurement_id"] == "M-014"].iloc[0]
        context = replace(
            self.model.build_context(row),
            mass_fractions=np.asarray([1.0, 0.0, 0.0]),
        )
        fractions, _ = self.model._species_fractions(
            context,
            context.initial_dissolved_ppm,
            np.ones(3),
            context.redox_shift_initial,
        )
        np.testing.assert_array_equal(fractions, np.asarray([1.0, 0.0, 0.0]))
        with self.assertRaisesRegex(ValueError, "at least one"):
            self.model._species_fractions(
                replace(context, mass_fractions=np.zeros(3)),
                context.initial_dissolved_ppm,
                np.ones(3),
                context.redox_shift_initial,
            )


if __name__ == "__main__":
    unittest.main()
