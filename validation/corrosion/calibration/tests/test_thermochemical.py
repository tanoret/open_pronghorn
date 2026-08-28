from __future__ import annotations

import json
import math
import os
import stat
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from msr_corrosion_bv.calibrate import load_model_from_json
from msr_corrosion_bv.mstdb import MSTDBPair, ThermochimicaRunner
from msr_corrosion_bv.thermochemical import MSTDBThermochemicalCorrosionModel


ROOT = Path(__file__).resolve().parents[1]


class ThermochemicalIdentityRegressionTests(unittest.TestCase):
    @staticmethod
    def model_without_database() -> MSTDBThermochemicalCorrosionModel:
        return MSTDBThermochemicalCorrosionModel(
            object(),
            SimpleNamespace(fluoride=object(), chloride=object()),
        )

    def test_nonfuel_fe_buffer_identity_does_not_require_database_lookup(self) -> None:
        model = self.model_without_database()
        for salt_class in ("flinak", "flibe", "fluoroborate", "chloride"):
            affinity = model.reaction_log_K_over_Q(
                "Fe",
                {
                    "material_class": "stainless_316h",
                    "salt_class": salt_class,
                    "redox_class": "purified_baseline",
                },
                923.15,
            )
            self.assertEqual(affinity, 0.0)

    def test_inventory_geometry_is_mandatory_without_database_access(self) -> None:
        model = self.model_without_database()
        with self.assertRaisesRegex(ValueError, "inventory_source_material"):
            model._inventory_source_material({})
        with self.assertRaisesRegex(ValueError, "area_to_salt_mass_cm2_g"):
            model._area_to_salt_mass({})


class ThermochimicaAdapterTests(unittest.TestCase):
    def test_missing_executable_is_explicit(self) -> None:
        with self.assertRaises(FileNotFoundError):
            ThermochimicaRunner("definitely-not-a-real-thermochimica-executable")

    def test_input_script_adapter_reads_json_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir_text:
            temp_dir = Path(temp_dir_text)
            executable = temp_dir / "InputScriptMode"
            executable.write_text(
                "#!/bin/sh\n"
                "set -eu\n"
                "test -f \"$1\"\n"
                "grep -q 'temperature       = 923.15' \"$1\"\n"
                "grep -q 'mass(24)' \"$1\"\n"
                "mkdir -p outputs\n"
                "printf '{\"1\": {\"status\": \"ok\"}}\\n' > outputs/thermoout.json\n",
                encoding="utf-8",
            )
            executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
            database = temp_dir / "database.dat"
            database.write_text("placeholder", encoding="utf-8")
            result = ThermochimicaRunner(executable).run_elements(
                database,
                923.15,
                {24: 1.0, 26: 2.0},
            )
            self.assertEqual(result["1"]["status"], "ok")


@unittest.skipUnless(os.environ.get("MSTDB_TC_DIR"), "set MSTDB_TC_DIR for thermochemical integration tests")
class ThermochemicalModelIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        base = load_model_from_json(ROOT / "data" / "parameters.json")
        mstdb = MSTDBPair.from_calibration_directory(os.environ["MSTDB_TC_DIR"])
        params = json.loads((ROOT / "results" / "advanced" / "thermochemical_parameters.json").read_text(encoding="utf-8"))
        cls.model = MSTDBThermochemicalCorrosionModel(base, mstdb, params)

    @staticmethod
    def base_row(**updates: object) -> dict[str, object]:
        row: dict[str, object] = {
            "material_class": "hastelloy_n",
            "salt_class": "fluoride_fuel",
            "redox_class": "purified_baseline",
            "temperature_K": 923.15,
            "flow_factor": 1.0,
            "delta_T_C": 100.0,
            "temperature_cold_C": 550.0,
            "position_class": "hot_leg",
            "time_years": 1.0,
            "experiment_family": "ORNL fluoride loops",
            "source_id": "unit-test",
            "inventory_source_material": "hastelloy_n",
            "area_to_salt_mass_cm2_g": 1.0,
            "inventory_coupling_factor": 1.0,
        }
        row.update(updates)
        return row

    def test_default_parameter_constructor_is_operational(self) -> None:
        base = load_model_from_json(ROOT / "data" / "parameters.json")
        mstdb = MSTDBPair.from_calibration_directory(os.environ["MSTDB_TC_DIR"])
        model = MSTDBThermochemicalCorrosionModel(base, mstdb)
        rate = model.corrosion_rate_um_y(self.base_row())
        self.assertTrue(math.isfinite(rate))
        self.assertGreater(rate, 0.0)

    def test_species_fractions_are_normalized_and_cr_selective(self) -> None:
        fractions = self.model.species_flux_fractions(self.base_row())
        self.assertAlmostEqual(sum(fractions.values()), 1.0, places=12)
        self.assertGreater(fractions["Cr"], fractions["Fe"])
        self.assertGreater(fractions["Cr"], fractions["Ni"])

    def test_explicit_redox_buffers_have_expected_order(self) -> None:
        row = self.base_row()
        baseline = self.model.corrosion_rate_um_y(row)
        oxidizing = self.model.corrosion_rate_um_y(row, redox_override="oxidizing_fef2")
        reducing = self.model.corrosion_rate_um_y(
            self.base_row(salt_class="flibe", redox_class="reducing_be")
        )
        flibe_baseline = self.model.corrosion_rate_um_y(self.base_row(salt_class="flibe"))
        self.assertGreater(oxidizing, baseline)
        self.assertLess(reducing, flibe_baseline)

    def test_small_positive_flow_uses_the_continuous_mass_transfer_law(self) -> None:
        rate = self.model.dissolution_front_rate_um_y(self.base_row(flow_factor=1.0e-4))
        self.assertTrue(math.isfinite(rate))
        self.assertGreater(rate, 0.0)
        self.assertNotEqual(rate, 1.0e-9)
        for invalid in (0.0, -1.0, math.inf, math.nan):
            with self.assertRaisesRegex(ValueError, "required|finite and positive"):
                self.model.dissolution_front_rate_um_y(self.base_row(flow_factor=invalid))

    def test_nonfuel_fe_buffer_identity_has_zero_affinity(self) -> None:
        for salt_class in ("flinak", "flibe", "fluoroborate", "chloride"):
            affinity = self.model.reaction_log_K_over_Q(
                "Fe",
                self.base_row(salt_class=salt_class),
                923.15,
            )
            self.assertEqual(affinity, 0.0)

    def test_cold_leg_saturation_drives_capture(self) -> None:
        hot = self.model.saturation_activity("Cr", "fluoride_fuel", 977.15)
        cold = self.model.saturation_activity("Cr", "fluoride_fuel", 811.15)
        self.assertLess(cold, hot)
        self.assertGreater(self.model.cold_capture_fraction("Cr", "fluoride_fuel", 977.15, 811.15), 0.0)

    def test_flinak_cold_leg_deposit_is_fe_rich(self) -> None:
        row = self.base_row(
            material_class="stainless_316h",
            salt_class="flinak",
            temperature_K=813.15,
            temperature_cold_C=540.0,
            delta_T_C=110.0,
            position_class="cold_leg",
            time_years=1000.0 / (365.25 * 24.0),
        )
        fractions = self.model.deposition_species_fractions(row)
        self.assertGreater(fractions["Fe"], fractions["Cr"])
        self.assertGreater(fractions["Fe"], fractions["Ni"])

    def test_species_inventory_is_finite_and_nonnegative(self) -> None:
        inventory = self.model.species_inventory_ppm(self.base_row(time_years=3.0))
        self.assertEqual(set(inventory), {"Cr", "Fe", "Ni"})
        for value in inventory.values():
            self.assertTrue(math.isfinite(value))
            self.assertGreaterEqual(value, 0.0)

    def test_zero_exposure_has_zero_extents_and_finite_rates(self) -> None:
        row = self.base_row(time_years=0.0)
        self.assertTrue(math.isfinite(self.model.dissolution_front_rate_um_y(row)))
        self.assertTrue(math.isfinite(self.model.corrosion_rate_um_y(row)))
        self.assertTrue(math.isfinite(self.model.cr_diffusion_cm2_s(row)))
        self.assertTrue(math.isfinite(self.model.redox_acceleration_ratio(row)))
        self.assertEqual(self.model.corrosion_depth_um(row), 0.0)
        self.assertEqual(self.model.mass_loss_mg_cm2(row), 0.0)
        self.assertEqual(self.model.mass_gain_mg_cm2(row), 0.0)
        self.assertEqual(self.model.igc_depth_um(row), 0.0)
        self.assertEqual(self.model.species_inventory_ppm(row), {"Cr": 0.0, "Fe": 0.0, "Ni": 0.0})

    def test_cold_leg_inventory_and_diagnostics_use_the_normalized_hot_temperature(self) -> None:
        common = {
            "material_class": "stainless_316h",
            "inventory_source_material": "stainless_316h",
            "salt_class": "flinak",
            "redox_class": "purified_baseline",
            "temperature_cold_C": 540.0,
            "delta_T_C": 110.0,
            "time_years": 1000.0 / (365.25 * 24.0),
        }
        cold_leg = self.base_row(
            **common, temperature_K=813.15, position_class="cold_leg"
        )
        explicit_pair = self.base_row(
            **common, temperature_K=923.15, position_class="hot_leg"
        )
        cold_inventory = self.model.species_inventory_ppm(cold_leg)
        pair_inventory = self.model.species_inventory_ppm(explicit_pair)
        for element in ("Cr", "Fe", "Ni"):
            self.assertAlmostEqual(cold_inventory[element], pair_inventory[element], places=13)
        cold_diagnostic = self.model.diagnostic_row(cold_leg)
        pair_diagnostic = self.model.diagnostic_row(explicit_pair)
        self.assertEqual(cold_diagnostic["temperature_K"], 923.15)
        for key in cold_diagnostic:
            if key.startswith(("ln_K_over_Q_", "dissolution_fraction_", "saturation_activity_")):
                self.assertAlmostEqual(cold_diagnostic[key], pair_diagnostic[key], places=13)

    def test_cold_leg_predict_response_uses_hot_governing_temperature_for_every_kind(self) -> None:
        common = {
            "material_class": "stainless_316h",
            "inventory_source_material": "stainless_316h",
            "salt_class": "flinak",
            "redox_class": "purified_baseline",
            "temperature_cold_C": 540.0,
            "delta_T_C": 110.0,
            "time_years": 1000.0 / (365.25 * 24.0),
        }
        cold_leg = self.base_row(
            **common, temperature_K=813.15, position_class="cold_leg"
        )
        explicit_pair = self.base_row(
            **common, temperature_K=923.15, position_class="hot_leg"
        )
        response_kinds = (
            "corrosion_rate_um_y",
            "corrosion_depth_um",
            "mass_loss_mg_cm2",
            "mass_gain_mg_cm2",
            "salt_cr_ppm",
            "salt_fe_decrease_ppm",
            "cr_diffusion_cm2_s",
            "redox_acceleration_ratio",
            "redox_acceleration_qualitative",
            "igc_depth_um",
        )
        for response_kind in response_kinds:
            with self.subTest(response_kind=response_kind):
                cold_prediction = self.model.predict_response(
                    {**cold_leg, "response_kind": response_kind}
                )
                explicit_prediction = self.model.predict_response(
                    {**explicit_pair, "response_kind": response_kind}
                )
                self.assertTrue(math.isfinite(cold_prediction))
                self.assertTrue(
                    math.isclose(
                        cold_prediction,
                        explicit_prediction,
                        rel_tol=2.0e-13,
                        abs_tol=1.0e-30,
                    ),
                    msg=(
                        f"{response_kind}: cold-leg={cold_prediction!r}, "
                        f"explicit-pair={explicit_prediction!r}"
                    ),
                )

    def test_advanced_geometry_labels_and_physical_boundaries_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "position_class"):
            self.model.mass_gain_mg_cm2(self.base_row(position_class="cold-leg-typo"))
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            self.model.mass_gain_mg_cm2(
                self.base_row(delta_T_C=0.0, temperature_cold_C=550.0)
            )
        for field, invalid in (
            ("material_class", "graphite"),
            ("material_class", "unknown-alloy"),
            ("salt_class", "no_salt"),
            ("salt_class", "unknown-salt"),
            ("redox_class", "unknown-redox"),
        ):
            with self.subTest(field=field, invalid=invalid), self.assertRaisesRegex(
                ValueError, field
            ):
                self.model.dissolution_front_rate_um_y(self.base_row(**{field: invalid}))
        for temperature in (0.0, -1.0, math.inf, math.nan):
            with self.subTest(temperature=temperature), self.assertRaisesRegex(
                ValueError, "temperature_K"
            ):
                self.model.dissolution_front_rate_um_y(
                    self.base_row(temperature_K=temperature)
                )
        for time_years in (-1.0, math.inf, math.nan):
            with self.subTest(time_years=time_years), self.assertRaisesRegex(
                ValueError, "time_years"
            ):
                self.model.corrosion_depth_um(self.base_row(time_years=time_years))

    def test_inventory_physics_rejects_missing_explicit_context(self) -> None:
        with self.assertRaisesRegex(ValueError, "inventory_source_material"):
            self.model.species_inventory_ppm(
                {key: value for key, value in self.base_row().items() if key != "inventory_source_material"}
            )
        with self.assertRaisesRegex(ValueError, "area_to_salt_mass_cm2_g"):
            self.model.species_inventory_ppm(
                {key: value for key, value in self.base_row().items() if key != "area_to_salt_mass_cm2_g"}
            )
        with self.assertRaisesRegex(ValueError, "inventory_coupling_factor"):
            self.model.species_inventory_ppm(
                {
                    key: value
                    for key, value in self.base_row().items()
                    if key != "inventory_coupling_factor"
                }
            )

    def test_inventory_coupling_factor_scales_the_pre_capacity_source(self) -> None:
        full = self.model.species_inventory_ppm(
            self.base_row(time_years=1.0e-4, inventory_coupling_factor=1.0)
        )
        half = self.model.species_inventory_ppm(
            self.base_row(time_years=1.0e-4, inventory_coupling_factor=0.5)
        )
        capacity = math.exp(float(self.model.params["log_inventory_capacity_ppm"]))
        full_total = sum(full.values())
        half_total = sum(half.values())
        raw_full = -capacity * math.log1p(-full_total / capacity)
        expected_half = capacity * (-math.expm1(-0.5 * raw_full / capacity))
        self.assertAlmostEqual(half_total, expected_half, places=12)
        self.assertLess(half_total, full_total)
        with self.assertRaisesRegex(ValueError, r"\[0, 1\]"):
            self.model.species_inventory_ppm(
                self.base_row(inventory_coupling_factor=1.01)
            )

    def test_mstdb_physics_is_invariant_to_reporting_metadata_renaming(self) -> None:
        original = self.base_row(response_kind="salt_cr_ppm")
        renamed = dict(original)
        renamed.update(
            {
                "response_kind": "renamed-observable",
                "source_id": "renamed-source",
                "experiment_family": "renamed-family",
            }
        )
        self.assertEqual(
            self.model.corrosion_rate_um_y(original),
            self.model.corrosion_rate_um_y(renamed),
        )
        self.assertEqual(
            self.model.species_inventory_ppm(original),
            self.model.species_inventory_ppm(renamed),
        )


if __name__ == "__main__":
    unittest.main()
