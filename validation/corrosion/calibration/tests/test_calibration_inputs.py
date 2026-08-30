from __future__ import annotations

import json
import math
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class CalibrationInputProvenanceTests(unittest.TestCase):
    def test_current_pr_corrosion_target_corrections_are_present(self) -> None:
        targets = pd.read_csv(ROOT / "data" / "targets.csv").set_index("measurement_id")
        for measurement_id in ("M-035", "M-036", "M-037"):
            self.assertEqual(float(targets.loc[measurement_id, "temperature_K"]), 973.15)
        self.assertEqual(str(targets.loc["M-041", "fit_role"]), "direct")
        self.assertEqual(float(targets.loc["M-041", "target_mid"]), 19.7)
        self.assertAlmostEqual(float(targets.loc["M-041", "time_h"]), 1079.0)

    def test_shared_effective_parameters_match_current_branch(self) -> None:
        parameters = json.loads(
            (ROOT / "data" / "parameters.json").read_text(encoding="utf-8")
        )
        self.assertEqual(len(parameters), 61)
        self.assertAlmostEqual(parameters["Ea_corr_kJ_mol"], 54.83553514725915)

    def test_advanced_context_is_complete_and_explicit(self) -> None:
        targets = pd.read_csv(ROOT / "data" / "targets.csv")
        context = pd.read_csv(ROOT / "results" / "advanced" / "advanced_case_context.csv")
        required = {
            "initial_dissolved_Cr_ppm",
            "initial_dissolved_Fe_ppm",
            "initial_dissolved_Ni_ppm",
            "explicit_inventory_scale",
            "transient_redox",
            "stress_interfacial_activation",
            "fluoride_impurity_interfacial_activation",
            "chloride_salt",
        }
        self.assertTrue(required.issubset(context.columns))
        self.assertEqual(set(targets["measurement_id"]), set(context["measurement_id"]))
        self.assertFalse(context.isna().any().any())
        self.assertEqual(set(context["context_revision"]), {"explicit-geometry-v2"})
        self.assertEqual(set(context["inventory_scale"]), {"msre", "loop"})
        self.assertEqual(set(context["deposition_closure"]), {"fuel", "flinak"})
        self.assertTrue(context["area_to_salt_mass_cm2_g"].gt(0.0).all())
        self.assertTrue(context["inventory_coupling_factor"].between(0.0, 1.0, inclusive="right").all())
        self.assertTrue(context["explicit_inventory_scale"].eq(1.0).all())

        thermochemical = json.loads(
            (ROOT / "results" / "advanced" / "thermochemical_parameters.json").read_text(
                encoding="utf-8"
            )
        )
        former_floor = math.exp(float(thermochemical["log_product_floor_ppm"]))
        for element in ("Cr", "Fe", "Ni"):
            self.assertTrue(
                context[f"initial_dissolved_{element}_ppm"].eq(former_floor).all(),
                element,
            )

        joined = context.merge(
            targets[["measurement_id", "redox_class", "salt_class"]],
            on="measurement_id",
            validate="one_to_one",
        )
        self.assertTrue(
            joined["transient_redox"].eq(
                joined["redox_class"].isin(
                    {"oxidizing_fef2", "reducing_be", "impure_moisture"}
                )
            ).all()
        )
        self.assertTrue(
            joined["stress_interfacial_activation"].eq(joined["redox_class"].eq("stressed")).all()
        )
        self.assertTrue(
            joined["fluoride_impurity_interfacial_activation"].eq(
                joined["redox_class"].eq("impure_moisture")
            ).all()
        )
        self.assertTrue(joined["chloride_salt"].eq(joined["salt_class"].eq("chloride")).all())

    def test_provenance_marks_conditional_sensitivity_as_nonindependent(self) -> None:
        metrics = json.loads(
            (ROOT / "results" / "advanced" / "thermochemical_comparison_metrics.json").read_text(
                encoding="utf-8"
            )
        )
        sensitivity = metrics["modern_family_conditional_sensitivity"]
        self.assertFalse(sensitivity["independent_holdout"])
        self.assertFalse(sensitivity["shared_effective_model_refit_on_training_only"])


if __name__ == "__main__":
    unittest.main()
