from __future__ import annotations

import math
import unittest

from msr_corrosion_bv.thermochemical_fit import _residual_weighting, relation_aware_log_error
from msr_corrosion_bv.dynamic_network_fit import range_aware_log_error


class RelationAwareResidualTests(unittest.TestCase):
    def test_range_is_zero_inside_and_nearest_bound_outside(self) -> None:
        row = {
            "fit_role": "range",
            "target_low": 10.0,
            "target_mid": 20.0,
            "target_high": 40.0,
        }
        self.assertEqual(relation_aware_log_error(row, 25.0), 0.0)
        self.assertAlmostEqual(relation_aware_log_error(row, 5.0), math.log(0.5))
        self.assertAlmostEqual(relation_aware_log_error(row, 80.0), math.log(2.0))

    def test_one_sided_bounds_only_penalize_violations(self) -> None:
        upper = {"fit_role": "upper", "target_mid": 10.0, "target_high": 10.0}
        lower = {"fit_role": "lower", "target_low": 10.0, "target_mid": 10.0}
        self.assertEqual(relation_aware_log_error(upper, 5.0), 0.0)
        self.assertAlmostEqual(relation_aware_log_error(upper, 20.0), math.log(2.0))
        self.assertEqual(relation_aware_log_error(lower, 20.0), 0.0)
        self.assertAlmostEqual(relation_aware_log_error(lower, 5.0), math.log(0.5))
        for evaluator in (relation_aware_log_error, range_aware_log_error):
            self.assertEqual(evaluator(upper, 5.0), 0.0)
            self.assertAlmostEqual(evaluator(upper, 20.0), math.log(2.0))
            self.assertEqual(evaluator(lower, 20.0), 0.0)
            self.assertAlmostEqual(evaluator(lower, 5.0), math.log(0.5))

    def test_residual_weighting_distinguishes_missing_from_explicit_zero(self) -> None:
        self.assertEqual(_residual_weighting({}), (0.75, 1.0))
        self.assertEqual(_residual_weighting({"default_sigma_ln": 0.1}), (0.15, 1.0))
        for field in ("default_sigma_ln", "quality_weight"):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, field):
                _residual_weighting({field: 0.0})


if __name__ == "__main__":
    unittest.main()
