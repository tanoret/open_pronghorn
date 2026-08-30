from __future__ import annotations

import math
import hashlib
import os
import tempfile
import unittest
from pathlib import Path

from msr_corrosion_bv.mstdb import (
    ChemSageDatabase,
    GibbsInterval,
    MSTDBPair,
    R_GAS,
    magnetic_gibbs_J_mol,
)


class GibbsPolynomialTests(unittest.TestCase):
    def test_polynomial_and_extra_terms(self) -> None:
        interval = GibbsInterval(
            upper_temperature_K=2000.0,
            coefficients=(1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
            additional_terms=((7.0, 99.0), (8.0, 0.5)),
        )
        temperature = 900.0
        expected = (
            1.0
            + 2.0 * temperature
            + 3.0 * temperature * math.log(temperature)
            + 4.0 * temperature**2
            + 5.0 * temperature**3
            + 6.0 / temperature
            + 7.0 * math.log(temperature)
            + 8.0 * temperature**0.5
        )
        self.assertAlmostEqual(interval.evaluate(temperature), expected, places=7)

    def test_interval_selection_and_parser(self) -> None:
        text = """System X\n1\nX_S1(s)\n4 2 1\n500 1 2 3 4 5 6\n1 7 99\n2000 10 20 30 40 50 60\n0\n"""
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "MSTDB-TC_V0.0_Test_No_Func.dat"
            path.write_text(text, encoding="utf-8")
            database = ChemSageDatabase(path)
            self.assertEqual(database.n_elements, 1)
            self.assertEqual(len(database.records), 1)
            self.assertEqual(database.metadata()["version"], "0.0")
            self.assertNotIn("path", database.metadata())
            lower = database.standard_gibbs_J_mol("X_S1(s)", 400.0)
            upper = database.standard_gibbs_J_mol("X_S1(s)", 800.0)
            expected_lower = 1 + 2 * 400 + 3 * 400 * math.log(400) + 4 * 400**2 + 5 * 400**3 + 6 / 400 + 7 * math.log(400)
            expected_upper = 10 + 20 * 800 + 30 * 800 * math.log(800) + 40 * 800**2 + 50 * 800**3 + 60 / 800
            self.assertAlmostEqual(lower, expected_lower, places=6)
            self.assertAlmostEqual(upper, expected_upper, places=6)
            at_boundary = database.standard_gibbs_J_mol("X_S1(s)", 2000.0)
            self.assertTrue(math.isfinite(at_boundary))
            with self.assertRaisesRegex(ValueError, "allow_extrapolation"):
                database.standard_gibbs_J_mol("X_S1(s)", 2000.1)
            self.assertTrue(
                math.isfinite(
                    database.standard_gibbs_J_mol(
                        "X_S1(s)", 2000.1, allow_extrapolation=True
                    )
                )
            )

            compatibility_database = ChemSageDatabase(path, allow_extrapolation=True)
            self.assertTrue(
                math.isfinite(compatibility_database.standard_gibbs_J_mol("X_S1(s)", 2000.1))
            )
            self.assertTrue(compatibility_database.metadata()["allow_extrapolation"])

    def test_pair_selection_requires_a_common_version(self) -> None:
        text = "System X\n1\nX_S1(s)\n1 1 1\n2000 0 0 0 0 0 0\n"
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "MSTDB-TC_V3.1_Fluorides_No_Func.dat").write_text(text, encoding="utf-8")
            (root / "MSTDB-TC_V3.1_Chlorides_No_Func.dat").write_text(text, encoding="utf-8")
            (root / "MSTDB-TC_V4.0_Fluorides_No_Func.dat").write_text(text, encoding="utf-8")
            pair = MSTDBPair.from_directory(root)
            self.assertEqual(pair.fluoride.metadata()["version"], "3.1")
            self.assertEqual(pair.chloride.metadata()["version"], "3.1")

    def test_pair_selection_rejects_mixed_versions(self) -> None:
        text = "System X\n1\nX_S1(s)\n1 1 1\n2000 0 0 0 0 0 0\n"
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "MSTDB-TC_V4.0_Fluorides_No_Func.dat").write_text(text, encoding="utf-8")
            (root / "MSTDB-TC_V3.1_Chlorides_No_Func.dat").write_text(text, encoding="utf-8")
            with self.assertRaises(ValueError):
                MSTDBPair.from_directory(root)

    def test_provenance_bound_loader_rejects_wrong_version_hash_and_duplicates(self) -> None:
        text = "System X\n1\nX_S1(s)\n1 1 1\n2000 0 0 0 0 0 0\n"
        expected_hash = hashlib.sha256(text.encode("utf-8")).hexdigest()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            fluoride = root / "MSTDB-TC_V3.1_Fluorides_No_Func.dat"
            chloride = root / "MSTDB-TC_V3.1_Chlorides_No_Func.dat"
            fluoride.write_text(text, encoding="utf-8")
            chloride.write_text(text, encoding="utf-8")
            pair = MSTDBPair.from_directory(
                root,
                expected_version="3.1",
                expected_fluoride_sha256=expected_hash,
                expected_chloride_sha256=expected_hash,
            )
            self.assertEqual(pair.fluoride.sha256, expected_hash)
            with self.assertRaisesRegex(ValueError, "Required MSTDB-TC V4.0"):
                MSTDBPair.from_directory(root, expected_version="4.0")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                MSTDBPair.from_directory(
                    root,
                    expected_version="3.1",
                    expected_fluoride_sha256="0" * 64,
                )
            duplicate_dir = root / "duplicate"
            duplicate_dir.mkdir()
            (duplicate_dir / fluoride.name).write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Ambiguous duplicate"):
                MSTDBPair.from_directory(root, expected_version="3.1")

    def test_sgte_magnetic_term_matches_independent_expression(self) -> None:
        temperature = 923.15
        critical = 1043.0
        moment = 2.22
        structure_factor = 1.0
        p = 0.4
        tau = temperature / critical
        invpmone = 1.0 / p - 1.0
        denominator = 518.0 / 1125.0 + (11692.0 / 15975.0) * invpmone
        g_tau = 1.0 - (
            79.0 / (140.0 * p * tau)
            + (474.0 / 497.0) * invpmone * (tau**3 / 6.0 + tau**9 / 135.0 + tau**15 / 600.0)
        ) / denominator
        expected = R_GAS * temperature * math.log(moment + 1.0) * g_tau
        self.assertAlmostEqual(
            magnetic_gibbs_J_mol(temperature, critical, moment, structure_factor, p),
            expected,
            places=10,
        )


@unittest.skipUnless(os.environ.get("MSTDB_TC_DIR"), "set MSTDB_TC_DIR for database integration tests")
class MSTDBIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.mstdb = MSTDBPair.from_calibration_directory(os.environ["MSTDB_TC_DIR"])

    def test_attached_v31_regression_values(self) -> None:
        fluoride = self.mstdb.fluoride
        chloride = self.mstdb.chloride
        self.assertEqual(len(fluoride.records), 677)
        self.assertEqual(len(chloride.records), 649)
        self.assertAlmostEqual(
            fluoride.standard_gibbs_J_mol("CrF2_L1(liq)", 923.15),
            -881852.8502324164,
            places=5,
        )
        self.assertAlmostEqual(
            chloride.standard_gibbs_J_mol("CrCl2_L1(liq)", 923.15),
            -528859.0468304644,
            places=5,
        )

    def test_reaction_direction_sanity(self) -> None:
        database = self.mstdb.fluoride
        temperature = 923.15
        cr_u = database.equilibrium_log_constant(
            {"CrF2_L1(liq)": 1, "UF3_L1(liq)": 2, "Cr_S1(s)": -1, "UF4_L1(liq)": -2},
            temperature,
        )
        fe_u = database.equilibrium_log_constant(
            {"FeF2_L1(liq)": 1, "UF3_L1(liq)": 2, "Fe_bcc(s)": -1, "UF4_L1(liq)": -2},
            temperature,
        )
        ni_u = database.equilibrium_log_constant(
            {"NiF2_L1(liq)": 1, "UF3_L1(liq)": 2, "Ni_fcc(s)": -1, "UF4_L1(liq)": -2},
            temperature,
        )
        self.assertGreater(cr_u, fe_u)
        self.assertGreater(fe_u, ni_u)
        self.assertGreater(
            database.equilibrium_log_constant(
                {"CrF2_L1(liq)": 1, "Fe_bcc(s)": 1, "Cr_S1(s)": -1, "FeF2_L1(liq)": -1},
                temperature,
            ),
            0.0,
        )
        self.assertLess(
            database.equilibrium_log_constant(
                {"CrF2_L1(liq)": 1, "Be_S1(s)": 1, "Cr_S1(s)": -1, "BeF2_L1(liq)": -1},
                temperature,
            ),
            0.0,
        )


if __name__ == "__main__":
    unittest.main()
