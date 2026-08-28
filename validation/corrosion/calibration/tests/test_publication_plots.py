from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pandas as pd

from msr_corrosion_bv.publication_plots import (
    apply_publication_style,
    figure_all_constraint_margins,
    figure_validation_ratios,
)


ROOT = Path(__file__).resolve().parents[1]


class PublicationPlotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.comparison = pd.read_csv(ROOT / "results" / "advanced" / "advanced_validation_summary.csv")

    def test_validation_plots_export_png_svg_and_pdf(self) -> None:
        apply_publication_style()
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir)
            generated = []
            generated.extend(figure_validation_ratios(self.comparison, output))
            generated.extend(figure_all_constraint_margins(self.comparison, output))
            self.assertEqual(len(generated), 6)
            for path in generated:
                self.assertTrue(path.exists(), path)
                self.assertGreater(path.stat().st_size, 1000, path)
            self.assertEqual({path.suffix for path in generated}, {".png", ".svg", ".pdf"})


if __name__ == "__main__":
    unittest.main()
