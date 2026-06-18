from TestHarness.validation import ValidationCase

import csv
import numpy as np


class TestCase(ValidationCase):
    """Validate the modeled solvated-electron decay against an Iwamatsu et al. PCCP 2026
    transient-absorption trace for a given Cr(II) or Cr(III) concentration in molten LiCl-KCl. The
    measured absorbance is proportional to [e_sol] up to an unknown scale factor, so the model trace
    is rescaled so that its maximum matches the measured maximum (the scale-free convention used by
    the standalone MoltenSaltRadiolysis harness) before computing the residual.

    The measurement file and the model output file are passed in via parameters so the same script
    serves every concentration trace."""

    @staticmethod
    def _read_csv(path, x_key, y_key):
        x, y = [], []
        with open(path) as f:
            for row in csv.DictReader(f):
                x.append(float(row[x_key]))
                y.append(float(row[y_key]))
        return np.array(x), np.array(y)

    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addRequiredParam("validation_data_file",
                                "Digitized absorbance trace (time [ns], absorbance)")
        params.addRequiredParam("validation_model_file",
                                "Model output CSV (time [s], c_e_sol)")
        params.addParam("validation_upper_bound", 0.30,
                        "Maximum acceptable scale-free relative RMSE")
        return params

    def initialize(self):
        model_t, model_e = self._read_csv(self.getParam("validation_model_file"), "time", "c_e_sol")
        exp_t, exp_a = self._read_csv(self.getParam("validation_data_file"), "time", "absorbance")
        exp_t = exp_t * 1.0e-9  # ns -> s

        model_on_exp = np.interp(exp_t, model_t, model_e)
        model_scaled = model_on_exp * (np.nanmax(exp_a) / model_on_exp.max())

        self.exp_t_ns = exp_t * 1.0e9
        self.exp_a = exp_a
        self.model_scaled = model_scaled
        self.rel_rmse = float(
            np.sqrt(np.mean((model_scaled - exp_a) ** 2)) / np.nanmax(exp_a)
        )

    def testValidation(self):
        upper = float(self.getParam("validation_upper_bound"))
        self.addVectorData(
            "absorbance",
            (self.exp_t_ns, "Time", "ns"),
            (self.exp_a, "Measured scale-free absorbance", "-"),
        )
        self.addVectorData(
            "model_scaled",
            (self.exp_t_ns, "Time", "ns"),
            (self.model_scaled, "Model e_sol rescaled to measurement", "-"),
        )
        self.addScalarData(
            key="scale_free_rel_rmse",
            value=self.rel_rmse,
            description="Scale-free relative RMSE of modeled e_sol decay vs measured absorbance",
            units="-",
            bounds=(0.0, upper),
        )
