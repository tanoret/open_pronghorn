from TestHarness.validation import ValidationCase

import csv
import numpy as np


class TestCase(ValidationCase):
    """Validate the modeled radiolytic F2 yield against the Davis et al. NSE 2022 Table III
    measurement for FLiBe-UF4. The effective G value is recovered from the linear F2 accumulation,
    G = (dC_F2/dt) / (dose_rate / (100 eV * N_A)), and compared to the measured value read from
    data/G_values_table_III.csv."""

    # Conversion constants (match MSR::gToSource)
    EV_TO_J = 1.602176634e-19
    N_A = 6.02214076e23

    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addParam("dose_rate", 1.0e6, "Volumetric dose rate used in the input [J/m^3/s]")
        params.addParam("validation_upper_bound", 0.20,
                        "Maximum acceptable relative error on the recovered F2 G value")
        return params

    def initialize(self):
        # Modeled F2 accumulation -> production rate -> effective G value.
        t, c = [], []
        with open("f2_yield_out.csv") as f:
            for row in csv.DictReader(f):
                t.append(float(row["time"]))
                c.append(float(row["c_F2_diss"]))
        t = np.array(t)
        c = np.array(c)
        dCdt = float(np.polyfit(t, c, 1)[0])
        dose_rate = float(self.getParam("dose_rate"))
        factor = 1.0 / (100.0 * self.EV_TO_J) / self.N_A
        self.g_model = dCdt / (dose_rate * factor)

        # Measured FLiBe-UF4 yield from Davis Table III.
        self.g_exp = None
        with open("data/G_values_table_III.csv") as f:
            for row in csv.reader(f):
                if not row or row[0].startswith("#") or row[0] == "salt":
                    continue
                if row[0] == "FLiBe_UF4":
                    self.g_exp = float(row[2])
        if self.g_exp is None:
            raise RuntimeError("Could not find the FLiBe_UF4 G value in Table III")

        self.rel_err = abs(self.g_model - self.g_exp) / self.g_exp

    def testValidation(self):
        upper = float(self.getParam("validation_upper_bound"))
        self.addScalarData(
            key="g_f2_model",
            value=self.g_model,
            description="Effective radiolytic F2 G value recovered from the modeled production rate",
            units="molecules/100eV",
        )
        self.addScalarData(
            key="g_f2_relative_error",
            value=self.rel_err,
            description="Relative error of the modeled F2 G value vs the Davis 2022 FLiBe-UF4 measurement",
            units="-",
            bounds=(0.0, upper),
        )
