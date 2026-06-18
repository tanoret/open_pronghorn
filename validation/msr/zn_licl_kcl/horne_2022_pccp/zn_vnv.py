from TestHarness.validation import ValidationCase

import csv
import numpy as np


class TestCase(ValidationCase):
    """Validate the modeled pseudo-first-order rate of e_sol capture by Zn(II) against the digitized
    Fig. 4B observation of Iwamatsu/Horne et al. PCCP 2022. With [Zn_II] >> [e_sol] the solvated
    electron decays as e(t) = e0 exp(-k_obs t); k_obs is extracted from the modeled decay and
    compared to the measured value at 400 C and 9.41 mM Zn(II)."""

    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addParam("validation_upper_bound", 0.15,
                        "Maximum acceptable relative error on the pseudo-first-order rate constant")
        return params

    def initialize(self):
        # Modeled e_sol(t)
        t, e = [], []
        with open("zn_pulse_out.csv") as f:
            for row in csv.DictReader(f):
                t.append(float(row["time"]))
                e.append(float(row["c_e_sol"]))
        t = np.array(t)
        e = np.array(e)
        # Pseudo-first-order fit of ln(e_sol) vs t while e_sol is well above zero.
        mask = (t > 0) & (e > 1e-2 * e[0])
        self.k_obs_model = float(-np.polyfit(t[mask], np.log(e[mask]), 1)[0])

        # Digitized measurement: 400 C, 9.41 mM, k_obs reported in units of 1e8 /s.
        self.k_obs_exp = None
        with open("data/vision_fig4B_Zn2_pseudo1st_order.csv") as f:
            for row in csv.reader(f):
                if not row or row[0].startswith("#") or row[0] == "T_C":
                    continue
                if abs(float(row[0]) - 400.0) < 1e-9 and abs(float(row[1]) - 9.41) < 1e-6:
                    self.k_obs_exp = float(row[2]) * 1.0e8
        if self.k_obs_exp is None:
            raise RuntimeError("Could not find the 400 C / 9.41 mM Fig. 4B observation")

        self.rel_err = abs(self.k_obs_model - self.k_obs_exp) / self.k_obs_exp

    def testValidation(self):
        upper = float(self.getParam("validation_upper_bound"))
        self.addScalarData(
            key="k_obs_model",
            value=self.k_obs_model,
            description="Modeled pseudo-first-order rate of e_sol capture by 9.41 mM Zn(II) at 400 C",
            units="1/s",
        )
        self.addScalarData(
            key="k_obs_relative_error",
            value=self.rel_err,
            description="Relative error of the modeled pseudo-first-order rate vs the Fig. 4B measurement",
            units="-",
            bounds=(0.0, upper),
        )
