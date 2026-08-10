"""
Arrhenius fit for chromium bulk diffusivity in INOR-8 / Hastelloy N.

Data source:
J. H. DeVan, M.S. thesis, University of Tennessee (1960),
Table XVI: INOR-8 Loop 1248, 500-h Cr-51 radiotracer exposure.

The model fitted is

    D_B(T) = D0 * exp(-Q / (R*T))

where:
    D_B : chromium bulk diffusivity [cm^2/s]
    T   : absolute temperature [K]
    D0  : pre-exponential factor [cm^2/s]
    Q   : activation energy [J/mol]
    R   : universal gas constant [J/(mol K)]
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# 1. Experimental data from DeVan Table XVI
# ---------------------------------------------------------------------------
# The listed log10(D) values are used directly because they are the clearest
# tabulated representation of the original diffusion coefficients.
specimen = np.arange(1, 16)

temperature_c = np.array(
    [832, 870, 867, 821, 849, 826, 801, 786, 765, 691, 700, 716, 728, 746, 807],
    dtype=float,
)

log10_diffusivity = np.array(
    [
        -12.566,
        -12.364,
        -12.390,
        -12.638,
        -12.478,
        -12.604,
        -12.710,
        -12.798,
        -12.915,
        -13.405,
        -13.320,
        -13.210,
        -13.158,
        -13.091,
        -12.658,
    ],
    dtype=float,
)

temperature_k = temperature_c + 273.15
diffusivity_cm2_s = 10.0 ** log10_diffusivity


# ---------------------------------------------------------------------------
# 2. Fit ln(D) = ln(D0) - Q/(R*T)
# ---------------------------------------------------------------------------
R = 8.31446261815324  # J/(mol K)

x_fit = 1.0 / temperature_k
y_fit = np.log(diffusivity_cm2_s)

slope, intercept = np.polyfit(x_fit, y_fit, 1)

D0 = np.exp(intercept)
Q_j_mol = -slope * R
Q_kj_mol = Q_j_mol / 1000.0

y_pred = intercept + slope * x_fit
ss_res = np.sum((y_fit - y_pred) ** 2)
ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
r_squared = 1.0 - ss_res / ss_tot

# ---------------------------------------------------------------------------
# Published chromium diffusivity correlations for SUS316 / SUS316L
# ---------------------------------------------------------------------------
# Source:
# Mizouchi et al., "Low Temperature Grain Boundary Diffusion of Chromium
# in SUS316 and 316L Stainless Steels," Materials Transactions (2004).
#
# All pre-exponential factors below have been converted from m^2/s to cm^2/s.
#
# General form:
#
#     D(T) = D0 * exp(-Q / (R*T))
#
# where:
#     D  : diffusivity [cm^2/s]
#     T  : absolute temperature [K]
#     D0 : pre-exponential factor [cm^2/s]
#     Q  : activation energy [J/mol]


# Bulk / lattice chromium diffusivity in solution-treated SUS316
D0_bulk_sus316_cm2_s = 1.13e-3
Q_bulk_sus316_j_mol = 234.0e3


# Grain-boundary chromium diffusivity in solution-treated SUS316
D0_gb_sus316_cm2_s = 1.02e-1
Q_gb_sus316_j_mol = 87.5e3


# Grain-boundary chromium diffusivity in solution-treated SUS316L
D0_gb_sus316l_cm2_s = 2.84e-1
Q_gb_sus316l_j_mol = 90.6e3


def diffusivity_arrhenius(
    temperature_k,
    D0_cm2_s,
    Q_j_mol,
):
    """
    Evaluate an Arrhenius diffusivity correlation.

    Parameters
    ----------
    temperature_k : float or numpy.ndarray
        Absolute temperature [K].
    D0_cm2_s : float
        Pre-exponential factor [cm^2/s].
    Q_j_mol : float
        Activation energy [J/mol].

    Returns
    -------
    float or numpy.ndarray
        Diffusivity [cm^2/s].
    """
    temperature_k = np.asarray(temperature_k, dtype=float)

    if np.any(temperature_k <= 0.0):
        raise ValueError("Temperature must be greater than 0 K.")

    return D0_cm2_s * np.exp(
        -Q_j_mol / (R * temperature_k)
    )


def D_bulk_sus316(temperature_k):
    """Bulk chromium diffusivity in solution-treated SUS316 [cm^2/s]."""
    return diffusivity_arrhenius(
        temperature_k,
        D0_bulk_sus316_cm2_s,
        Q_bulk_sus316_j_mol,
    )


def D_gb_sus316(temperature_k):
    """Grain-boundary chromium diffusivity in solution-treated SUS316 [cm^2/s]."""
    return diffusivity_arrhenius(
        temperature_k,
        D0_gb_sus316_cm2_s,
        Q_gb_sus316_j_mol,
    )


def D_gb_sus316l(temperature_k):
    """Grain-boundary chromium diffusivity in solution-treated SUS316L [cm^2/s]."""
    return diffusivity_arrhenius(
        temperature_k,
        D0_gb_sus316l_cm2_s,
        Q_gb_sus316l_j_mol,
    )


# ---------------------------------------------------------------------------
# 3. Print a reproducible summary
# ---------------------------------------------------------------------------
print("INOR-8 / Hastelloy N chromium diffusivity fit")
print("------------------------------------------------")
print(f"Number of data points : {len(temperature_k)}")
print(f"Temperature range     : {temperature_c.min():.0f} to "
      f"{temperature_c.max():.0f} degC")
print(f"D0                    : {D0:.6e} cm^2/s")
print(f"Q                     : {Q_kj_mol:.3f} kJ/mol")
print(f"R^2                   : {r_squared:.6f}")
print()
print("Correlation:")
print(f"D_B(T) = {D0:.6e} * exp[-{Q_j_mol:.3f} / (R*T)] cm^2/s")


# ---------------------------------------------------------------------------
# 4. Make the Arrhenius plot
# ---------------------------------------------------------------------------
arrhenius_x = 1000.0 / temperature_k

temperature_line_k = np.linspace(temperature_k.min(), temperature_k.max(), 500)
arrhenius_line_x = 1000.0 / temperature_line_k
diffusivity_line = D0 * np.exp(-Q_j_mol / (R * temperature_line_k))
log10_line = np.log10(diffusivity_line)

fig, ax = plt.subplots(figsize=(8.2, 6.1))

ax.scatter(
    arrhenius_x,
    log10_diffusivity,
    s=58,
    marker="o",
    label="DeVan Table XVI data",
    zorder=3,
)

ax.plot(
    arrhenius_line_x,
    log10_line,
    linewidth=2.2,
    label="Arrhenius fit",
    zorder=2,
)

ax.set_xlabel(r"$1000/T$ (K$^{-1}$)", fontsize=13)
ax.set_ylabel(r"$\log_{10}\!\left[D_B\;(\mathrm{cm^2\,s^{-1}})\right]$", fontsize=13)
ax.set_title(
    "Chromium Bulk Diffusivity in INOR-8 / Hastelloy N",
    fontsize=15,
    pad=14,
)

ax.grid(True, which="major", alpha=0.28)
ax.minorticks_on()
ax.tick_params(axis="both", which="major", labelsize=11, direction="in")
ax.tick_params(axis="both", which="minor", direction="in")
ax.legend(frameon=True, fontsize=10, loc="best")

fit_text = (
    r"$D_B(T)=D_0\exp\!\left(-Q/RT\right)$"
    "\n"
    rf"$D_0={D0:.3e}\ \mathrm{{cm^2\,s^{{-1}}}}$"
    "\n"
    rf"$Q={Q_kj_mol:.1f}\ \mathrm{{kJ\,mol^{{-1}}}}$"
    "\n"
    rf"$R^2={r_squared:.4f}$"
)

ax.text(
    0.04,
    0.05,
    fit_text,
    transform=ax.transAxes,
    fontsize=11,
    va="bottom",
    bbox={"boxstyle": "round", "alpha": 0.85},
)

# Add a temperature scale across the top.
def arrhenius_to_celsius(value):
    value = np.asarray(value, dtype=float)
    return 1000.0 / value - 273.15


def celsius_to_arrhenius(value):
    value = np.asarray(value, dtype=float)
    return 1000.0 / (value + 273.15)


top_axis = ax.secondary_xaxis(
    "top",
    functions=(arrhenius_to_celsius, celsius_to_arrhenius),
)
top_axis.set_xlabel("Temperature (degC)", fontsize=12, labelpad=9)
top_axis.tick_params(labelsize=10, direction="in")

fig.tight_layout()

output_directory = Path(__file__).resolve().parent
png_path = output_directory / "inor8_chromium_arrhenius_fit.png"
pdf_path = output_directory / "inor8_chromium_arrhenius_fit.pdf"

fig.savefig(png_path, dpi=400, bbox_inches="tight")
fig.savefig(pdf_path, bbox_inches="tight")
plt.close(fig)

print()
print(f"Saved plot: {png_path}")
print(f"Saved plot: {pdf_path}")

"""
Arrhenius fit for chromium bulk diffusivity in INOR-8 / Hastelloy N.

Data source:
J. H. DeVan, M.S. thesis, University of Tennessee (1960),
Table XVI: INOR-8 Loop 1248, 500-h Cr-51 radiotracer exposure.

The model fitted is

    D_B(T) = D0 * exp(-Q / (R*T))

where:
    D_B : chromium bulk diffusivity [cm^2/s]
    T   : absolute temperature [K]
    D0  : pre-exponential factor [cm^2/s]
    Q   : activation energy [J/mol]
    R   : universal gas constant [J/(mol K)]
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# 1. Experimental data from DeVan Table XVI
# ---------------------------------------------------------------------------
# The listed log10(D) values are used directly because they are the clearest
# tabulated representation of the original diffusion coefficients.
specimen = np.arange(1, 16)

temperature_c = np.array(
    [832, 870, 867, 821, 849, 826, 801, 786, 765, 691, 700, 716, 728, 746, 807],
    dtype=float,
)

log10_diffusivity = np.array(
    [
        -12.566,
        -12.364,
        -12.390,
        -12.638,
        -12.478,
        -12.604,
        -12.710,
        -12.798,
        -12.915,
        -13.405,
        -13.320,
        -13.210,
        -13.158,
        -13.091,
        -12.658,
    ],
    dtype=float,
)

temperature_k = temperature_c + 273.15
diffusivity_cm2_s = 10.0 ** log10_diffusivity


# ---------------------------------------------------------------------------
# 2. Fit ln(D) = ln(D0) - Q/(R*T)
# ---------------------------------------------------------------------------
R = 8.31446261815324  # J/(mol K)

x_fit = 1.0 / temperature_k
y_fit = np.log(diffusivity_cm2_s)

slope, intercept = np.polyfit(x_fit, y_fit, 1)

D0 = np.exp(intercept)
Q_j_mol = -slope * R
Q_kj_mol = Q_j_mol / 1000.0

y_pred = intercept + slope * x_fit
ss_res = np.sum((y_fit - y_pred) ** 2)
ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
r_squared = 1.0 - ss_res / ss_tot

# ---------------------------------------------------------------------------
# Published chromium diffusivity correlations for SUS316 / SUS316L
# ---------------------------------------------------------------------------
# Source:
# Mizouchi et al., "Low Temperature Grain Boundary Diffusion of Chromium
# in SUS316 and 316L Stainless Steels," Materials Transactions (2004).
#
# All pre-exponential factors below have been converted from m^2/s to cm^2/s.
#
# General form:
#
#     D(T) = D0 * exp(-Q / (R*T))
#
# where:
#     D  : diffusivity [cm^2/s]
#     T  : absolute temperature [K]
#     D0 : pre-exponential factor [cm^2/s]
#     Q  : activation energy [J/mol]


# Bulk / lattice chromium diffusivity in solution-treated SUS316
D0_bulk_sus316_cm2_s = 1.13e-3
Q_bulk_sus316_j_mol = 234.0e3


# Grain-boundary chromium diffusivity in solution-treated SUS316
D0_gb_sus316_cm2_s = 1.02e-1
Q_gb_sus316_j_mol = 87.5e3


# Grain-boundary chromium diffusivity in solution-treated SUS316L
D0_gb_sus316l_cm2_s = 2.84e-1
Q_gb_sus316l_j_mol = 90.6e3


def diffusivity_arrhenius(
    temperature_k,
    D0_cm2_s,
    Q_j_mol,
):
    """
    Evaluate an Arrhenius diffusivity correlation.

    Parameters
    ----------
    temperature_k : float or numpy.ndarray
        Absolute temperature [K].
    D0_cm2_s : float
        Pre-exponential factor [cm^2/s].
    Q_j_mol : float
        Activation energy [J/mol].

    Returns
    -------
    float or numpy.ndarray
        Diffusivity [cm^2/s].
    """
    temperature_k = np.asarray(temperature_k, dtype=float)

    if np.any(temperature_k <= 0.0):
        raise ValueError("Temperature must be greater than 0 K.")

    return D0_cm2_s * np.exp(
        -Q_j_mol / (R * temperature_k)
    )


def D_bulk_sus316(temperature_k):
    """Bulk chromium diffusivity in solution-treated SUS316 [cm^2/s]."""
    return diffusivity_arrhenius(
        temperature_k,
        D0_bulk_sus316_cm2_s,
        Q_bulk_sus316_j_mol,
    )


def D_gb_sus316(temperature_k):
    """Grain-boundary chromium diffusivity in solution-treated SUS316 [cm^2/s]."""
    return diffusivity_arrhenius(
        temperature_k,
        D0_gb_sus316_cm2_s,
        Q_gb_sus316_j_mol,
    )


def D_gb_sus316l(temperature_k):
    """Grain-boundary chromium diffusivity in solution-treated SUS316L [cm^2/s]."""
    return diffusivity_arrhenius(
        temperature_k,
        D0_gb_sus316l_cm2_s,
        Q_gb_sus316l_j_mol,
    )


# ---------------------------------------------------------------------------
# 3. Print a reproducible summary
# ---------------------------------------------------------------------------
print("INOR-8 / Hastelloy N chromium diffusivity fit")
print("------------------------------------------------")
print(f"Number of data points : {len(temperature_k)}")
print(f"Temperature range     : {temperature_c.min():.0f} to "
      f"{temperature_c.max():.0f} degC")
print(f"D0                    : {D0:.6e} cm^2/s")
print(f"Q                     : {Q_kj_mol:.3f} kJ/mol")
print(f"R^2                   : {r_squared:.6f}")
print()
print("Correlation:")
print(f"D_B(T) = {D0:.6e} * exp[-{Q_j_mol:.3f} / (R*T)] cm^2/s")


# ---------------------------------------------------------------------------
# 4. Make the Arrhenius plot
# ---------------------------------------------------------------------------
arrhenius_x = 1000.0 / temperature_k

temperature_line_k = np.linspace(temperature_k.min(), temperature_k.max(), 500)
arrhenius_line_x = 1000.0 / temperature_line_k
diffusivity_line = D0 * np.exp(-Q_j_mol / (R * temperature_line_k))
log10_line = np.log10(diffusivity_line)

fig, ax = plt.subplots(figsize=(8.2, 6.1))

ax.scatter(
    arrhenius_x,
    log10_diffusivity,
    s=58,
    marker="o",
    label="DeVan Table XVI data",
    zorder=3,
)

ax.plot(
    arrhenius_line_x,
    log10_line,
    linewidth=2.2,
    label="Arrhenius fit",
    zorder=2,
)

ax.set_xlabel(r"$1000/T$ (K$^{-1}$)", fontsize=13)
ax.set_ylabel(r"$\log_{10}\!\left[D_B\;(\mathrm{cm^2\,s^{-1}})\right]$", fontsize=13)
ax.set_title(
    "Chromium Bulk Diffusivity in INOR-8 / Hastelloy N",
    fontsize=15,
    pad=14,
)

ax.grid(True, which="major", alpha=0.28)
ax.minorticks_on()
ax.tick_params(axis="both", which="major", labelsize=11, direction="in")
ax.tick_params(axis="both", which="minor", direction="in")
ax.legend(frameon=True, fontsize=10, loc="best")

fit_text = (
    r"$D_B(T)=D_0\exp\!\left(-Q/RT\right)$"
    "\n"
    rf"$D_0={D0:.3e}\ \mathrm{{cm^2\,s^{{-1}}}}$"
    "\n"
    rf"$Q={Q_kj_mol:.1f}\ \mathrm{{kJ\,mol^{{-1}}}}$"
    "\n"
    rf"$R^2={r_squared:.4f}$"
)

ax.text(
    0.04,
    0.05,
    fit_text,
    transform=ax.transAxes,
    fontsize=11,
    va="bottom",
    bbox={"boxstyle": "round", "alpha": 0.85},
)

# Add a temperature scale across the top.
def arrhenius_to_celsius(value):
    value = np.asarray(value, dtype=float)
    return 1000.0 / value - 273.15


def celsius_to_arrhenius(value):
    value = np.asarray(value, dtype=float)
    return 1000.0 / (value + 273.15)


top_axis = ax.secondary_xaxis(
    "top",
    functions=(arrhenius_to_celsius, celsius_to_arrhenius),
)
top_axis.set_xlabel("Temperature (degC)", fontsize=12, labelpad=9)
top_axis.tick_params(labelsize=10, direction="in")

fig.tight_layout()

output_directory = Path(__file__).resolve().parent
png_path = output_directory / "inor8_chromium_arrhenius_fit.png"
pdf_path = output_directory / "inor8_chromium_arrhenius_fit.pdf"

fig.savefig(png_path, dpi=400, bbox_inches="tight")
fig.savefig(pdf_path, bbox_inches="tight")
plt.close(fig)

print()
print(f"Saved plot: {png_path}")
print(f"Saved plot: {pdf_path}")

# ---------------------------------------------------------------------------
# 4. Compare bulk chromium diffusivity correlations
# ---------------------------------------------------------------------------

# Common temperature range for comparing the two alloys.
temperature_plot_c = np.linspace(600.0, 900.0, 500)
temperature_plot_k = temperature_plot_c + 273.15
arrhenius_plot_x = 1000.0 / temperature_plot_k

# INOR-8 / Hastelloy N bulk diffusivity from the DeVan fit.
D_bulk_inor8_values = D0 * np.exp(
    -Q_j_mol / (R * temperature_plot_k)
)

# SUS316 bulk diffusivity from Mizouchi et al.
D_bulk_sus316_values = D_bulk_sus316(temperature_plot_k)

# Convert diffusivities to log10 values for the Arrhenius plot.
log10_D_bulk_inor8 = np.log10(D_bulk_inor8_values)
log10_D_bulk_sus316 = np.log10(D_bulk_sus316_values)

# Experimental DeVan data coordinates.
arrhenius_data_x = 1000.0 / temperature_k


fig, ax = plt.subplots(figsize=(8.6, 6.4))

# DeVan experimental measurements.
ax.scatter(
    arrhenius_data_x,
    log10_diffusivity,
    s=58,
    marker="o",
    label="INOR-8: DeVan Table XVI data",
    zorder=3,
)

# INOR-8 fitted bulk-diffusion correlation.
ax.plot(
    arrhenius_plot_x,
    log10_D_bulk_inor8,
    linewidth=2.3,
    label="INOR-8 / Hastelloy N bulk fit",
    zorder=2,
)

# SUS316 published bulk-diffusion correlation.
ax.plot(
    arrhenius_plot_x,
    log10_D_bulk_sus316,
    linewidth=2.3,
    label="SUS316 bulk correlation",
    zorder=2,
)


ax.set_xlabel(
    r"$1000/T$ (K$^{-1}$)",
    fontsize=13,
)

ax.set_ylabel(
    r"$\log_{10}\!\left[D_B\;(\mathrm{cm^2\,s^{-1}})\right]$",
    fontsize=13,
)

ax.set_title(
    "Chromium Bulk Diffusivity in INOR-8 and SUS316",
    fontsize=15,
    pad=14,
)

ax.grid(
    True,
    which="major",
    alpha=0.28,
)

ax.minorticks_on()

ax.tick_params(
    axis="both",
    which="major",
    labelsize=11,
    direction="in",
)

ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
)

ax.legend(
    frameon=True,
    fontsize=10,
    loc="best",
)


# Display both Arrhenius expressions on the plot.
fit_text = (
    r"$\mathrm{INOR\!-\!8:}$"
    "\n"
    rf"$D_B={D0:.3e}"
    rf"\exp\!\left(-{Q_kj_mol:.1f}\times10^3/RT\right)$"
    "\n"
    r"$\mathrm{SUS316:}$"
    "\n"
    rf"$D_B={D0_bulk_sus316_cm2_s:.3e}"
    rf"\exp\!\left(-{Q_bulk_sus316_j_mol / 1000.0:.1f}"
    r"\times10^3/RT\right)$"
    "\n"
    r"$D_B$ in $\mathrm{cm^2\,s^{-1}}$"
)

ax.text(
    0.04,
    0.05,
    fit_text,
    transform=ax.transAxes,
    fontsize=10.5,
    va="bottom",
    bbox={
        "boxstyle": "round",
        "alpha": 0.85,
    },
)


# ---------------------------------------------------------------------------
# Add a temperature scale across the top
# ---------------------------------------------------------------------------
def arrhenius_to_celsius(value):
    value = np.asarray(value, dtype=float)
    return 1000.0 / value - 273.15


def celsius_to_arrhenius(value):
    value = np.asarray(value, dtype=float)
    return 1000.0 / (value + 273.15)


top_axis = ax.secondary_xaxis(
    "top",
    functions=(
        arrhenius_to_celsius,
        celsius_to_arrhenius,
    ),
)

top_axis.set_xlabel(
    "Temperature (degC)",
    fontsize=12,
    labelpad=9,
)

top_axis.tick_params(
    labelsize=10,
    direction="in",
)


fig.tight_layout()


# ---------------------------------------------------------------------------
# Save the comparison plot
# ---------------------------------------------------------------------------
output_directory = Path(__file__).resolve().parent

png_path = (
    output_directory
    / "chromium_bulk_diffusivity_inor8_sus316.png"
)

pdf_path = (
    output_directory
    / "chromium_bulk_diffusivity_inor8_sus316.pdf"
)

fig.savefig(
    png_path,
    dpi=400,
    bbox_inches="tight",
)

fig.savefig(
    pdf_path,
    bbox_inches="tight",
)

plt.close(fig)

print()
print(f"Saved plot: {png_path}")
print(f"Saved plot: {pdf_path}")