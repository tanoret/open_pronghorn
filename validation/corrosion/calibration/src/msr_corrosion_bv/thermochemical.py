"""MSTDB-TC-grounded corrosion and species-inventory model.

The implementation is split into data, equilibrium, transport, and fitting modules
so each physical layer can be reviewed and tested independently.
"""

from .thermochemical_data import (
    ALLOY_ELEMENT_MASS_FRACTIONS,
    CHLORIDE_SPECIES,
    ELEMENT_MOLAR_MASS_G_MOL,
    FLUORIDE_SPECIES,
    HALIDE_MOLAR_MASS_G_MOL,
    MEAN_SALT_MOLAR_MASS_G_MOL,
    THERMOCHEMICAL_PARAMETER_NAMES,
    THERMOCHEMICAL_PARAMETER_SPECS,
    ThermochemicalParameterSpec,
    thermochemical_initial_vector,
    thermochemical_lower_bounds,
    thermochemical_upper_bounds,
    thermochemical_vector_to_params,
)
from .thermochemical_fit import (
    compare_thermochemical_models,
    fit_thermochemical_model,
    relation_aware_log_error,
    save_thermochemical_outputs,
    thermochemical_physics_residuals,
    thermochemical_prior_residuals,
    thermochemical_residuals_for_targets,
    validate_held_out_measurements,
)
from .thermochemical_model import MSTDBThermochemicalCorrosionModel

__all__ = [
    "ALLOY_ELEMENT_MASS_FRACTIONS",
    "CHLORIDE_SPECIES",
    "ELEMENT_MOLAR_MASS_G_MOL",
    "FLUORIDE_SPECIES",
    "HALIDE_MOLAR_MASS_G_MOL",
    "MEAN_SALT_MOLAR_MASS_G_MOL",
    "MSTDBThermochemicalCorrosionModel",
    "THERMOCHEMICAL_PARAMETER_NAMES",
    "THERMOCHEMICAL_PARAMETER_SPECS",
    "ThermochemicalParameterSpec",
    "compare_thermochemical_models",
    "fit_thermochemical_model",
    "relation_aware_log_error",
    "save_thermochemical_outputs",
    "thermochemical_initial_vector",
    "thermochemical_lower_bounds",
    "thermochemical_physics_residuals",
    "thermochemical_prior_residuals",
    "thermochemical_residuals_for_targets",
    "thermochemical_upper_bounds",
    "thermochemical_vector_to_params",
    "validate_held_out_measurements",
]
