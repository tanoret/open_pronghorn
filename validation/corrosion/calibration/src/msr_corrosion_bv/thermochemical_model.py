"""Concrete MSTDB-TC thermochemical corrosion model."""

from __future__ import annotations

from .thermochemical_equilibrium import ThermochemicalEquilibriumMixin
from .thermochemical_transport import ThermochemicalTransportMixin


class MSTDBThermochemicalCorrosionModel(
    ThermochemicalEquilibriumMixin,
    ThermochemicalTransportMixin,
):
    """Species-resolved corrosion model driven by MSTDB-TC reaction affinities."""
