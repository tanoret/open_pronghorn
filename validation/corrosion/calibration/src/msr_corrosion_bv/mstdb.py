"""MSTDB-TC/ChemSage thermodynamic data access for corrosion modeling.

The parser reads function-expanded ChemSage ``*_No_Func.dat`` files and evaluates
standard molar Gibbs energies with the same polynomial form used by
Thermochimica.  It intentionally does not reproduce the full Gibbs-energy
minimizer or the modified-quasichemical solution model; those are available
through :class:`ThermochimicaRunner` when a native Thermochimica executable is
installed.  The internal backend is sufficient for standard-state reaction
energies, Nernst activity corrections, pure-solid/liquid saturation ratios, and
reproducible unit tests without redistributing MSTDB-TC.
"""

from __future__ import annotations

import functools
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

R_GAS = 8.31446261815324
_SUPPORTED_EQUATION_TYPES = {1, 4, 10, 13, 16}
CALIBRATED_MSTDB_VERSION = "3.1"
CALIBRATED_MSTDB_FLUORIDE_SHA256 = (
    "a7ca53e6061c3aa16d1932e35f908fa61b328b88625e8361e63881bfc8592000"
)
CALIBRATED_MSTDB_CHLORIDE_SHA256 = (
    "460e8d7148cc76c597064e167b4afa98ea762e59a666f9c42962f3d2d013e5c1"
)


def _as_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


@dataclass(frozen=True)
class GibbsInterval:
    upper_temperature_K: float
    coefficients: tuple[float, float, float, float, float, float]
    additional_terms: tuple[tuple[float, float], ...] = ()

    def evaluate(self, temperature_K: float) -> float:
        """Return the standard molar Gibbs energy in J/mol."""
        if temperature_K <= 0.0 or not math.isfinite(temperature_K):
            raise ValueError("temperature_K must be finite and positive")
        A, B, C, D, E, F = self.coefficients
        T = float(temperature_K)
        value = A + B * T + C * T * math.log(T) + D * T**2 + E * T**3 + F / T
        for coefficient, exponent in self.additional_terms:
            value += coefficient * math.log(T) if exponent == 99.0 else coefficient * T**exponent
        return value


@dataclass(frozen=True)
class ChemSageSpecies:
    name: str
    equation_type: int
    stoichiometry: tuple[float, ...]
    intervals: tuple[GibbsInterval, ...]
    magnetic_parameters: tuple[float, ...] = ()
    source_line: int = 0

    def standard_gibbs_J_mol(
        self,
        temperature_K: float,
        include_magnetic: bool = True,
        allow_extrapolation: bool = False,
    ) -> float:
        interval = next(
            (item for item in self.intervals if temperature_K <= item.upper_temperature_K),
            None,
        )
        if interval is None:
            if not allow_extrapolation:
                raise ValueError(
                    f"temperature_K={temperature_K} exceeds the last standard-state interval "
                    f"({self.intervals[-1].upper_temperature_K} K) for species {self.name!r}; "
                    "set allow_extrapolation=True only with a documented basis"
                )
            interval = self.intervals[-1]
        value = interval.evaluate(temperature_K)
        if include_magnetic and len(self.magnetic_parameters) >= 4:
            value += magnetic_gibbs_J_mol(temperature_K, *self.magnetic_parameters[:4])
        return value


def magnetic_gibbs_J_mol(
    temperature_K: float,
    critical_temperature_K: float,
    average_magnetic_moment: float,
    structure_factor: float,
    p: float,
) -> float:
    """Evaluate the SGTE magnetic contribution used by Thermochimica.

    Adapted from Thermochimica ``src/gem/CompGibbsMagnetic.f90`` at commit
    ``0c35c8d7d1cf2084b4e2ca5d6608f7dcdf60adad``; see ``THIRD_PARTY_NOTICES.md``.
    """
    if critical_temperature_K == 0.0 or p == 0.0:
        return 0.0
    Tcritical = float(critical_temperature_K)
    B = float(average_magnetic_moment)
    if Tcritical < 0.0:
        Tcritical = -Tcritical * structure_factor
        B = -B * structure_factor
    if Tcritical <= 0.0 or B + 1.0 <= 0.0:
        return 0.0
    tau = temperature_K / Tcritical
    invpmone = 1.0 / p - 1.0
    denominator = 518.0 / 1125.0 + (11692.0 / 15975.0) * invpmone
    if tau > 1.0:
        g_tau = -(tau**-5 / 10.0 + tau**-15 / 315.0 + tau**-25 / 1500.0) / denominator
    else:
        g_tau = 1.0 - (
            79.0 / (140.0 * p * tau)
            + (474.0 / 497.0)
            * invpmone
            * (tau**3 / 6.0 + tau**9 / 135.0 + tau**15 / 600.0)
        ) / denominator
    return R_GAS * temperature_K * math.log(B + 1.0) * g_tau


class ChemSageDatabase:
    """Function-expanded ChemSage database with exact standard-state evaluation."""

    def __init__(self, path: str | Path, allow_extrapolation: bool = False) -> None:
        self.path = Path(path).expanduser().resolve()
        self.allow_extrapolation = bool(allow_extrapolation)
        if not self.path.is_file():
            raise FileNotFoundError(self.path)
        self.sha256 = hashlib.sha256(self.path.read_bytes()).hexdigest()
        self.system_line, self.n_elements, records = self._parse(self.path)
        self.records = tuple(records)
        by_name: dict[str, list[ChemSageSpecies]] = {}
        for record in self.records:
            by_name.setdefault(record.name, []).append(record)
        self._by_name = {name: tuple(items) for name, items in by_name.items()}

    @staticmethod
    def _parse(path: Path) -> tuple[str, int, list[ChemSageSpecies]]:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        if len(lines) < 3 or not lines[0].lstrip().startswith("System"):
            raise ValueError(f"{path} does not look like a ChemSage data file")
        header_tokens = lines[1].split()
        if not header_tokens:
            raise ValueError(f"Unable to read ChemSage header in {path}")
        n_elements = int(header_tokens[0])
        records: list[ChemSageSpecies] = []
        index = 2
        while index < len(lines) - 1:
            raw_name = lines[index]
            name = raw_name[:26].strip()
            if not name or not any(char.isalpha() for char in name):
                index += 1
                continue
            equation_tokens = lines[index + 1].split()
            try:
                equation_type = int(equation_tokens[0])
                n_intervals = int(equation_tokens[1])
            except (ValueError, IndexError):
                index += 1
                continue
            if equation_type not in _SUPPORTED_EQUATION_TYPES or not 1 <= n_intervals <= 50:
                index += 1
                continue

            try:
                stoichiometry = [_as_float(item) for item in equation_tokens[2:]]
                cursor = index + 2
                while len(stoichiometry) < n_elements:
                    stoichiometry.extend(_as_float(item) for item in lines[cursor].split())
                    cursor += 1
                stoichiometry = stoichiometry[:n_elements]

                intervals: list[GibbsInterval] = []
                for _ in range(n_intervals):
                    values: list[float] = []
                    while len(values) < 7:
                        values.extend(_as_float(item) for item in lines[cursor].split())
                        cursor += 1
                    upper_temperature = values[0]
                    coefficients = tuple(values[1:7])
                    extra_terms: list[tuple[float, float]] = []
                    if equation_type in {4, 16}:
                        extra_tokens: list[str] = []
                        while not extra_tokens:
                            extra_tokens.extend(lines[cursor].split())
                            cursor += 1
                        n_extra = int(float(extra_tokens[0]))
                        while len(extra_tokens) < 1 + 2 * n_extra:
                            extra_tokens.extend(lines[cursor].split())
                            cursor += 1
                        for extra_index in range(n_extra):
                            coefficient = _as_float(extra_tokens[1 + 2 * extra_index])
                            exponent = _as_float(extra_tokens[2 + 2 * extra_index])
                            extra_terms.append((coefficient, exponent))
                    intervals.append(
                        GibbsInterval(
                            upper_temperature_K=upper_temperature,
                            coefficients=coefficients,  # type: ignore[arg-type]
                            additional_terms=tuple(extra_terms),
                        )
                    )

                magnetic: tuple[float, ...] = ()
                if equation_type in {13, 16}:
                    magnetic_values = tuple(_as_float(item) for item in lines[cursor].split())
                    if len(magnetic_values) >= 2:
                        magnetic = magnetic_values
                        cursor += 1
            except (ValueError, IndexError):
                index += 1
                continue

            records.append(
                ChemSageSpecies(
                    name=name,
                    equation_type=equation_type,
                    stoichiometry=tuple(stoichiometry),
                    intervals=tuple(intervals),
                    magnetic_parameters=magnetic,
                    source_line=index + 1,
                )
            )
            index = cursor
        if not records:
            raise ValueError(f"No Gibbs-energy species blocks were parsed from {path}")
        return lines[0].strip(), n_elements, records

    def occurrences(self, name: str) -> tuple[ChemSageSpecies, ...]:
        return self._by_name.get(name, ())

    def species(self, name: str, occurrence: int = -1) -> ChemSageSpecies:
        records = self.occurrences(name)
        if not records:
            available = ", ".join(sorted(key for key in self._by_name if key.startswith(name[:3]))[:8])
            raise KeyError(f"Species {name!r} not found in {self.path.name}; nearby names: {available}")
        return records[occurrence]

    @functools.lru_cache(maxsize=8192)
    def standard_gibbs_J_mol(
        self,
        name: str,
        temperature_K: float,
        occurrence: int = -1,
        allow_extrapolation: bool | None = None,
    ) -> float:
        extrapolate = self.allow_extrapolation if allow_extrapolation is None else allow_extrapolation
        return self.species(name, occurrence).standard_gibbs_J_mol(
            temperature_K, allow_extrapolation=extrapolate
        )

    def reaction_gibbs_J_mol(self, stoichiometry: Mapping[str, float], temperature_K: float) -> float:
        """Return ``sum(nu_i * G_i)`` with products carrying positive coefficients."""
        return sum(
            float(coefficient) * self.standard_gibbs_J_mol(name, temperature_K)
            for name, coefficient in stoichiometry.items()
        )

    def equilibrium_log_constant(self, stoichiometry: Mapping[str, float], temperature_K: float) -> float:
        return -self.reaction_gibbs_J_mol(stoichiometry, temperature_K) / (R_GAS * temperature_K)

    def metadata(self) -> dict[str, object]:
        version_match = re.search(r"V(\d+(?:\.\d+)*)", self.path.name, flags=re.IGNORECASE)
        return {
            "filename": self.path.name,
            "version": version_match.group(1) if version_match else "unknown",
            "sha256": self.sha256,
            "system": self.system_line,
            "n_elements": self.n_elements,
            "n_species_records": len(self.records),
            "allow_extrapolation": self.allow_extrapolation,
        }


@dataclass(frozen=True)
class MSTDBPair:
    fluoride: ChemSageDatabase
    chloride: ChemSageDatabase

    @classmethod
    def from_directory(
        cls,
        directory: str | Path,
        allow_extrapolation: bool = False,
        *,
        expected_version: str | None = None,
        expected_fluoride_sha256: str | None = None,
        expected_chloride_sha256: str | None = None,
    ) -> "MSTDBPair":
        directory = Path(directory).expanduser().resolve()
        fluoride_candidates = list(directory.rglob("*Fluorides_No_Func.dat"))
        chloride_candidates = list(directory.rglob("*Chlorides_No_Func.dat"))
        if not fluoride_candidates or not chloride_candidates:
            raise FileNotFoundError(
                f"Expected MSTDB-TC *Fluorides_No_Func.dat and *Chlorides_No_Func.dat below {directory}"
            )

        def version(path: Path) -> tuple[int, ...]:
            match = re.search(r"V(\d+(?:\.\d+)*)", path.name, flags=re.IGNORECASE)
            return tuple(int(part) for part in match.group(1).split(".")) if match else ()

        fluoride_by_version: dict[tuple[int, ...], list[Path]] = {}
        chloride_by_version: dict[tuple[int, ...], list[Path]] = {}
        for path in fluoride_candidates:
            if version(path):
                fluoride_by_version.setdefault(version(path), []).append(path)
        for path in chloride_candidates:
            if version(path):
                chloride_by_version.setdefault(version(path), []).append(path)
        common_versions = set(fluoride_by_version).intersection(chloride_by_version)
        if not common_versions:
            raise ValueError(
                "Fluoride and chloride MSTDB-TC files must come from the same database version; "
                f"found fluoride versions {sorted(fluoride_by_version)} and chloride versions {sorted(chloride_by_version)}"
            )
        if expected_version is None:
            selected_version = max(common_versions)
        else:
            selected_version = tuple(int(part) for part in expected_version.split("."))
            if selected_version not in common_versions:
                raise ValueError(
                    f"Required MSTDB-TC V{expected_version} fluoride/chloride pair was not found; "
                    f"available common versions are {sorted(common_versions)}"
                )
        if len(fluoride_by_version[selected_version]) != 1 or len(
            chloride_by_version[selected_version]
        ) != 1:
            raise ValueError(
                f"Ambiguous duplicate MSTDB-TC V{'.'.join(map(str, selected_version))} files: "
                f"fluoride={fluoride_by_version[selected_version]}, "
                f"chloride={chloride_by_version[selected_version]}"
            )
        pair = cls(
            fluoride=ChemSageDatabase(
                fluoride_by_version[selected_version][0],
                allow_extrapolation=allow_extrapolation,
            ),
            chloride=ChemSageDatabase(
                chloride_by_version[selected_version][0],
                allow_extrapolation=allow_extrapolation,
            ),
        )
        for label, actual, expected in (
            ("fluoride", pair.fluoride.sha256, expected_fluoride_sha256),
            ("chloride", pair.chloride.sha256, expected_chloride_sha256),
        ):
            if expected is not None and actual != expected:
                raise ValueError(
                    f"MSTDB-TC {label} SHA-256 mismatch: expected {expected}, got {actual}"
                )
        return pair

    @classmethod
    def from_calibration_directory(
        cls, directory: str | Path, *, allow_extrapolation: bool = False
    ) -> "MSTDBPair":
        """Load the exact V3.1 bytes used by the checked advanced calibration."""
        return cls.from_directory(
            directory,
            allow_extrapolation=allow_extrapolation,
            expected_version=CALIBRATED_MSTDB_VERSION,
            expected_fluoride_sha256=CALIBRATED_MSTDB_FLUORIDE_SHA256,
            expected_chloride_sha256=CALIBRATED_MSTDB_CHLORIDE_SHA256,
        )

    def metadata(self) -> dict[str, object]:
        return {"fluoride": self.fluoride.metadata(), "chloride": self.chloride.metadata()}


class ThermochimicaRunner:
    """Optional native Thermochimica ``InputScriptMode`` adapter.

    The corrosion package does not vendor Thermochimica or MSTDB-TC.  This adapter
    writes a standard input script, launches a user-provided executable, and reads
    ``thermoout.json``.  It is intended for full SUBQ equilibrium checks, while the
    internal parser provides deterministic standard-state calculations.
    """

    def __init__(self, executable: str | Path = "InputScriptMode") -> None:
        candidate = str(executable)
        resolved = shutil.which(candidate) if not Path(candidate).is_file() else str(Path(candidate).resolve())
        if resolved is None:
            raise FileNotFoundError(
                f"Thermochimica executable {candidate!r} was not found. Build the public ORNL-CEES/thermochimica "
                "repository and pass the path to bin/InputScriptMode."
            )
        self.executable = resolved

    def run_elements(
        self,
        database_path: str | Path,
        temperature_K: float,
        element_moles_by_atomic_number: Mapping[int, float],
        *,
        pressure_atm: float = 1.0,
        timeout_s: float = 120.0,
    ) -> dict[str, object]:
        with tempfile.TemporaryDirectory(prefix="thermochimica-") as temp_dir_text:
            temp_dir = Path(temp_dir_text)
            input_path = temp_dir / "corrosion.ti"
            lines = [
                f"pressure          = {pressure_atm:.16g}",
                f"temperature       = {temperature_K:.16g}",
            ]
            for atomic_number, moles in sorted(element_moles_by_atomic_number.items()):
                if moles > 0.0:
                    lines.append(f"mass({int(atomic_number)})           = {float(moles):.16g}")
            lines.extend(
                [
                    "temperature unit  = K",
                    "pressure unit     = atm",
                    "mass unit         = moles",
                    f"data file         = {Path(database_path).expanduser().resolve()}",
                    "print mode        = 0",
                    "debug mode        = .FALSE.",
                    "reinit            = .TRUE.",
                    "write json        = .TRUE.",
                ]
            )
            input_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            completed = subprocess.run(
                [self.executable, str(input_path)],
                cwd=temp_dir,
                text=True,
                capture_output=True,
                timeout=timeout_s,
                check=False,
            )
            if completed.returncode != 0:
                raise RuntimeError(
                    f"Thermochimica exited with status {completed.returncode}:\n{completed.stdout}\n{completed.stderr}"
                )
            output_candidates = [temp_dir / "outputs" / "thermoout.json", temp_dir / "thermoout.json"]
            output_path = next((path for path in output_candidates if path.is_file()), None)
            if output_path is None:
                raise RuntimeError(
                    "Thermochimica completed but no thermoout.json was produced. "
                    f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
                )
            return json.loads(output_path.read_text(encoding="utf-8"))


def resolve_mstdb_directory(explicit: str | Path | None = None) -> Path:
    """Resolve MSTDB-TC from an explicit path or ``MSTDB_TC_DIR``."""
    value = explicit if explicit is not None else os.environ.get("MSTDB_TC_DIR")
    if value is None:
        raise FileNotFoundError(
            "MSTDB-TC is not bundled. Pass --mstdb-dir or set MSTDB_TC_DIR to a directory containing "
            "the fluoride and chloride *_No_Func.dat files."
        )
    return Path(value).expanduser().resolve()
