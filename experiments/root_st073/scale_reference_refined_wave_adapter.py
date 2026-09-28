"""Callable adapter for the selected refined-grid wave-amplitude candidate.

The parent is the frozen corrected reference from
``scale_reference_trust_nonlinear_fit.json``.  This adapter adds exactly one
original mode-1 wave amplitude increment and the selected 90 pressure-gradient
increments from ``scale_reference_refined_wave.json``.  Existing parent
pressure is retained as the baseline; it is not reconstructed or added again.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scale_reference_candidate as candidate_module  # noqa: E402
from scale_generator_pressure_projection import _pressure_patch_design  # noqa: E402
from supported_fourier_analytic_jets import basis_jets  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402


REFINED_REPORT_PATH = ROOT / "scale_reference_refined_wave.json"
ACTUAL_REPLAY_PATH = ROOT / "scale_reference_trust_nonlinear_actual.json"
GEOMETRY_PATH = ROOT / "full_wave_frozen_cache.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _points(value: Any) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim == 1 and array.shape == (3,):
        array = array[None, :]
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    if not np.all(np.isfinite(array)):
        raise ValueError("points must be finite")
    return array


def _decode_pressure_block(block: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    block = np.asarray(block, dtype=float)
    if block.shape != (45,):
        raise ValueError(f"pressure block must have shape (45,), got {block.shape}")
    return (
        block[:9].astype(complex),
        block[9:27:2] + 1j * block[10:27:2],
        block[27:45:2] + 1j * block[28:45:2],
    )


class RefinedWaveReference:
    """Frozen parent plus the selected original-wave amplitude correction."""

    def __init__(self, parent: Any, refined_report: dict, actual_report: dict,
                 geometry: dict, source_hashes: dict[str, str]):
        selected = refined_report["selected"]
        if not bool(selected.get("accepted", False)):
            raise ValueError("Refined wave report does not contain an accepted simultaneous improvement")
        self.reference = parent
        self.base_field = parent.base_field
        self.parent = parent
        self.tau0 = float(parent.tau0)
        self.h = float(parent.h)
        self.nu = float(parent.nu)
        self.amplitude = float(selected["wave_amplitude"])
        self.amplitude_shift = float(selected["amplitude_shift"])
        self.pressure_coefficients = np.asarray(selected["pressure_coefficients"], dtype=float)
        if self.pressure_coefficients.shape != (90,):
            raise ValueError("Selected pressure increment must contain 90 coefficients")
        self.refined_report = refined_report
        self.actual_report = actual_report
        self.geometry = geometry
        self.source_hashes = dict(source_hashes)
        self.wave_center = np.asarray(geometry["center"], dtype=float)
        self.wave_widths = np.asarray(geometry["widths"], dtype=float)
        self.wave_degree = int(geometry["degree"])
        self.wave_mode = int(geometry["mode"])
        self.wave_carrier = np.asarray(geometry["carrier"], dtype=float)
        packed = np.asarray(geometry["coefficients_original"], dtype=float)
        if packed.shape != (27, 2):
            raise ValueError(f"Original wave coefficients must have shape (27,2), got {packed.shape}")
        self.wave_coefficients = packed[:, 0] + 1j * packed[:, 1]
        pressure_geometry = refined_report["pressure_projection"]["geometry"]
        self.pressure_centers = (
            np.asarray(pressure_geometry["first_center"], dtype=float),
            np.asarray(pressure_geometry["second_center"], dtype=float),
        )
        self.pressure_widths = (
            np.asarray(pressure_geometry["first_widths"], dtype=float),
            np.asarray(pressure_geometry["second_widths"], dtype=float),
        )
        self._pressure_modes = (
            _decode_pressure_block(self.pressure_coefficients[:45]),
            _decode_pressure_block(self.pressure_coefficients[45:]),
        )

    def _wave_values(self, points: np.ndarray) -> np.ndarray:
        values = basis_data(
            points, self.wave_center, self.wave_widths,
            self.wave_mode, self.wave_degree, self.wave_carrier,
        )[0]
        return np.einsum("niq,q->ni", values, self.wave_coefficients).real

    def wave_jets(self, points: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the original wave value, Cartesian Jacobian, and Laplacian."""

        points_array = _points(points)
        values, gradients, viscous, _, _ = basis_jets(
            points_array, self.wave_center, self.wave_widths,
            self.wave_mode, self.wave_degree, self.wave_carrier, self.nu,
        )
        value = np.einsum("niq,q->ni", values, self.wave_coefficients).real
        jacobian = np.einsum("nijq,q->nij", gradients, self.wave_coefficients).real
        laplacian = -np.einsum("niq,q->ni", viscous, self.wave_coefficients).real / self.nu
        return value, jacobian, laplacian

    def pressure_increment(self, points: Any) -> np.ndarray:
        points_array = _points(points)
        result = np.zeros(len(points_array), dtype=float)
        for patch in range(2):
            for mode in (0, 1, 2):
                scalars = basis_data(
                    points_array,
                    self.pressure_centers[patch],
                    self.pressure_widths[patch],
                    mode,
                    self.wave_degree,
                    np.zeros(2, dtype=float) if mode == 0 else mode * self.wave_carrier,
                )[1]
                result += np.einsum("nq,q->n", scalars, self._pressure_modes[patch][mode]).real
        return result

    def pressure_gradient_increment(self, points: Any) -> np.ndarray:
        points_array = _points(points)
        design = np.column_stack((
            _pressure_patch_design(points_array, self.pressure_centers[0], self.pressure_widths[0], self.wave_carrier),
            _pressure_patch_design(points_array, self.pressure_centers[1], self.pressure_widths[1], self.wave_carrier),
        ))
        return (design @ self.pressure_coefficients).reshape(-1, 3)

    def fields(self, points: Any, tau: float | None = None) -> tuple[np.ndarray, np.ndarray]:
        points_array = _points(points)
        if tau is not None and not np.isclose(float(tau), self.tau0, rtol=0.0, atol=2.0e-15):
            raise ValueError("Refined wave reference is frozen at tau0")
        parent_velocity, parent_pressure = self.parent.fields(points_array, self.tau0)
        wave_increment = self.amplitude_shift * self._wave_values(points_array)
        pressure_increment = self.pressure_increment(points_array)
        return (
            np.asarray(parent_velocity, dtype=float) + wave_increment,
            np.asarray(parent_pressure, dtype=float) + pressure_increment,
        )

    def velocity(self, points: Any) -> np.ndarray:
        return self.fields(points, self.tau0)[0]

    def pressure(self, points: Any) -> np.ndarray:
        return self.fields(points, self.tau0)[1]


def load_reference(report_path: Path | str = REFINED_REPORT_PATH) -> RefinedWaveReference:
    report_path = Path(report_path)
    refined_report = json.loads(report_path.read_text(encoding="utf-8"))
    if refined_report.get("status") != "completed":
        raise ValueError("Refined wave report is not completed")
    refined_sources = refined_report["sources"]
    replay_name = refined_sources["actual_replay"]["path"]
    actual_path = ROOT / replay_name
    actual_report = json.loads(actual_path.read_text(encoding="utf-8"))
    if actual_report.get("status") != "completed":
        raise ValueError("Parent actual replay is not completed")
    geometry = json.loads(GEOMETRY_PATH.read_text(encoding="utf-8"))["inputs"]["wave"]
    actual_sources = actual_report["sources"]
    expected = {
        "refined_report": (report_path, refined_sources.get("actual_replay", {}).get("sha256")),
        "actual_replay": (actual_path, refined_sources.get("actual_replay", {}).get("sha256")),
        "step": (ROOT / actual_sources["step"]["path"], actual_sources["step"]["sha256"]),
        "candidate_loader": (ROOT / actual_sources["candidate_loader"]["path"], actual_sources["candidate_loader"]["sha256"]),
        "geometry": (GEOMETRY_PATH, refined_sources["geometry"]["sha256"]),
    }
    source_hashes = {}
    for name, (path, expected_hash) in expected.items():
        if name == "refined_report":
            source_hashes[name] = _sha256(path)
            continue
        if not path.exists():
            raise FileNotFoundError(f"Missing {name} source: {path}")
        actual_hash = _sha256(path)
        if expected_hash and actual_hash != expected_hash:
            raise ValueError(f"{name} hash mismatch: expected {expected_hash}, got {actual_hash}")
        source_hashes[name] = actual_hash
    step_path = ROOT / actual_sources["step"]["path"]
    parent = candidate_module.load_reference(step=step_path)
    return RefinedWaveReference(parent, refined_report, actual_report, geometry, source_hashes)


load = load_reference


if __name__ == "__main__":
    reference = load_reference()
    probe = np.array([[0.001, 0.0, 0.0]], dtype=float)
    velocity, pressure = reference.fields(probe, reference.tau0)
    print(json.dumps({
        "status": "loaded",
        "tau0": reference.tau0,
        "amplitude": reference.amplitude,
        "amplitude_shift": reference.amplitude_shift,
        "velocity_shape": list(velocity.shape),
        "pressure_shape": list(pressure.shape),
        "source_hashes": reference.source_hashes,
        "pressure_is_increment": True,
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
    }, sort_keys=True))
