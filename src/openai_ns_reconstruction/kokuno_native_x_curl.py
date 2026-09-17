"""Native-X localized complete-curl oscillatory correction.

This module makes one deliberately narrow structural change to the retained
Agent-2 oscillatory surrogate: localization is moved from a fixed Cartesian box
to Kokuno's public native similarity coordinate

    X = (x^2 + y^2) / (2 q),
    q - z^2 q^(2h) = 1 - t.

The corrected 2026-09-09 reconstruction requires oscillatory localization to be
performed at the vector-potential level and the *complete* curl to retain the
cutoff-gradient remainder.  We preserve that contract exactly.  The particular
compact annulus used here,

    X0 < X < 2 X0,   X0 = 4 / Lambda,

is an autonomous, predeclared reconstruction choice chosen for a first transfer
probe; it is not claimed to be a Kokuno/OpenAI parameter.  The affine Cartesian
phase, polarization, frequency and numerical amplitude likewise remain the
existing autonomous Agent-2 surrogate choices.

With a constant phase covector ``n`` and transverse direction ``t``, define

    d = (n x t) / |n|^2,
    A = -a B(X) d sin(psi),
    psi = n . (x-center) + omega*time + phase.

The velocity is evaluated as the analytic complete curl

    curl A = a [B t cos(psi) - (grad B x d) sin(psi)].

Since the field is returned as a curl of a smooth compactly supported potential,
its divergence vanishes identically.  The inner support radius is strictly
positive in X, so the implementation is also axis safe without cylindrical
1/r divisions.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from math import pi
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-native-x-complete-curl-v1"

ArrayLike = Any
VelocityCallable = Callable[[ArrayLike, ArrayLike, ArrayLike, ArrayLike], np.ndarray]


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _triple(values: tuple[float, float, float], name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != (3,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain exactly three finite values")
    return array


def _compact_c4(value: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(1-s^2)^5_+`` and its first derivative in ``s``."""
    inside = np.abs(value) < 1.0
    base = np.where(inside, 1.0 - value * value, 0.0)
    window = base**5
    derivative = np.where(inside, -10.0 * value * base**4, 0.0)
    return window, derivative


@dataclass(frozen=True)
class KokunoNativeXCompleteCurlCorrection:
    """One exact-curl wave with a compact window in source-native ``X``.

    Only the coordinate used for localization is source native.  The support
    multipliers ``(1,2)`` and all wave parameters are autonomous reconstruction
    choices and are intentionally frozen in the first K2-OSC-006 screen.
    """

    amplitude: float = 0.125
    h: float = 0.005
    Lambda: float = 10.0
    wave_vector: tuple[float, float, float] = (5.0, -3.0, 4.0)
    polarization: tuple[float, float, float] = (1.0, 2.0, 1.0)
    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    omega: float = 2.0
    phase: float = 0.0

    _coordinates: KokunoNativeSimilarityCoordinates = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        amplitude = float(self.amplitude)
        h = float(self.h)
        Lambda = float(self.Lambda)
        omega = float(self.omega)
        phase = float(self.phase)
        if not np.isfinite(amplitude) or abs(amplitude) > 2.0:
            raise ValueError("amplitude must be finite and lie in [-2, 2]")
        if not np.isfinite(h) or not (0.0 < h < 1.0e-2):
            raise ValueError("h must satisfy the public source bound 0<h<1e-2")
        if not np.isfinite(Lambda) or not (4.0 <= Lambda <= 32.0):
            raise ValueError("Lambda must be finite and lie in the autonomous bound [4,32]")
        if not np.isfinite(omega) or abs(omega) > 16.0 * pi:
            raise ValueError("omega must be finite and lie in [-16*pi,16*pi]")
        if not np.isfinite(phase) or abs(phase) > pi:
            raise ValueError("phase must be finite and lie in [-pi,pi]")

        n = _triple(self.wave_vector, "wave_vector")
        if np.max(np.abs(n)) > 32.0 or float(np.dot(n, n)) < 1.0e-12:
            raise ValueError("wave_vector components must lie in [-32,32] and be nonzero")
        p = _triple(self.polarization, "polarization")
        transverse = p - n * (float(np.dot(p, n)) / float(np.dot(n, n)))
        if float(np.linalg.norm(transverse)) < 1.0e-10:
            raise ValueError("polarization must have a transverse component")
        _triple(self.center, "center")

        object.__setattr__(self, "amplitude", amplitude)
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "Lambda", Lambda)
        object.__setattr__(self, "omega", omega)
        object.__setattr__(self, "phase", phase)
        object.__setattr__(self, "_coordinates", KokunoNativeSimilarityCoordinates(h=h))

    @property
    def X0(self) -> float:
        """Source reference-continuation start ``X0=4/Lambda``."""
        return 4.0 / self.Lambda

    @property
    def X_support(self) -> tuple[float, float]:
        """Autonomous first annular probe: ``(X0, 2 X0)``."""
        return self.X0, 2.0 * self.X0

    @property
    def transverse_polarization(self) -> np.ndarray:
        n = _triple(self.wave_vector, "wave_vector")
        p = _triple(self.polarization, "polarization")
        transverse = p - n * (float(np.dot(p, n)) / float(np.dot(n, n)))
        return transverse / np.linalg.norm(transverse)

    @property
    def potential_direction(self) -> np.ndarray:
        n = _triple(self.wave_vector, "wave_vector")
        transverse = self.transverse_polarization
        return np.cross(n, transverse) / float(np.dot(n, n))

    @staticmethod
    def _finite_array(value: Any, name: str) -> np.ndarray:
        array = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{name} must contain only finite values")
        return array

    def _broadcast(self, x: Any, y: Any, z: Any, time: Any):
        return np.broadcast_arrays(
            self._finite_array(x, "x"),
            self._finite_array(y, "y"),
            self._finite_array(z, "z"),
            self._finite_array(time, "time"),
        )

    def envelope_data(
        self, x: Any, y: Any, z: Any, time: Any
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return native ``X``, compact envelope ``B`` and Cartesian ``grad B``."""
        x, y, z, time = self._broadcast(x, y, z, time)
        native = self._coordinates.evaluate(x, y, z, time)
        q = np.asarray(native["q"], dtype=float)
        X = np.asarray(native["X"], dtype=float)

        lower, upper = self.X_support
        midpoint = 0.5 * (lower + upper)
        half_width = 0.5 * (upper - lower)
        s = (X - midpoint) / half_width
        envelope, dwindow_ds = _compact_c4(s)
        dB_dX = dwindow_ds / half_width

        grad_X = np.stack(
            (
                x / q,
                y / q,
                np.asarray(native["X_z"], dtype=float),
            ),
            axis=-1,
        )
        grad_envelope = dB_dX[..., None] * grad_X
        return X, envelope, grad_envelope

    def _phase(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, time: np.ndarray) -> np.ndarray:
        center = _triple(self.center, "center")
        shifted = np.stack((x, y, z), axis=-1) - center
        n = _triple(self.wave_vector, "wave_vector")
        return np.einsum("...i,i->...", shifted, n) + self.omega * time + self.phase

    def vector_potential(self, x: Any, y: Any, z: Any, time: Any) -> np.ndarray:
        """Evaluate the source-coordinate-localized vector potential."""
        x, y, z, time = self._broadcast(x, y, z, time)
        _, envelope, _ = self.envelope_data(x, y, z, time)
        psi = self._phase(x, y, z, time)
        direction = self.potential_direction
        return (
            -self.amplitude
            * envelope[..., None]
            * np.sin(psi)[..., None]
            * direction
        )

    def velocity(self, x: Any, y: Any, z: Any, time: Any) -> np.ndarray:
        """Evaluate the analytic complete curl ``curl A`` in ``(...,3)`` layout."""
        x, y, z, time = self._broadcast(x, y, z, time)
        _, envelope, grad_envelope = self.envelope_data(x, y, z, time)
        psi = self._phase(x, y, z, time)
        transverse = self.transverse_polarization
        direction = self.potential_direction
        leading = envelope[..., None] * np.cos(psi)[..., None] * transverse
        cutoff_remainder = -np.cross(grad_envelope, direction) * np.sin(psi)[..., None]
        return self.amplitude * (leading + cutoff_remainder)

    __call__ = velocity

    def at_points(self, points: ArrayLike, time: ArrayLike) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
            raise ValueError("points must have shape (n,3) and contain finite values")
        return self.velocity(points[:, 0], points[:, 1], points[:, 2], time)

    def compose_velocity(self, base_velocity: VelocityCallable) -> VelocityCallable:
        """Return ``base_velocity + u_osc`` without mutating either component."""
        if not callable(base_velocity):
            raise TypeError("base_velocity must be callable")

        def composed(x: Any, y: Any, z: Any, time: Any) -> np.ndarray:
            base = np.asarray(base_velocity(x, y, z, time), dtype=float)
            correction = self.velocity(x, y, z, time)
            if base.shape != correction.shape:
                raise ValueError(
                    "base velocity must use the same (...,3) layout as the correction"
                )
            return base + correction

        return composed

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "source_derived": [
                    "native similarity X=(x^2+y^2)/(2q)",
                    "reference start X0=4/Lambda",
                    "potential-level localization",
                    "complete curl including cutoff-gradient remainder",
                ],
            },
            "autonomous_choices": {
                "X_support": [self.X0, 2.0 * self.X0],
                "support_rule": "predeclared first probe X0<X<2*X0; not source exact",
                "amplitude": self.amplitude,
                "wave_vector": list(self.wave_vector),
                "polarization": list(self.polarization),
                "center": list(self.center),
                "omega": self.omega,
                "phase": self.phase,
                "phase_frame": "retained affine Cartesian surrogate",
            },
            "parameters": {"h": self.h, "Lambda": self.Lambda},
            "contract": {
                "vector_potential": "A=-a*B(X)*d*sin(psi)",
                "potential_direction": "d=(n cross t)/|n|^2",
                "velocity": "curl(A)=a*(B*t*cos(psi)-(grad B cross d)*sin(psi))",
                "axis_handling": "support excludes X=0; no cylindrical 1/r division",
            },
            "truth_boundary": {
                "source_native_localization_coordinate": True,
                "complete_curl_executable": True,
                "source_exact_support": False,
                "source_exact_phase_frame": False,
                "parameter_selection_performed": False,
                "pde_validated": False,
                "paper_exact": False,
                "openai_field_identified": False,
            },
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target
