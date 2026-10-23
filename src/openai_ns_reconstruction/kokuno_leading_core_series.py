"""Executable Kokuno leading-core velocity built from the existing nonlinear series.

This module is deliberately an adapter, not a second implementation of the radial
profile equations.  It reuses :class:`PaperCoreSeries`, whose finite nonlinear
near-axis recurrence implements the repository's source-aligned Eqs. (4.9)/(4.13),
and binds it to the Kokuno-native ``q/X/eta`` convention and provenance.

The resulting field is a genuine callable three-dimensional *core* candidate on
``Lambda*X <= 4.1``.  It is not the source's completed outer/heat-matched leading
field, is not compactly supported, and is not an independently validated NS
solution.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_leading_axis_profiles import KokunoLeadingAxisProfile
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates
from .paper_core_reference import PaperCoreReference
from .paper_core_series import PaperCoreSeries


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
CORRECTED_PDF = "released_ns_reader_corrected_20260909.pdf"
SCHEMA = "kokuno-leading-core-series-candidate-v1"
CORE_Y_MAX = 4.1

_SOURCE_MAP = {
    "coordinates": "tau=1-t=q(1-eta^2), z=q^(1/2-h) eta, X=(x^2+y^2)/(2q)",
    "q_equation": "q-z^2*q^(2h)=1-t",
    "cartesian_velocity": (
        "u1=x*v0/(2q)-y*q^(-1-h)*F; "
        "u2=y*v0/(2q)+x*q^(-1-h)*F; u3=q^(-1/2-h)*U"
    ),
    "axis_U": "U(0,eta)=4*eta+j0",
    "pressure_radial_identity": "Pi_X=F^2",
    "radial_series": "finite nonlinear triangular X-series for source Eqs. (4.9)/(4.13)",
    "core_domain": "Lambda*X<=4.1",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "native_similarity_coordinates": True,
    "finite_nonlinear_core_series_executable": True,
    "velocity_api_compatible": True,
    "core_domain_only": True,
    "outer_heat_profile_matched": False,
    "global_leading_profile_reconstructed": False,
    "complete_kokuno_composite_velocity": False,
    "independent_full_momentum_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_array(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _bounded_integer(value: Any, name: str, lower: int, upper: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer")
    result = int(value)
    if result < lower or result > upper:
        raise ValueError(f"{name} must lie in [{lower}, {upper}]")
    return result


@dataclass(frozen=True)
class KokunoLeadingCoreSeriesCandidate:
    """Replayable finite nonlinear Kokuno leading-core candidate.

    The defaults use the already-screened stable inner-series choice
    ``sigma=.5, degree=14, eta_nodes=257``.  Every finite parameter is an
    autonomous reconstruction choice within the source/repository contracts;
    none is presented as a recovered hidden OpenAI parameter.
    """

    h: float = 0.005
    j0: float = 0.02
    sigma: float = 0.5
    Lambda: float = 10.0
    C: float = 2.0
    pressure_scale: float = 1.0
    maxdegree: int = 14
    eta_nodes: int = 257
    quadrature_points: int = 32

    _reference: PaperCoreReference = field(init=False, repr=False, compare=False)
    _series: PaperCoreSeries = field(init=False, repr=False, compare=False)
    _coordinates: KokunoNativeSimilarityCoordinates = field(
        init=False, repr=False, compare=False
    )
    _axis_profile: KokunoLeadingAxisProfile = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        maxdegree = _bounded_integer(self.maxdegree, "maxdegree", 0, 16)
        eta_nodes = _bounded_integer(self.eta_nodes, "eta_nodes", 3, 513)
        quadrature_points = _bounded_integer(
            self.quadrature_points, "quadrature_points", 8, 128
        )
        reference = PaperCoreReference(
            h=float(self.h),
            j=float(self.j0),
            sigma=float(self.sigma),
            Lambda=float(self.Lambda),
            C=float(self.C),
            pressure_scale=float(self.pressure_scale),
        )
        coordinates = KokunoNativeSimilarityCoordinates(h=reference.h)
        axis_profile = KokunoLeadingAxisProfile(h=reference.h, j0=reference.j)
        series = PaperCoreSeries(
            reference=reference, maxdegree=maxdegree, eta_nodes=eta_nodes
        )
        object.__setattr__(self, "h", reference.h)
        object.__setattr__(self, "j0", reference.j)
        object.__setattr__(self, "sigma", reference.sigma)
        object.__setattr__(self, "Lambda", reference.Lambda)
        object.__setattr__(self, "C", reference.C)
        object.__setattr__(self, "pressure_scale", reference.pressure_scale)
        object.__setattr__(self, "maxdegree", maxdegree)
        object.__setattr__(self, "eta_nodes", eta_nodes)
        object.__setattr__(self, "quadrature_points", quadrature_points)
        object.__setattr__(self, "_reference", reference)
        object.__setattr__(self, "_coordinates", coordinates)
        object.__setattr__(self, "_axis_profile", axis_profile)
        object.__setattr__(self, "_series", series)

    @property
    def reference(self) -> PaperCoreReference:
        return self._reference

    @property
    def series(self) -> PaperCoreSeries:
        return self._series

    @property
    def coordinates_model(self) -> KokunoNativeSimilarityCoordinates:
        return self._coordinates

    @property
    def axis_profile(self) -> KokunoLeadingAxisProfile:
        return self._axis_profile

    @property
    def max_core_X(self) -> float:
        return CORE_Y_MAX / self.Lambda

    def coordinates(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        values = self._coordinates.evaluate(x, y, z, t)
        if np.any(self.Lambda * values["X"] > CORE_Y_MAX + 2.0e-12):
            raise ValueError("Kokuno leading core candidate requires Lambda*X <= 4.1")
        return values

    def profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return low-dimensional profile values and their X/eta partials.

        X derivatives are analytic derivatives of the retained polynomial.
        Eta derivatives are the existing Chebyshev spectral derivatives used by
        ``PaperCoreSeries``; they are finite numerical approximations, not an AD
        proof of the infinite source fixed point.
        """

        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        return {
            "F": self._series.F(X_array, eta_array),
            "F_X": self._series.F_radial_derivative(X_array, eta_array),
            "F_eta": self._series.F_eta(X_array, eta_array),
            "U": self._series.U(X_array, eta_array),
            "U_X": self._series.U_radial_derivative(X_array, eta_array),
            "U_eta": self._series.dU_deta(X_array, eta_array),
            "Pi": self._series.Pi(X_array, eta_array),
            "Pi_X": self._series.Pi_radial_derivative(X_array, eta_array),
            "Pi_eta": self._series.Pi_eta(X_array, eta_array),
        }

    @staticmethod
    def _time_array(t: Any) -> np.ndarray:
        values = _finite_array(t, "t")
        if np.any((values < 0.0) | (values >= 1.0)):
            raise ValueError("time must satisfy 0 <= t < 1")
        return values

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Evaluate Cartesian ``[...,3]`` velocity on the source-native core."""

        x_array, y_array, z_array, t_array = np.broadcast_arrays(
            _finite_array(x, "x"),
            _finite_array(y, "y"),
            _finite_array(z, "z"),
            self._time_array(t),
        )
        self.coordinates(x_array, y_array, z_array, t_array)
        shape = x_array.shape
        result = np.empty(shape + (3,), dtype=float)
        flat_result = result.reshape((-1, 3))
        for index, (xx, yy, zz, tt) in enumerate(
            zip(
                x_array.reshape(-1),
                y_array.reshape(-1),
                z_array.reshape(-1),
                t_array.reshape(-1),
                strict=True,
            )
        ):
            flat_result[index] = self._series.velocity(
                float(xx),
                float(yy),
                float(zz),
                float(tt),
                quadrature_points=self.quadrature_points,
            )
        if not np.all(np.isfinite(result)):
            raise RuntimeError("Kokuno leading core velocity became non-finite")
        return result

    __call__ = velocity

    def at_points(self, points_xyz: Any, t: Any) -> np.ndarray:
        points = _finite_array(points_xyz, "points_xyz")
        if points.ndim < 2 or points.shape[-1] != 3:
            raise ValueError("points_xyz must have shape (...,3)")
        return self.velocity(points[..., 0], points[..., 1], points[..., 2], t)

    def pressure(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        x_array, y_array, z_array, t_array = np.broadcast_arrays(
            _finite_array(x, "x"),
            _finite_array(y, "y"),
            _finite_array(z, "z"),
            self._time_array(t),
        )
        self.coordinates(x_array, y_array, z_array, t_array)
        result = np.empty(x_array.shape, dtype=float)
        flat = result.reshape(-1)
        for index, (xx, yy, zz, tt) in enumerate(
            zip(
                x_array.reshape(-1),
                y_array.reshape(-1),
                z_array.reshape(-1),
                t_array.reshape(-1),
                strict=True,
            )
        ):
            flat[index] = self._series.pressure(
                float(xx), float(yy), float(zz), float(tt)
            )
        return result

    def grid(
        self,
        x: Any,
        y: Any,
        z: Any,
        times: Any,
    ) -> np.ndarray:
        """Return ``(time,x,y,z,component)`` samples without extrapolation."""

        x_axis = _finite_array(x, "x")
        y_axis = _finite_array(y, "y")
        z_axis = _finite_array(z, "z")
        time_axis = self._time_array(times)
        for name, axis in (("x", x_axis), ("y", y_axis), ("z", z_axis), ("times", time_axis)):
            if axis.ndim != 1 or axis.size == 0:
                raise ValueError(f"{name} must be a nonempty one-dimensional axis")
        Xg, Yg, Zg = np.meshgrid(x_axis, y_axis, z_axis, indexing="ij")
        output = np.empty(
            (time_axis.size, x_axis.size, y_axis.size, z_axis.size, 3), dtype=float
        )
        for index, time in enumerate(time_axis):
            output[index] = self.velocity(Xg, Yg, Zg, float(time))
        return output

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "corrected_pdf": CORRECTED_PDF,
                "formula_evidence": "pinned_public_workbench_commit",
                "radial_implementation": "reuse_existing_PaperCoreSeries_no_reimplementation",
            },
            "parameters": {
                "h": self.h,
                "j0": self.j0,
                "sigma": self.sigma,
                "Lambda": self.Lambda,
                "C": self.C,
                "pressure_scale": self.pressure_scale,
                "maxdegree": self.maxdegree,
                "eta_nodes": self.eta_nodes,
                "quadrature_points": self.quadrature_points,
                "origin": "autonomous_finite_reconstruction_choices_not_hidden_parameters",
                "stable_inner_screen": "sigma=.5,maxdegree=14,eta_nodes=257",
            },
            "formula_map": dict(_SOURCE_MAP),
            "dependencies": {
                "native_coordinates_sha256": self._coordinates.sha256,
                "leading_axis_profile_sha256": self._axis_profile.sha256,
                "paper_core_series_metadata": self._series.metadata(),
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoLeadingCoreSeriesCandidate":
        if not isinstance(payload, dict):
            raise ValueError("candidate payload must be a JSON object")
        expected = {"schema", "source", "parameters", "formula_map", "dependencies", "truth_boundary"}
        if set(payload) - {"sha256"} != expected:
            raise ValueError("candidate payload schema keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported Kokuno leading-core schema")
        if payload["formula_map"] != _SOURCE_MAP:
            raise ValueError("Kokuno source formula map changed")
        if payload["truth_boundary"] != _TRUTH_BOUNDARY:
            raise ValueError("Kokuno truth-boundary metadata changed")

        source = payload["source"]
        expected_source = cls().to_payload()["source"]
        if source != expected_source:
            raise ValueError("Kokuno source/provenance metadata changed")
        parameters = payload["parameters"]
        expected_parameter_keys = {
            "h", "j0", "sigma", "Lambda", "C", "pressure_scale", "maxdegree",
            "eta_nodes", "quadrature_points", "origin", "stable_inner_screen"
        }
        if set(parameters) != expected_parameter_keys:
            raise ValueError("Kokuno candidate parameter metadata changed")
        if parameters["origin"] != "autonomous_finite_reconstruction_choices_not_hidden_parameters":
            raise ValueError("Kokuno candidate parameter provenance changed")
        if parameters["stable_inner_screen"] != "sigma=.5,maxdegree=14,eta_nodes=257":
            raise ValueError("Kokuno stable-inner-screen provenance changed")

        candidate = cls(
            h=float(parameters["h"]),
            j0=float(parameters["j0"]),
            sigma=float(parameters["sigma"]),
            Lambda=float(parameters["Lambda"]),
            C=float(parameters["C"]),
            pressure_scale=float(parameters["pressure_scale"]),
            maxdegree=int(parameters["maxdegree"]),
            eta_nodes=int(parameters["eta_nodes"]),
            quadrature_points=int(parameters["quadrature_points"]),
        )
        if payload["dependencies"] != candidate.to_payload()["dependencies"]:
            raise ValueError("Kokuno leading-core dependency identity changed")
        if "sha256" in payload and payload["sha256"] != candidate.sha256:
            raise ValueError("Kokuno leading-core candidate SHA mismatch")
        return candidate

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoLeadingCoreSeriesCandidate":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


__all__ = ["KokunoLeadingCoreSeriesCandidate"]
