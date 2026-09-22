"""Same-identity auxiliary-T2 lift of the frozen current-I4 A2 oscillation.

This Kokuno Agent-2 adapter closes one interface seam requested by Agent 3's
RF30 auxiliary-Haar bridge.  It does *not* invent a new velocity field.  The
provider reuses the exact sign-resolved public-z complete-curl runtime already
bound into A2 PR #1080 and lifts its two signed covariance contributions onto an
explicit auxiliary two-torus by independent angular phase translations.

For one fixed-Q current-I4 preflight, source-chart coordinates are mapped back to
the physical coordinates of the frozen candidate by

    r = Q^(1/2) R,   z = Q^D Z,   t = 1 - Q T.

The sigma=+ and sigma=- sign-resolved complete-curl contributions are evaluated
at independent angles theta_+=2*pi*y1 and theta_-=2*pi*y2, rotated back to
cylindrical components, summed over beta, and normalized by Q^A.  On the torus
diagonal y1=y2 this exactly replays the frozen candidate's cylindrical
oscillatory field at the corresponding physical angle.

The *existence of an auxiliary T^2 Haar variable* is source structure from the
corrected 2026-09-09 Kokuno reconstruction.  The present identification of y1/y2
with independent translations of the two repository sign channels is an
explicit repository-autonomous candidate lift; it is not recovered source-exact
auxiliary mode data.  Slow R/Z derivatives are recomputed from this raw lifted
wave with a fixed centered FD4 stencil.  They are numerical approximations, not
paper-exact derivatives.  No pre-averaged covariance, RF30 defect, forcing,
pressure, held-out residual, or correction coefficient enters this provider.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_COMPOSITE_HEAD,
    AGENT2_COMPOSITE_PR,
    AGENT2_COMPOSITE_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
)
from .kokuno_current_i4_rf30_fixedq_preflight import (
    CurrentI4FixedQObservablePreflight,
    materialize_current_i4_rf30_fixedq_preflight,
)
from .kokuno_public_z_pullback_velocity import default_field as _default_oscillatory_field
from .kokuno_rf30_auxiliary_t2_haar_bridge import (
    AuxiliaryT2ProviderMetadata,
    AuxiliaryT2WaveSamples,
)

TASK = "K2-OSC-111"
SCHEMA = "kokuno-a2-current-i4-auxiliary-t2-provider-v1"
PARENT_AGENT3_PR = 1190
PARENT_AGENT3_HEAD = "21b0ade966fda3b7344273c7bff2efe7e7b9a93e"
PARENT_AGENT2_PR = 1080
PARENT_AGENT2_HEAD = "c40d8ddecd2971544a6e07dab093436b423cf326"
PUBLIC_Z_SOURCE_BLOB = "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"

# Fixed numerical derivative realization.  This is not a scientific tolerance
# and is deliberately not caller configurable.
FD4_RELATIVE_STEP = 2.0 ** -12
CHART_REPLAY_RTOL = 2.0e-12


class CurrentI4AuxiliaryT2ProviderError(RuntimeError):
    """Raised when the exact candidate/provider identity contract is violated."""


def _canonical_json(payload: Any) -> str:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def _sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _git_blob_sha1(path: str | Path) -> str:
    data = Path(path).read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _finite_scalar(value: Any, label: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise CurrentI4AuxiliaryT2ProviderError(f"{label} must be finite")
    return out


def _finite_vector(value: Any, label: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim != 1 or out.size < 1 or not np.all(np.isfinite(out)):
        raise CurrentI4AuxiliaryT2ProviderError(f"{label} must be one finite vector")
    return out


def _runtime_payload_sha(field: Any) -> str:
    payload = field.to_payload()
    if not isinstance(payload, dict):
        raise CurrentI4AuxiliaryT2ProviderError("oscillatory runtime payload is not a mapping")
    return _sha256(payload)


def _runtime_source_blob(field: Any) -> str:
    path = inspect.getsourcefile(type(field))
    if path is None:
        raise CurrentI4AuxiliaryT2ProviderError("oscillatory runtime source is unavailable")
    return _git_blob_sha1(path)


def _provider_source_blob() -> str:
    path = inspect.getsourcefile(KokunoCurrentI4AuxiliaryT2Provider)
    if path is None:
        raise CurrentI4AuxiliaryT2ProviderError("provider source is unavailable")
    return _git_blob_sha1(path)


@dataclass(frozen=True)
class _BoundChart:
    Q: float
    A: float
    D: float
    Z: float
    T: float
    radius_R: tuple[float, ...]


class KokunoCurrentI4AuxiliaryT2Provider:
    """Raw auxiliary-T2 wave provider bound to one exact current-I4 preflight."""

    def __init__(
        self,
        backend: ExactCurrentI4NonlinearBackend,
        preflight: CurrentI4FixedQObservablePreflight,
    ) -> None:
        if not isinstance(backend, ExactCurrentI4NonlinearBackend):
            raise TypeError("backend must be ExactCurrentI4NonlinearBackend")
        if not isinstance(preflight, CurrentI4FixedQObservablePreflight):
            raise TypeError("preflight must be CurrentI4FixedQObservablePreflight")
        if backend.composite_source_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
            raise CurrentI4AuxiliaryT2ProviderError("A2 #1080 composite source blob drifted")
        if preflight.candidate_semantic_sha256 != backend.composite_semantic_sha256:
            raise CurrentI4AuxiliaryT2ProviderError("candidate semantic identity drifted")
        if preflight.oscillatory_runtime_sha256 != backend.oscillatory_runtime_sha256:
            raise CurrentI4AuxiliaryT2ProviderError("oscillatory runtime identity drifted")
        leading = getattr(backend.composite_field, "leading_backend", None)
        leading_semantic = str(getattr(leading, "semantic_sha256", ""))
        if preflight.leading_semantic_sha256 != leading_semantic:
            raise CurrentI4AuxiliaryT2ProviderError("leading semantic identity drifted")
        if not preflight.all_strict_I4:
            raise CurrentI4AuxiliaryT2ProviderError("provider requires one strict current-I4 preflight")

        Q = _finite_scalar(preflight.Q, "Q")
        A = _finite_scalar(preflight.A, "A")
        D = _finite_scalar(preflight.D, "D")
        Z = _finite_scalar(preflight.Z, "Z")
        T = _finite_scalar(preflight.T, "T")
        radius_R = tuple(float(v) for v in _finite_vector(preflight.radius_R, "radius_R"))
        if not (0.0 < Q <= 1.0 and A > 0.0 and D > 0.0):
            raise CurrentI4AuxiliaryT2ProviderError("fixed-Q chart exponents are invalid")

        field = _default_oscillatory_field()
        if _runtime_source_blob(field) != PUBLIC_Z_SOURCE_BLOB:
            raise CurrentI4AuxiliaryT2ProviderError("public-z oscillatory runtime source blob drifted")
        runtime_sha = _runtime_payload_sha(field)
        if runtime_sha != backend.oscillatory_runtime_sha256:
            raise CurrentI4AuxiliaryT2ProviderError(
                "public-z runtime payload is not the oscillatory identity bound by A2 #1080"
            )

        self._backend = backend
        self._preflight = preflight
        self._field = field
        self._chart = _BoundChart(Q=Q, A=A, D=D, Z=Z, T=T, radius_R=radius_R)
        self._runtime_sha = runtime_sha

    @property
    def chart(self) -> _BoundChart:
        return self._chart

    def semantic_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent3": {"pr": PARENT_AGENT3_PR, "head": PARENT_AGENT3_HEAD},
            "exact_current_i4_candidate": {
                "agent2_pr": AGENT2_COMPOSITE_PR,
                "agent2_head": AGENT2_COMPOSITE_HEAD,
                "candidate_semantic_sha256": self._preflight.candidate_semantic_sha256,
                "oscillatory_runtime_sha256": self._preflight.oscillatory_runtime_sha256,
                "leading_semantic_sha256": self._preflight.leading_semantic_sha256,
                "public_z_source_blob_sha1": PUBLIC_Z_SOURCE_BLOB,
            },
            "fixed_q_chart": {
                "Q": self.chart.Q,
                "A": self.chart.A,
                "D": self.chart.D,
                "Z": self.chart.Z,
                "T": self.chart.T,
                "radius_R": list(self.chart.radius_R),
                "map_to_physical": "r=sqrt(Q)R; z=Q^D Z; t=1-QT",
            },
            "auxiliary_t2_lift": {
                "coordinates": "y1,y2 in [0,1)",
                "sigma_plus_angle": "theta_plus=2*pi*y1",
                "sigma_minus_angle": "theta_minus=2*pi*y2",
                "diagonal_replay": "y1=y2 replays the actual sign-summed candidate at that physical angle",
                "classification": "repository_autonomous_candidate_lift_of_actual_sign_resolved_complete_curl",
                "source_exact_claimed": False,
            },
            "derivatives": {
                "operator": "centered_fd4_on_raw_lifted_wave_in_source_R_Z",
                "relative_step": FD4_RELATIVE_STEP,
                "caller_configurable": False,
                "source_exact_claimed": False,
            },
            "source_provenance": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_READER_HEAD,
                "date": SOURCE_READER_DATE,
                "path": SOURCE_READER_PATH,
                "blob_sha1": SOURCE_READER_BLOB,
                "scope": "auxiliary-T2 Haar structure plus localized complete-curl organization",
            },
            "truth_boundary": truth_boundary(),
        }

    @property
    def semantic_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    @property
    def metadata(self) -> AuxiliaryT2ProviderMetadata:
        return AuxiliaryT2ProviderMetadata(
            candidate_semantic_sha256=self._preflight.candidate_semantic_sha256,
            oscillatory_runtime_sha256=self._preflight.oscillatory_runtime_sha256,
            leading_semantic_sha256=self._preflight.leading_semantic_sha256,
            provider_semantic_sha256=self.semantic_sha256,
            source_blob_sha1=_provider_source_blob(),
            provider_kind="repository_candidate",
            source_auxiliary_t2_field_materialized=True,
            normalized_haar_measure_total_mass_one=True,
            recomputed_from_actual_candidate=True,
            surrogate_or_preaveraged_covariance_used=False,
            heldout_data_used=False,
            residual_as_forcing_used=False,
        )

    def _validate_bound_chart(self, R: Any, Z: Any, T: Any) -> np.ndarray:
        radius = _finite_vector(R, "R")
        expected = np.asarray(self.chart.radius_R, dtype=float)
        if radius.shape != expected.shape or not np.allclose(
            radius, expected, rtol=CHART_REPLAY_RTOL, atol=0.0
        ):
            raise CurrentI4AuxiliaryT2ProviderError("R drifted from the bound fixed-Q preflight")
        z = _finite_scalar(Z, "Z")
        t = _finite_scalar(T, "T")
        if not math.isclose(z, self.chart.Z, rel_tol=CHART_REPLAY_RTOL, abs_tol=0.0):
            raise CurrentI4AuxiliaryT2ProviderError("Z drifted from the bound fixed-Q preflight")
        if not math.isclose(t, self.chart.T, rel_tol=CHART_REPLAY_RTOL, abs_tol=0.0):
            raise CurrentI4AuxiliaryT2ProviderError("T drifted from the bound fixed-Q preflight")
        return radius

    @staticmethod
    def _torus_grid(y1: Any, y2: Any) -> tuple[np.ndarray, np.ndarray]:
        a = np.asarray(y1, dtype=float)
        b = np.asarray(y2, dtype=float)
        if a.ndim != 2 or b.ndim != 2 or a.shape != b.shape or a.shape[0] != a.shape[1]:
            raise CurrentI4AuxiliaryT2ProviderError("y1,y2 must be one matching square torus grid")
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            raise CurrentI4AuxiliaryT2ProviderError("auxiliary torus coordinates must be finite")
        if np.any((a < 0.0) | (a >= 1.0)) or np.any((b < 0.0) | (b >= 1.0)):
            raise CurrentI4AuxiliaryT2ProviderError("auxiliary torus coordinates must lie in [0,1)")
        return a, b

    def _sign_cylindrical(
        self,
        R: np.ndarray,
        Z: float,
        T: float,
        theta: np.ndarray,
        sign_index: int,
    ) -> np.ndarray:
        q = self.chart.Q
        r = math.sqrt(q) * np.asarray(R, dtype=float)
        z_phys = (q ** self.chart.D) * float(Z)
        t_phys = 1.0 - q * float(T)
        theta = np.asarray(theta, dtype=float)
        rr = r[:, None, None]
        tt = theta[None, :, :]
        x = rr * np.cos(tt)
        y = rr * np.sin(tt)
        z = np.full_like(x, z_phys)
        time = np.full_like(x, t_phys)
        out = self._field.evaluate(x, y, z, time)
        values = np.asarray(out["velocity_cartesian_by_beta_sign"], dtype=float)
        expected_prefix = x.shape
        if values.shape[:3] != expected_prefix or values.shape[-2:] != (2, 3):
            raise CurrentI4AuxiliaryT2ProviderError("sign-resolved runtime returned an unexpected shape")
        selected = values[..., sign_index, :]
        c = np.cos(tt)[..., None]
        s = np.sin(tt)[..., None]
        ur_beta = selected[..., 0] * c - 0.0 + selected[..., 1] * s
        ut_beta = -selected[..., 0] * s + selected[..., 1] * c
        uz_beta = selected[..., 2]
        return np.stack(
            (
                np.sum(ur_beta, axis=-1),
                np.sum(ut_beta, axis=-1),
                np.sum(uz_beta, axis=-1),
            ),
            axis=-1,
        )

    def _raw_wave(
        self,
        R: np.ndarray,
        Z: float,
        T: float,
        y1: np.ndarray,
        y2: np.ndarray,
    ) -> np.ndarray:
        plus = self._sign_cylindrical(R, Z, T, 2.0 * math.pi * y1, 0)
        minus = self._sign_cylindrical(R, Z, T, 2.0 * math.pi * y2, 1)
        wave = (self.chart.Q ** self.chart.A) * (plus + minus)
        if not np.all(np.isfinite(wave)):
            raise CurrentI4AuxiliaryT2ProviderError("raw auxiliary-T2 wave became non-finite")
        return wave

    def _fd4_derivatives(
        self,
        R: np.ndarray,
        Z: float,
        T: float,
        y1: np.ndarray,
        y2: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        hR = FD4_RELATIVE_STEP * np.maximum(1.0, np.abs(R))
        if np.any(R - 2.0 * hR <= 0.0):
            raise CurrentI4AuxiliaryT2ProviderError("FD4 radial stencil would cross R=0")
        rp2 = self._raw_wave(R + 2.0 * hR, Z, T, y1, y2)
        rp1 = self._raw_wave(R + hR, Z, T, y1, y2)
        rm1 = self._raw_wave(R - hR, Z, T, y1, y2)
        rm2 = self._raw_wave(R - 2.0 * hR, Z, T, y1, y2)
        dR = (-rp2 + 8.0 * rp1 - 8.0 * rm1 + rm2) / (
            12.0 * hR[:, None, None, None]
        )

        hZ = FD4_RELATIVE_STEP * max(1.0, abs(float(Z)))
        zp2 = self._raw_wave(R, Z + 2.0 * hZ, T, y1, y2)
        zp1 = self._raw_wave(R, Z + hZ, T, y1, y2)
        zm1 = self._raw_wave(R, Z - hZ, T, y1, y2)
        zm2 = self._raw_wave(R, Z - 2.0 * hZ, T, y1, y2)
        dZ = (-zp2 + 8.0 * zp1 - 8.0 * zm1 + zm2) / (12.0 * hZ)
        if not np.all(np.isfinite(dR)) or not np.all(np.isfinite(dZ)):
            raise CurrentI4AuxiliaryT2ProviderError("raw auxiliary-T2 derivatives became non-finite")
        return dR, dZ

    def sample_auxiliary_wave(
        self,
        *,
        R: Any,
        Z: float,
        T: float,
        y1: Any,
        y2: Any,
    ) -> AuxiliaryT2WaveSamples:
        radius = self._validate_bound_chart(R, Z, T)
        yy1, yy2 = self._torus_grid(y1, y2)
        wave = self._raw_wave(radius, float(Z), float(T), yy1, yy2)
        dR, dZ = self._fd4_derivatives(radius, float(Z), float(T), yy1, yy2)
        return AuxiliaryT2WaveSamples(
            w_r=wave[..., 0],
            w_theta=wave[..., 1],
            w_z=wave[..., 2],
            dR_w_r=dR[..., 0],
            dR_w_theta=dR[..., 1],
            dR_w_z=dR[..., 2],
            dZ_w_r=dZ[..., 0],
            dZ_w_theta=dZ[..., 1],
            dZ_w_z=dZ[..., 2],
        )


def bind_current_i4_auxiliary_t2_provider(
    backend: ExactCurrentI4NonlinearBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> KokunoCurrentI4AuxiliaryT2Provider:
    """Bind one provider to the exact fixed-Q preflight used by Agent 3."""
    preflight = materialize_current_i4_rf30_fixedq_preflight(backend, radius, z, t)
    return KokunoCurrentI4AuxiliaryT2Provider(backend, preflight)


def truth_boundary() -> dict[str, Any]:
    params = tuple(inspect.signature(bind_current_i4_auxiliary_t2_provider).parameters)
    forbidden = {
        "residual", "defect", "forcing", "pressure", "target", "gain", "damping",
        "threshold", "viscosity", "nu", "covariance", "correction", "fd_step",
        "phase_gain", "amplitude",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_composite_pr": PARENT_AGENT2_PR,
        "agent2_composite_head": PARENT_AGENT2_HEAD,
        "public_parameters": params,
        "forbidden_scientific_controls_exposed": bool(forbidden.intersection(params)),
        "same_identity_current_i4_candidate_bound": True,
        "actual_sign_resolved_complete_curl_runtime_reused": True,
        "raw_auxiliary_t2_wave_materialized": True,
        "normalized_haar_domain_total_mass_one": True,
        "preaveraged_covariance_supplied": False,
        "physical_theta_mean_relabelled_as_auxiliary_haar": False,
        "auxiliary_sign_phase_lift_repository_autonomous": True,
        "source_exact_auxiliary_t2_mode_family_recovered": False,
        "source_exact_auxiliary_phase_assignment_recovered": False,
        "slow_R_Z_derivatives_recomputed_from_raw_wave": True,
        "slow_derivative_realization": "centered_fd4_repository_fixed",
        "source_exact_slow_derivatives_claimed": False,
        "agent3_rf30_defect_materialized_here": False,
        "agent3_mean_correction_materialized_here": False,
        "forcing_or_pressure_added": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


__all__ = [
    "CurrentI4AuxiliaryT2ProviderError",
    "FD4_RELATIVE_STEP",
    "KokunoCurrentI4AuxiliaryT2Provider",
    "bind_current_i4_auxiliary_t2_provider",
    "truth_boundary",
]
