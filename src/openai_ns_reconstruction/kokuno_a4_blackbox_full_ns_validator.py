"""Implementation-distinct black-box full-NS validator for Kokuno candidates.

This Agent-4 module is the fail-closed bridge from the current strict-inner
transport precursor to the eventual complete candidate.  It consumes only
candidate-facing public surfaces

    velocity(x,y,z,t) -> (...,3)
    pressure(x,y,z,t) -> (...,)
    forcing(x,y,z,t) -> (...,3)

plus a fixed physical-contract identity.  It never consumes training loss,
construction tensors, Agent-1 analytic derivatives, Agent-2 production
Jacobians/Laplacians, or an Agent-3 correction defect.

The independent operator is a Cartesian fourth-order five-point stencil, with
endpoint-aware fourth-order time differentiation.  This is deliberately
different from the Richardson transport oracle in parent Agent-4 PR #889.

The current strict-inner Kokuno artifact is *not* eligible for this validator:
it still lacks the global/outer leading join, matched pressure, preregistered
restricted forcing, and a real correction velocity.  The module therefore
provides a complete-residual evaluator and admission contract without
inventing any missing scientific surface.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from typing import Any, Mapping, Protocol

import numpy as np

TASK = "KOKUNO-A4-BLACKBOX-FULL-NS-VALIDATOR-081"
SCHEMA = "kokuno-a4-blackbox-full-ns-validator-v1"
PARENT_A4_PR = 889
PARENT_A4_HEAD = "574fd5f7d67e5fbf905ef91a483a7e9e8124698d"
LATEST_A1_PR = 893
LATEST_A1_HEAD = "d89cc4c21c1f08a124776e38495e46ce09e36a78"
LATEST_A2_PR = 894
LATEST_A2_HEAD = "a467bc175a801ed69c197fc9e38ca68ac8c0b86b"
LATEST_A3_PR = 895
LATEST_A3_HEAD = "d305b64edcdde4f16dd8a48908f53fa23c2b9baf"

NU = 0.01
EVALUATION_BOX = ((-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0))
TIME_INTERVAL = (0.25, 0.75)
VALIDATION_TIMES = (0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75)
HELDOUT_SEED = 914027
HELDOUT_POINTS = 4096
DERIVATIVE_STEPS = (2.0e-2, 1.0e-2, 5.0e-3)
RESIDUAL_SCALE = 1.0
MOMENTUM_GATE = 1.0e-3
DIVERGENCE_GATE = 1.0e-5
BOUNDARY_VELOCITY_GATE = 1.0e-10
BOUNDARY_PRESSURE_GATE = 1.0e-10
ENERGY_REFERENCE_TIME = 0.25
ENERGY_REFERENCE = 1.0
ENERGY_REFERENCE_TOLERANCE = 1.0e-3
ENERGY_MIN = 0.1
ENERGY_MAX = 10.0
RESOLUTION_RELATIVE_CHANGE_GATE = 5.0e-2
RESOLUTION_RATIO_GATE = 1.5
RESOLUTION_FLOOR = 5.0e-8
NONTRIVIAL_VELOCITY_RMS_FLOOR = 1.0e-8
NU_PERTURBATION_FRACTION = 1.0e-3


class CandidateAdmissionError(RuntimeError):
    """Raised when a candidate does not expose the frozen complete-NS contract."""


class BlackBoxFullNSCandidate(Protocol):
    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...
    def pressure(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...
    def forcing(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...
    def validation_metadata(self) -> Mapping[str, Any]: ...


def _sha256_payload(payload: Mapping[str, Any]) -> str:
    raw = json.dumps(
        dict(payload), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


@dataclass(frozen=True)
class FixedPhysicalContract:
    contract_id: str
    viscosity: float
    physical_domain: str
    evaluation_box: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    time_interval: tuple[float, float]
    support_condition: str
    residual_scale: float
    matched_pressure_identity_sha256: str
    restricted_forcing_family_id: str
    restricted_forcing_parameters_sha256: str
    restricted_forcing_preregistration_sha256: str
    physical_contract_sha256: str

    def identity_payload(self) -> dict[str, Any]:
        return {
            "contract_id": self.contract_id,
            "viscosity": float(self.viscosity),
            "physical_domain": self.physical_domain,
            "evaluation_box": [list(v) for v in self.evaluation_box],
            "time_interval": list(self.time_interval),
            "support_condition": self.support_condition,
            "residual_scale": float(self.residual_scale),
            "matched_pressure_identity_sha256": self.matched_pressure_identity_sha256,
            "restricted_forcing_family_id": self.restricted_forcing_family_id,
            "restricted_forcing_parameters_sha256": self.restricted_forcing_parameters_sha256,
            "restricted_forcing_preregistration_sha256": self.restricted_forcing_preregistration_sha256,
        }

    def __post_init__(self) -> None:
        if float(self.viscosity) != NU:
            raise ValueError("viscosity drifted from frozen CR001 value")
        if self.evaluation_box != EVALUATION_BOX:
            raise ValueError("evaluation box drifted from frozen CR001 value")
        if self.time_interval != TIME_INTERVAL:
            raise ValueError("time interval drifted from frozen CR001 value")
        if float(self.residual_scale) != RESIDUAL_SCALE:
            raise ValueError("residual normalization scale drifted from frozen value")
        text = (
            self.contract_id,
            self.physical_domain,
            self.support_condition,
            self.matched_pressure_identity_sha256,
            self.restricted_forcing_family_id,
            self.restricted_forcing_parameters_sha256,
            self.restricted_forcing_preregistration_sha256,
            self.physical_contract_sha256,
        )
        if any(not str(v) for v in text):
            raise ValueError("physical-contract identity fields must be nonempty")
        expected = _sha256_payload(self.identity_payload())
        if self.physical_contract_sha256 != expected:
            raise ValueError("physical_contract_sha256 does not match structured identity")


def build_fixed_physical_contract(
    *,
    contract_id: str,
    matched_pressure_identity_sha256: str,
    restricted_forcing_parameters_sha256: str,
    restricted_forcing_preregistration_sha256: str,
    restricted_forcing_family_id: str = "restricted_two_parameter_family",
) -> FixedPhysicalContract:
    payload = {
        "contract_id": contract_id,
        "viscosity": NU,
        "physical_domain": "R^3",
        "evaluation_box": [list(v) for v in EVALUATION_BOX],
        "time_interval": list(TIME_INTERVAL),
        "support_condition": "r < 2 and abs(z) < 2; smooth zero extension",
        "residual_scale": RESIDUAL_SCALE,
        "matched_pressure_identity_sha256": matched_pressure_identity_sha256,
        "restricted_forcing_family_id": restricted_forcing_family_id,
        "restricted_forcing_parameters_sha256": restricted_forcing_parameters_sha256,
        "restricted_forcing_preregistration_sha256": restricted_forcing_preregistration_sha256,
    }
    return FixedPhysicalContract(
        contract_id=contract_id,
        viscosity=NU,
        physical_domain="R^3",
        evaluation_box=EVALUATION_BOX,
        time_interval=TIME_INTERVAL,
        support_condition="r < 2 and abs(z) < 2; smooth zero extension",
        residual_scale=RESIDUAL_SCALE,
        matched_pressure_identity_sha256=matched_pressure_identity_sha256,
        restricted_forcing_family_id=restricted_forcing_family_id,
        restricted_forcing_parameters_sha256=restricted_forcing_parameters_sha256,
        restricted_forcing_preregistration_sha256=restricted_forcing_preregistration_sha256,
        physical_contract_sha256=_sha256_payload(payload),
    )


@dataclass(frozen=True)
class FullNSLevel:
    derivative_step: float
    velocity: np.ndarray
    velocity_dt: np.ndarray
    jacobian: np.ndarray
    advection: np.ndarray
    pressure_gradient: np.ndarray
    laplacian: np.ndarray
    forcing: np.ndarray
    residual: np.ndarray
    divergence: np.ndarray
    normalized_sampled_max: float
    normalized_sampled_rms: float


@dataclass(frozen=True)
class BoundarySupportMetrics:
    velocity_max: float
    pressure_max: float
    passed: bool


@dataclass(frozen=True)
class EnergyMetrics:
    by_time: dict[str, float]
    reference_energy: float | None
    passed: bool
    reason: str | None


@dataclass(frozen=True)
class FullNSValidationReceipt:
    schema: str
    task: str
    candidate_id: str
    candidate_sha256: str
    stage: str
    physical_contract_sha256: str
    derivative_steps: tuple[float, float, float]
    momentum_sampled_max_by_level: tuple[float, float, float]
    momentum_sampled_rms_by_level: tuple[float, float, float]
    momentum_volume_l2_by_level: tuple[float | None, float | None, float | None]
    divergence_sampled_max_by_level: tuple[float, float, float]
    divergence_sampled_rms_by_level: tuple[float, float, float]
    divergence_volume_l2_by_level: tuple[float | None, float | None, float | None]
    resolution_relative_change: float
    resolution_ratio: float
    boundary_support: BoundarySupportMetrics
    energy: EnergyMetrics
    nu_minus_relative_rms_change: float
    nu_plus_relative_rms_change: float
    gate_failures: tuple[str, ...]
    heldout_normalized_ns_residual_assessed: bool
    quadrature_ladder_assessed: bool
    pde_validated: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _broadcast_xyzt(x: Any, y: Any, z: Any, t: Any) -> tuple[np.ndarray, ...]:
    arrays = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    if not all(np.all(np.isfinite(v)) for v in arrays):
        raise ValueError("x,y,z,t must be finite and broadcastable")
    return tuple(arrays)


def _vector_value(provider: Any, x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray, name: str) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise ValueError(f"{name} must return shape {expected}, got {value.shape}")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} returned non-finite values")
    return value


def _scalar_value(provider: Any, x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray, name: str) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    if value.shape != x.shape:
        raise ValueError(f"{name} must return shape {x.shape}, got {value.shape}")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} returned non-finite values")
    return value


def vector_rms(value: Any) -> float:
    arr = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def vector_sampled_max(value: Any) -> float:
    arr = np.asarray(value, dtype=float)
    return float(np.max(np.linalg.norm(arr, axis=-1)))


def _shift(coords: tuple[np.ndarray, ...], axis: int, offset: float) -> tuple[np.ndarray, ...]:
    out = [np.array(v, copy=True) for v in coords]
    out[axis] += float(offset)
    return tuple(out)


def _fd4_first_spatial(provider: Any, coords: tuple[np.ndarray, ...], axis: int, step: float, *, vector: bool) -> np.ndarray:
    if axis not in (0, 1, 2):
        raise ValueError("spatial axis must be 0,1,2")
    h = float(step)
    eval_fn = _vector_value if vector else _scalar_value
    name = "vector provider" if vector else "scalar provider"
    fm2 = eval_fn(provider, *_shift(coords, axis, -2.0 * h), name)
    fm1 = eval_fn(provider, *_shift(coords, axis, -h), name)
    fp1 = eval_fn(provider, *_shift(coords, axis, h), name)
    fp2 = eval_fn(provider, *_shift(coords, axis, 2.0 * h), name)
    return (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h)


def _fd4_second_spatial(provider: Any, coords: tuple[np.ndarray, ...], axis: int, step: float) -> np.ndarray:
    h = float(step)
    fm2 = _vector_value(provider, *_shift(coords, axis, -2.0 * h), "velocity")
    fm1 = _vector_value(provider, *_shift(coords, axis, -h), "velocity")
    f0 = _vector_value(provider, *coords, "velocity")
    fp1 = _vector_value(provider, *_shift(coords, axis, h), "velocity")
    fp2 = _vector_value(provider, *_shift(coords, axis, 2.0 * h), "velocity")
    return (-fp2 + 16.0 * fp1 - 30.0 * f0 + 16.0 * fm1 - fm2) / (12.0 * h * h)


def _fd4_time_first(provider: Any, coords: tuple[np.ndarray, ...], step: float) -> np.ndarray:
    h = float(step)
    t = coords[3]
    out = np.empty(t.shape + (3,), dtype=float)
    low = t - 2.0 * h < TIME_INTERVAL[0]
    high = t + 2.0 * h > TIME_INTERVAL[1]
    if np.any(low & high):
        raise ValueError("time step too large for frozen time interval")
    interior = ~(low | high)

    def evaluate(mask: np.ndarray, offsets: tuple[float, ...], coeffs: tuple[float, ...]) -> None:
        if not np.any(mask):
            return
        sub = tuple(v[mask] for v in coords)
        acc = np.zeros((int(np.sum(mask)), 3), dtype=float)
        for offset, coeff in zip(offsets, coeffs):
            shifted = _shift(sub, 3, offset * h)
            acc += coeff * _vector_value(provider, *shifted, "velocity")
        out[mask] = acc / (12.0 * h)

    evaluate(interior, (-2.0, -1.0, 1.0, 2.0), (1.0, -8.0, 8.0, -1.0))
    evaluate(low, (0.0, 1.0, 2.0, 3.0, 4.0), (-25.0, 48.0, -36.0, 16.0, -3.0))
    evaluate(high, (0.0, -1.0, -2.0, -3.0, -4.0), (25.0, -48.0, 36.0, -16.0, 3.0))
    return out


def full_ns_level(candidate: BlackBoxFullNSCandidate, x: Any, y: Any, z: Any, t: Any, *, derivative_step: float, nu: float = NU) -> FullNSLevel:
    h = float(derivative_step)
    if h not in DERIVATIVE_STEPS:
        raise ValueError("derivative step must be one of the frozen levels")
    if float(nu) != NU:
        raise ValueError("scientific evaluation viscosity is frozen at 0.01")
    coords = _broadcast_xyzt(x, y, z, t)
    velocity = _vector_value(candidate.velocity, *coords, "velocity")
    forcing = _vector_value(candidate.forcing, *coords, "forcing")
    velocity_dt = _fd4_time_first(candidate.velocity, coords, h)
    spatial = [_fd4_first_spatial(candidate.velocity, coords, axis, h, vector=True) for axis in range(3)]
    jacobian = np.stack(spatial, axis=-1)
    advection = np.einsum("...j,...ij->...i", velocity, jacobian)
    pressure_gradient = np.stack(
        [_fd4_first_spatial(candidate.pressure, coords, axis, h, vector=False) for axis in range(3)],
        axis=-1,
    )
    laplacian = sum(
        (_fd4_second_spatial(candidate.velocity, coords, axis, h) for axis in range(3)),
        start=np.zeros_like(velocity),
    )
    residual = velocity_dt + advection + pressure_gradient - NU * laplacian - forcing
    divergence = np.trace(jacobian, axis1=-2, axis2=-1)
    return FullNSLevel(
        derivative_step=h,
        velocity=velocity,
        velocity_dt=velocity_dt,
        jacobian=jacobian,
        advection=advection,
        pressure_gradient=pressure_gradient,
        laplacian=laplacian,
        forcing=forcing,
        residual=residual,
        divergence=divergence,
        normalized_sampled_max=vector_sampled_max(residual) / RESIDUAL_SCALE,
        normalized_sampled_rms=vector_rms(residual) / RESIDUAL_SCALE,
    )


def three_resolution_levels(candidate: BlackBoxFullNSCandidate, x: Any, y: Any, z: Any, t: Any) -> tuple[FullNSLevel, FullNSLevel, FullNSLevel]:
    return tuple(full_ns_level(candidate, x, y, z, t, derivative_step=h) for h in DERIVATIVE_STEPS)  # type: ignore[return-value]


def _volume_l2(vector_or_scalar: np.ndarray, weights: np.ndarray | None) -> float | None:
    if weights is None:
        return None
    w = np.asarray(weights, dtype=float)
    if np.any(~np.isfinite(w)) or np.any(w < 0.0):
        raise ValueError("volume weights must be finite and nonnegative")
    if vector_or_scalar.shape[: w.ndim] != w.shape:
        raise ValueError("volume weights shape must match sample shape")
    if vector_or_scalar.ndim == w.ndim + 1:
        mag2 = np.sum(vector_or_scalar * vector_or_scalar, axis=-1)
    else:
        mag2 = vector_or_scalar * vector_or_scalar
    return float(np.sqrt(np.sum(w * mag2)))


def _resolution_metrics(levels: tuple[FullNSLevel, FullNSLevel, FullNSLevel]) -> tuple[float, float]:
    coarse, medium, fine = levels
    c2m = vector_rms(coarse.residual - medium.residual)
    m2f = vector_rms(medium.residual - fine.residual)
    relative = m2f / max(vector_rms(fine.residual), 1.0e-30)
    ratio = c2m / max(m2f, 1.0e-300)
    return relative, ratio


def candidate_admission_failures(candidate: Any, contract: FixedPhysicalContract, *, expected_stage: str) -> list[str]:
    failures: list[str] = []
    for name in ("velocity", "pressure", "forcing", "validation_metadata"):
        if not callable(getattr(candidate, name, None)):
            failures.append(f"missing_public_surface:{name}")
    if failures:
        return failures
    metadata = dict(candidate.validation_metadata())
    required_text = ("candidate_id", "candidate_sha256", "physical_contract_sha256", "stage")
    for key in required_text:
        if not str(metadata.get(key, "")):
            failures.append(f"missing_metadata:{key}")
    if metadata.get("physical_contract_sha256") != contract.physical_contract_sha256:
        failures.append("physical_contract_identity_mismatch")
    if metadata.get("stage") != expected_stage:
        failures.append("stage_mismatch")
    if not bool(metadata.get("global_leading_velocity_materialized", False)):
        failures.append("global_leading_velocity_missing")
    if not bool(metadata.get("outer_join_materialized", False)):
        failures.append("outer_join_missing")
    if not bool(metadata.get("matched_pressure_included", False)):
        failures.append("matched_pressure_missing")
    if not bool(metadata.get("restricted_forcing_included", False)):
        failures.append("restricted_forcing_missing")
    if not bool(metadata.get("restricted_forcing_preregistered", False)):
        failures.append("restricted_forcing_not_preregistered")
    if bool(metadata.get("residual_defined_forcing", True)):
        failures.append("residual_defined_forcing_forbidden")
    if not bool(metadata.get("whole_domain_support_materialized", False)):
        failures.append("whole_domain_support_missing")
    if not bool(metadata.get("leading_included", False)):
        failures.append("leading_component_missing")
    if expected_stage == "leading_only":
        if bool(metadata.get("oscillatory_included", False)) or bool(metadata.get("correction_included", False)):
            failures.append("leading_only_stage_contaminated")
    elif expected_stage == "leading_plus_oscillatory":
        if not bool(metadata.get("oscillatory_included", False)):
            failures.append("oscillatory_component_missing")
        if bool(metadata.get("correction_included", False)):
            failures.append("pre_correction_stage_contaminated")
    elif expected_stage == "after_correction":
        if not bool(metadata.get("oscillatory_included", False)):
            failures.append("oscillatory_component_missing")
        if not bool(metadata.get("correction_included", False)):
            failures.append("correction_component_missing")
    else:
        failures.append("unsupported_stage")
    return failures


def frozen_boundary_support_metrics(candidate: BlackBoxFullNSCandidate) -> BoundarySupportMetrics:
    points = np.asarray(
        [
            (2.0, 0.0, 0.0), (-2.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, -2.0, 0.0),
            (0.0, 0.0, 2.0), (0.0, 0.0, -2.0), (2.0, 0.0, 2.0), (-2.0, 0.0, -2.0),
            (0.0, 2.0, 2.0), (0.0, -2.0, -2.0),
        ],
        dtype=float,
    )
    xyz = np.repeat(points, len(VALIDATION_TIMES), axis=0)
    times = np.tile(np.asarray(VALIDATION_TIMES, dtype=float), len(points))
    velocity = np.asarray(candidate.velocity(xyz[:, 0], xyz[:, 1], xyz[:, 2], times), dtype=float)
    pressure = np.asarray(candidate.pressure(xyz[:, 0], xyz[:, 1], xyz[:, 2], times), dtype=float)
    velocity_max = vector_sampled_max(velocity)
    pressure_max = float(np.max(np.abs(pressure)))
    return BoundarySupportMetrics(
        velocity_max=velocity_max,
        pressure_max=pressure_max,
        passed=bool(velocity_max <= BOUNDARY_VELOCITY_GATE and pressure_max <= BOUNDARY_PRESSURE_GATE),
    )


def _energy_metrics(candidate: BlackBoxFullNSCandidate, x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray, weights: np.ndarray | None) -> EnergyMetrics:
    if weights is None:
        return EnergyMetrics(by_time={}, reference_energy=None, passed=False, reason="volume_weights_missing")
    by_time: dict[str, float] = {}
    for time_value in sorted(float(v) for v in np.unique(t)):
        mask = np.isclose(t, time_value, rtol=0.0, atol=1.0e-14)
        if not np.any(mask):
            continue
        u = np.asarray(candidate.velocity(x[mask], y[mask], z[mask], t[mask]), dtype=float)
        w = np.asarray(weights, dtype=float)[mask]
        by_time[f"{time_value:.10g}"] = float(0.5 * np.sum(w * np.sum(u * u, axis=-1)))
    reference = by_time.get(f"{ENERGY_REFERENCE_TIME:.10g}")
    if reference is None:
        return EnergyMetrics(by_time=by_time, reference_energy=None, passed=False, reason="reference_time_missing")
    if abs(reference - ENERGY_REFERENCE) > ENERGY_REFERENCE_TOLERANCE:
        return EnergyMetrics(by_time=by_time, reference_energy=reference, passed=False, reason="reference_energy_gate")
    if any(value < ENERGY_MIN or value > ENERGY_MAX for value in by_time.values()):
        return EnergyMetrics(by_time=by_time, reference_energy=reference, passed=False, reason="validation_time_energy_gate")
    return EnergyMetrics(by_time=by_time, reference_energy=reference, passed=True, reason=None)


def _nu_sensitivity(level: FullNSLevel) -> tuple[float, float]:
    base = vector_rms(level.residual)
    delta = NU_PERTURBATION_FRACTION * NU * level.laplacian
    minus = vector_rms(level.residual + delta)
    plus = vector_rms(level.residual - delta)
    denom = max(base, 1.0e-30)
    return abs(minus - base) / denom, abs(plus - base) / denom


def validate_complete_candidate(
    candidate: BlackBoxFullNSCandidate,
    contract: FixedPhysicalContract,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    expected_stage: str,
    volume_weights: Any | None = None,
    quadrature_ladder_assessed: bool = False,
) -> FullNSValidationReceipt:
    admission = candidate_admission_failures(candidate, contract, expected_stage=expected_stage)
    if admission:
        raise CandidateAdmissionError(";".join(admission))
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    weights = None if volume_weights is None else np.broadcast_to(np.asarray(volume_weights, dtype=float), xb.shape)
    levels = three_resolution_levels(candidate, xb, yb, zb, tb)
    relative_change, ratio = _resolution_metrics(levels)
    boundary = frozen_boundary_support_metrics(candidate)
    energy = _energy_metrics(candidate, xb, yb, zb, tb, weights)
    fine = levels[-1]
    momentum_l2 = tuple(_volume_l2(level.residual / RESIDUAL_SCALE, weights) for level in levels)
    divergence_l2 = tuple(_volume_l2(level.divergence, weights) for level in levels)
    failures: list[str] = []
    if fine.normalized_sampled_max > MOMENTUM_GATE:
        failures.append("momentum_sampled_max")
    if momentum_l2[-1] is None:
        failures.append("momentum_volume_l2_missing")
    elif momentum_l2[-1] > MOMENTUM_GATE:
        failures.append("momentum_volume_l2")
    div_max = float(np.max(np.abs(fine.divergence)))
    div_rms = float(np.sqrt(np.mean(fine.divergence * fine.divergence)))
    if div_max > DIVERGENCE_GATE:
        failures.append("divergence_sampled_max")
    if divergence_l2[-1] is None:
        failures.append("divergence_volume_l2_missing")
    elif divergence_l2[-1] > DIVERGENCE_GATE:
        failures.append("divergence_volume_l2")
    m2f = vector_rms(levels[-2].residual - fine.residual)
    if relative_change > RESOLUTION_RELATIVE_CHANGE_GATE:
        failures.append("residual_resolution_relative_change")
    if m2f > RESOLUTION_FLOOR and ratio < RESOLUTION_RATIO_GATE:
        failures.append("residual_resolution_ratio")
    if vector_rms(fine.velocity) < NONTRIVIAL_VELOCITY_RMS_FLOOR:
        failures.append("velocity_nontriviality")
    if not boundary.passed:
        failures.append("boundary_support")
    if not energy.passed:
        failures.append("energy_nontriviality")
    if not quadrature_ladder_assessed:
        failures.append("quadrature_ladder_unassessed")
    metadata = dict(candidate.validation_metadata())
    nu_minus, nu_plus = _nu_sensitivity(fine)
    heldout_assessed = bool(not admission)
    pde_validated = bool(heldout_assessed and quadrature_ladder_assessed and not failures)
    return FullNSValidationReceipt(
        schema=SCHEMA,
        task=TASK,
        candidate_id=str(metadata["candidate_id"]),
        candidate_sha256=str(metadata["candidate_sha256"]),
        stage=expected_stage,
        physical_contract_sha256=contract.physical_contract_sha256,
        derivative_steps=DERIVATIVE_STEPS,
        momentum_sampled_max_by_level=tuple(level.normalized_sampled_max for level in levels),
        momentum_sampled_rms_by_level=tuple(level.normalized_sampled_rms for level in levels),
        momentum_volume_l2_by_level=momentum_l2,
        divergence_sampled_max_by_level=tuple(float(np.max(np.abs(level.divergence))) for level in levels),
        divergence_sampled_rms_by_level=tuple(float(np.sqrt(np.mean(level.divergence * level.divergence))) for level in levels),
        divergence_volume_l2_by_level=divergence_l2,
        resolution_relative_change=relative_change,
        resolution_ratio=ratio,
        boundary_support=boundary,
        energy=energy,
        nu_minus_relative_rms_change=nu_minus,
        nu_plus_relative_rms_change=nu_plus,
        gate_failures=tuple(failures),
        heldout_normalized_ns_residual_assessed=heldout_assessed,
        quadrature_ladder_assessed=bool(quadrature_ladder_assessed),
        pde_validated=pde_validated,
    )


def frozen_heldout_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the frozen 4096-point x six-time off-grid validation cloud.

    Weights are uniform Monte-Carlo volume weights per time.  They support a
    reproducible held-out volume-L2 estimate, but the separate canonical
    24/48/96 quadrature ladder still has to be assessed before PDE promotion.
    """
    rng = np.random.default_rng(HELDOUT_SEED)
    xyz0 = rng.uniform(-2.0, 2.0, size=(HELDOUT_POINTS, 3))
    xyz = np.tile(xyz0, (len(VALIDATION_TIMES), 1))
    t = np.repeat(np.asarray(VALIDATION_TIMES, dtype=float), HELDOUT_POINTS)
    volume = 4.0 ** 3
    weights = np.full(t.shape, volume / HELDOUT_POINTS, dtype=float)
    return xyz[:, 0], xyz[:, 1], xyz[:, 2], t, weights


def frozen_axis_near_probes() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    radii = np.asarray((0.0, 1.0e-12, 1.0e-9, 1.0e-6), dtype=float)
    angles = np.asarray((0.0, 0.371, 1.117), dtype=float)
    rows: list[tuple[float, float, float, float]] = []
    for time_value in VALIDATION_TIMES:
        for radius in radii:
            for angle in angles:
                rows.append((radius * math.cos(angle), radius * math.sin(angle), 0.13, time_value))
    arr = np.asarray(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def current_strict_inner_ineligibility() -> dict[str, Any]:
    """Machine-visible boundary for the currently available A2/A3 surfaces."""
    return {
        "latest_agent1_pr": LATEST_A1_PR,
        "latest_agent1_head": LATEST_A1_HEAD,
        "latest_agent2_pr": LATEST_A2_PR,
        "latest_agent2_head": LATEST_A2_HEAD,
        "latest_agent3_pr": LATEST_A3_PR,
        "latest_agent3_head": LATEST_A3_HEAD,
        "global_leading_velocity_materialized": False,
        "outer_join_materialized": False,
        "matched_pressure_included": False,
        "restricted_forcing_included": False,
        "agent3_correction_velocity_materialized": False,
        "complete_blackbox_candidate_available": False,
        "leading_only_ns_residual_assessed": False,
        "leading_plus_oscillatory_ns_residual_assessed": False,
        "after_correction_ns_residual_assessed": False,
        "pde_validated": False,
    }


def public_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_a4_pr": PARENT_A4_PR,
        "parent_a4_head": PARENT_A4_HEAD,
        "equation": "u_t + (u dot grad)u + grad(p) - nu Laplacian(u) - f",
        "public_surfaces_only": ["velocity(x,y,z,t)", "pressure(x,y,z,t)", "forcing(x,y,z,t)"],
        "independent_operator": "Cartesian FD4 five-point; endpoint-aware FD4 time derivative",
        "heldout_seed": HELDOUT_SEED,
        "heldout_points": HELDOUT_POINTS,
        "validation_times": list(VALIDATION_TIMES),
        "derivative_steps": list(DERIVATIVE_STEPS),
        "nu": NU,
        "residual_scale": RESIDUAL_SCALE,
        "frozen_gates": {
            "momentum_max": MOMENTUM_GATE,
            "momentum_volume_l2": MOMENTUM_GATE,
            "divergence_max": DIVERGENCE_GATE,
            "divergence_volume_l2": DIVERGENCE_GATE,
            "boundary_velocity_max": BOUNDARY_VELOCITY_GATE,
            "boundary_pressure_max": BOUNDARY_PRESSURE_GATE,
        },
        "staged_validation": ["leading_only", "leading_plus_oscillatory", "after_correction"],
        "current_ineligibility": current_strict_inner_ineligibility(),
        "truth_boundary": {
            "training_loss_is_validation": False,
            "construction_derivative_path_used": False,
            "residual_defined_free_forcing_allowed": False,
            "low_resolution_only_can_pass": False,
            "quadrature_ladder_required_for_pde_promotion": True,
            "pde_validated": False,
        },
    }
