"""Typed normalized-Haar bridge for Kokuno RF30 auxiliary-T2 covariance.

The corrected Kokuno RF30 bar is normalized Haar average on an auxiliary
two-torus, not the physical azimuthal mean.  PR #1161 froze the exact current-I4
fixed-Q kinematics while deliberately refusing to relabel the physical theta
average.  This module closes only the *mean-operation interface*:

    source auxiliary-T2 wave samples
      -> W^{rr}, W^{zr}, W^{theta theta}, W^{z theta}, W^{zz}
      -> normalized Haar averages
      -> d_R bar(W^{rr}), d_Z bar(W^{zr})

The quadratic products and slow derivatives are formed here from raw wave and
raw slow-derivative samples.  A provider cannot hand in pre-averaged covariance,
RF30 P/J defects, forcing, residuals, or correction coefficients.

No repository candidate is authorized merely by satisfying this Python
Protocol.  A future Agent-2 auxiliary-T2 provider must be checksum-pinned in a
separate increment before this bridge can be promoted to repository-candidate
RF30 evidence.  Until then the executable bridge is mechanics/provenance
infrastructure only.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
import math
from typing import Any, Protocol, runtime_checkable

import numpy as np

from .kokuno_current_i4_nonlinear_mean_attribution import ExactCurrentI4NonlinearBackend
from .kokuno_current_i4_rf30_fixedq_preflight import (
    CurrentI4FixedQObservablePreflight,
    materialize_current_i4_rf30_fixedq_preflight,
)

TASK = "KOKUNO-A3-AUXT2-HAAR-RF30-BRIDGE-127"
SCHEMA = "kokuno-a3-auxt2-haar-rf30-bridge-v1"
PARENT_AGENT3_PR = 1161
PARENT_AGENT3_HEAD = "cafab2098f00c3dc2ed5668e1b12bdb38a1f31e8"
PARENT_PREFLIGHT_SOURCE_BLOB = "b2fac4318a40342d6c08170bf7093c31685c53ab"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_HAAR_CONVENTION = "normalized Haar measure of total mass one on auxiliary T^2"
SOURCE_DERIVATIVE_RULE = (
    "torus-translation invariance: slow derivatives commute with normalized Haar mean"
)

TORUS_ORDERS = (12, 24, 48)
# Deliberately absent until an implementation-distinct Agent-2 source-T2 provider exists.
PINNED_REPOSITORY_PROVIDER_BLOB: str | None = None

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


class AuxiliaryT2HaarBridgeError(RuntimeError):
    """Raised when the source auxiliary-T2/Haar contract is violated."""


def _canonical_sha(payload: dict[str, object]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _valid_hex(value: str, length: int) -> bool:
    return len(value) == length and all(ch in "0123456789abcdef" for ch in value)


def _relative_difference(new: np.ndarray, old: np.ndarray) -> float:
    n = np.asarray(new, dtype=float)
    o = np.asarray(old, dtype=float)
    denom = float(np.linalg.norm(n.ravel()))
    numer = float(np.linalg.norm((n - o).ravel()))
    if denom == 0.0:
        return 0.0 if numer == 0.0 else math.inf
    return numer / denom


@dataclass(frozen=True)
class AuxiliaryT2ProviderMetadata:
    candidate_semantic_sha256: str
    oscillatory_runtime_sha256: str
    leading_semantic_sha256: str
    provider_semantic_sha256: str
    source_blob_sha1: str
    provider_kind: str
    source_auxiliary_t2_field_materialized: bool
    normalized_haar_measure_total_mass_one: bool
    recomputed_from_actual_candidate: bool
    surrogate_or_preaveraged_covariance_used: bool
    heldout_data_used: bool
    residual_as_forcing_used: bool


@dataclass(frozen=True)
class AuxiliaryT2WaveSamples:
    w_r: np.ndarray
    w_theta: np.ndarray
    w_z: np.ndarray
    dR_w_r: np.ndarray
    dR_w_theta: np.ndarray
    dR_w_z: np.ndarray
    dZ_w_r: np.ndarray
    dZ_w_theta: np.ndarray
    dZ_w_z: np.ndarray


@runtime_checkable
class SourceAuxiliaryT2WaveProvider(Protocol):
    metadata: AuxiliaryT2ProviderMetadata

    def sample_auxiliary_wave(
        self,
        *,
        R: np.ndarray,
        Z: float,
        T: float,
        y1: np.ndarray,
        y2: np.ndarray,
    ) -> AuxiliaryT2WaveSamples:
        """Return raw cylindrical wave and slow derivatives on the source T2 grid."""


@dataclass(frozen=True)
class AuxiliaryT2HaarCovarianceReceipt:
    candidate_semantic_sha256: str
    oscillatory_runtime_sha256: str
    leading_semantic_sha256: str
    fixed_q_preflight_sha256: str
    provider_semantic_sha256: str
    provider_source_blob_sha1: str
    provider_kind: str
    Q: float
    radius_R: tuple[float, ...]
    Z: float
    T: float
    torus_orders: tuple[int, ...]
    successive_covariance_relative_differences: tuple[float, ...]
    mean_W_rr_auxiliary_haar: tuple[float, ...]
    mean_W_zr_auxiliary_haar: tuple[float, ...]
    mean_W_thetatheta_auxiliary_haar: tuple[float, ...]
    mean_W_ztheta_auxiliary_haar: tuple[float, ...]
    mean_W_zz_auxiliary_haar: tuple[float, ...]
    dR_mean_W_rr_auxiliary_haar: tuple[float, ...]
    dZ_mean_W_zr_auxiliary_haar: tuple[float, ...]
    mean_w_r_auxiliary_haar: tuple[float, ...]
    mean_w_theta_auxiliary_haar: tuple[float, ...]
    mean_w_z_auxiliary_haar: tuple[float, ...]
    normalized_source_auxiliary_t2_haar_mean_used: bool
    covariance_formed_in_bridge_from_raw_wave_samples: bool
    covariance_derivatives_formed_in_bridge_from_raw_wave_derivatives: bool
    provider_checksum_pinned_for_repository_candidate: bool
    rf30_covariance_bridge_materialized: bool
    rf30_repository_candidate_state_authorized: bool
    receipt_sha256: str

    def to_dict(self) -> dict[str, object]:
        return {
            **asdict(self),
            "schema": SCHEMA,
            "task": TASK,
            "source_provenance": source_provenance(),
            "truth_boundary": truth_boundary(),
        }


def source_provenance() -> dict[str, object]:
    return {
        "repository": SOURCE_READER_REPO,
        "commit": SOURCE_READER_HEAD,
        "date": SOURCE_READER_DATE,
        "path": SOURCE_READER_PATH,
        "blob_sha1": SOURCE_READER_BLOB,
        "haar_convention": SOURCE_HAAR_CONVENTION,
        "derivative_rule": SOURCE_DERIVATIVE_RULE,
        "classification": "public_structure_provenance_not_independent_validation",
    }


def _validate_metadata(
    metadata: AuxiliaryT2ProviderMetadata,
    preflight: CurrentI4FixedQObservablePreflight,
) -> None:
    if not isinstance(metadata, AuxiliaryT2ProviderMetadata):
        raise TypeError("provider.metadata must be AuxiliaryT2ProviderMetadata")
    identities = {
        "candidate_semantic_sha256": preflight.candidate_semantic_sha256,
        "oscillatory_runtime_sha256": preflight.oscillatory_runtime_sha256,
        "leading_semantic_sha256": preflight.leading_semantic_sha256,
    }
    for name, expected in identities.items():
        actual = getattr(metadata, name)
        if actual != expected:
            raise AuxiliaryT2HaarBridgeError(f"{name} drifted from exact #1161 preflight")
    if not _valid_hex(metadata.provider_semantic_sha256, 64):
        raise AuxiliaryT2HaarBridgeError("provider_semantic_sha256 must be lowercase sha256")
    if not _valid_hex(metadata.source_blob_sha1, 40):
        raise AuxiliaryT2HaarBridgeError("source_blob_sha1 must be lowercase git blob sha1")
    if metadata.provider_kind not in {"mechanics_fixture", "repository_candidate"}:
        raise AuxiliaryT2HaarBridgeError("provider_kind must be mechanics_fixture or repository_candidate")
    if not metadata.source_auxiliary_t2_field_materialized:
        raise AuxiliaryT2HaarBridgeError("provider must expose the source auxiliary-T2 field")
    if not metadata.normalized_haar_measure_total_mass_one:
        raise AuxiliaryT2HaarBridgeError("provider must use the source normalized Haar convention")
    if metadata.surrogate_or_preaveraged_covariance_used:
        raise AuxiliaryT2HaarBridgeError("surrogate/pre-averaged covariance is forbidden")
    if metadata.heldout_data_used:
        raise AuxiliaryT2HaarBridgeError("held-out data cannot construct the correction state")
    if metadata.residual_as_forcing_used:
        raise AuxiliaryT2HaarBridgeError("residual-as-forcing shortcut is forbidden")
    if metadata.provider_kind == "repository_candidate" and not metadata.recomputed_from_actual_candidate:
        raise AuxiliaryT2HaarBridgeError(
            "repository candidate provider must recompute from the actual candidate"
        )


def _as_sample_array(value: Any, shape: tuple[int, int, int], label: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != shape:
        raise AuxiliaryT2HaarBridgeError(
            f"{label} has shape {out.shape}; expected {shape}"
        )
    if not np.all(np.isfinite(out)):
        raise AuxiliaryT2HaarBridgeError(f"{label} contains non-finite values")
    return out


def _haar_covariance_for_order(
    provider: SourceAuxiliaryT2WaveProvider,
    preflight: CurrentI4FixedQObservablePreflight,
    order: int,
) -> tuple[np.ndarray, np.ndarray]:
    R = np.asarray(preflight.radius_R, dtype=float)
    y = np.arange(order, dtype=float) / float(order)
    y1, y2 = np.meshgrid(y, y, indexing="ij")
    # Broadcasted provider target is (n_R, order, order).  Equal weights implement
    # normalized Haar measure with total mass one on the discrete torus.
    sample = provider.sample_auxiliary_wave(
        R=R,
        Z=float(preflight.Z),
        T=float(preflight.T),
        y1=y1,
        y2=y2,
    )
    if not isinstance(sample, AuxiliaryT2WaveSamples):
        raise TypeError("sample_auxiliary_wave must return AuxiliaryT2WaveSamples")
    shape = (R.size, order, order)
    wr = _as_sample_array(sample.w_r, shape, "w_r")
    wt = _as_sample_array(sample.w_theta, shape, "w_theta")
    wz = _as_sample_array(sample.w_z, shape, "w_z")
    dR_wr = _as_sample_array(sample.dR_w_r, shape, "dR_w_r")
    _as_sample_array(sample.dR_w_theta, shape, "dR_w_theta")
    _as_sample_array(sample.dR_w_z, shape, "dR_w_z")
    dZ_wr = _as_sample_array(sample.dZ_w_r, shape, "dZ_w_r")
    _as_sample_array(sample.dZ_w_theta, shape, "dZ_w_theta")
    dZ_wz = _as_sample_array(sample.dZ_w_z, shape, "dZ_w_z")

    mean = lambda value: np.mean(value, axis=(1, 2))
    covariance = np.stack(
        (
            mean(wr * wr),
            mean(wz * wr),
            mean(wt * wt),
            mean(wz * wt),
            mean(wz * wz),
            mean(2.0 * wr * dR_wr),
            mean(dZ_wz * wr + wz * dZ_wr),
        ),
        axis=-1,
    )
    wave_mean = np.stack((mean(wr), mean(wt), mean(wz)), axis=-1)
    if not np.all(np.isfinite(covariance)) or not np.all(np.isfinite(wave_mean)):
        raise AuxiliaryT2HaarBridgeError("normalized Haar observables became non-finite")
    return covariance, wave_mean


def materialize_current_i4_rf30_auxiliary_haar_bridge(
    backend: ExactCurrentI4NonlinearBackend,
    provider: SourceAuxiliaryT2WaveProvider,
    radius: Any,
    z: Any,
    t: Any,
) -> AuxiliaryT2HaarCovarianceReceipt:
    """Compute RF30 covariance observables from raw source-T2 wave samples.

    This function always replays #1161 first, so the auxiliary provider is bound
    to the same candidate identity and fixed-Q chart.  It deliberately does not
    construct P, J_theta, J_z or apply a correction.
    """
    preflight = materialize_current_i4_rf30_fixedq_preflight(
        backend, radius, z, t
    )
    metadata = getattr(provider, "metadata", None)
    _validate_metadata(metadata, preflight)
    if not callable(getattr(provider, "sample_auxiliary_wave", None)):
        raise AuxiliaryT2HaarBridgeError("provider lacks sample_auxiliary_wave")

    covariances: list[np.ndarray] = []
    wave_means: list[np.ndarray] = []
    for order in TORUS_ORDERS:
        covariance, wave_mean = _haar_covariance_for_order(
            provider, preflight, int(order)
        )
        covariances.append(covariance)
        wave_means.append(wave_mean)
    convergence = tuple(
        _relative_difference(covariances[i], covariances[i - 1])
        for i in range(1, len(covariances))
    )
    finest = covariances[-1]
    finest_mean = wave_means[-1]

    pinned = (
        metadata.provider_kind == "repository_candidate"
        and PINNED_REPOSITORY_PROVIDER_BLOB is not None
        and metadata.source_blob_sha1 == PINNED_REPOSITORY_PROVIDER_BLOB
        and metadata.recomputed_from_actual_candidate
    )
    # A typed provider alone is not scientific provenance.  Authorization remains
    # false until a future increment pins the exact implementation-distinct A2 blob.
    authorized = bool(
        pinned
        and preflight.all_strict_I4
        and metadata.source_auxiliary_t2_field_materialized
        and metadata.normalized_haar_measure_total_mass_one
    )

    payload: dict[str, object] = {
        "schema": SCHEMA,
        "task": TASK,
        "candidate_semantic_sha256": preflight.candidate_semantic_sha256,
        "oscillatory_runtime_sha256": preflight.oscillatory_runtime_sha256,
        "leading_semantic_sha256": preflight.leading_semantic_sha256,
        "fixed_q_preflight_sha256": preflight.preflight_sha256,
        "provider_semantic_sha256": metadata.provider_semantic_sha256,
        "provider_source_blob_sha1": metadata.source_blob_sha1,
        "provider_kind": metadata.provider_kind,
        "Q": preflight.Q,
        "radius_R": list(preflight.radius_R),
        "Z": preflight.Z,
        "T": preflight.T,
        "torus_orders": list(TORUS_ORDERS),
        "convergence": list(convergence),
        "covariance": finest.tolist(),
        "wave_mean": finest_mean.tolist(),
        "normalized_source_auxiliary_t2_haar_mean_used": True,
        "covariance_formed_in_bridge_from_raw_wave_samples": True,
        "covariance_derivatives_formed_in_bridge_from_raw_wave_derivatives": True,
        "provider_checksum_pinned_for_repository_candidate": pinned,
        "rf30_repository_candidate_state_authorized": authorized,
    }
    receipt_sha = _canonical_sha(payload)
    return AuxiliaryT2HaarCovarianceReceipt(
        candidate_semantic_sha256=preflight.candidate_semantic_sha256,
        oscillatory_runtime_sha256=preflight.oscillatory_runtime_sha256,
        leading_semantic_sha256=preflight.leading_semantic_sha256,
        fixed_q_preflight_sha256=preflight.preflight_sha256,
        provider_semantic_sha256=metadata.provider_semantic_sha256,
        provider_source_blob_sha1=metadata.source_blob_sha1,
        provider_kind=metadata.provider_kind,
        Q=float(preflight.Q),
        radius_R=tuple(float(v) for v in preflight.radius_R),
        Z=float(preflight.Z),
        T=float(preflight.T),
        torus_orders=tuple(int(v) for v in TORUS_ORDERS),
        successive_covariance_relative_differences=convergence,
        mean_W_rr_auxiliary_haar=tuple(float(v) for v in finest[:, 0]),
        mean_W_zr_auxiliary_haar=tuple(float(v) for v in finest[:, 1]),
        mean_W_thetatheta_auxiliary_haar=tuple(float(v) for v in finest[:, 2]),
        mean_W_ztheta_auxiliary_haar=tuple(float(v) for v in finest[:, 3]),
        mean_W_zz_auxiliary_haar=tuple(float(v) for v in finest[:, 4]),
        dR_mean_W_rr_auxiliary_haar=tuple(float(v) for v in finest[:, 5]),
        dZ_mean_W_zr_auxiliary_haar=tuple(float(v) for v in finest[:, 6]),
        mean_w_r_auxiliary_haar=tuple(float(v) for v in finest_mean[:, 0]),
        mean_w_theta_auxiliary_haar=tuple(float(v) for v in finest_mean[:, 1]),
        mean_w_z_auxiliary_haar=tuple(float(v) for v in finest_mean[:, 2]),
        normalized_source_auxiliary_t2_haar_mean_used=True,
        covariance_formed_in_bridge_from_raw_wave_samples=True,
        covariance_derivatives_formed_in_bridge_from_raw_wave_derivatives=True,
        provider_checksum_pinned_for_repository_candidate=pinned,
        rf30_covariance_bridge_materialized=True,
        rf30_repository_candidate_state_authorized=authorized,
        receipt_sha256=receipt_sha,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_i4_rf30_auxiliary_haar_bridge)
    forbidden = {
        "residual", "defect", "forcing", "pressure", "target", "gain", "damping",
        "threshold", "viscosity", "nu", "correction", "covariance", "P",
        "J_theta", "J_z", "heldout",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_preflight_source_blob": PARENT_PREFLIGHT_SOURCE_BLOB,
        "source_reader_head": SOURCE_READER_HEAD,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(
            forbidden.intersection(signature.parameters)
        ),
        "normalized_source_auxiliary_t2_haar_operator_executable": True,
        "haar_measure_total_mass_one": True,
        "covariance_formed_in_bridge_from_raw_wave_samples": True,
        "covariance_derivatives_formed_in_bridge_from_raw_wave_derivatives": True,
        "preaveraged_covariance_input_exposed": False,
        "physical_theta_mean_promoted_to_source_auxiliary_haar": False,
        "agent2_curl_or_jacobian_reimplemented": False,
        "repository_provider_blob_pinned": PINNED_REPOSITORY_PROVIDER_BLOB is not None,
        "rf30_repository_candidate_state_authorized": False,
        "rf30_defect_materialized": False,
        "rf31_five_row_system_materialized_from_this_bridge": False,
        "correction_velocity_materialized": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proof_claimed": False,
        "pde_validated": False,
    }
