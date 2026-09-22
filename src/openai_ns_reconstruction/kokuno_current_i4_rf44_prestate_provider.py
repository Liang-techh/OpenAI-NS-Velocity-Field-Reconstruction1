"""Current-I4 RF44 pre-update state from the actual complete-curl wave provider.

This is one narrow Kokuno Agent-2 handoff surface.  It stacks on the exact
pinned-provider correction frontier from A3 PR #1210, but it does not implement
or solve any mean correction.  Instead it extends the exact A2 PR #1198 raw
auxiliary-T2 provider with the *pre-update* state required by the existing A3
RF44--RF49 operator.

The source-side facts used here are deliberately limited:

* physical velocity is Q^(-A) times the fixed-Q normalized velocity;
* RF30/RF44 use normalized Haar covariance on the auxiliary T^2;
* before the mean update beta=gamma=v=0 and RF44 reduces to RF30;
* RF44 retains the actual normalized radial background b.

The repository realization reuses #1198 unchanged for the raw complete-curl
wave and its numerical slow derivatives.  The #1198 y1/y2 -> sign-channel
phase lift remains repository-autonomous and is not promoted to source-exact
mode data.  The radial background below is reconstructed from the exact bound
current-I4 leading Cartesian velocity at the same physical fixed-Q probes and
scaled by Q^A.  No post-update state, correction increment, pressure, forcing,
held-out residual, gain criterion, or scientific threshold is produced here.
"""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_current_i4_auxiliary_t2_provider import (
    KokunoCurrentI4AuxiliaryT2Provider,
)
from .kokuno_current_i4_nonlinear_mean_attribution import (
    ExactCurrentI4NonlinearBackend,
)
from .kokuno_current_i4_rf30_fixedq_preflight import (
    CurrentI4FixedQObservablePreflight,
    materialize_current_i4_rf30_fixedq_preflight,
)
from .kokuno_current_i4_rf30_typed_defect import (
    CurrentI4RF30TypedDefectReceipt,
    materialize_current_i4_rf30_typed_defect,
)
from .kokuno_rf34_rf39_compact_mean_correction import RF34RF39CorrectionReceipt
from .kokuno_rf44_rf49_postupdate_recompute import RF44MeanState

TASK = "K2-OSC-113"
SCHEMA = "kokuno-a2-current-i4-rf44-prestate-provider-v1"

PARENT_AGENT3_PR = 1210
PARENT_AGENT3_HEAD = "e75940c33127b0725ddc703eaae97590158f2a56"
PARENT_AGENT2_RAW_PROVIDER_PR = 1198
PARENT_AGENT2_RAW_PROVIDER_HEAD = "5765c2b3bab7482df514efa87f9b6bba48e04b7f"
PARENT_AGENT2_RAW_PROVIDER_BLOB = "96168ac6583ad5e71aa044bd558c6feeadf051f4"

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_FORMULAS = ("RF30", "RF44")
SOURCE_PHYSICAL_VELOCITY_SCALING = "u_phys=Q^(-A)u_*"

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


class CurrentI4RF44PreStateProviderError(RuntimeError):
    """Raised when the pre-update state leaves its exact candidate/chart identity."""


def _canonical_sha256(payload: Any) -> str:
    raw = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _git_blob_sha1(path: str | Path) -> str:
    data = Path(path).read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _finite_vector(value: Any, label: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim != 1 or out.size < 3 or not np.all(np.isfinite(out)):
        raise CurrentI4RF44PreStateProviderError(
            f"{label} must be one finite vector with at least three entries"
        )
    return out


def _current_candidate(current: Any) -> Any:
    candidate = getattr(current, "candidate", None)
    if candidate is None:
        raise CurrentI4RF44PreStateProviderError("current correction receipt lacks candidate identity")
    for name in ("candidate_id", "candidate_sha256"):
        if not getattr(candidate, name, None):
            raise CurrentI4RF44PreStateProviderError(
                f"current correction candidate lacks {name}"
            )
    return candidate


def source_provenance() -> dict[str, object]:
    return {
        "repository": SOURCE_REPOSITORY,
        "commit": SOURCE_READER_HEAD,
        "date": SOURCE_READER_DATE,
        "path": SOURCE_READER_PATH,
        "blob_sha1": SOURCE_READER_BLOB,
        "source_formulas": SOURCE_FORMULAS,
        "physical_velocity_scaling": SOURCE_PHYSICAL_VELOCITY_SCALING,
        "rf30_preupdate_mean_tuple": "beta=gamma=v=0; covariance W retained",
        "rf44_preupdate_reduction": (
            "with beta=gamma=v=0, tstar_beta=epsilon_laplacian_beta=0, "
            "dR_rr_bundle=dR Wrr and dZ_zr_bundle=dZ Wzr, RF44 reduces to RF30"
        ),
        "raw_wave_provider": {
            "pr": PARENT_AGENT2_RAW_PROVIDER_PR,
            "head": PARENT_AGENT2_RAW_PROVIDER_HEAD,
            "blob_sha1": PARENT_AGENT2_RAW_PROVIDER_BLOB,
            "classification": (
                "actual_complete_curl_runtime_with_repository_autonomous_auxiliary_phase_lift"
            ),
        },
        "classification": "public_structure_provenance_not_independent_validation",
    }


class KokunoCurrentI4RF44PreStateProvider:
    """Delegate exact #1198 raw waves and materialize only the RF44 pre-state."""

    def __init__(
        self,
        backend: ExactCurrentI4NonlinearBackend,
        preflight: CurrentI4FixedQObservablePreflight,
    ) -> None:
        if not isinstance(backend, ExactCurrentI4NonlinearBackend):
            raise TypeError("backend must be ExactCurrentI4NonlinearBackend")
        if not isinstance(preflight, CurrentI4FixedQObservablePreflight):
            raise TypeError("preflight must be CurrentI4FixedQObservablePreflight")
        if not preflight.all_strict_I4:
            raise CurrentI4RF44PreStateProviderError(
                "RF44 pre-state provider requires one strict current-I4 preflight"
            )

        raw = KokunoCurrentI4AuxiliaryT2Provider(backend, preflight)
        if raw.metadata.source_blob_sha1 != PARENT_AGENT2_RAW_PROVIDER_BLOB:
            raise CurrentI4RF44PreStateProviderError(
                "exact A2 #1198 raw-provider implementation blob drifted"
            )
        if raw.metadata.provider_semantic_sha256 != raw.semantic_sha256:
            raise CurrentI4RF44PreStateProviderError(
                "raw-provider semantic identity is internally inconsistent"
            )
        self._backend = backend
        self._preflight = preflight
        self._raw = raw

    @property
    def metadata(self):
        """Preserve the exact #1198 metadata used by the A3 RF30 checksum pin."""
        return self._raw.metadata

    @property
    def semantic_sha256(self) -> str:
        """Preserve the exact #1198 raw-wave semantic identity for RF30."""
        return self._raw.semantic_sha256

    @property
    def chart(self):
        return self._raw.chart

    def sample_auxiliary_wave(self, **kwargs):
        """Delegate the exact #1198 raw complete-curl wave surface unchanged."""
        return self._raw.sample_auxiliary_wave(**kwargs)

    def _physical_probe_coordinates(self) -> tuple[np.ndarray, float, float]:
        p = self._preflight
        radius = _finite_vector(p.radius_physical, "radius_physical")
        z_phys = float((p.Q ** p.D) * p.Z)
        t_phys = float(1.0 - p.Q * p.T)
        if not math.isfinite(z_phys) or not math.isfinite(t_phys):
            raise CurrentI4RF44PreStateProviderError(
                "fixed-Q chart mapped to non-finite physical coordinates"
            )
        return radius, z_phys, t_phys

    def _normalized_leading_radial_background(self) -> np.ndarray:
        """Return source-chart b=Q^A u_r^lead on the exact bound physical rings."""
        leading = getattr(self._backend.composite_field, "leading_backend", None)
        if leading is None or not callable(getattr(leading, "velocity", None)):
            raise CurrentI4RF44PreStateProviderError(
                "bound current-I4 leading backend lacks Cartesian velocity"
            )
        if str(getattr(leading, "semantic_sha256", "")) != self._preflight.leading_semantic_sha256:
            raise CurrentI4RF44PreStateProviderError(
                "leading semantic identity drifted before RF44 pre-state materialization"
            )

        radius, z_phys, t_phys = self._physical_probe_coordinates()
        zeros = np.zeros_like(radius)
        velocity = np.asarray(
            leading.velocity(
                radius,
                zeros,
                np.full_like(radius, z_phys),
                np.full_like(radius, t_phys),
            ),
            dtype=float,
        )
        if velocity.shape != (radius.size, 3) or not np.all(np.isfinite(velocity)):
            raise CurrentI4RF44PreStateProviderError(
                "bound leading velocity returned an invalid radial-line batch"
            )

        # The physical probes lie at theta=0, so Cartesian u_x is u_r.  The
        # corrected fixed-Q scaling is u_phys=Q^(-A)u_*, hence b=u_{*,r}=Q^A u_r.
        b = (self._preflight.Q ** self._preflight.A) * velocity[:, 0]
        if not np.all(np.isfinite(b)):
            raise CurrentI4RF44PreStateProviderError(
                "normalized current-I4 radial background became non-finite"
            )
        return b

    def _typed_pre_mean_state(self) -> CurrentI4RF30TypedDefectReceipt:
        radius, z_phys, t_phys = self._physical_probe_coordinates()
        typed = materialize_current_i4_rf30_typed_defect(
            self._backend,
            self,
            radius,
            z_phys,
            t_phys,
        )
        if not isinstance(typed, CurrentI4RF30TypedDefectReceipt):
            raise CurrentI4RF44PreStateProviderError(
                "RF30 typed operator did not return its typed receipt"
            )
        if typed.fixed_q_preflight_sha256 != self._preflight.preflight_sha256:
            raise CurrentI4RF44PreStateProviderError(
                "RF30 typed state detached from the exact bound fixed-Q preflight"
            )
        if typed.provider_semantic_sha256 != self.metadata.provider_semantic_sha256:
            raise CurrentI4RF44PreStateProviderError(
                "RF30 typed state detached from the exact #1198 provider identity"
            )
        state = typed.state
        if not all(
            (
                state.actual_candidate_recomputed,
                state.oscillatory_covariance_recomputed,
                state.normalized_haar_mean_used,
                state.source_fixed_q_chart_used,
            )
        ):
            raise CurrentI4RF44PreStateProviderError(
                "RF30 typed state lacks actual-candidate/Haar/fixed-Q provenance"
            )
        if (
            state.surrogate_defect_used
            or state.residual_as_forcing_shortcut_used
            or state.heldout_samples_used_to_construct_state
        ):
            raise CurrentI4RF44PreStateProviderError(
                "RF30 typed state used surrogate/held-out/residual-as-forcing evidence"
            )
        return typed

    @staticmethod
    def _require_handoff_identity(
        current: Any,
        correction: RF34RF39CorrectionReceipt,
        typed: CurrentI4RF30TypedDefectReceipt,
    ) -> None:
        if not isinstance(correction, RF34RF39CorrectionReceipt):
            raise TypeError("correction must be RF34RF39CorrectionReceipt")
        candidate = _current_candidate(current)
        identities = (
            (candidate.candidate_id, candidate.candidate_sha256),
            (correction.candidate.candidate_id, correction.candidate.candidate_sha256),
            (typed.candidate.candidate_id, typed.candidate.candidate_sha256),
        )
        if len(set(identities)) != 1:
            raise CurrentI4RF44PreStateProviderError(
                "RF44 pre-state candidate identity drifted across current/correction/RF30"
            )
        current_chart = (
            str(getattr(current, "source_chart_id", "")),
            str(getattr(current, "source_chart_sha256", "")),
        )
        correction_chart = (correction.source_chart_id, correction.source_chart_sha256)
        typed_chart = (typed.source_chart_id, typed.source_chart_sha256)
        if not current_chart[0] or not current_chart[1]:
            raise CurrentI4RF44PreStateProviderError(
                "current correction receipt lacks fixed-Q source-chart identity"
            )
        if not (current_chart == correction_chart == typed_chart):
            raise CurrentI4RF44PreStateProviderError(
                "RF44 pre-state fixed-Q source-chart identity drifted"
            )

    def rf44_pre_update_state(
        self,
        *,
        current: Any,
        correction: RF34RF39CorrectionReceipt,
    ) -> RF44MeanState:
        """Build the actual current-I4 RF44 state *before* any mean correction."""
        typed = self._typed_pre_mean_state()
        self._require_handoff_identity(current, correction, typed)
        state = typed.state
        R = _finite_vector(state.radius_R, "RF30 radius_R")
        b = self._normalized_leading_radial_background()
        if b.shape != R.shape:
            raise CurrentI4RF44PreStateProviderError(
                "leading radial background and RF30 radial grid sizes differ"
            )
        zeros = tuple(0.0 for _ in range(R.size))

        # RF30 is exactly the pre-mean specialization of RF44: beta=gamma=v=0.
        # The two derivative bundles therefore reduce to derivatives of the Haar
        # wave covariance already materialized by the typed RF30 state.
        return RF44MeanState(
            source_candidate_id=typed.candidate.candidate_id,
            source_candidate_sha256=typed.candidate.candidate_sha256,
            source_chart_id=typed.source_chart_id,
            source_chart_sha256=typed.source_chart_sha256,
            radius_R=tuple(float(v) for v in R),
            base_b=tuple(float(v) for v in b),
            base_G=tuple(float(v) for v in state.base_G),
            base_V=tuple(float(v) for v in state.base_V),
            beta=zeros,
            gamma=zeros,
            v=zeros,
            mean_W_rr=tuple(float(v) for v in state.mean_W_rr),
            mean_W_zr=tuple(float(v) for v in state.mean_W_zr),
            mean_W_thetatheta=tuple(float(v) for v in state.mean_W_thetatheta),
            mean_W_ztheta=tuple(float(v) for v in state.mean_W_ztheta),
            mean_W_zz=tuple(float(v) for v in state.mean_W_zz),
            tstar_beta=zeros,
            dR_rr_bundle=tuple(float(v) for v in state.dR_mean_W_rr),
            dZ_zr_bundle=tuple(float(v) for v in state.dZ_mean_W_zr),
            epsilon_laplacian_beta=zeros,
            actual_candidate_recomputed=True,
            normalized_haar_mean_used=True,
            source_fixed_q_chart_used=True,
            operator_terms_recomputed_from_state=True,
            correction_applied=False,
            heldout_samples_used_to_construct_state=False,
            surrogate_defect_used=False,
            residual_as_forcing_shortcut_used=False,
        )

    def rf44_prestate_semantic_payload(self) -> dict[str, object]:
        b = self._normalized_leading_radial_background()
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent3": {"pr": PARENT_AGENT3_PR, "head": PARENT_AGENT3_HEAD},
            "raw_provider": {
                "pr": PARENT_AGENT2_RAW_PROVIDER_PR,
                "head": PARENT_AGENT2_RAW_PROVIDER_HEAD,
                "blob_sha1": PARENT_AGENT2_RAW_PROVIDER_BLOB,
                "provider_semantic_sha256": self.semantic_sha256,
            },
            "fixed_q_preflight_sha256": self._preflight.preflight_sha256,
            "leading_semantic_sha256": self._preflight.leading_semantic_sha256,
            "normalized_radial_background_b": [float(v) for v in b],
            "pre_mean_tuple": {"beta": 0.0, "gamma": 0.0, "v": 0.0},
            "source_provenance": source_provenance(),
            "truth_boundary": truth_boundary(),
        }

    @property
    def rf44_prestate_semantic_sha256(self) -> str:
        return _canonical_sha256(self.rf44_prestate_semantic_payload())

    @property
    def rf44_prestate_source_blob_sha1(self) -> str:
        path = inspect.getsourcefile(type(self))
        if path is None:
            raise CurrentI4RF44PreStateProviderError(
                "RF44 pre-state provider source file is unavailable"
            )
        return _git_blob_sha1(path)


def bind_current_i4_rf44_prestate_provider(
    backend: ExactCurrentI4NonlinearBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> KokunoCurrentI4RF44PreStateProvider:
    """Bind the exact current-I4 fixed-Q identity before RF44 state evaluation."""
    preflight = materialize_current_i4_rf30_fixedq_preflight(backend, radius, z, t)
    return KokunoCurrentI4RF44PreStateProvider(backend, preflight)


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(bind_current_i4_rf44_prestate_provider)
    forbidden = {
        "residual",
        "defect",
        "forcing",
        "pressure",
        "target",
        "gain",
        "damping",
        "threshold",
        "viscosity",
        "nu",
        "coefficient",
        "amplitude",
        "phase",
        "heldout",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent2_raw_provider_pr": PARENT_AGENT2_RAW_PROVIDER_PR,
        "parent_agent2_raw_provider_head": PARENT_AGENT2_RAW_PROVIDER_HEAD,
        "parent_agent2_raw_provider_blob": PARENT_AGENT2_RAW_PROVIDER_BLOB,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_blob": SOURCE_READER_BLOB,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(
            forbidden.intersection(signature.parameters)
        ),
        "exact_a2_1198_raw_provider_delegated_unchanged": True,
        "actual_complete_curl_wave_reused": True,
        "normalized_source_auxiliary_t2_haar_state_reused": True,
        "current_i4_leading_radial_background_recomputed": True,
        "physical_to_normalized_velocity_scaling_Q_to_A_applied": True,
        "rf30_preupdate_tuple_beta_gamma_v_zero": True,
        "rf44_preupdate_state_materialized": True,
        "rf44_preupdate_reduces_exactly_to_rf30_operator": True,
        "rf44_prestate_has_separate_semantic_identity": True,
        "raw_provider_blob_pin_relabelled_as_rf44_prestate_pin": False,
        "auxiliary_sign_phase_lift_repository_autonomous": True,
        "source_exact_auxiliary_phase_assignment_recovered": False,
        "source_exact_slow_derivatives_claimed": False,
        "source_exact_rf44_prestate_claimed": False,
        "rf44_postupdate_state_materialized": False,
        "rf44_correction_increment_materialized": False,
        "agent3_mean_correction_solved_here": False,
        "cartesian_correction_velocity_materialized": False,
        "rf44_rf49_repository_remainder_recomputed": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "residual_defined_free_forcing_allowed": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proof_claimed": False,
        "pde_validated": False,
    }
