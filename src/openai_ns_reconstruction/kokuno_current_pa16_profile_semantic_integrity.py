"""Authenticate the #959 joined-profile configuration before #961 consumption.

CR002 scope only.  Agent-3 #961 pins the exact #959 source blob and carries the
joined profile's reported ``semantic_sha256``, but its runtime binder only checks
that the report field is a 64-character string.  #959 defines that semantic
identity more strongly as

    sha256(canonical_json(KokunoPA16CurrentJoinedProfile.configuration()))

This guard replays that definition from a typed backend before delegating to the
existing #961 application binder.  It also pins the authenticated identity/report
objects during delegation so a second report fetch cannot silently swap identity.

This is repository-autonomous provenance machinery.  It changes no profile,
velocity, pressure, forcing, PA.16 coefficient, validation sample, or scientific
threshold, and it does not turn source-coordinate repair application into a
Cartesian correction velocity or Navier--Stokes evidence.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import re
from typing import Any, Mapping, Protocol

from . import kokuno_current_pa16_profile_application as app


TASK = "CR002-KOKUNO-JOINED-PROFILE-SEMANTIC-INTEGRITY-095"
SCHEMA = "cr002-kokuno-joined-profile-semantic-integrity-v1"
PARENT_AGENT3_PR = 961
PARENT_AGENT3_HEAD = "670f7b71d72735d441dcf48f2391c271ee5c163c"
PARENT_AGENT3_SOURCE_PATH = (
    "src/openai_ns_reconstruction/kokuno_current_pa16_profile_application.py"
)
PARENT_AGENT3_SOURCE_BLOB = "cb3dcd25e0c55a7935daee1b77ffb1741320f83b"
UPSTREAM_AGENT1_PR = 959
UPSTREAM_AGENT1_HEAD = app.UPSTREAM_AGENT1_HEAD
UPSTREAM_AGENT1_SOURCE_PATH = app.UPSTREAM_AGENT1_SOURCE_PATH
UPSTREAM_AGENT1_SOURCE_BLOB = app.UPSTREAM_AGENT1_SOURCE_BLOB
UPSTREAM_AGENT1_SCHEMA = app.UPSTREAM_AGENT1_SCHEMA
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class CurrentPA16ProfileSemanticIntegrityError(RuntimeError):
    """Raised when the runtime joined-profile semantic identity is not authentic."""


class CurrentPA16ProfileSemanticIntegrityBackend(
    app.CurrentPA16ProfileApplicationBackend, Protocol
):
    def joined_profile_configuration(self) -> Mapping[str, Any]: ...


@dataclass(frozen=True)
class CurrentPA16ProfileSemanticIntegrityReceipt:
    parent_pr_number: int
    parent_exact_head: str
    parent_receipt_sha256: str
    joined_profile_pr_number: int
    joined_profile_exact_head: str
    joined_profile_source_blob: str
    joined_profile_schema: str
    joined_profile_configuration_sha256: str
    joined_profile_reported_semantic_sha256: str
    configuration_stable_across_authentication: bool
    parent_application_status: str
    parent_selected_route_ready: bool
    receipt_sha256: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            **asdict(self),
            "truth_boundary": truth_boundary(),
        }


def _canonical_json(value: Any) -> str:
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise CurrentPA16ProfileSemanticIntegrityError(
            "joined-profile configuration/report must be canonical JSON data"
        ) from exc


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise CurrentPA16ProfileSemanticIntegrityError(f"{name} must be a mapping")
    return value


def _snapshot_mapping(value: Any, name: str) -> tuple[dict[str, Any], str]:
    mapping = _mapping(value, name)
    text = _canonical_json(mapping)
    snapshot = json.loads(text)
    if not isinstance(snapshot, dict):
        raise CurrentPA16ProfileSemanticIntegrityError(f"{name} must encode a JSON object")
    return snapshot, text


def _require_sha256(value: Any, name: str) -> str:
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        raise CurrentPA16ProfileSemanticIntegrityError(
            f"{name} must be a lowercase hexadecimal SHA-256"
        )
    return value


class _PinnedJoinedProfileBackend:
    """Delegate all parent methods while freezing the authenticated #959 identity/report."""

    def __init__(
        self,
        backend: CurrentPA16ProfileSemanticIntegrityBackend,
        *,
        identity: app.JoinedProfileIdentity,
        report_snapshot: Mapping[str, Any],
        configuration_snapshot: Mapping[str, Any],
    ) -> None:
        self._backend = backend
        self._identity = identity
        self._report_json = _canonical_json(report_snapshot)
        self._configuration_json = _canonical_json(configuration_snapshot)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._backend, name)

    def joined_profile_identity(self) -> app.JoinedProfileIdentity:
        return self._identity

    def joined_profile_report(self) -> Mapping[str, Any]:
        return json.loads(self._report_json)

    def joined_profile_configuration(self) -> Mapping[str, Any]:
        return json.loads(self._configuration_json)


def _authenticate_joined_profile(
    backend: CurrentPA16ProfileSemanticIntegrityBackend,
) -> tuple[app.JoinedProfileIdentity, dict[str, Any], dict[str, Any], str]:
    configuration_before, configuration_before_json = _snapshot_mapping(
        backend.joined_profile_configuration(), "joined_profile_configuration"
    )
    if configuration_before.get("schema") != UPSTREAM_AGENT1_SCHEMA:
        raise CurrentPA16ProfileSemanticIntegrityError(
            "joined-profile configuration schema is not exact #959"
        )

    report_snapshot, _ = _snapshot_mapping(
        backend.joined_profile_report(), "joined_profile_report"
    )
    reported_sha = _require_sha256(
        report_snapshot.get("semantic_sha256"), "joined_profile_report.semantic_sha256"
    )

    # Read the configuration again around the report fetch.  This is not a claim
    # of process-wide immutability; it closes the concrete authentication-time
    # swap seam and the pinned wrapper prevents #961 from refetching the report.
    configuration_after, configuration_after_json = _snapshot_mapping(
        backend.joined_profile_configuration(), "joined_profile_configuration"
    )
    if configuration_before_json != configuration_after_json:
        raise CurrentPA16ProfileSemanticIntegrityError(
            "joined-profile configuration changed during semantic authentication"
        )

    computed_sha = _sha256_text(configuration_before_json)
    if reported_sha != computed_sha:
        raise CurrentPA16ProfileSemanticIntegrityError(
            "joined-profile reported semantic SHA does not match canonical configuration"
        )

    identity = backend.joined_profile_identity()
    if not isinstance(identity, app.JoinedProfileIdentity):
        raise CurrentPA16ProfileSemanticIntegrityError(
            "joined-profile identity must use the typed #961 identity object"
        )
    expected = app.JoinedProfileIdentity(
        pr_number=UPSTREAM_AGENT1_PR,
        exact_head=UPSTREAM_AGENT1_HEAD,
        source_path=UPSTREAM_AGENT1_SOURCE_PATH,
        source_blob=UPSTREAM_AGENT1_SOURCE_BLOB,
        report_schema=UPSTREAM_AGENT1_SCHEMA,
    )
    if identity != expected:
        raise CurrentPA16ProfileSemanticIntegrityError(
            "joined-profile identity is not the exact pinned #959 artifact"
        )

    return identity, report_snapshot, configuration_before, computed_sha


def materialize_semantic_verified_current_pa16_profile_application_receipt(
    backend: CurrentPA16ProfileSemanticIntegrityBackend,
) -> CurrentPA16ProfileSemanticIntegrityReceipt:
    """Replay #959 semantic identity, then consume the existing #961 binder.

    No residual, defect, discrepancy, coefficient, pressure, forcing, viscosity,
    threshold, sample, derivative step, or scientific tuning argument is exposed.
    """
    required = (
        "joined_profile_configuration",
        "joined_profile_identity",
        "joined_profile_report",
        "joined_profile_values",
        "joined_profile_replay_with_coefficients",
    )
    missing = [name for name in required if not callable(getattr(backend, name, None))]
    if missing:
        raise CurrentPA16ProfileSemanticIntegrityError(
            "backend is missing semantic-integrity methods: " + ", ".join(missing)
        )

    identity, report_snapshot, configuration_snapshot, computed_sha = (
        _authenticate_joined_profile(backend)
    )
    pinned = _PinnedJoinedProfileBackend(
        backend,
        identity=identity,
        report_snapshot=report_snapshot,
        configuration_snapshot=configuration_snapshot,
    )
    parent = app.materialize_current_pa16_profile_application_receipt(pinned)
    if not isinstance(parent, app.CurrentPA16ProfileApplicationReceipt):
        raise CurrentPA16ProfileSemanticIntegrityError(
            "#961 did not return the typed profile-application receipt"
        )
    if parent.joined_profile_semantic_sha256 != computed_sha:
        raise CurrentPA16ProfileSemanticIntegrityError(
            "#961 receipt semantic SHA disagrees with authenticated #959 configuration"
        )
    parent_receipt_sha = _require_sha256(parent.receipt_sha256, "#961 receipt_sha256")

    unsigned = {
        "schema": SCHEMA,
        "parent_pr_number": PARENT_AGENT3_PR,
        "parent_exact_head": PARENT_AGENT3_HEAD,
        "parent_receipt_sha256": parent_receipt_sha,
        "joined_profile_pr_number": identity.pr_number,
        "joined_profile_exact_head": identity.exact_head,
        "joined_profile_source_blob": identity.source_blob,
        "joined_profile_schema": identity.report_schema,
        "joined_profile_configuration_sha256": computed_sha,
        "joined_profile_reported_semantic_sha256": computed_sha,
        "configuration_stable_across_authentication": True,
        "parent_application_status": str(parent.status),
        "parent_selected_route_ready": bool(parent.selected_route_ready),
    }
    receipt_sha = _sha256_text(_canonical_json(unsigned))
    return CurrentPA16ProfileSemanticIntegrityReceipt(
        parent_pr_number=PARENT_AGENT3_PR,
        parent_exact_head=PARENT_AGENT3_HEAD,
        parent_receipt_sha256=parent_receipt_sha,
        joined_profile_pr_number=identity.pr_number,
        joined_profile_exact_head=identity.exact_head,
        joined_profile_source_blob=identity.source_blob,
        joined_profile_schema=identity.report_schema,
        joined_profile_configuration_sha256=computed_sha,
        joined_profile_reported_semantic_sha256=computed_sha,
        configuration_stable_across_authentication=True,
        parent_application_status=str(parent.status),
        parent_selected_route_ready=bool(parent.selected_route_ready),
        receipt_sha256=receipt_sha,
    )


def truth_boundary() -> dict[str, bool]:
    return {
        "joined_profile_configuration_semantic_sha_replayed": True,
        "joined_profile_report_pinned_during_parent_delegation": True,
        "parent_961_application_receipt_consumed": True,
        "parent_961_entrypoint_itself_modified": False,
        "runtime_evaluator_process_immutability_proved": False,
        "source_coordinate_repair_application_is_cartesian_velocity": False,
        "source_T_sh_lower_bound_verified": False,
        "source_global_inner_to_outer_join_admitted": False,
        "global_cartesian_leading_velocity_materialized": False,
        "current_real_ns_correction_velocity_materialized": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "discrepancy_from_complete_ns_defect": False,
        "authorized_for_gain_gated_ns_stage": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
