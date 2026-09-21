"""Bind the authenticated current PA.16 joined profile to Agent-1 #965 Cartesian leading.

Agent-3 #964 authenticates the exact current joined-profile configuration and the
#961 repair-application receipt.  Agent-1 #965 independently consumes that joined
profile in the governed Kokuno PA.10 Cartesian leading-velocity map through the
current source-coordinate exit ``X_h``.

This module closes only the cross-branch provenance seam: it proves that the
joined configuration authenticated by #964 is exactly the nested ``joined``
configuration carried by the exact #965 Cartesian-leading artifact, and replays
#965's reported semantic SHA from the report fields that define it.

No Cartesian map, PA.16 repair, oscillatory curl, pressure, forcing, residual,
correction velocity, validation sample, scientific threshold, or derivative
parameter is implemented or tunable here.  In particular, this is not a global
velocity completion beyond ``X_h`` and is not Navier--Stokes residual evidence.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import re
from typing import Any, Mapping, Protocol

from . import kokuno_current_pa16_profile_semantic_integrity as parent_guard


TASK = "KOKUNO-A3-CURRENT-CARTESIAN-LEADING-BIND-096"
SCHEMA = "kokuno-a3-current-pa16-cartesian-leading-binding-v1"
PARENT_CR002_PR = 964
PARENT_CR002_HEAD = "6002a3bb0ba24a97cd9f609206a4aa0a276a2486"
PARENT_CR002_SOURCE_PATH = (
    "src/openai_ns_reconstruction/kokuno_current_pa16_profile_semantic_integrity.py"
)
PARENT_CR002_SOURCE_BLOB = "a8e125839f767f665736cb8e5e81ac924887fcb9"
UPSTREAM_AGENT1_PR = 965
UPSTREAM_AGENT1_HEAD = "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"
UPSTREAM_AGENT1_SOURCE_PATH = (
    "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_leading_velocity.py"
)
UPSTREAM_AGENT1_SOURCE_BLOB = "ca8b80b0451be1a8f31deaf620f542d8c2e92c0c"
UPSTREAM_AGENT1_SCHEMA = "kokuno-pa16-current-cartesian-leading-velocity-v1"
KOKUNO_PUBLIC_SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

_EXPECTED_AGENT1_TRUTH = {
    "public_reconstruction_source": True,
    "governed_source_coordinate_inverse_reused": True,
    "current_joined_F_U_E_consumed": True,
    "current_joined_incompressibility_primitive_materialized": True,
    "current_joined_radial_v0_materialized": True,
    "current_joined_cartesian_spacetime_leading_velocity_through_Xh_materialized": True,
    "axis_regular_cartesian_formula_executable": True,
    "velocity_interface_vectorized": True,
    "velocity_configuration_serializable": True,
    "eta_derivative_is_repository_numerical_realization": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_prepared_appendixA_pressure_stress_materialized": False,
    "source_admitted_kappa0_materialized": False,
    "source_global_inner_to_outer_join_admitted": False,
    "outer_global_leading_velocity_materialized": False,
    "velocity_beyond_Xh_materialized": False,
    "unified_global_cartesian_velocity_export_ready": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


class CurrentPA16CartesianLeadingBindingError(RuntimeError):
    """Raised when #964 and exact #965 cannot be bound to one immutable lineage."""


@dataclass(frozen=True)
class CartesianLeadingIdentity:
    pr_number: int
    exact_head: str
    source_path: str
    source_blob: str
    report_schema: str


class CurrentPA16CartesianLeadingBindingBackend(
    parent_guard.CurrentPA16ProfileSemanticIntegrityBackend, Protocol
):
    def cartesian_leading_identity(self) -> CartesianLeadingIdentity: ...

    def cartesian_leading_configuration(self) -> Mapping[str, Any]: ...

    def cartesian_leading_report(self) -> Mapping[str, Any]: ...


@dataclass(frozen=True)
class CurrentPA16CartesianLeadingBindingReceipt:
    parent_cr002_pr: int
    parent_cr002_exact_head: str
    parent_receipt_sha256: str
    agent1_pr: int
    agent1_exact_head: str
    agent1_source_blob: str
    agent1_report_schema: str
    joined_configuration_sha256: str
    joined_reported_semantic_sha256: str
    cartesian_configuration_sha256: str
    cartesian_reported_semantic_sha256: str
    configuration_and_report_stable_across_binding: bool
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
        raise CurrentPA16CartesianLeadingBindingError(
            "Cartesian-leading configuration/report must be canonical JSON data"
        ) from exc


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _snapshot_mapping(value: Any, name: str) -> tuple[dict[str, Any], str]:
    if not isinstance(value, Mapping):
        raise CurrentPA16CartesianLeadingBindingError(f"{name} must be a mapping")
    text = _canonical_json(value)
    snapshot = json.loads(text)
    if not isinstance(snapshot, dict):
        raise CurrentPA16CartesianLeadingBindingError(f"{name} must encode a JSON object")
    return snapshot, text


def _require_sha256(value: Any, name: str) -> str:
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        raise CurrentPA16CartesianLeadingBindingError(
            f"{name} must be a lowercase hexadecimal SHA-256"
        )
    return value


def _expected_identity() -> CartesianLeadingIdentity:
    return CartesianLeadingIdentity(
        pr_number=UPSTREAM_AGENT1_PR,
        exact_head=UPSTREAM_AGENT1_HEAD,
        source_path=UPSTREAM_AGENT1_SOURCE_PATH,
        source_blob=UPSTREAM_AGENT1_SOURCE_BLOB,
        report_schema=UPSTREAM_AGENT1_SCHEMA,
    )


def _replay_agent1_semantic(report: Mapping[str, Any]) -> str:
    source = report.get("source")
    if not isinstance(source, Mapping):
        raise CurrentPA16CartesianLeadingBindingError("#965 report.source must be a mapping")
    source_commit = source.get("commit")
    if source_commit != KOKUNO_PUBLIC_SOURCE_COMMIT:
        raise CurrentPA16CartesianLeadingBindingError(
            "#965 report is not pinned to the corrected public Kokuno source commit"
        )
    source_formulas = report.get("source_formulas")
    numerical_realization = report.get("numerical_realization")
    configuration = report.get("configuration")
    if not isinstance(source_formulas, Mapping):
        raise CurrentPA16CartesianLeadingBindingError(
            "#965 report.source_formulas must be a mapping"
        )
    if not isinstance(numerical_realization, Mapping):
        raise CurrentPA16CartesianLeadingBindingError(
            "#965 report.numerical_realization must be a mapping"
        )
    if not isinstance(configuration, Mapping):
        raise CurrentPA16CartesianLeadingBindingError(
            "#965 report.configuration must be a mapping"
        )
    payload = {
        "schema": UPSTREAM_AGENT1_SCHEMA,
        "source_commit": source_commit,
        "source_formulas": source_formulas,
        "numerical_realization": numerical_realization,
        "configuration": configuration,
    }
    return _sha256_text(_canonical_json(payload))


def _authenticate_agent1_cartesian(
    backend: CurrentPA16CartesianLeadingBindingBackend,
    parent: parent_guard.CurrentPA16ProfileSemanticIntegrityReceipt,
) -> tuple[CartesianLeadingIdentity, dict[str, Any], dict[str, Any], str, str]:
    identity = backend.cartesian_leading_identity()
    if not isinstance(identity, CartesianLeadingIdentity):
        raise CurrentPA16CartesianLeadingBindingError(
            "Cartesian-leading identity must use the typed Agent-3 identity object"
        )
    if identity != _expected_identity():
        raise CurrentPA16CartesianLeadingBindingError(
            "Cartesian-leading identity is not exact Agent-1 #965"
        )

    configuration_before, configuration_before_json = _snapshot_mapping(
        backend.cartesian_leading_configuration(), "cartesian_leading_configuration"
    )
    report_before, report_before_json = _snapshot_mapping(
        backend.cartesian_leading_report(), "cartesian_leading_report"
    )
    configuration_after, configuration_after_json = _snapshot_mapping(
        backend.cartesian_leading_configuration(), "cartesian_leading_configuration"
    )
    report_after, report_after_json = _snapshot_mapping(
        backend.cartesian_leading_report(), "cartesian_leading_report"
    )
    if configuration_before_json != configuration_after_json:
        raise CurrentPA16CartesianLeadingBindingError(
            "Cartesian-leading configuration changed during binding"
        )
    if report_before_json != report_after_json:
        raise CurrentPA16CartesianLeadingBindingError(
            "Cartesian-leading report changed during binding"
        )

    if configuration_before.get("schema") != UPSTREAM_AGENT1_SCHEMA:
        raise CurrentPA16CartesianLeadingBindingError(
            "Cartesian-leading configuration schema is not exact #965"
        )
    if report_before.get("configuration") != configuration_before:
        raise CurrentPA16CartesianLeadingBindingError(
            "#965 report.configuration disagrees with the runtime configuration"
        )

    joined_configuration = configuration_before.get("joined")
    if not isinstance(joined_configuration, Mapping):
        raise CurrentPA16CartesianLeadingBindingError(
            "#965 configuration.joined must be a mapping"
        )
    joined_sha = _sha256_text(_canonical_json(joined_configuration))
    if joined_sha != parent.joined_profile_configuration_sha256:
        raise CurrentPA16CartesianLeadingBindingError(
            "#965 nested joined configuration is not the joined profile authenticated by #964"
        )
    parent_joined_semantic = _require_sha256(
        report_before.get("parent_joined_semantic_sha256"),
        "#965 report.parent_joined_semantic_sha256",
    )
    if parent_joined_semantic != parent.joined_profile_reported_semantic_sha256:
        raise CurrentPA16CartesianLeadingBindingError(
            "#965 parent joined semantic SHA disagrees with #964 authentication"
        )
    if parent_joined_semantic != joined_sha:
        raise CurrentPA16CartesianLeadingBindingError(
            "#965 nested joined configuration and reported parent semantic SHA disagree"
        )

    truth = report_before.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise CurrentPA16CartesianLeadingBindingError("#965 truth_boundary must be a mapping")
    for key, expected in _EXPECTED_AGENT1_TRUTH.items():
        if truth.get(key) is not expected:
            raise CurrentPA16CartesianLeadingBindingError(
                f"#965 truth boundary mismatch for {key}"
            )

    reported_semantic = _require_sha256(
        report_before.get("semantic_sha256"), "#965 report.semantic_sha256"
    )
    replayed_semantic = _replay_agent1_semantic(report_before)
    if reported_semantic != replayed_semantic:
        raise CurrentPA16CartesianLeadingBindingError(
            "#965 reported semantic SHA does not match replayed report semantics"
        )

    configuration_sha = _sha256_text(configuration_before_json)
    return (
        identity,
        configuration_before,
        report_before,
        configuration_sha,
        replayed_semantic,
    )


def materialize_semantic_verified_current_pa16_cartesian_leading_receipt(
    backend: CurrentPA16CartesianLeadingBindingBackend,
) -> CurrentPA16CartesianLeadingBindingReceipt:
    """Bind #964's authenticated repair lineage to exact #965 through ``X_h``.

    The only public argument is a typed backend.  No residual, defect,
    discrepancy, coefficient, correction, pressure, forcing, gain, damping,
    viscosity, derivative step, validation sample, stage budget, or scientific
    threshold is accepted.
    """
    required = (
        "cartesian_leading_identity",
        "cartesian_leading_configuration",
        "cartesian_leading_report",
    )
    missing = [name for name in required if not callable(getattr(backend, name, None))]
    if missing:
        raise CurrentPA16CartesianLeadingBindingError(
            "backend is missing Cartesian-leading binding methods: " + ", ".join(missing)
        )

    parent = parent_guard.materialize_semantic_verified_current_pa16_profile_application_receipt(
        backend
    )
    if not isinstance(parent, parent_guard.CurrentPA16ProfileSemanticIntegrityReceipt):
        raise CurrentPA16CartesianLeadingBindingError(
            "#964 did not return the typed semantic-integrity receipt"
        )
    parent_receipt_sha = _require_sha256(parent.receipt_sha256, "#964 receipt_sha256")

    identity, configuration, report, configuration_sha, semantic_sha = (
        _authenticate_agent1_cartesian(backend, parent)
    )
    joined_sha = _sha256_text(_canonical_json(configuration["joined"]))
    joined_reported_semantic = _require_sha256(
        report["parent_joined_semantic_sha256"],
        "#965 report.parent_joined_semantic_sha256",
    )

    unsigned = {
        "schema": SCHEMA,
        "parent_cr002_pr": PARENT_CR002_PR,
        "parent_cr002_exact_head": PARENT_CR002_HEAD,
        "parent_receipt_sha256": parent_receipt_sha,
        "agent1_pr": identity.pr_number,
        "agent1_exact_head": identity.exact_head,
        "agent1_source_blob": identity.source_blob,
        "agent1_report_schema": identity.report_schema,
        "joined_configuration_sha256": joined_sha,
        "joined_reported_semantic_sha256": joined_reported_semantic,
        "cartesian_configuration_sha256": configuration_sha,
        "cartesian_reported_semantic_sha256": semantic_sha,
        "configuration_and_report_stable_across_binding": True,
    }
    receipt_sha = _sha256_text(_canonical_json(unsigned))
    return CurrentPA16CartesianLeadingBindingReceipt(
        parent_cr002_pr=PARENT_CR002_PR,
        parent_cr002_exact_head=PARENT_CR002_HEAD,
        parent_receipt_sha256=parent_receipt_sha,
        agent1_pr=identity.pr_number,
        agent1_exact_head=identity.exact_head,
        agent1_source_blob=identity.source_blob,
        agent1_report_schema=identity.report_schema,
        joined_configuration_sha256=joined_sha,
        joined_reported_semantic_sha256=joined_reported_semantic,
        cartesian_configuration_sha256=configuration_sha,
        cartesian_reported_semantic_sha256=semantic_sha,
        configuration_and_report_stable_across_binding=True,
        receipt_sha256=receipt_sha,
    )


def truth_boundary() -> dict[str, bool]:
    return {
        "parent_964_semantic_integrity_receipt_consumed": True,
        "exact_agent1_965_identity_pinned": True,
        "authenticated_joined_profile_bound_to_cartesian_leading_through_Xh": True,
        "agent1_965_semantic_sha_replayed": True,
        "current_joined_cartesian_spacetime_leading_velocity_through_Xh_materialized": True,
        "outer_global_leading_velocity_materialized": False,
        "velocity_beyond_Xh_materialized": False,
        "unified_global_cartesian_velocity_export_ready": False,
        "current_real_ns_correction_velocity_materialized": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "discrepancy_from_complete_ns_defect": False,
        "authorized_for_gain_gated_ns_stage": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
