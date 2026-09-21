from __future__ import annotations

import hashlib
import inspect
import json

import pytest

import openai_ns_reconstruction.kokuno_current_pa16_cartesian_leading_binding as binding
import openai_ns_reconstruction.kokuno_current_pa16_profile_semantic_integrity as parent_guard


def _canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(value) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _joined_configuration(*, quadrature_points: int = 96) -> dict:
    return {
        "schema": "kokuno-pa16-current-joined-profile-v1",
        "certificate": {"schema": "focused-current-pa16-certificate-v1"},
        "quadrature_points": quadrature_points,
        "absolute_tolerance": 2.0e-11,
    }


def _cartesian_configuration(*, joined: dict | None = None, eta_fd_step: float = 2.0e-5) -> dict:
    return {
        "schema": binding.UPSTREAM_AGENT1_SCHEMA,
        "joined": joined or _joined_configuration(),
        "geometry": {"schema": "focused-governed-pa10-cartesian-geometry-v1"},
        "moment_quadrature_order": 48,
        "eta_fd_step": eta_fd_step,
    }


def _cartesian_report(configuration: dict, *, truth_updates: dict | None = None) -> dict:
    truth = dict(binding._EXPECTED_AGENT1_TRUTH)
    if truth_updates:
        truth.update(truth_updates)
    report = {
        "source": {
            "repository": "KokunoYumeto/navier-stokes-3d-axi",
            "commit": binding.KOKUNO_PUBLIC_SOURCE_COMMIT,
            "path": "reconstruction",
            "corrected_release": "corrected",
            "corrected_release_date": "2026-09-09",
        },
        "source_formulas": {
            "coordinates": "focused exact-head test payload",
            "cartesian_velocity": "focused exact-head test payload",
        },
        "numerical_realization": {
            "coordinate_inverse": "reuse governed inverse",
            "M_eta_over_X": "frozen FD4 realization",
        },
        "configuration": json.loads(_canonical(configuration)),
        "parent_joined_semantic_sha256": _sha(configuration["joined"]),
        "governed_geometry_field_sha256": "a" * 64,
        "machine_checks": {"axis_velocity_finite": True},
        "domain": {"X_i": 110.0, "X_h": 120.0},
        "truth_boundary": truth,
    }
    semantic_payload = {
        "schema": binding.UPSTREAM_AGENT1_SCHEMA,
        "source_commit": report["source"]["commit"],
        "source_formulas": report["source_formulas"],
        "numerical_realization": report["numerical_realization"],
        "configuration": report["configuration"],
    }
    report["semantic_sha256"] = _sha(semantic_payload)
    return report


def _identity() -> binding.CartesianLeadingIdentity:
    return binding.CartesianLeadingIdentity(
        pr_number=binding.UPSTREAM_AGENT1_PR,
        exact_head=binding.UPSTREAM_AGENT1_HEAD,
        source_path=binding.UPSTREAM_AGENT1_SOURCE_PATH,
        source_blob=binding.UPSTREAM_AGENT1_SOURCE_BLOB,
        report_schema=binding.UPSTREAM_AGENT1_SCHEMA,
    )


def _parent_receipt(joined: dict) -> parent_guard.CurrentPA16ProfileSemanticIntegrityReceipt:
    semantic = _sha(joined)
    return parent_guard.CurrentPA16ProfileSemanticIntegrityReceipt(
        parent_pr_number=961,
        parent_exact_head="670f7b71d72735d441dcf48f2391c271ee5c163c",
        parent_receipt_sha256="c" * 64,
        joined_profile_pr_number=959,
        joined_profile_exact_head="1" * 40,
        joined_profile_source_blob="2" * 40,
        joined_profile_schema="kokuno-pa16-current-joined-profile-v1",
        joined_profile_configuration_sha256=semantic,
        joined_profile_reported_semantic_sha256=semantic,
        configuration_stable_across_authentication=True,
        parent_application_status="ready",
        parent_selected_route_ready=True,
        receipt_sha256="d" * 64,
    )


class FakeBackend:
    def __init__(
        self,
        *,
        configurations: list[dict] | None = None,
        reports: list[dict] | None = None,
        identity: binding.CartesianLeadingIdentity | None = None,
    ) -> None:
        base = _cartesian_configuration()
        self.configurations = configurations or [base]
        self.reports = reports or [_cartesian_report(self.configurations[0])]
        self.identity = identity or _identity()
        self.configuration_calls = 0
        self.report_calls = 0

    def cartesian_leading_identity(self):
        return self.identity

    def cartesian_leading_configuration(self):
        index = min(self.configuration_calls, len(self.configurations) - 1)
        self.configuration_calls += 1
        return json.loads(_canonical(self.configurations[index]))

    def cartesian_leading_report(self):
        index = min(self.report_calls, len(self.reports) - 1)
        self.report_calls += 1
        return json.loads(_canonical(self.reports[index]))


def _patch_parent(monkeypatch, joined: dict):
    receipt = _parent_receipt(joined)
    monkeypatch.setattr(
        parent_guard,
        "materialize_semantic_verified_current_pa16_profile_application_receipt",
        lambda _backend: receipt,
    )
    return receipt


def test_exact_joined_configuration_is_bound_to_exact_965_cartesian_report(monkeypatch):
    joined = _joined_configuration()
    configuration = _cartesian_configuration(joined=joined)
    report = _cartesian_report(configuration)
    backend = FakeBackend(configurations=[configuration], reports=[report])
    parent = _patch_parent(monkeypatch, joined)

    receipt = binding.materialize_semantic_verified_current_pa16_cartesian_leading_receipt(
        backend
    )

    assert receipt.parent_receipt_sha256 == parent.receipt_sha256
    assert receipt.agent1_exact_head == binding.UPSTREAM_AGENT1_HEAD
    assert receipt.agent1_source_blob == binding.UPSTREAM_AGENT1_SOURCE_BLOB
    assert receipt.joined_configuration_sha256 == _sha(joined)
    assert receipt.joined_reported_semantic_sha256 == _sha(joined)
    assert receipt.cartesian_configuration_sha256 == _sha(configuration)
    assert receipt.cartesian_reported_semantic_sha256 == report["semantic_sha256"]
    assert receipt.configuration_and_report_stable_across_binding is True
    truth = receipt.to_dict()["truth_boundary"]
    assert truth["authenticated_joined_profile_bound_to_cartesian_leading_through_Xh"] is True
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["heldout_normalized_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert backend.configuration_calls == 2
    assert backend.report_calls == 2


def test_nested_joined_configuration_mismatch_fails_closed(monkeypatch):
    authenticated_joined = _joined_configuration(quadrature_points=96)
    different_joined = _joined_configuration(quadrature_points=128)
    configuration = _cartesian_configuration(joined=different_joined)
    backend = FakeBackend(
        configurations=[configuration], reports=[_cartesian_report(configuration)]
    )
    _patch_parent(monkeypatch, authenticated_joined)

    with pytest.raises(
        binding.CurrentPA16CartesianLeadingBindingError,
        match="nested joined configuration is not",
    ):
        binding.materialize_semantic_verified_current_pa16_cartesian_leading_receipt(
            backend
        )


def test_parent_joined_semantic_mismatch_fails_closed(monkeypatch):
    joined = _joined_configuration()
    configuration = _cartesian_configuration(joined=joined)
    report = _cartesian_report(configuration)
    report["parent_joined_semantic_sha256"] = "e" * 64
    # Recompute #965 semantic: parent_joined_semantic_sha256 is deliberately not
    # part of #965's semantic payload, so this catches the separate lineage seam.
    backend = FakeBackend(configurations=[configuration], reports=[report])
    _patch_parent(monkeypatch, joined)

    with pytest.raises(
        binding.CurrentPA16CartesianLeadingBindingError,
        match="parent joined semantic SHA disagrees",
    ):
        binding.materialize_semantic_verified_current_pa16_cartesian_leading_receipt(
            backend
        )


def test_promoted_scientific_truth_boundary_is_rejected(monkeypatch):
    joined = _joined_configuration()
    configuration = _cartesian_configuration(joined=joined)
    report = _cartesian_report(
        configuration,
        truth_updates={"matched_global_pressure_materialized": True},
    )
    backend = FakeBackend(configurations=[configuration], reports=[report])
    _patch_parent(monkeypatch, joined)

    with pytest.raises(
        binding.CurrentPA16CartesianLeadingBindingError,
        match="truth boundary mismatch for matched_global_pressure_materialized",
    ):
        binding.materialize_semantic_verified_current_pa16_cartesian_leading_receipt(
            backend
        )


def test_wrong_agent1_identity_is_rejected_before_cross_lineage_acceptance(monkeypatch):
    joined = _joined_configuration()
    configuration = _cartesian_configuration(joined=joined)
    bad_identity = binding.CartesianLeadingIdentity(
        pr_number=binding.UPSTREAM_AGENT1_PR,
        exact_head="0" * 40,
        source_path=binding.UPSTREAM_AGENT1_SOURCE_PATH,
        source_blob=binding.UPSTREAM_AGENT1_SOURCE_BLOB,
        report_schema=binding.UPSTREAM_AGENT1_SCHEMA,
    )
    backend = FakeBackend(
        configurations=[configuration],
        reports=[_cartesian_report(configuration)],
        identity=bad_identity,
    )
    _patch_parent(monkeypatch, joined)

    with pytest.raises(
        binding.CurrentPA16CartesianLeadingBindingError,
        match="not exact Agent-1 #965",
    ):
        binding.materialize_semantic_verified_current_pa16_cartesian_leading_receipt(
            backend
        )


def test_configuration_toctou_mutation_is_rejected(monkeypatch):
    joined = _joined_configuration()
    first = _cartesian_configuration(joined=joined, eta_fd_step=2.0e-5)
    second = _cartesian_configuration(joined=joined, eta_fd_step=3.0e-5)
    backend = FakeBackend(
        configurations=[first, second],
        reports=[_cartesian_report(first), _cartesian_report(first)],
    )
    _patch_parent(monkeypatch, joined)

    with pytest.raises(
        binding.CurrentPA16CartesianLeadingBindingError,
        match="configuration changed during binding",
    ):
        binding.materialize_semantic_verified_current_pa16_cartesian_leading_receipt(
            backend
        )


def test_report_toctou_mutation_is_rejected(monkeypatch):
    joined = _joined_configuration()
    configuration = _cartesian_configuration(joined=joined)
    first = _cartesian_report(configuration)
    second = _cartesian_report(configuration)
    second["machine_checks"]["axis_velocity_finite"] = False
    backend = FakeBackend(configurations=[configuration], reports=[first, second])
    _patch_parent(monkeypatch, joined)

    with pytest.raises(
        binding.CurrentPA16CartesianLeadingBindingError,
        match="report changed during binding",
    ):
        binding.materialize_semantic_verified_current_pa16_cartesian_leading_receipt(
            backend
        )


def test_report_semantic_replay_rejects_tampering(monkeypatch):
    joined = _joined_configuration()
    configuration = _cartesian_configuration(joined=joined)
    report = _cartesian_report(configuration)
    report["source_formulas"]["coordinates"] = "tampered after semantic hash"
    backend = FakeBackend(configurations=[configuration], reports=[report])
    _patch_parent(monkeypatch, joined)

    with pytest.raises(
        binding.CurrentPA16CartesianLeadingBindingError,
        match="reported semantic SHA does not match replayed",
    ):
        binding.materialize_semantic_verified_current_pa16_cartesian_leading_receipt(
            backend
        )


def test_public_entrypoint_exposes_no_scientific_knobs():
    signature = inspect.signature(
        binding.materialize_semantic_verified_current_pa16_cartesian_leading_receipt
    )
    assert tuple(signature.parameters) == ("backend",)


def test_exact_upstream_source_identity_is_frozen():
    assert binding.PARENT_CR002_PR == 964
    assert binding.PARENT_CR002_HEAD == "6002a3bb0ba24a97cd9f609206a4aa0a276a2486"
    assert binding.PARENT_CR002_SOURCE_BLOB == "a8e125839f767f665736cb8e5e81ac924887fcb9"
    assert binding.UPSTREAM_AGENT1_PR == 965
    assert binding.UPSTREAM_AGENT1_HEAD == "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"
    assert binding.UPSTREAM_AGENT1_SOURCE_BLOB == "ca8b80b0451be1a8f31deaf620f542d8c2e92c0c"
