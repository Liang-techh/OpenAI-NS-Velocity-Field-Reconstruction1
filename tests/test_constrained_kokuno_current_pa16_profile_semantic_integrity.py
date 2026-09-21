from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import openai_ns_reconstruction.kokuno_current_pa16_profile_application as app
import openai_ns_reconstruction.kokuno_current_pa16_profile_semantic_integrity as guard


def _canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(value) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _configuration(*, quadrature_points: int = 96) -> dict:
    return {
        "schema": app.UPSTREAM_AGENT1_SCHEMA,
        "certificate": {
            "schema": "focused-test-certificate",
            "selected_numerical_B0": 1.25,
            "selected_T_sh": 2.0,
        },
        "quadrature_points": quadrature_points,
        "absolute_tolerance": 2.0e-11,
    }


def _identity() -> app.JoinedProfileIdentity:
    return app.JoinedProfileIdentity(
        pr_number=app.UPSTREAM_AGENT1_PR,
        exact_head=app.UPSTREAM_AGENT1_HEAD,
        source_path=app.UPSTREAM_AGENT1_SOURCE_PATH,
        source_blob=app.UPSTREAM_AGENT1_SOURCE_BLOB,
        report_schema=app.UPSTREAM_AGENT1_SCHEMA,
    )


def _parent_receipt(semantic_sha: str) -> app.CurrentPA16ProfileApplicationReceipt:
    return app.CurrentPA16ProfileApplicationReceipt(
        upstream=_identity(),
        parent_repair_receipt_sha256="1" * 64,
        parent_tsh_semantic_sha256="2" * 64,
        joined_profile_semantic_sha256=semantic_sha,
        selected_route_ready=False,
        selected_T_sh=2.0,
        log_X_R=10.0,
        log_X_h=5.0,
        status="blocked_by_current_tsh_geometry",
        blocked_reason="focused semantic-integrity test",
        samples=(),
        max_application_abs=0.0,
        max_application_rel=0.0,
        max_zero_coefficient_effect_l2=0.0,
        any_nontrivial_repair_effect=False,
        receipt_sha256="c" * 64,
    )


class FakeBackend:
    def __init__(
        self,
        *,
        configurations: list[dict] | None = None,
        reported_sha: str | None = None,
        drift_report_on_second_fetch: bool = False,
    ) -> None:
        self.configurations = configurations or [_configuration()]
        self.configuration_calls = 0
        self.report_calls = 0
        self.reported_sha = reported_sha or _sha(self.configurations[0])
        self.drift_report_on_second_fetch = drift_report_on_second_fetch

    def joined_profile_configuration(self):
        index = min(self.configuration_calls, len(self.configurations) - 1)
        self.configuration_calls += 1
        return json.loads(_canonical(self.configurations[index]))

    def joined_profile_identity(self):
        return _identity()

    def joined_profile_report(self):
        self.report_calls += 1
        semantic = self.reported_sha
        if self.drift_report_on_second_fetch and self.report_calls > 1:
            semantic = "d" * 64
        return {
            "semantic_sha256": semantic,
            "source": {"repository": "focused-test"},
            "truth_boundary": {"pde_validated": False},
        }

    # Existing #961 surface; focused tests monkeypatch the parent materializer.
    def joined_profile_values(self, eta, X):
        raise AssertionError("profile evaluation is outside this focused guard test")

    def joined_profile_replay_with_coefficients(self, eta, X, coefficients):
        raise AssertionError("profile replay is outside this focused guard test")


def test_replays_configuration_sha_and_pins_report_during_parent_delegation(monkeypatch):
    backend = FakeBackend(drift_report_on_second_fetch=True)
    expected_sha = _sha(_configuration())

    def fake_parent(pinned_backend):
        # A second call by #961 must hit the pinned snapshot, not the backend's
        # deliberately drifting second report.
        assert pinned_backend.joined_profile_report()["semantic_sha256"] == expected_sha
        return _parent_receipt(expected_sha)

    monkeypatch.setattr(app, "materialize_current_pa16_profile_application_receipt", fake_parent)
    receipt = guard.materialize_semantic_verified_current_pa16_profile_application_receipt(
        backend
    )

    assert receipt.joined_profile_configuration_sha256 == expected_sha
    assert receipt.joined_profile_reported_semantic_sha256 == expected_sha
    assert receipt.configuration_stable_across_authentication is True
    assert receipt.parent_receipt_sha256 == "c" * 64
    assert backend.report_calls == 1
    assert backend.configuration_calls == 2
    assert receipt.to_dict()["truth_boundary"]["pde_validated"] is False
    assert receipt.to_dict()["truth_boundary"]["parent_961_entrypoint_itself_modified"] is False


def test_stale_report_sha_after_configuration_mutation_fails_before_parent(monkeypatch):
    old = _configuration(quadrature_points=96)
    new = _configuration(quadrature_points=128)
    backend = FakeBackend(configurations=[new], reported_sha=_sha(old))
    called = False

    def fake_parent(_backend):
        nonlocal called
        called = True
        return _parent_receipt(_sha(new))

    monkeypatch.setattr(app, "materialize_current_pa16_profile_application_receipt", fake_parent)
    with pytest.raises(
        guard.CurrentPA16ProfileSemanticIntegrityError,
        match="does not match canonical configuration",
    ):
        guard.materialize_semantic_verified_current_pa16_profile_application_receipt(backend)
    assert called is False


def test_configuration_swap_during_authentication_fails_closed(monkeypatch):
    first = _configuration(quadrature_points=96)
    second = _configuration(quadrature_points=128)
    backend = FakeBackend(configurations=[first, second], reported_sha=_sha(first))
    monkeypatch.setattr(
        app,
        "materialize_current_pa16_profile_application_receipt",
        lambda _backend: pytest.fail("parent must not run after configuration drift"),
    )
    with pytest.raises(
        guard.CurrentPA16ProfileSemanticIntegrityError,
        match="changed during semantic authentication",
    ):
        guard.materialize_semantic_verified_current_pa16_profile_application_receipt(backend)


def test_parent_961_receipt_must_carry_the_authenticated_semantic_sha(monkeypatch):
    backend = FakeBackend()
    monkeypatch.setattr(
        app,
        "materialize_current_pa16_profile_application_receipt",
        lambda _backend: _parent_receipt("0" * 64),
    )
    with pytest.raises(
        guard.CurrentPA16ProfileSemanticIntegrityError,
        match="#961 receipt semantic SHA disagrees",
    ):
        guard.materialize_semantic_verified_current_pa16_profile_application_receipt(backend)


def test_configuration_schema_is_exact_959(monkeypatch):
    bad = _configuration()
    bad["schema"] = "lookalike-joined-profile-v0"
    backend = FakeBackend(configurations=[bad], reported_sha=_sha(bad))
    monkeypatch.setattr(
        app,
        "materialize_current_pa16_profile_application_receipt",
        lambda _backend: pytest.fail("parent must not run for wrong schema"),
    )
    with pytest.raises(
        guard.CurrentPA16ProfileSemanticIntegrityError,
        match="configuration schema is not exact #959",
    ):
        guard.materialize_semantic_verified_current_pa16_profile_application_receipt(backend)


def test_scope_contract_preserves_canonical_cr001_and_delivery_boundaries():
    root = Path(__file__).resolve().parents[1]
    scope = json.loads(
        (root / "configs/kokuno_current_pa16_profile_semantic_integrity_scope.json").read_text()
    )
    constraints = json.loads((root / "configs/constraints.json").read_text())
    frozen = scope["canonical_cr001_binding"]

    assert constraints["nu"] == frozen["nu"] == 0.01
    assert constraints["domain"]["physical"] == frozen["physical_domain"] == "R^3"
    assert constraints["domain"]["evaluation_box"] == frozen["evaluation_box"]
    assert constraints["domain"]["support"] == frozen["support"]
    assert constraints["domain"]["time_interval"] == frozen["time_interval"]
    assert constraints["forcing"]["mode"] == frozen["forcing_mode"]
    assert constraints["validation"]["seed"] == frozen["validation_seed"] == 914027
    assert constraints["validation"]["held_out_points"] == frozen["held_out_points"] == 4096
    assert constraints["validation"]["derivative_steps"] == frozen["derivative_steps"]
    assert constraints["validation"]["quadrature_orders_per_axis"] == frozen[
        "quadrature_orders_per_axis"
    ]
    thresholds = constraints["validation"]["thresholds"]
    assert thresholds["pde_residual_max"] == frozen["pde_residual_max"] == 1e-3
    assert thresholds["pde_residual_L2"] == frozen["pde_residual_L2"] == 1e-3
    assert thresholds["divergence_max"] == frozen["divergence_max"] == 1e-5
    assert thresholds["divergence_L2"] == frozen["divergence_L2"] == 1e-5
    assert frozen["free_residual_defined_forcing_allowed"] is False
    assert frozen["collapsed_candidate_allowed"] is False
    assert frozen["posthoc_threshold_relaxation_allowed"] is False

    delivery = scope["delivery_state_boundary"]
    assert delivery["canonical_eq45_velocity_export_ready"] is True
    assert delivery["kokuno_global_cartesian_velocity_materialized"] is False
    assert delivery["visual_correspondence_verified"] is False
    assert delivery["pde_validated"] is False
    assert delivery["paper_exact"] is False
    assert delivery["openai_field_identified"] is False
