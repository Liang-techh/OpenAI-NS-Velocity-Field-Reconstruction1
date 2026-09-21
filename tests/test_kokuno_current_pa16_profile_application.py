from __future__ import annotations

from dataclasses import replace
import math

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_current_pa16_profile_application as app
from openai_ns_reconstruction.kokuno_current_pa16_repair_receipt import (
    CurrentPA16RepairReceipt,
    PA16RepairSample,
    UpstreamPA16Identity,
)


COEFFICIENTS = (0.20, -0.10, 0.05, 0.03, -0.02)
TSH_SHA = "a" * 64
PROFILE_SHA = "b" * 64


def _upstream_pa16_identity() -> UpstreamPA16Identity:
    return UpstreamPA16Identity(
        pr_number=953,
        exact_head="bf16b88db92d23923b3e049c00e000a5bfa09992",
        parent_pr_number=947,
        parent_exact_head="2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b",
        tsh_source_path="src/openai_ns_reconstruction/kokuno_pa16_current_tsh_certificate.py",
        tsh_source_blob="5fd610b417ccaa1c4673ea64985296b94dfe5946",
        join_source_path="src/openai_ns_reconstruction/kokuno_inner_join_exit.py",
        join_source_blob="baa283c9592e19653e20f5c24371230d680704f5",
    )


def _repair_sample() -> PA16RepairSample:
    incoming = (0.25, -0.10, 0.08, 0.12, -0.04)
    pre = (0.20, -0.08, 0.05, 0.10, -0.03)
    achieved = (-0.19, 0.075, -0.048, -0.096, 0.028)
    exit_error = tuple(float(a + b) for a, b in zip(pre, achieved))
    eta = 0.25
    pre_t = (
        pre[0],
        pre[2] - 4.0 * eta * pre[1],
        pre[1],
        pre[3] - 8.0 * eta * pre[0],
        pre[4],
    )
    exit_t = (
        exit_error[0],
        exit_error[2] - 4.0 * eta * exit_error[1],
        exit_error[1],
        exit_error[3] - 8.0 * eta * exit_error[0],
        exit_error[4],
    )
    pre_norm = float(np.linalg.norm(pre))
    exit_norm = float(np.linalg.norm(exit_error))
    pre_t_norm = float(np.linalg.norm(pre_t))
    exit_t_norm = float(np.linalg.norm(exit_t))
    return PA16RepairSample(
        eta=eta,
        incoming_scaled_discrepancy=incoming,
        pre_repair_scaled_discrepancy=pre,
        repair_achieved_scaled_moments=achieved,
        exit_scaled_discrepancy=exit_error,
        pa16_transformed_pre_discrepancy=pre_t,
        pa16_transformed_exit_discrepancy=exit_t,
        coefficient_vector=COEFFICIENTS,
        incoming_alignment_abs_max=0.0,
        repair_algebra_abs_max=0.0,
        target_alignment_abs_max=0.0,
        solver_residual_alignment_abs_max=0.0,
        pre_raw_l2_norm=pre_norm,
        exit_raw_l2_norm=exit_norm,
        pre_transformed_l2_norm=pre_t_norm,
        exit_transformed_l2_norm=exit_t_norm,
        raw_contraction_ratio=exit_norm / pre_norm,
        transformed_contraction_ratio=exit_t_norm / pre_t_norm,
        coefficient_l2_norm=float(np.linalg.norm(COEFFICIENTS)),
        coefficient_max_abs=max(abs(value) for value in COEFFICIENTS),
        coefficient_nontrivial=True,
        max_abs_exit_discrepancy=max(abs(value) for value in exit_error),
        scaled_jacobian_condition=3.0,
    )


def _repair_receipt(*, ready: bool = True) -> CurrentPA16RepairReceipt:
    samples = (_repair_sample(),) if ready else ()
    return CurrentPA16RepairReceipt(
        upstream=_upstream_pa16_identity(),
        xi_handoff_sha256="1" * 64,
        xi_upstream_semantic_sha256="2" * 64,
        tsh_semantic_sha256=TSH_SHA,
        selected_route_ready=ready,
        selected_T_sh=2.0,
        separation_geometry_feasible=ready,
        geometry_margin=0.25 if ready else -0.25,
        status=(
            "current_source_moment_repair_verified"
            if ready
            else "blocked_by_current_tsh_geometry"
        ),
        blocked_reason=None if ready else "selected numerical T_sh route is not executable",
        samples=samples,
        max_incoming_alignment_abs=0.0,
        max_repair_algebra_abs=0.0,
        max_exit_raw_l2_norm=max((sample.exit_raw_l2_norm for sample in samples), default=0.0),
        max_exit_transformed_l2_norm=max(
            (sample.exit_transformed_l2_norm for sample in samples), default=0.0
        ),
        any_nontrivial_correction=any(sample.coefficient_nontrivial for sample in samples),
        receipt_sha256="3" * 64,
    )


def _truth_boundary() -> dict[str, bool]:
    return {
        "current_lineage_prefix_0_to_Xi_bound": True,
        "current_candidate_side_join_Xi_to_Xh_materialized": True,
        "current_candidate_side_PA16_coefficients_consumed": True,
        "physical_X_F_U_E_profile_executable_through_Xh": True,
        "source_T_sh_lower_bound_verified": False,
        "source_global_inner_to_outer_join_admitted": False,
        "outer_global_leading_velocity_materialized": False,
        "cartesian_velocity_materialized": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "heldout_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


class FakeBackend:
    def __init__(
        self,
        *,
        ready: bool = True,
        parent_tsh_sha: str = TSH_SHA,
        identity: app.JoinedProfileIdentity | None = None,
        actual_coefficient_vector=COEFFICIENTS,
        zero_has_same_effect: bool = False,
    ):
        self.ready = ready
        self.parent_tsh_sha = parent_tsh_sha
        self.identity = identity or app.JoinedProfileIdentity(
            pr_number=app.UPSTREAM_AGENT1_PR,
            exact_head=app.UPSTREAM_AGENT1_HEAD,
            source_path=app.UPSTREAM_AGENT1_SOURCE_PATH,
            source_blob=app.UPSTREAM_AGENT1_SOURCE_BLOB,
            report_schema=app.UPSTREAM_AGENT1_SCHEMA,
        )
        self.actual_coefficients = tuple(float(v) for v in actual_coefficient_vector)
        self.zero_has_same_effect = zero_has_same_effect
        self.profile_calls = 0
        self.replay_calls = 0

    # Parent #958/#957 methods are deliberately dummies in focused tests because
    # materialize_checksum_verified_current_pa16_repair_receipt is monkeypatched.
    def current_xi_bridge_receipt(self):
        raise AssertionError("parent repair path should be monkeypatched in focused test")

    def upstream_identity(self):
        raise AssertionError("parent repair path should be monkeypatched in focused test")

    def tsh_report(self):
        raise AssertionError("parent repair path should be monkeypatched in focused test")

    def pa16_solve_at_eta(self, eta: float):
        raise AssertionError("parent repair path should be monkeypatched in focused test")

    def joined_profile_identity(self):
        return self.identity

    def joined_profile_report(self):
        return {
            "source": {
                "repository": app.SOURCE_READER_REPO,
                "commit": app.SOURCE_READER_HEAD,
                "path": app.SOURCE_READER_PATH,
            },
            "parent_certificate_semantic_sha256": self.parent_tsh_sha,
            "geometry": {
                "route_ready": self.ready,
                "selected_T_sh": 2.0,
                "log_X_R": 10.0,
                "log_X_h": 5.0,
            },
            "truth_boundary": _truth_boundary(),
            "semantic_sha256": PROFILE_SHA,
        }

    @staticmethod
    def _values(eta: float, X: tuple[float, ...], coefficients) -> dict[str, np.ndarray]:
        X_array = np.asarray(X, dtype=float)
        c = np.asarray(coefficients, dtype=float)
        # Five linearly independent coefficient signatures over the repair probes.
        logx = np.log(X_array) - 10.0
        basis = np.stack(
            [
                np.ones_like(logx),
                logx + 5.5,
                (logx + 5.5) ** 2,
                np.sin(logx + 6.0),
                np.cos(logx + 6.0),
            ],
            axis=-1,
        )
        effect = basis @ c
        F = 1.0 + 0.05 * eta + 0.02 * effect
        U = 0.3 + eta + 0.04 * effect
        E = 0.7 + 0.03 * eta - 0.01 * effect
        return {"F": F, "U": U, "E": E}

    def joined_profile_values(self, eta: float, X: tuple[float, ...]):
        self.profile_calls += 1
        return self._values(eta, X, self.actual_coefficients)

    def joined_profile_replay_with_coefficients(self, eta, X, coefficients):
        self.replay_calls += 1
        effective = coefficients
        if self.zero_has_same_effect and not any(float(v) != 0.0 for v in coefficients):
            effective = self.actual_coefficients
        return self._values(float(eta), X, effective)


def _patch_parent(monkeypatch, receipt: CurrentPA16RepairReceipt):
    monkeypatch.setattr(
        app,
        "materialize_checksum_verified_current_pa16_repair_receipt",
        lambda backend: receipt,
    )


def test_ready_profile_is_bound_to_exact_repair_coefficients(monkeypatch):
    repair = _repair_receipt(ready=True)
    _patch_parent(monkeypatch, repair)
    backend = FakeBackend()
    receipt = app.materialize_current_pa16_profile_application_receipt(backend)

    assert receipt.status == "current_source_moment_repair_applied_in_joined_profile"
    assert receipt.selected_route_ready is True
    assert len(receipt.samples) == 1
    assert receipt.samples[0].coefficient_vector == COEFFICIENTS
    assert receipt.samples[0].application_abs_max == pytest.approx(0.0)
    assert receipt.samples[0].application_rel_max == pytest.approx(0.0)
    assert receipt.samples[0].zero_coefficient_effect_l2 > 0.0
    assert receipt.samples[0].nontrivial_repair_effect_observed is True
    assert receipt.any_nontrivial_repair_effect is True
    assert backend.profile_calls == 1
    assert backend.replay_calls == 2

    truth = receipt.to_dict()["truth_boundary"]
    assert truth["source_coordinate_repair_application_verified"] is True
    assert truth["current_real_ns_correction_velocity_materialized"] is False
    assert truth["discrepancy_from_complete_ns_defect"] is False
    assert truth["heldout_normalized_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_joined_profile_must_share_exact_tsh_certificate_identity(monkeypatch):
    _patch_parent(monkeypatch, _repair_receipt(ready=True))
    backend = FakeBackend(parent_tsh_sha="f" * 64)
    with pytest.raises(app.CurrentPA16ProfileApplicationError, match="same #953 T_sh certificate"):
        app.materialize_current_pa16_profile_application_receipt(backend)


def test_wrong_applied_coefficients_fail_closed(monkeypatch):
    _patch_parent(monkeypatch, _repair_receipt(ready=True))
    backend = FakeBackend(actual_coefficient_vector=(0.18, -0.10, 0.05, 0.03, -0.02))
    with pytest.raises(app.CurrentPA16ProfileApplicationError, match="does not reproduce"):
        app.materialize_current_pa16_profile_application_receipt(backend)


def test_nontrivial_receipt_cannot_be_ignored_by_joined_profile(monkeypatch):
    _patch_parent(monkeypatch, _repair_receipt(ready=True))
    backend = FakeBackend(zero_has_same_effect=True)
    with pytest.raises(app.CurrentPA16ProfileApplicationError, match="no observable effect"):
        app.materialize_current_pa16_profile_application_receipt(backend)


def test_infeasible_tsh_route_remains_blocked_and_never_evaluates_profile(monkeypatch):
    repair = _repair_receipt(ready=False)
    _patch_parent(monkeypatch, repair)
    backend = FakeBackend(ready=False)
    receipt = app.materialize_current_pa16_profile_application_receipt(backend)
    assert receipt.status == "blocked_by_current_tsh_geometry"
    assert receipt.samples == ()
    assert receipt.any_nontrivial_repair_effect is False
    assert backend.profile_calls == 0
    assert backend.replay_calls == 0


def test_joined_profile_exact_identity_is_fail_closed(monkeypatch):
    _patch_parent(monkeypatch, _repair_receipt(ready=True))
    forged = app.JoinedProfileIdentity(
        pr_number=app.UPSTREAM_AGENT1_PR,
        exact_head="0" * 40,
        source_path=app.UPSTREAM_AGENT1_SOURCE_PATH,
        source_blob=app.UPSTREAM_AGENT1_SOURCE_BLOB,
        report_schema=app.UPSTREAM_AGENT1_SCHEMA,
    )
    backend = FakeBackend(identity=forged)
    with pytest.raises(app.CurrentPA16ProfileApplicationError, match="pinned Agent-1 #959"):
        app.materialize_current_pa16_profile_application_receipt(backend)
