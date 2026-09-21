from __future__ import annotations

from dataclasses import replace
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_pa16_repair_receipt import (
    CurrentPA16RepairReceiptError,
    UpstreamPA16Identity,
    materialize_current_pa16_repair_receipt,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_current_xi_moment_discrepancy_bridge import (
    CurrentXiMomentBridgeReceipt,
    CurrentXiMomentSample,
    UpstreamA1Identity,
)


XI_UPSTREAM = UpstreamA1Identity(
    pr_number=947,
    exact_head="2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b",
    source_path="src/openai_ns_reconstruction/kokuno_pa10_actual_xi_prefix_moments.py",
    source_blob="92dce65c9a8793e06497a7e7347c832861032238",
    report_schema="kokuno-pa10-actual-xi-prefix-moments-v1",
)
REPAIR_UPSTREAM = UpstreamPA16Identity(
    pr_number=953,
    exact_head="bf16b88db92d23923b3e049c00e000a5bfa09992",
    parent_pr_number=947,
    parent_exact_head="2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b",
    tsh_source_path="src/openai_ns_reconstruction/kokuno_pa16_current_tsh_certificate.py",
    tsh_source_blob="5fd610b417ccaa1c4673ea64985296b94dfe5946",
    join_source_path="src/openai_ns_reconstruction/kokuno_inner_join_exit.py",
    join_source_blob="baa283c9592e19653e20f5c24371230d680704f5",
)


def _bridge() -> CurrentXiMomentBridgeReceipt:
    row = (0.2, -0.1, 0.05, 0.3, -0.04)
    eta = 0.25
    transformed = (0.2, 0.15, -0.1, -0.1, -0.04)
    sample = CurrentXiMomentSample(
        eta=eta,
        ell_i=0.02,
        G_i=1.01,
        incoming_scaled_discrepancy=row,
        pa16_row_transformed_discrepancy=transformed,
        raw_l2_norm=float(np.linalg.norm(row)),
        pa16_transformed_l2_norm=float(np.linalg.norm(transformed)),
        correction_needed=True,
    )
    return CurrentXiMomentBridgeReceipt(
        upstream=XI_UPSTREAM,
        upstream_semantic_sha256="1" * 64,
        upstream_report_sha256="2" * 64,
        outer_schedule_sha256="3" * 64,
        samples=(sample,),
        max_raw_l2_norm=sample.raw_l2_norm,
        max_pa16_transformed_l2_norm=sample.pa16_transformed_l2_norm,
        any_correction_needed=True,
        handoff_sha256="4" * 64,
    )


def _tsh_report(*, ready: bool = True) -> dict:
    return {
        "source": {
            "repository": "KokunoYumeto/yang-mills-interacting-workbench",
            "commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
            "path": "navier-stokes/navier_stokes_workbench.tex",
        },
        "current_parent_semantic_sha256": "1" * 64,
        "semantic_sha256": "5" * 64,
        "selected_route_ready": ready,
        "geometry": {
            "selected_T_sh": 120.0,
            "geometry_margin": 5.0 if ready else -5.0,
            "separation_geometry_feasible": ready,
        },
        "truth_boundary": {
            "current_lineage_Xi_data_bound": True,
            "PA10_T_sh_formula_executable": True,
            "selected_T_sh_lower_bound_numerically_instantiated": True,
            "source_T_sh_lower_bound_verified": False,
            "five_moment_repair_applied": False,
        },
    }


def _solve() -> dict:
    incoming = np.array([0.2, -0.1, 0.05, 0.3, -0.04])
    pre = np.array([0.24, -0.08, 0.07, 0.25, -0.03])
    achieved = -pre + np.array([1e-13, -2e-13, 0.0, 1e-13, 0.0])
    target = -pre
    residual = achieved - target
    exit_error = pre + achieved
    return {
        "eta": 0.25,
        "incoming_scaled_discrepancy": incoming.tolist(),
        "pre_repair_scaled_discrepancy": pre.tolist(),
        "repair_achieved_scaled_moments": achieved.tolist(),
        "repair_target_scaled_moments": target.tolist(),
        "repair_residual_scaled_moments": residual.tolist(),
        "exit_scaled_discrepancy": exit_error.tolist(),
        "coefficients": [0.03, -0.02, 0.01, 0.005, -0.004],
        "max_abs_exit_discrepancy": float(np.max(np.abs(exit_error))),
        "scaled_jacobian_condition": 6.0,
    }


class Backend:
    def __init__(self, *, ready: bool = True, solve: dict | None = None):
        self.bridge = _bridge()
        self.report = _tsh_report(ready=ready)
        self.solve = _solve() if solve is None else solve
        self.solve_calls = 0

    def current_xi_bridge_receipt(self):
        return self.bridge

    def upstream_identity(self):
        return REPAIR_UPSTREAM

    def tsh_report(self):
        return self.report

    def pa16_solve_at_eta(self, eta: float):
        self.solve_calls += 1
        return self.solve


def test_current_pa16_repair_receipt_closes_and_stays_scoped():
    backend = Backend()
    receipt = materialize_current_pa16_repair_receipt(backend)
    assert receipt.status == "current_source_moment_repair_verified"
    assert receipt.selected_route_ready is True
    assert backend.solve_calls == 1
    assert receipt.max_incoming_alignment_abs == 0.0
    assert receipt.max_repair_algebra_abs < 5e-12
    assert receipt.any_nontrivial_correction is True
    assert receipt.samples[0].raw_contraction_ratio < 1e-10
    assert receipt.samples[0].transformed_contraction_ratio < 1e-10
    assert receipt.receipt_sha256 and len(receipt.receipt_sha256) == 64
    truth = receipt.to_dict()["truth_boundary"]
    assert truth["current_pa16_repair_receipt_materialized_when_route_ready"] is True
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["discrepancy_from_complete_ns_defect"] is False
    assert truth["authorized_for_gain_gated_ns_stage"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["pde_validated"] is False


def test_current_pa16_repair_fails_on_incoming_discrepancy_laundering():
    solve = _solve()
    solve["incoming_scaled_discrepancy"][0] += 1e-3
    with pytest.raises(CurrentPA16RepairReceiptError, match="not the #949 handoff"):
        materialize_current_pa16_repair_receipt(Backend(solve=solve))


def test_current_pa16_repair_fails_on_algebra_drift():
    solve = _solve()
    solve["exit_scaled_discrepancy"][1] += 1e-4
    solve["max_abs_exit_discrepancy"] = max(abs(v) for v in solve["exit_scaled_discrepancy"])
    with pytest.raises(CurrentPA16RepairReceiptError, match="does not reproduce exit"):
        materialize_current_pa16_repair_receipt(Backend(solve=solve))


def test_current_pa16_repair_fails_on_target_drift():
    solve = _solve()
    solve["repair_target_scaled_moments"][2] += 1e-4
    with pytest.raises(CurrentPA16RepairReceiptError, match="target is not -pre_repair"):
        materialize_current_pa16_repair_receipt(Backend(solve=solve))


def test_current_pa16_repair_fails_on_exact_upstream_identity_drift():
    backend = Backend()
    bad = replace(REPAIR_UPSTREAM, exact_head="0" * 40)
    backend.upstream_identity = lambda: bad
    with pytest.raises(CurrentPA16RepairReceiptError, match="not pinned Agent-1 #953"):
        materialize_current_pa16_repair_receipt(backend)


def test_current_pa16_repair_fails_on_parent_semantic_drift():
    backend = Backend()
    backend.report["current_parent_semantic_sha256"] = "9" * 64
    with pytest.raises(CurrentPA16RepairReceiptError, match="parent semantic SHA"):
        materialize_current_pa16_repair_receipt(backend)


def test_infeasible_current_tsh_route_is_blocked_without_solver_fallback():
    backend = Backend(ready=False)
    receipt = materialize_current_pa16_repair_receipt(backend)
    assert receipt.status == "blocked_by_current_tsh_geometry"
    assert receipt.samples == ()
    assert backend.solve_calls == 0
    assert receipt.blocked_reason is not None
    assert receipt.to_dict()["truth_boundary"]["source_T_sh_lower_bound_verified"] is False


def test_public_api_exposes_no_scientific_knobs():
    signature = inspect.signature(materialize_current_pa16_repair_receipt)
    assert list(signature.parameters) == ["backend"]
    forbidden = {
        "residual",
        "defect",
        "discrepancy",
        "T_sh",
        "correction",
        "pressure",
        "forcing",
        "gain",
        "damping",
        "viscosity",
        "tolerance",
        "threshold",
        "scientific_threshold",
    }
    assert forbidden.isdisjoint(signature.parameters)
    truth = truth_boundary()
    assert truth["existing_agent1_pa16_solver_reused"] is True
    assert truth["real_candidate_finite_correction_cycle_run"] is False
