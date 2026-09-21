from __future__ import annotations

from dataclasses import replace
import hashlib
import json

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_current_pa16_bridge_integrity as guard
from openai_ns_reconstruction.kokuno_current_pa16_repair_receipt import (
    CurrentPA16RepairReceiptError,
)
from openai_ns_reconstruction.kokuno_current_xi_moment_discrepancy_bridge import (
    CurrentXiMomentBridgeReceipt,
    CurrentXiMomentSample,
    UpstreamA1Identity,
    SCHEMA as XI_BRIDGE_SCHEMA,
)


XI_UPSTREAM = UpstreamA1Identity(
    pr_number=947,
    exact_head="2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b",
    source_path="src/openai_ns_reconstruction/kokuno_pa10_actual_xi_prefix_moments.py",
    source_blob="92dce65c9a8793e06497a7e7347c832861032238",
    report_schema="kokuno-pa10-actual-xi-prefix-moments-v1",
)


def _sha256_json(value) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _bridge() -> CurrentXiMomentBridgeReceipt:
    row = (0.2, -0.1, 0.05, 0.3, -0.04)
    transformed = (0.2, 0.15, -0.1, -0.1, -0.04)
    sample = CurrentXiMomentSample(
        eta=0.25,
        ell_i=0.02,
        G_i=1.01,
        incoming_scaled_discrepancy=row,
        pa16_row_transformed_discrepancy=transformed,
        raw_l2_norm=float(np.linalg.norm(row)),
        pa16_transformed_l2_norm=float(np.linalg.norm(transformed)),
        correction_needed=True,
    )
    unsigned = {
        "schema": XI_BRIDGE_SCHEMA,
        "upstream": {
            "pr_number": XI_UPSTREAM.pr_number,
            "exact_head": XI_UPSTREAM.exact_head,
            "source_path": XI_UPSTREAM.source_path,
            "source_blob": XI_UPSTREAM.source_blob,
            "report_schema": XI_UPSTREAM.report_schema,
        },
        "upstream_semantic_sha256": "1" * 64,
        "upstream_report_sha256": "2" * 64,
        "outer_schedule_sha256": "3" * 64,
        "samples": [
            {
                "eta": sample.eta,
                "ell_i": sample.ell_i,
                "G_i": sample.G_i,
                "incoming_scaled_discrepancy": sample.incoming_scaled_discrepancy,
                "pa16_row_transformed_discrepancy": sample.pa16_row_transformed_discrepancy,
                "raw_l2_norm": sample.raw_l2_norm,
                "pa16_transformed_l2_norm": sample.pa16_transformed_l2_norm,
                "correction_needed": sample.correction_needed,
            }
        ],
    }
    return CurrentXiMomentBridgeReceipt(
        upstream=XI_UPSTREAM,
        upstream_semantic_sha256="1" * 64,
        upstream_report_sha256="2" * 64,
        outer_schedule_sha256="3" * 64,
        samples=(sample,),
        max_raw_l2_norm=sample.raw_l2_norm,
        max_pa16_transformed_l2_norm=sample.pa16_transformed_l2_norm,
        any_correction_needed=True,
        handoff_sha256=_sha256_json(unsigned),
    )


def test_valid_source_bridge_handoff_checksum_replays():
    receipt = _bridge()
    assert guard.replay_current_xi_handoff_sha256(receipt) == receipt.handoff_sha256
    assert guard.validate_current_xi_bridge_handoff(receipt) is receipt


def test_stale_source_bridge_checksum_fails_closed_after_payload_mutation():
    receipt = _bridge()
    sample = receipt.samples[0]
    forged_sample = replace(
        sample,
        incoming_scaled_discrepancy=(
            sample.incoming_scaled_discrepancy[0] + 1.0e-3,
            *sample.incoming_scaled_discrepancy[1:],
        ),
    )
    forged = replace(receipt, samples=(forged_sample,))
    assert forged.handoff_sha256 == receipt.handoff_sha256
    with pytest.raises(CurrentPA16RepairReceiptError, match="handoff SHA does not replay"):
        guard.validate_current_xi_bridge_handoff(forged)


def test_unhashed_summary_drift_also_fails_closed():
    receipt = _bridge()
    forged = replace(receipt, max_raw_l2_norm=receipt.max_raw_l2_norm + 1.0e-3)
    with pytest.raises(CurrentPA16RepairReceiptError, match="max_raw_l2_norm"):
        guard.validate_current_xi_bridge_handoff(forged)


def test_guarded_entrypoint_pins_validated_bridge_before_parent_delegation(monkeypatch):
    first = _bridge()
    stale = replace(
        first,
        samples=(replace(first.samples[0], G_i=first.samples[0].G_i + 0.1),),
    )

    class SwappingBackend:
        def __init__(self):
            self.calls = 0

        def current_xi_bridge_receipt(self):
            self.calls += 1
            return first if self.calls == 1 else stale

        def upstream_identity(self):
            return object()

        def tsh_report(self):
            return {}

        def pa16_solve_at_eta(self, eta: float):
            return {}

    seen = {}

    def fake_parent(backend):
        seen["bridge"] = backend.current_xi_bridge_receipt()
        return "sentinel"

    monkeypatch.setattr(guard, "_materialize_parent", fake_parent)
    backend = SwappingBackend()
    assert guard.materialize_checksum_verified_current_pa16_repair_receipt(backend) == "sentinel"
    assert backend.calls == 1
    assert seen["bridge"] is first


def test_truth_boundary_keeps_scientific_promotions_false():
    truth = guard.truth_boundary()
    assert truth["source_bridge_checksum_payload_replayed"] is True
    assert truth["validated_bridge_pinned_across_parent_delegation"] is True
    assert truth["parent_957_entrypoint_itself_modified"] is False
    assert truth["source_moment_receipt_promoted_to_complete_ns_defect"] is False
    assert truth["authorized_for_gain_gated_ns_stage"] is False
    assert truth["scientific_threshold_modified"] is False
    assert truth["pde_validated"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
