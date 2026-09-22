from __future__ import annotations

import copy
import unittest

from openai_ns_reconstruction.st052_stable_tip_basis_growth_decision import (
    CONTROL_HEAD,
    CONTROL_RUN,
    EXPECTED_CANDIDATE_IDENTITY,
    EXPECTED_VELOCITY_IDENTITY,
    canonical_sha256,
    make_decision,
    validate_source_receipt,
)


def source_receipt(*, tapered: bool = True):
    upper, lower, central = (0.7, 0.8, 1.0) if tapered else (1.1, 0.9, 1.0)
    worst = max(upper, lower)
    measurement = {
        "metrics": {
            "tip_radial_rms_upper": upper,
            "tip_radial_rms_lower": lower,
            "central_radial_rms": central,
            "worst_tip_radial_rms": worst,
            "worst_tip_to_central_ratio": worst / central,
            "upper_lower_tip_relative_mismatch": abs(upper - lower) / worst,
            "both_tips_radially_tapered_proxy": worst < central,
        }
    }
    r = {
        "schema": "st052-stable-identity-bound-tip-taper-proxy/v1",
        "task_id": "CR003-ST052M-STABLE-TIP-TAPER-PROXY-141",
        "prereg_issue": 1194,
        "stable_candidate_identity": {
            "candidate_id": "ST052-M-linear-temporal-child-v1",
            "candidate_semantic_identity_sha256": EXPECTED_CANDIDATE_IDENTITY,
            "velocity_semantic_identity_sha256": EXPECTED_VELOCITY_IDENTITY,
        },
        "public_context": {"public_numerical_target": None},
        "measurement": measurement,
        "measurement_sha256": canonical_sha256(measurement),
        "truth_boundary": {
            "candidate_velocity_changed": False,
            "basis_dimension_changed": False,
            "coefficient_selected": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "visualization_ready": False,
            "source_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    r["receipt_sha256"] = canonical_sha256(r)
    return r


def run(*, status="queued", conclusion=None):
    return {"id": CONTROL_RUN, "head_sha": CONTROL_HEAD, "status": status, "conclusion": conclusion}


class TipDecisionTests(unittest.TestCase):
    def test_tapered_closes_only_coarse_growth_trigger(self):
        d = make_decision(source_receipt(tapered=True), run())
        self.assertEqual(d["routing"]["classification"], "coarse_blunt_tip_basis_trigger_closed")
        self.assertFalse(d["routing"]["new_tip_basis_authorized"])
        self.assertEqual(d["truth_boundary"]["direct_visualization_fingerprint_improvement"], 0.0)

    def test_deficit_blocks_while_existing_control_receipt_unresolved(self):
        d = make_decision(source_receipt(tapered=False), run())
        self.assertEqual(d["routing"]["classification"], "tip_taper_deficit_existing_control_receipt_blocked")
        self.assertFalse(d["routing"]["new_tip_basis_authorized"])

    def test_deficit_routes_to_existing_control_when_exact_run_succeeds(self):
        d = make_decision(source_receipt(tapered=False), run(status="completed", conclusion="success"))
        self.assertEqual(d["routing"]["classification"], "tip_taper_deficit_existing_control_receipt_ready")
        self.assertIn("pr1078", d["routing"]["next_minimal_action"])

    def test_tampered_metric_fails_closed(self):
        r = source_receipt()
        r["measurement"]["metrics"]["central_radial_rms"] = 2.0
        with self.assertRaises(ValueError):
            validate_source_receipt(r)

    def test_semantic_identity_drift_fails_closed(self):
        r = source_receipt()
        r["stable_candidate_identity"]["velocity_semantic_identity_sha256"] = "0" * 64
        r["receipt_sha256"] = canonical_sha256({k: v for k, v in r.items() if k != "receipt_sha256"})
        with self.assertRaises(ValueError):
            validate_source_receipt(r)

    def test_truth_promotion_fails_closed(self):
        r = source_receipt()
        r["truth_boundary"]["pde_validated"] = True
        r["receipt_sha256"] = canonical_sha256({k: v for k, v in r.items() if k != "receipt_sha256"})
        with self.assertRaises(ValueError):
            validate_source_receipt(r)

    def test_control_run_identity_drift_fails_closed(self):
        with self.assertRaises(ValueError):
            make_decision(source_receipt(False), {"id": CONTROL_RUN, "head_sha": "0" * 40, "status": "queued", "conclusion": None})


if __name__ == "__main__":
    unittest.main()
