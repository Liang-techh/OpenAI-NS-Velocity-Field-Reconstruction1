from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent7_st052m_axial_evidence_admissibility as audit


def good_run():
    return {
        "id": audit.A9_RUN,
        "head_sha": audit.A9_HEAD,
        "status": "completed",
        "conclusion": "failure",
    }


def good_receipt():
    return {
        "schema": audit.A9_SCHEMA,
        "task_id": audit.A9_TASK_ID,
        "candidate_identity": {
            "candidate_id": "ST052-M-linear-temporal-child-v1",
            "candidate_sha256": audit.REPLAY_CANDIDATE_SHA256,
            "velocity_identity_sha256": audit.REPLAY_VELOCITY_SHA256,
        },
        "protocol": {
            "time": 0.5,
            "grid_resolution": 25,
            "box": [-2.0, 2.0],
            "cartesian_curl_order": 2,
            "integration_rule": "tensor_product_trapezoid",
            "source_numeric_targets_used": False,
            "renderer_or_camera_used": False,
            "pixel_loss_used": False,
        },
        "public_observable": {"numerical_target": None},
        "measurement": {
            "metrics": {
                "full_axial_rms": audit.EXPECTED_AXIAL_RMS,
                "full_radial_rms": audit.EXPECTED_RADIAL_RMS,
                "full_aspect_ratio": audit.EXPECTED_ASPECT,
                "axially_elongated_proxy": True,
            }
        },
        "measurement_sha256": audit.MEASUREMENT_SHA256,
        "truth_boundary": {
            "axial_stretching_correspondence_verified": False,
            "visualization_ready": False,
            "pde_validated": False,
        },
    }


class AdmissibilityTests(unittest.TestCase):
    def test_failed_but_measured_receipt_is_quarantined(self):
        report = audit.classify(good_run(), good_receipt())
        decision = report["decision"]
        self.assertEqual(decision["status"], "blocked_identity_admissibility")
        self.assertEqual(decision["action"], "repair_semantic_identity_before_basis_decision")
        self.assertFalse(decision["measurement_admitted_for_basis_routing"])
        self.assertFalse(decision["axial_proxy_no_deficit_conclusion_admitted"])
        self.assertFalse(decision["existing_control_selection_authorized"])
        self.assertFalse(decision["targeted_basis_preflight_authorized"])
        self.assertIsNone(decision["recommended_existing_control"])
        self.assertEqual(report["descriptive_unadmitted_measurement"]["full_aspect_ratio"], audit.EXPECTED_ASPECT)
        self.assertTrue(report["descriptive_unadmitted_measurement"]["axially_elongated_proxy"])
        self.assertEqual(report["truth"]["direct_visualization_fingerprint_improvement"], 0.0)

    def test_successful_run_is_not_reclassified_by_failure_quarantine(self):
        run = good_run()
        run["conclusion"] = "success"
        with self.assertRaisesRegex(ValueError, "completed/failed"):
            audit.classify(run, good_receipt())

    def test_old_admitted_identity_cannot_be_substituted_for_failed_replay(self):
        receipt = good_receipt()
        receipt["candidate_identity"]["candidate_sha256"] = audit.ADMITTED_CANDIDATE_SHA256
        with self.assertRaisesRegex(ValueError, "replay candidate"):
            audit.classify(good_run(), receipt)

    def test_velocity_identity_drift_is_fail_closed(self):
        receipt = good_receipt()
        receipt["candidate_identity"]["velocity_identity_sha256"] = audit.ADMITTED_VELOCITY_SHA256
        with self.assertRaisesRegex(ValueError, "replay velocity"):
            audit.classify(good_run(), receipt)

    def test_public_numeric_target_laundering_is_rejected(self):
        receipt = good_receipt()
        receipt["public_observable"]["numerical_target"] = 1.4
        with self.assertRaisesRegex(ValueError, "numerical target"):
            audit.classify(good_run(), receipt)

    def test_measurement_digest_drift_is_rejected(self):
        receipt = good_receipt()
        receipt["measurement_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "measurement digest"):
            audit.classify(good_run(), receipt)

    def test_metric_mutation_is_rejected(self):
        receipt = good_receipt()
        receipt["measurement"]["metrics"]["full_aspect_ratio"] += 1e-3
        with self.assertRaisesRegex(ValueError, "aspect drift"):
            audit.classify(good_run(), receipt)

    def test_truth_promotion_is_rejected(self):
        receipt = good_receipt()
        receipt["truth_boundary"]["visualization_ready"] = True
        with self.assertRaisesRegex(ValueError, "visual/PDE"):
            audit.classify(good_run(), receipt)

    def test_report_keeps_identity_mismatch_explicit(self):
        report = audit.classify(good_run(), good_receipt())
        identity = report["identity_admissibility"]
        self.assertFalse(identity["candidate_identity_matches"])
        self.assertFalse(identity["velocity_identity_matches"])
        self.assertFalse(identity["same_identity_gate_passed"])
        self.assertEqual(identity["classification"], "identity_drifted_measurement_observed_not_admitted")
        self.assertFalse(report["decision"]["hash_replacement_alone_is_sufficient"])


if __name__ == "__main__":
    unittest.main()
