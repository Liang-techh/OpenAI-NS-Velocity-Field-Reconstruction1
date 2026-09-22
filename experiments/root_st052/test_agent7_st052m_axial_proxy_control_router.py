from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent7_st052m_axial_proxy_control_router as router


def axial_receipt(*, aspect=1.2, elongated=True):
    radial = 2.0
    axial = aspect * radial
    return {
        "schema": router.A9_SCHEMA,
        "task_id": router.A9_TASK_ID,
        "candidate_identity": {
            "candidate_id": router.EXPECTED_CANDIDATE_ID,
            "candidate_sha256": router.EXPECTED_CANDIDATE_SHA256,
            "velocity_identity_sha256": router.EXPECTED_VELOCITY_IDENTITY_SHA256,
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
        "measurement": {"metrics": {
            "full_axial_rms": axial,
            "full_radial_rms": radial,
            "full_aspect_ratio": aspect,
            "axially_elongated_proxy": elongated,
        }},
        "measurement_sha256": "a" * 64,
        "truth_boundary": {
            "axial_stretching_correspondence_verified": False,
            "visualization_ready": False,
            "pde_validated": False,
        },
    }


def control_receipt(*, derivative=0.4, dominant="axial_turnover", unresponsive=False, geometry=True, stability=True):
    derivs = {channel: 0.01 for channel in router.CHANNELS}
    if dominant is not None:
        derivs[dominant] = derivative
    return {
        "task_id": router.A7_TASK_ID,
        "prereg_issue": router.A7_PREREG_ISSUE,
        "source_parent": {"pr": router.A7_PARENT_PR, "head": router.A7_PARENT_HEAD},
        "frozen_protocol": {
            "channels": list(router.CHANNELS),
            "defect_coordinates": [router.COORDINATE],
            "new_basis_dimension": 0,
            "coefficient_selected": None,
            "public_openai_numeric_target": None,
        },
        "defect_coordinate_controllability": {
            "normalized_column_rank": 5,
            "normalized_condition": 4.0,
            "compressed_routing_geometry_gate_passes": geometry,
            "all_derivative_stability_gates_pass": stability,
            "coordinates": {router.COORDINATE: {
                "fine_raw_derivative_by_channel": derivs,
                "fine_raw_row_l2": 0.0 if unresponsive else math.sqrt(sum(v*v for v in derivs.values())),
                "structurally_unresponsive_at_numeric_floor": unresponsive,
                "dominant_aligned_existing_channel": dominant,
            }},
        },
        "decision": {
            "candidate_mutation_authorized_by_this_audit": False,
            "direct_visualization_fingerprint_improvement": 0.0,
        },
    }


class RouterTests(unittest.TestCase):
    def test_elongated_proxy_does_not_authorize_change(self):
        report = router.route(axial_receipt(aspect=1.2, elongated=True), control_receipt())
        d = report["decision"]
        self.assertEqual(d["action"], "no_axial_proxy_deficit_identified")
        self.assertIsNone(d["recommended_existing_control"])
        self.assertFalse(d["bounded_existing_control_child_test_justified"])
        self.assertFalse(d["targeted_axial_tip_basis_preflight_justified"])
        self.assertFalse(d["candidate_mutation_authorized_by_this_router"])
        self.assertEqual(report["truth"]["direct_visualization_fingerprint_improvement"], 0.0)

    def test_deficit_routes_positive_dominant_existing_control_upward(self):
        report = router.route(
            axial_receipt(aspect=0.8, elongated=False),
            control_receipt(derivative=0.45, dominant="axial_turnover"),
        )
        d = report["decision"]
        self.assertEqual(d["action"], "test_existing_control_first")
        self.assertEqual(d["recommended_existing_control"], "axial_turnover")
        self.assertEqual(d["recommended_control_direction_to_increase_axial_aspect"], "increase")
        self.assertTrue(d["bounded_existing_control_child_test_justified"])
        self.assertFalse(d["targeted_axial_tip_basis_preflight_justified"])

    def test_deficit_routes_negative_dominant_existing_control_downward(self):
        report = router.route(
            axial_receipt(aspect=0.8, elongated=False),
            control_receipt(derivative=-0.2, dominant="radial_shape"),
        )
        d = report["decision"]
        self.assertEqual(d["recommended_existing_control"], "radial_shape")
        self.assertEqual(d["recommended_control_direction_to_increase_axial_aspect"], "decrease")

    def test_unresponsive_row_justifies_targeted_basis_preflight_only(self):
        report = router.route(
            axial_receipt(aspect=0.8, elongated=False),
            control_receipt(derivative=0.0, dominant=None, unresponsive=True),
        )
        d = report["decision"]
        self.assertEqual(d["action"], "targeted_axial_tip_basis_preflight_only")
        self.assertTrue(d["targeted_axial_tip_basis_preflight_justified"])
        self.assertFalse(d["basis_growth_authorized_for_candidate_promotion"])
        self.assertFalse(d["candidate_mutation_authorized_by_this_router"])

    def test_unstable_geometry_routes_reparameterization_not_basis(self):
        report = router.route(
            axial_receipt(aspect=0.8, elongated=False),
            control_receipt(geometry=False, stability=True),
        )
        d = report["decision"]
        self.assertEqual(d["action"], "reparameterize_existing_controls_before_basis_growth")
        self.assertFalse(d["bounded_existing_control_child_test_justified"])
        self.assertFalse(d["targeted_axial_tip_basis_preflight_justified"])

    def test_identity_and_public_target_drift_fail_closed(self):
        a = axial_receipt()
        a["candidate_identity"]["candidate_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "identity drift"):
            router.route(a, control_receipt())
        b = axial_receipt()
        b["public_observable"]["numerical_target"] = 1.5
        with self.assertRaisesRegex(ValueError, "numerical target"):
            router.route(b, control_receipt())

    def test_blocked_report_never_routes(self):
        report = router.blocked_report({"a9": "queued", "a7": "queued"})
        self.assertEqual(report["decision"]["status"], "blocked_upstream")
        self.assertIsNone(report["decision"]["recommended_existing_control"])
        self.assertFalse(report["decision"]["targeted_axial_tip_basis_preflight_justified"])
        self.assertEqual(report["truth"]["direct_visualization_fingerprint_improvement"], 0.0)


if __name__ == "__main__":
    unittest.main()
