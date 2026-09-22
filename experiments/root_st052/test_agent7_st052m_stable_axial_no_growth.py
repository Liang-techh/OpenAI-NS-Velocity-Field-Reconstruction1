from __future__ import annotations

import copy
import unittest

from experiments.root_st052 import agent7_st052m_stable_axial_no_growth as audit


class StableAxialNoGrowthTests(unittest.TestCase):
    def make_receipt(self):
        measurement = {
            "box": [-2.0, 2.0],
            "grid_resolution": 25,
            "metrics": {
                "axially_elongated_proxy": True,
                "full_aspect_ratio": 1.0813103938217585,
                "full_axial_rms": 1.0237804347017856,
                "full_enstrophy_trapezoid": 22.16496417748613,
                "full_radial_rms": 0.9467960731269397,
            },
            "proxy_interpretation": (
                "aspect>1 means only that this candidate's enstrophy-weighted axial extent "
                "exceeds its radial extent under this autonomous cube diagnostic; it is not "
                "an OpenAI numerical target or a correspondence verdict"
            ),
            "spacing": 0.16666666666666674,
            "time": 0.5,
        }
        self.assertEqual(audit._canonical_sha256(measurement), audit.EXPECTED_MEASUREMENT_SHA256)
        return {
            "task_id": audit.A9_TASK_ID,
            "schema": audit.A9_SCHEMA,
            "stable_candidate_identity": {
                "candidate_id": audit.EXPECTED_CANDIDATE_ID,
                "candidate_semantic_identity_sha256": audit.EXPECTED_CANDIDATE_SEMANTIC_ID,
                "velocity_semantic_identity_sha256": audit.EXPECTED_VELOCITY_SEMANTIC_ID,
            },
            "measurement": measurement,
            "measurement_sha256": audit.EXPECTED_MEASUREMENT_SHA256,
            "protocol": {
                "box": [-2.0, 2.0],
                "cartesian_curl_order": 2,
                "grid_resolution": 25,
                "integration_rule": "tensor_product_trapezoid",
                "pixel_loss_used": False,
                "renderer_or_camera_used": False,
                "source_numeric_targets_used": False,
                "time": 0.5,
            },
            "public_observable": {
                "id": "axial_stretching_trajectories",
                "numerical_target": None,
                "source_contract_admitted_by_this_increment": False,
            },
            "truth_boundary": {
                "stable_semantic_identity_bound": True,
                "axial_vorticity_aspect_measured": True,
                "historical_cr_a9_104_receipt_reused": False,
                "candidate_changed": False,
                "velocity_coefficients_changed": False,
                "pressure_changed": False,
                "forcing_changed": False,
                "scientific_threshold_changed": False,
                "public_source_numeric_targets_used": False,
                "renderer_or_camera_used": False,
                "pixel_loss_used": False,
                "axial_stretching_correspondence_verified": False,
                "visual_correspondence_verified": False,
                "source_correspondence_verified": False,
                "visualization_ready": False,
                "velocity_export_ready": False,
                "pde_validated": False,
                "paper_exact": False,
                "openai_field_identified": False,
                "blowup_proved": False,
            },
        }

    def test_exact_stable_receipt_closes_axial_only_growth_trigger(self):
        report = audit.build_report(self.make_receipt())
        self.assertEqual(report["task_id"], audit.TASK_ID)
        self.assertEqual(
            report["decision"]["status"],
            "no_candidate_side_axial_extent_deficit_established",
        )
        self.assertFalse(report["decision"]["axial_specific_basis_growth_authorized"])
        self.assertFalse(report["decision"]["axial_turnover_adjustment_authorized"])
        self.assertFalse(report["decision"]["candidate_mutation_authorized"])
        self.assertIsNone(report["decision"]["coefficient_selected"])
        self.assertTrue(
            report["decision"]["distinct_stable_identity_bound_defect_required_for_future_axial_mutation"]
        )
        self.assertEqual(report["truth"]["direct_visualization_fingerprint_improvement"], 0.0)
        self.assertFalse(report["truth"]["pde_validated"])
        self.assertFalse(report["truth"]["visual_correspondence_verified"])

    def test_rejects_candidate_semantic_identity_drift(self):
        receipt = self.make_receipt()
        receipt["stable_candidate_identity"]["candidate_semantic_identity_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "candidate semantic identity drift"):
            audit.build_report(receipt)

    def test_rejects_velocity_semantic_identity_drift(self):
        receipt = self.make_receipt()
        receipt["stable_candidate_identity"]["velocity_semantic_identity_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "velocity semantic identity drift"):
            audit.build_report(receipt)

    def test_rejects_measurement_tamper_even_with_frozen_sha_field(self):
        receipt = self.make_receipt()
        receipt["measurement"]["metrics"]["full_axial_rms"] += 1.0e-6
        with self.assertRaisesRegex(ValueError, "measurement content/hash mismatch"):
            audit.build_report(receipt)

    def test_rejects_public_numeric_axial_target(self):
        receipt = self.make_receipt()
        receipt["public_observable"]["numerical_target"] = 1.2
        with self.assertRaisesRegex(ValueError, "public axial numerical target"):
            audit.build_report(receipt)

    def test_rejects_visual_truth_promotion(self):
        receipt = self.make_receipt()
        receipt["truth_boundary"]["visual_correspondence_verified"] = True
        with self.assertRaisesRegex(ValueError, "forbidden truth promotion"):
            audit.build_report(receipt)

    def test_rejects_stale_historical_receipt_reuse(self):
        receipt = self.make_receipt()
        receipt["truth_boundary"]["historical_cr_a9_104_receipt_reused"] = True
        with self.assertRaisesRegex(ValueError, "historical stale receipt reused"):
            audit.build_report(receipt)

    def test_rejects_protocol_drift(self):
        receipt = self.make_receipt()
        receipt["protocol"]["grid_resolution"] = 33
        with self.assertRaisesRegex(ValueError, "grid resolution drift"):
            audit.build_report(receipt)


if __name__ == "__main__":
    unittest.main()
