from __future__ import annotations

import importlib.util
import json
import pathlib
import unittest

MODULE_PATH = pathlib.Path(__file__).with_name("agent7_st052m_public_observable_coverage.py")
spec = importlib.util.spec_from_file_location("audit", MODULE_PATH)
audit = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(audit)


def fixture() -> dict:
    return {
        "candidate_identity": {
            "candidate_id": audit.EXPECTED_CANDIDATE_ID,
            "candidate_sha256": audit.EXPECTED_CANDIDATE_SHA256,
            "velocity_identity_sha256": audit.EXPECTED_VELOCITY_IDENTITY_SHA256,
        },
        "morphology_receipt": {
            "candidate_identity_bound": True,
            "cylindrical_morphology_diagnostic_ready": True,
            "measurement_sha256": audit.EXPECTED_MEASUREMENT_SHA256,
            "candidate_identity": {
                "candidate_id": audit.EXPECTED_CANDIDATE_ID,
                "candidate_sha256": audit.EXPECTED_CANDIDATE_SHA256,
                "velocity_identity_sha256": audit.EXPECTED_VELOCITY_IDENTITY_SHA256,
            },
            "measurements": {
                "radial_variation": [
                    {
                        "all_rings_inward": True,
                        "all_rings_swirl_nonzero": True,
                        "angular_rate_span_across_radii": 0.2,
                        "circulation_speed_span_across_radii": 0.1,
                    },
                    {
                        "all_rings_inward": True,
                        "all_rings_swirl_nonzero": True,
                        "angular_rate_span_across_radii": 0.3,
                        "circulation_speed_span_across_radii": 0.15,
                    },
                ],
                "ring_rows": [
                    {"swirl_nonzero_fraction": 1.0, "inward_swirl_fraction": 1.0},
                    {"swirl_nonzero_fraction": 1.0, "inward_swirl_fraction": 1.0},
                ],
                "time_evolution_proxy": {
                    "speed_weighted_radius_decreased": True,
                    "peak_mean_speed_increased": True,
                },
            },
        },
    }


class CoverageAuditTests(unittest.TestCase):
    def test_five_covered_axial_unmeasured_routes_to_existing_coordinate(self) -> None:
        report = audit.audit_receipt(fixture())
        self.assertEqual(report["coverage_summary"]["covered_by_candidate_bound_proxy_count"], 5)
        self.assertEqual(
            report["coverage_summary"]["unresolved_or_unmeasured"],
            ["axial_stretching_elongation"],
        )
        self.assertEqual(report["routing"]["next_existing_agent7_coordinate"], "axial_vorticity_aspect_t050")
        self.assertFalse(report["routing"]["basis_growth_authorized"])
        self.assertEqual(report["truth"]["direct_visualization_fingerprint_improvement"], 0.0)

    def test_identity_drift_fails_closed(self) -> None:
        data = fixture()
        data["candidate_identity"]["candidate_id"] = "wrong"
        with self.assertRaises(RuntimeError):
            audit.audit_receipt(data)

    def test_missing_inward_proxy_is_reported_not_reinterpreted(self) -> None:
        data = fixture()
        data["morphology_receipt"]["measurements"]["radial_variation"][0]["all_rings_inward"] = False
        report = audit.audit_receipt(data)
        self.assertEqual(
            report["coverage"]["inward_spiraling_trajectories"]["status"],
            "not_covered_by_current_receipt",
        )
        self.assertIsNone(report["routing"]["next_existing_agent7_coordinate"])
        self.assertFalse(report["routing"]["candidate_mutation_authorized"])

    def test_zero_variation_is_not_called_public_match(self) -> None:
        data = fixture()
        data["morphology_receipt"]["measurements"]["radial_variation"][1]["angular_rate_span_across_radii"] = 0.0
        report = audit.audit_receipt(data)
        self.assertEqual(
            report["coverage"]["spatial_variation_in_angular_rotation"]["status"],
            "not_covered_by_current_receipt",
        )
        self.assertFalse(report["coverage_summary"]["visual_correspondence_established"])

    def test_explicit_axial_aspect_is_measured_but_not_promoted_to_match(self) -> None:
        data = fixture()
        data["morphology_receipt"]["measurements"]["vorticity_shape"] = {"full_aspect_ratio": 1.4}
        report = audit.audit_receipt(data)
        self.assertEqual(
            report["coverage"]["axial_stretching_elongation"]["status"],
            "measured_but_no_public_magnitude_target",
        )
        self.assertFalse(report["coverage_summary"]["visual_correspondence_established"])
        self.assertFalse(report["routing"]["basis_growth_authorized"])

    def test_receipt_json_roundtrip(self) -> None:
        report = audit.audit_receipt(json.loads(json.dumps(fixture())))
        self.assertEqual(tuple(report["coverage"]), audit.PUBLIC_OBSERVABLES)
        self.assertFalse(any(report["truth"][key] for key in (
            "candidate_velocity_changed",
            "basis_dimension_changed",
            "coefficient_selected",
            "public_image_numeric_target_used",
            "pde_validated",
            "openai_field_identified",
        )))


if __name__ == "__main__":
    unittest.main()
