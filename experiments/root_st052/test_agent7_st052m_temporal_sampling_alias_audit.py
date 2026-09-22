from __future__ import annotations

import unittest

import agent7_st052m_temporal_sampling_alias_audit as audit


class TemporalSamplingAliasAuditTest(unittest.TestCase):
    def test_frozen_time_skew_activation_is_aliased_on_existing_times(self) -> None:
        report = audit.audit_activation(audit.rep._activation_skew)
        self.assertTrue(report["temporal_sampling_alias_detected"])
        self.assertEqual(report["observation_activation"], {
            "0.250": 0.0,
            "0.500": 0.0,
            "0.750": 0.0,
        })
        self.assertAlmostEqual(report["offgrid_witness_activation"]["0.375"], -0.375)
        self.assertAlmostEqual(report["offgrid_witness_activation"]["0.625"], 0.375)

    def test_nonaliased_control_is_not_misclassified(self) -> None:
        report = audit.audit_activation(lambda t: float(t) - 0.25)
        self.assertFalse(report["temporal_sampling_alias_detected"])
        self.assertFalse(report["zero_at_all_observation_times"])

    def test_actual_probe_tangent_is_zero_on_observation_set_and_nonzero_offgrid(self) -> None:
        report = audit.build_report()
        self.assertTrue(report["temporal_sampling_alias_detected"])
        self.assertTrue(report["field_zero_at_observation_times"])
        self.assertTrue(report["field_nonzero_offgrid"])
        self.assertEqual(
            report["routing_consequence"],
            "exclude_time_skewed_probe_from_zero_derivative_capacity_ranking_until_offgrid_temporal_morphology_coordinate_exists",
        )

    def test_truth_boundary_remains_fail_closed(self) -> None:
        report = audit.build_report()
        truth = report["truth"]
        self.assertEqual(truth["direct_visualization_fingerprint_improvement"], 0.0)
        for key in (
            "candidate_velocity_changed",
            "basis_dimension_changed",
            "coefficient_selected",
            "pressure_or_force_changed",
            "held_out_pde_residual_evaluated",
            "visualization_ready",
            "visual_correspondence_verified",
            "pde_validated",
            "paper_exact",
            "openai_field_identified",
            "blowup_proved",
        ):
            self.assertFalse(truth[key], key)


if __name__ == "__main__":
    unittest.main()
