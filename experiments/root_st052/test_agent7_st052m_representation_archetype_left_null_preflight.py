from __future__ import annotations

import unittest

import numpy as np

import agent7_st052m_representation_archetype_left_null_preflight as target


class RepresentationArchetypeLeftNullPreflightTests(unittest.TestCase):
    def test_time_activations_are_frozen_and_endpoint_safe(self) -> None:
        self.assertEqual(target._activation_g1(0.25), 0.0)
        self.assertEqual(target._activation_g1(0.75), 1.0)
        self.assertEqual(target._activation_skew(0.25), 0.0)
        self.assertEqual(target._activation_skew(0.50), 0.0)
        self.assertEqual(target._activation_skew(0.75), 0.0)
        self.assertLess(target._activation_skew(0.375), 0.0)
        self.assertGreater(target._activation_skew(0.625), 0.0)

    def test_analysis_only_probe_geometry_passes_structural_guards(self) -> None:
        for name, tangent in target._probe_tangents().items():
            with self.subTest(name=name):
                report = target._structural_guard(tangent)
                self.assertTrue(report["passes"])
                self.assertLessEqual(
                    report["centered_divergence_max_abs"],
                    target.DIVERGENCE_FD_MAX,
                )
                self.assertLessEqual(
                    report["outside_support_max_abs_velocity"],
                    target.OUTSIDE_SUPPORT_MAX,
                )
                self.assertTrue(report["axis_values_finite"])
                self.assertLessEqual(
                    report["axis_transverse_max_abs_velocity"],
                    target.AXIS_TRANSVERSE_MAX,
                )

    def test_compact_toroidal_probe_is_purely_azimuthal(self) -> None:
        points = np.asarray(
            [[0.4, 0.2, 0.3], [-0.5, 0.3, -0.8], [0.7, -0.4, 1.0]],
            dtype=float,
        )
        velocity = target._toroidal_unit(points)
        radius = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)
        radial = (
            points[:, 0] * velocity[:, 0] + points[:, 1] * velocity[:, 1]
        ) / radius
        self.assertLessEqual(float(np.max(np.abs(radial))), 1.0e-14)
        self.assertLessEqual(float(np.max(np.abs(velocity[:, 2]))), 1.0e-14)
        self.assertGreater(float(np.linalg.norm(velocity)), 0.0)

    def test_poloidal_probe_is_regular_on_axis_and_nontrivial_off_axis(self) -> None:
        axis = np.asarray([[0.0, 0.0, -0.8], [0.0, 0.0, 0.8]], dtype=float)
        on_axis = target._poloidal_unit(axis)
        self.assertTrue(np.all(np.isfinite(on_axis)))
        self.assertLessEqual(float(np.max(np.abs(on_axis[:, :2]))), 1.0e-14)
        off_axis = target._poloidal_unit(
            np.asarray([[0.4, 0.0, 0.8], [0.4, 0.0, -0.8]], dtype=float)
        )
        self.assertGreater(float(np.linalg.norm(off_axis)), 0.0)

    def test_synthetic_left_null_probe_is_independent_and_eligible(self) -> None:
        existing = np.zeros((8, 5), dtype=float)
        existing[:5, :] = np.eye(5)
        projector = np.diag([1.0] * 5 + [0.0] * 3)
        null_projector = np.eye(8) - projector
        row_norms = np.ones(8, dtype=float)
        derivative = np.zeros(8, dtype=float)
        derivative[5] = 1.0
        structural = {"passes": True}

        report = target._score_probe(
            existing,
            projector,
            null_projector,
            row_norms,
            derivative,
            derivative,
            structural,
        )
        self.assertEqual(report["existing_rank"], 5)
        self.assertEqual(report["augmented_rank"], 6)
        self.assertEqual(report["rank_gain"], 1)
        self.assertAlmostEqual(report["left_null_novelty_fraction_l2"], 1.0)
        self.assertAlmostEqual(report["existing_span_fraction_l2"], 0.0)
        self.assertAlmostEqual(report["augmented_normalized_column_condition"], 1.0)
        self.assertAlmostEqual(report["max_abs_probe_existing_cosine"], 0.0)
        self.assertTrue(report["independent_preflight_eligible"])

    def test_synthetic_redundant_probe_is_not_basis_growth_evidence(self) -> None:
        existing = np.zeros((8, 5), dtype=float)
        existing[:5, :] = np.eye(5)
        projector = np.diag([1.0] * 5 + [0.0] * 3)
        null_projector = np.eye(8) - projector
        row_norms = np.ones(8, dtype=float)
        derivative = np.zeros(8, dtype=float)
        derivative[0] = 1.0

        report = target._score_probe(
            existing,
            projector,
            null_projector,
            row_norms,
            derivative,
            derivative,
            {"passes": True},
        )
        self.assertEqual(report["rank_gain"], 0)
        self.assertAlmostEqual(report["left_null_novelty_fraction_l2"], 0.0)
        self.assertFalse(report["independent_preflight_eligible"])

    def test_derivative_instability_fails_closed(self) -> None:
        existing = np.zeros((8, 5), dtype=float)
        existing[:5, :] = np.eye(5)
        projector = np.diag([1.0] * 5 + [0.0] * 3)
        null_projector = np.eye(8) - projector
        row_norms = np.ones(8, dtype=float)
        fine = np.zeros(8, dtype=float)
        fine[5] = 1.0
        coarse = -fine

        report = target._score_probe(
            existing,
            projector,
            null_projector,
            row_norms,
            coarse,
            fine,
            {"passes": True},
        )
        self.assertFalse(report["derivative_stable"])
        self.assertFalse(report["independent_preflight_eligible"])

    def test_truth_boundary_keeps_actual_velocity_unchanged(self) -> None:
        self.assertFalse(target.TRUTH["candidate_velocity_changed"])
        self.assertFalse(target.TRUTH["basis_dimension_changed"])
        self.assertFalse(target.TRUTH["coefficient_selected"])
        self.assertFalse(target.TRUTH["held_out_pde_residual_evaluated"])
        self.assertEqual(
            target.TRUTH["direct_visualization_fingerprint_improvement"], 0.0
        )
        self.assertFalse(target.TRUTH["pde_validated"])
        self.assertFalse(target.TRUTH["openai_field_identified"])


if __name__ == "__main__":
    unittest.main()
