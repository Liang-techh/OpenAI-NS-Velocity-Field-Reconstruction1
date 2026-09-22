from __future__ import annotations

import unittest

import numpy as np

import agent7_st052m_one_direction_obstruction_reduction as audit


class OneDirectionObstructionReductionTests(unittest.TestCase):
    @staticmethod
    def _five_coordinate_identity() -> np.ndarray:
        j = np.zeros((len(audit.DEFECT_COORDINATES), len(audit.CHANNELS)), dtype=float)
        for k in range(len(audit.CHANNELS)):
            j[k, k] = 1.0
        return j

    def test_independent_sixth_direction_removes_one_null_coordinate(self) -> None:
        j5 = self._five_coordinate_identity()
        d = np.zeros(len(audit.DEFECT_COORDINATES), dtype=float)
        d[5] = 1.0
        r = audit.obstruction_reduction_for_probe(
            j5, d, {"independent_preflight_eligible": True}
        )
        coordinate = audit.DEFECT_COORDINATES[5]
        row = r["coordinate_obstruction_reduction"][coordinate]
        self.assertEqual(r["existing_rank"], 5)
        self.assertEqual(r["augmented_rank"], 6)
        self.assertEqual(r["rank_gain"], 1)
        self.assertEqual(r["left_nullity_before"], 3)
        self.assertEqual(r["left_nullity_after"], 2)
        self.assertAlmostEqual(r["left_null_trace_reduction"], 1.0, places=12)
        self.assertAlmostEqual(row["unavoidable_residual_l2_before"], 1.0, places=12)
        self.assertAlmostEqual(row["unavoidable_residual_l2_after"], 0.0, places=12)
        self.assertAlmostEqual(row["absolute_obstruction_reduction_l2"], 1.0, places=12)
        self.assertAlmostEqual(row["fractional_obstruction_reduction"], 1.0, places=12)
        self.assertAlmostEqual(row["appended_probe_coefficient"], 1.0, places=12)
        self.assertTrue(r["at_least_two_uncontrolled_combinations_remain_if_rank5_to6"])
        self.assertTrue(r["eligible_for_descriptive_obstruction_ranking"])

    def test_coupled_sixth_direction_reports_unavoidable_spillover(self) -> None:
        j5 = self._five_coordinate_identity()
        d = np.zeros(len(audit.DEFECT_COORDINATES), dtype=float)
        d[5] = 1.0
        d[6] = 1.0
        r = audit.obstruction_reduction_for_probe(
            j5, d, {"independent_preflight_eligible": True}
        )
        target = audit.DEFECT_COORDINATES[5]
        coupled = audit.DEFECT_COORDINATES[6]
        row = r["coordinate_obstruction_reduction"][target]
        self.assertEqual(row["largest_off_target_projected_coordinate_after"], coupled)
        self.assertAlmostEqual(row["largest_off_target_projected_abs_after"], 0.5, places=12)
        self.assertAlmostEqual(row["off_target_projected_l2_after"], 0.5, places=12)
        self.assertGreater(row["absolute_obstruction_reduction_l2"], 0.0)
        self.assertLess(row["absolute_obstruction_reduction_l2"], 1.0)
        self.assertAlmostEqual(r["left_null_trace_reduction"], 1.0, places=12)

    def test_redundant_probe_does_not_reduce_obstruction(self) -> None:
        j5 = self._five_coordinate_identity()
        d = j5[:, 0].copy()
        r = audit.obstruction_reduction_for_probe(
            j5, d, {"independent_preflight_eligible": False}
        )
        self.assertEqual(r["rank_gain"], 0)
        self.assertAlmostEqual(r["left_null_trace_reduction"], 0.0, places=12)
        self.assertAlmostEqual(r["maximum_absolute_obstruction_reduction_l2"], 0.0, places=12)
        self.assertFalse(r["eligible_for_descriptive_obstruction_ranking"])
        self.assertFalse(r["candidate_promotion_authorized"])
        self.assertFalse(r["coefficient_selection_authorized"])

    def test_projector_reduction_is_coordinatewise_nonnegative_up_to_roundoff(self) -> None:
        rng = np.random.default_rng(9174137)
        j5 = rng.normal(size=(len(audit.DEFECT_COORDINATES), len(audit.CHANNELS)))
        d = rng.normal(size=len(audit.DEFECT_COORDINATES))
        r = audit.obstruction_reduction_for_probe(
            j5, d, {"independent_preflight_eligible": True}
        )
        for row in r["coordinate_obstruction_reduction"].values():
            self.assertGreaterEqual(row["absolute_obstruction_reduction_l2"], -1.0e-10)
        self.assertLessEqual(r["left_nullity_after"], r["left_nullity_before"])

    def test_zero_probe_fails_closed(self) -> None:
        j5 = self._five_coordinate_identity()
        with self.assertRaises(ValueError):
            audit.obstruction_reduction_for_probe(
                j5,
                np.zeros(len(audit.DEFECT_COORDINATES)),
                {"independent_preflight_eligible": True},
            )

    def test_nonfinite_input_fails_closed(self) -> None:
        j5 = self._five_coordinate_identity()
        d = np.zeros(len(audit.DEFECT_COORDINATES), dtype=float)
        d[5] = np.nan
        with self.assertRaises(ValueError):
            audit.obstruction_reduction_for_probe(
                j5, d, {"independent_preflight_eligible": True}
            )

    def test_truth_boundary_is_frozen(self) -> None:
        self.assertEqual(audit.TASK_ID, "CR003-ST052M-ONE-DIRECTION-OBSTRUCTION-REDUCTION-137")
        self.assertEqual(audit.PREREG_ISSUE, 1155)
        self.assertEqual(audit.STACK_BASE_PR, 1147)
        self.assertEqual(audit.STACK_BASE_HEAD, "f20a883b0b2721d940e9012651b737d52583349b")
        self.assertFalse(audit.TRUTH["candidate_velocity_changed"])
        self.assertFalse(audit.TRUTH["basis_dimension_changed"])
        self.assertFalse(audit.TRUTH["coefficient_selected"])
        self.assertFalse(audit.TRUTH["pressure_or_force_changed"])
        self.assertFalse(audit.TRUTH["held_out_pde_residual_evaluated"])
        self.assertFalse(audit.TRUTH["public_image_numeric_target_used"])
        self.assertEqual(audit.TRUTH["direct_visualization_fingerprint_improvement"], 0.0)
        self.assertFalse(audit.TRUTH["pde_validated"])
        self.assertFalse(audit.TRUTH["paper_exact"])
        self.assertFalse(audit.TRUTH["openai_field_identified"])


if __name__ == "__main__":
    unittest.main()
