from __future__ import annotations

import unittest

import numpy as np

import agent7_st052m_defect_space_projection as mod


def audit_from_matrix(matrix: np.ndarray) -> dict:
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (8, 5):
        raise ValueError("expected 8x5 matrix")
    return {
        "coordinates": {
            coordinate: {
                "fine_raw_derivative_by_channel": {
                    channel: float(matrix[i, j])
                    for j, channel in enumerate(mod.CHANNELS)
                }
            }
            for i, coordinate in enumerate(mod.DEFECT_COORDINATES)
        }
    }


def full_rank_fixture() -> np.ndarray:
    return np.asarray(
        [
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 1.0],
        ],
        dtype=float,
    )


class ProjectionAuditTests(unittest.TestCase):
    def test_full_rank_has_three_dimensional_left_null_space(self) -> None:
        report = mod.projection_audit_from_parent(audit_from_matrix(full_rank_fixture()))
        self.assertEqual(report["rank"], 5)
        self.assertEqual(report["left_nullity"], 3)
        self.assertAlmostEqual(report["projector_trace"], 5.0, places=10)
        self.assertAlmostEqual(report["left_null_projector_trace"], 3.0, places=10)
        self.assertLess(report["projector_idempotence_error_l2"], 1e-10)
        self.assertLess(report["left_null_projector_idempotence_error_l2"], 1e-10)
        self.assertLess(report["projector_null_orthogonality_error_l2"], 1e-10)
        self.assertFalse(report["rank_below_five_existing_control_degeneracy"])
        self.assertFalse(
            report["independent_eight_coordinate_control_possible_with_five_local_directions"]
        )

    def test_row_balancing_is_invariant_to_positive_per_coordinate_scale(self) -> None:
        matrix = full_rank_fixture()
        scales = np.asarray([1e-3, 2.0, 17.0, 0.2, 9.0, 3.0, 11.0, 0.7])
        a = mod.projection_audit_from_parent(audit_from_matrix(matrix))
        b = mod.projection_audit_from_parent(audit_from_matrix(matrix * scales[:, None]))
        self.assertTrue(
            np.allclose(
                np.asarray(a["row_balanced_jacobian"]),
                np.asarray(b["row_balanced_jacobian"]),
                rtol=1e-12,
                atol=1e-12,
            )
        )
        self.assertTrue(
            np.allclose(
                np.asarray(a["existing_control_projector"]),
                np.asarray(b["existing_control_projector"]),
                rtol=1e-11,
                atol=1e-11,
            )
        )

    def test_coordinate_projection_norms_are_consistent_with_leverage(self) -> None:
        report = mod.projection_audit_from_parent(audit_from_matrix(full_rank_fixture()))
        for row in report["coordinate_capacity"].values():
            leverage = row["leverage"]
            obstruction = row["obstruction_exposure"]
            self.assertGreaterEqual(leverage, -1e-12)
            self.assertLessEqual(leverage, 1.0 + 1e-12)
            self.assertAlmostEqual(leverage + obstruction, 1.0, places=10)
            self.assertAlmostEqual(
                row["best_attainable_projection_l2"] ** 2,
                leverage,
                places=10,
            )
            self.assertAlmostEqual(
                row["unavoidable_orthogonal_residual_l2"] ** 2,
                obstruction,
                places=10,
            )

    def test_rank_deficiency_is_reported_without_inventing_new_basis(self) -> None:
        matrix = full_rank_fixture()
        matrix[:, 4] = matrix[:, 3]
        report = mod.projection_audit_from_parent(audit_from_matrix(matrix))
        self.assertEqual(report["rank"], 4)
        self.assertEqual(report["left_nullity"], 4)
        self.assertTrue(report["rank_below_five_existing_control_degeneracy"])

    def test_unresponsive_coordinate_fails_closed(self) -> None:
        matrix = full_rank_fixture()
        matrix[6, :] = 0.0
        with self.assertRaisesRegex(RuntimeError, "structurally unresponsive"):
            mod.projection_audit_from_parent(audit_from_matrix(matrix))

    def test_nonfinite_jacobian_fails_closed(self) -> None:
        matrix = full_rank_fixture()
        matrix[2, 1] = np.nan
        with self.assertRaisesRegex(RuntimeError, "nonfinite"):
            mod.projection_audit_from_parent(audit_from_matrix(matrix))

    def test_null_basis_sign_canonicalization(self) -> None:
        basis = np.asarray(
            [
                [-0.9, 0.1],
                [0.1, -0.8],
                [0.2, 0.3],
            ],
            dtype=float,
        )
        out = mod._canonicalize_null_basis(basis)
        for j in range(out.shape[1]):
            pivot = int(np.argmax(np.abs(out[:, j])))
            self.assertGreaterEqual(out[pivot, j], 0.0)

    def test_truth_boundary_remains_fail_closed(self) -> None:
        self.assertFalse(mod.TRUTH["candidate_velocity_changed"])
        self.assertFalse(mod.TRUTH["basis_dimension_changed"])
        self.assertFalse(mod.TRUTH["coefficient_selected"])
        self.assertFalse(mod.TRUTH["held_out_pde_residual_evaluated"])
        self.assertFalse(mod.TRUTH["public_image_numeric_target_used"])
        self.assertEqual(mod.TRUTH["direct_visualization_fingerprint_improvement"], 0.0)
        self.assertFalse(mod.TRUTH["pde_validated"])
        self.assertFalse(mod.TRUTH["paper_exact"])
        self.assertFalse(mod.TRUTH["openai_field_identified"])


if __name__ == "__main__":
    unittest.main()
