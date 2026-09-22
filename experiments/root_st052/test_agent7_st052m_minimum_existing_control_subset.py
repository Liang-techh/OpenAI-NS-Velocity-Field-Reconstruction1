from __future__ import annotations

import unittest

import numpy as np

import agent7_st052m_minimum_existing_control_subset as audit_mod


class MinimumExistingControlSubsetTests(unittest.TestCase):
    def _audit(self, raw: np.ndarray, unstable: set[str] | None = None) -> dict:
        raw = np.asarray(raw, dtype=float)
        self.assertEqual(raw.shape, (len(audit_mod.parent.DEFECT_COORDINATES), len(audit_mod.CHANNELS)))
        norms = np.linalg.norm(raw, axis=0)
        normalized = np.zeros_like(raw)
        for index, norm in enumerate(norms):
            if norm > np.finfo(float).tiny:
                normalized[:, index] = raw[:, index] / norm
        unstable = unstable or set()
        coordinates = {}
        for row_index, coordinate in enumerate(audit_mod.parent.DEFECT_COORDINATES):
            coordinates[coordinate] = {
                "fine_raw_derivative_by_channel": {
                    channel: float(raw[row_index, column_index])
                    for column_index, channel in enumerate(audit_mod.CHANNELS)
                },
                "normalized_column_alignment_by_channel": {
                    channel: float(normalized[row_index, column_index])
                    for column_index, channel in enumerate(audit_mod.CHANNELS)
                },
            }
        return {
            "coordinates": coordinates,
            "derivative_stability": {
                channel: {"passes": channel not in unstable}
                for channel in audit_mod.CHANNELS
            },
        }

    @staticmethod
    def _five_required_matrix() -> np.ndarray:
        raw = np.zeros((8, 5), dtype=float)
        raw[:5, :] = np.eye(5)
        raw[5:, :] = 1.0
        return raw

    @staticmethod
    def _four_sufficient_matrix() -> np.ndarray:
        raw = np.zeros((8, 5), dtype=float)
        raw[:4, :4] = np.eye(4)
        raw[4:, :4] = 1.0
        raw[:, 4] = raw[:, 0]
        return raw

    def test_all_five_required_when_each_has_unique_coordinate_response(self) -> None:
        result = audit_mod.analyze_subsets(self._audit(self._five_required_matrix()))
        self.assertEqual(result["all_nonempty_subsets_evaluated"], 31)
        self.assertEqual(result["minimum_admissible_subset_size"], 5)
        self.assertEqual(result["minimum_admissible_subsets"], [list(audit_mod.CHANNELS)])
        self.assertTrue(result["all_five_channels_required_under_frozen_gates"])
        for omitted, row in result["leave_one_out_four_channel_audit"].items():
            self.assertFalse(row["passes_all_frozen_gates"], omitted)
            self.assertIn("coordinate_response", row["failed_gates"])

    def test_redundant_fifth_channel_yields_four_channel_minimum(self) -> None:
        result = audit_mod.analyze_subsets(self._audit(self._four_sufficient_matrix()))
        self.assertEqual(result["minimum_admissible_subset_size"], 4)
        self.assertIn(list(audit_mod.CHANNELS[:4]), result["minimum_admissible_subsets"])
        self.assertTrue(result["existing_five_channel_family_locally_overcomplete_for_these_coordinates"])
        self.assertFalse(result["all_five_channels_required_under_frozen_gates"])

    def test_unresponsive_coordinate_yields_no_admissible_subset(self) -> None:
        raw = self._five_required_matrix()
        raw[7, :] = 0.0
        result = audit_mod.analyze_subsets(self._audit(raw))
        self.assertIsNone(result["minimum_admissible_subset_size"])
        self.assertTrue(result["no_admissible_existing_subset"])
        self.assertTrue(all(
            "coordinate_response" in row["failed_gates"]
            for row in result["subsets_by_size"]["5"]
        ))

    def test_unstable_channel_cannot_be_used_to_satisfy_subset(self) -> None:
        result = audit_mod.analyze_subsets(
            self._audit(self._five_required_matrix(), unstable={"amplitude"})
        )
        self.assertIsNone(result["minimum_admissible_subset_size"])
        full = result["subsets_by_size"]["5"][0]
        self.assertIn("derivative_stability", full["failed_gates"])

    def test_subset_metrics_rejects_empty_or_unknown_subset(self) -> None:
        audit = self._audit(self._five_required_matrix())
        with self.assertRaises(ValueError):
            audit_mod.subset_metrics(audit, ())
        with self.assertRaises(ValueError):
            audit_mod.subset_metrics(audit, ("not_a_channel",))

    def test_truth_boundary_keeps_velocity_and_basis_unchanged(self) -> None:
        self.assertFalse(audit_mod.TRUTH["candidate_velocity_changed"])
        self.assertFalse(audit_mod.TRUTH["basis_dimension_changed"])
        self.assertFalse(audit_mod.TRUTH["held_out_pde_residual_evaluated"])
        self.assertEqual(audit_mod.TRUTH["direct_visualization_fingerprint_improvement"], 0.0)
        self.assertFalse(audit_mod.TRUTH["pde_validated"])
        self.assertFalse(audit_mod.TRUTH["openai_field_identified"])


if __name__ == "__main__":
    unittest.main()
