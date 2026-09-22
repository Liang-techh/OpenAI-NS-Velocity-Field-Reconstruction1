from __future__ import annotations

from copy import deepcopy
import unittest

from experiments.root_st052 import agent7_st052m_stable_identity_handoff as mod


STABLE_CANDIDATE = "1" * 64
STABLE_VELOCITY = "2" * 64


def run(status="queued", conclusion=None):
    return {
        "id": mod.STABLE_RUN,
        "head_sha": mod.STABLE_HEAD,
        "status": status,
        "conclusion": conclusion,
    }


def stable_receipt():
    return {
        "schema": mod.STABLE_SCHEMA,
        "task_id": mod.STABLE_TASK_ID,
        "candidate_id": mod.CANDIDATE_ID,
        "stable_identity": {
            "candidate_semantic_identity_sha256": STABLE_CANDIDATE,
            "velocity_semantic_identity_sha256": STABLE_VELOCITY,
        },
        "materialization_evidence": {
            "included_in_stable_candidate_identity": False,
            "included_in_stable_velocity_identity": False,
        },
        "migration_boundary": {"legacy_identity_reused_as_stable_identity": False},
        "truth_boundary": {
            "stable_semantic_identity_ready": True,
            "velocity_export_ready": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def axial_receipt(identity=STABLE_VELOCITY, aspect=1.08):
    return {
        "stable_velocity_semantic_identity_sha256": identity,
        "protocol": {"source_numeric_targets_used": False},
        "measurement": {"metrics": {"full_aspect_ratio": aspect}},
        "truth_boundary": {"visualization_ready": False, "pde_validated": False},
    }


class StableIdentityHandoffTests(unittest.TestCase):
    def test_queued_stable_identity_blocks_without_promoting_old_measurement(self):
        r = mod.classify(run())
        self.assertEqual(r["decision"]["status"], "blocked_stable_identity_ci")
        self.assertFalse(r["historical_failed_axial_evidence"]["admitted_for_basis_routing"])
        self.assertFalse(r["decision"]["handoff_to_existing_axial_router_authorized"])
        self.assertEqual(r["truth"]["direct_visualization_fingerprint_improvement"], 0.0)

    def test_completed_failure_stays_blocked(self):
        r = mod.classify(run("completed", "failure"))
        self.assertEqual(r["decision"]["status"], "blocked_stable_identity_ci")
        self.assertFalse(r["decision"]["stable_identity_admitted"])

    def test_success_requires_exact_stable_receipt(self):
        with self.assertRaisesRegex(ValueError, "requires its exact receipt"):
            mod.classify(run("completed", "success"))

    def test_stable_identity_success_requires_fresh_axial_rerun(self):
        r = mod.classify(run("completed", "success"), stable_receipt())
        self.assertEqual(
            r["decision"]["status"],
            "stable_identity_ready_fresh_axial_rerun_required",
        )
        self.assertTrue(r["decision"]["stable_identity_admitted"])
        self.assertFalse(r["decision"]["fresh_axial_measurement_admitted"])
        self.assertFalse(r["decision"]["candidate_mutation_authorized"])

    def test_fresh_axial_identity_mismatch_fails_closed(self):
        r = mod.classify(
            run("completed", "success"),
            stable_receipt(),
            axial_receipt(identity="3" * 64),
        )
        self.assertEqual(r["decision"]["status"], "stable_identity_measurement_mismatch")
        self.assertFalse(r["decision"]["handoff_to_existing_axial_router_authorized"])
        self.assertFalse(r["decision"]["targeted_basis_preflight_authorized"])

    def test_exact_stable_identity_match_only_authorizes_router_handoff(self):
        r = mod.classify(
            run("completed", "success"), stable_receipt(), axial_receipt()
        )
        self.assertEqual(
            r["decision"]["status"],
            "stable_identity_axial_measurement_ready_for_existing_router",
        )
        self.assertTrue(r["decision"]["fresh_axial_measurement_admitted"])
        self.assertTrue(r["decision"]["handoff_to_existing_axial_router_authorized"])
        self.assertFalse(r["decision"]["existing_control_selection_authorized"])
        self.assertFalse(r["decision"]["targeted_basis_preflight_authorized"])
        self.assertFalse(r["decision"]["candidate_mutation_authorized"])

    def test_legacy_velocity_hash_cannot_be_relabelled_stable(self):
        receipt = stable_receipt()
        receipt["stable_identity"]["velocity_semantic_identity_sha256"] = (
            mod.LEGACY_REPLAY_VELOCITY_SHA256
        )
        with self.assertRaisesRegex(ValueError, "legacy velocity hash"):
            mod.classify(run("completed", "success"), receipt)

    def test_materialization_hash_leak_is_rejected(self):
        receipt = stable_receipt()
        receipt["materialization_evidence"]["included_in_stable_velocity_identity"] = True
        with self.assertRaisesRegex(ValueError, "materialization evidence leaked"):
            mod.classify(run("completed", "success"), receipt)

    def test_fresh_measurement_cannot_use_openai_numeric_target(self):
        receipt = axial_receipt()
        receipt["protocol"]["source_numeric_targets_used"] = True
        with self.assertRaisesRegex(ValueError, "public numerical target"):
            mod.classify(run("completed", "success"), stable_receipt(), receipt)

    def test_run_head_is_exactly_pinned(self):
        bad = run()
        bad["head_sha"] = "0" * 40
        with self.assertRaisesRegex(ValueError, "run head"):
            mod.classify(bad)


if __name__ == "__main__":
    unittest.main()
