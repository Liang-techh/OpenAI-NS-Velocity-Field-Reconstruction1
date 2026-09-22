from __future__ import annotations

from copy import deepcopy
import hashlib
import unittest

from experiments.root_st052 import agent7_st052m_stable_identity_freshness as m


def h(label: str) -> str:
    return hashlib.sha256(label.encode()).hexdigest()


def pr_meta() -> dict:
    return {"number": m.STABLE_PR, "state": "open", "head": {"sha": m.LATEST_STABLE_HEAD}}


def compare_meta() -> dict:
    return {
        "status": "ahead",
        "ahead_by": m.EXPECTED_AHEAD_BY,
        "behind_by": 0,
        "base_commit": {"sha": m.OLD_STABLE_HEAD},
        "head_commit": {"sha": m.LATEST_STABLE_HEAD},
        "files": [
            {"filename": name, "status": "modified"}
            for name in sorted(m.EXPECTED_CHANGED_PATHS)
        ],
    }


def run_meta(status: str = "queued", conclusion=None) -> dict:
    return {
        "id": m.LATEST_STABLE_RUN,
        "head_sha": m.LATEST_STABLE_HEAD,
        "status": status,
        "conclusion": conclusion,
    }


def stable_receipt() -> dict:
    candidate_id = h("candidate")
    velocity_payload = {
        "candidate_semantic_identity_sha256": candidate_id,
        "constrained_callable_runtime": {
            "repository": m.CONSTRAINED_RUNTIME_REPOSITORY,
            "commit": m.CONSTRAINED_RUNTIME_IMPLEMENTATION_COMMIT,
        },
    }
    return {
        "schema": m.STABLE_SCHEMA,
        "task_id": m.STABLE_TASK_ID,
        "candidate_id": m.CANDIDATE_ID,
        "stable_identity": {
            "candidate_semantic_identity_sha256": candidate_id,
            "velocity_semantic_identity_sha256": h("velocity"),
            "velocity_semantic_payload": velocity_payload,
        },
        "truth_boundary": {
            "candidate_changed": False,
            "velocity_coefficients_changed": False,
            "constrained_callable_implementation_identity_bound": True,
            "velocity_export_ready": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "source_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


class FreshnessTests(unittest.TestCase):
    def test_queued_latest_run_blocks_without_promoting(self):
        r = m.build_receipt(
            upstream_pr=pr_meta(),
            upstream_compare=compare_meta(),
            latest_run=run_meta(),
        )
        self.assertEqual(r["decision"]["status"], "blocked_latest_stable_identity_ci")
        self.assertTrue(r["decision"]["fresh_axial_rerun_required"])
        self.assertFalse(r["historical_failed_axial_evidence"]["admitted_for_basis_routing"])
        self.assertEqual(r["truth"]["direct_visualization_fingerprint_improvement"], 0.0)

    def test_success_requires_latest_stable_receipt(self):
        with self.assertRaisesRegex(ValueError, "requires its receipt"):
            m.build_receipt(
                upstream_pr=pr_meta(),
                upstream_compare=compare_meta(),
                latest_run=run_meta("completed", "success"),
            )

    def test_success_still_requires_fresh_axial_rerun(self):
        r = m.build_receipt(
            upstream_pr=pr_meta(),
            upstream_compare=compare_meta(),
            latest_run=run_meta("completed", "success"),
            stable_receipt=stable_receipt(),
        )
        self.assertEqual(
            r["decision"]["status"],
            "latest_stable_identity_ready_fresh_axial_rerun_required",
        )
        self.assertFalse(r["decision"]["handoff_to_existing_axial_router_authorized"])
        self.assertEqual(
            r["latest_stable_identity"]["constrained_callable_implementation_commit"],
            m.CONSTRAINED_RUNTIME_IMPLEMENTATION_COMMIT,
        )

    def test_stale_upstream_head_fails_closed(self):
        p = pr_meta()
        p["head"]["sha"] = m.OLD_STABLE_HEAD
        with self.assertRaisesRegex(ValueError, "rolled again"):
            m.build_receipt(
                upstream_pr=p, upstream_compare=compare_meta(), latest_run=run_meta()
            )

    def test_unexpected_upstream_path_fails_closed(self):
        c = compare_meta()
        c["files"].append(
            {"filename": "src/openai_ns_reconstruction/velocity.py", "status": "modified"}
        )
        with self.assertRaisesRegex(ValueError, "unexpected upstream changed path"):
            m.build_receipt(
                upstream_pr=pr_meta(), upstream_compare=c, latest_run=run_meta()
            )

    def test_latest_callable_commit_drift_fails_closed(self):
        s = stable_receipt()
        s["stable_identity"]["velocity_semantic_payload"]["constrained_callable_runtime"]["commit"] = "1" * 40
        with self.assertRaisesRegex(ValueError, "implementation commit drift"):
            m.build_receipt(
                upstream_pr=pr_meta(),
                upstream_compare=compare_meta(),
                latest_run=run_meta("completed", "success"),
                stable_receipt=s,
            )

    def test_truth_promotion_fails_closed(self):
        s = stable_receipt()
        s["truth_boundary"]["visualization_ready"] = True
        with self.assertRaisesRegex(ValueError, "truth promotion forbidden"):
            m.build_receipt(
                upstream_pr=pr_meta(),
                upstream_compare=compare_meta(),
                latest_run=run_meta("completed", "success"),
                stable_receipt=s,
            )

    def test_unresolved_run_rejects_receipt(self):
        with self.assertRaisesRegex(ValueError, "cannot admit"):
            m.build_receipt(
                upstream_pr=pr_meta(),
                upstream_compare=compare_meta(),
                latest_run=run_meta(),
                stable_receipt=stable_receipt(),
            )


if __name__ == "__main__":
    unittest.main()
