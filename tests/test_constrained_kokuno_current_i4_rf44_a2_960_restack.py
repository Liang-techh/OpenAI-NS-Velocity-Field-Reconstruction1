from __future__ import annotations

import hashlib
import importlib.util
import inspect
from pathlib import Path

import numpy as np

import openai_ns_reconstruction.kokuno_current_i4_exact_backend_correction_firewall as fw
import openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution as mean
import openai_ns_reconstruction.kokuno_current_i4_rf44_prestate_provider as pre
import openai_ns_reconstruction.kokuno_oscillatory_batch_axis_safety as axis
import openai_ns_reconstruction.kokuno_oscillatory_batch_differentials as diff
import openai_ns_reconstruction.kokuno_public_oscillatory_vorticity_diagnostic as vort


EXPECTED = {
    "diff": "12df3bf6c949baeebf609b366f2973e7f515a157",
    "axis": "599baa190ec742d0e532e5517078534124e68d6b",
    "vort": "4ae525bad1e4b83c9dcabc1a97d3931562d099a5",
    "pre": "0d8bb6471a4b7825b1b4b428d3b2c38d207729fd",
    "firewall": "0152a85f7fbaf9c32c24c9d5c7d0d11166e3b157",
}


def _git_blob_sha1(module: object) -> str:
    source = inspect.getsourcefile(module)
    assert source is not None
    raw = Path(source).read_bytes()
    return hashlib.sha1(f"blob {len(raw)}\0".encode("ascii") + raw).hexdigest()


def test_exact_a2_960_runtime_and_rf44_prestate_share_one_lineage() -> None:
    assert _git_blob_sha1(diff) == EXPECTED["diff"]
    assert _git_blob_sha1(axis) == EXPECTED["axis"]
    assert _git_blob_sha1(vort) == EXPECTED["vort"]
    assert _git_blob_sha1(pre) == EXPECTED["pre"]
    assert _git_blob_sha1(fw) == EXPECTED["firewall"]

    assert diff.evaluate_oscillatory_batch_differentials.__module__ == (
        "openai_ns_reconstruction.kokuno_oscillatory_batch_differentials"
    )
    assert diff.evaluate_oscillatory_batch_differentials.__name__ == (
        "evaluate_oscillatory_batch_differentials"
    )
    assert mean.AGENT2_DIFFERENTIAL_SOURCE_BLOB == EXPECTED["diff"]
    assert mean.AGENT2_DIFFERENTIAL_HEAD == "6d2fb1f701a34f783dca15a267ae2ce0734ba741"

    contract = diff.public_contract()
    assert contract["fixed_spatial_step"] == 1.0e-3
    assert contract["caller_tunable_step"] is False
    assert contract["forbidden_inputs_present"] == []
    semantic = diff.batch_differentials_sha256()
    assert isinstance(semantic, str) and len(semantic) == 64

    result = diff.evaluate_oscillatory_batch_differentials(
        np.asarray(((0.0, 0.0, 0.0), (3.5, 0.0, 0.2)), dtype=float),
        np.asarray((0.51, 0.51), dtype=float),
    )
    assert np.array_equal(result.interior_mask, np.asarray((False, False)))
    assert np.array_equal(result.velocity, np.zeros((2, 3)))
    assert np.array_equal(result.velocity_jacobian, np.zeros((2, 3, 3)))
    assert np.array_equal(result.divergence, np.zeros(2))
    assert np.array_equal(result.vorticity, np.zeros((2, 3)))

    pre_truth = pre.truth_boundary()
    fw_truth = fw.truth_boundary()
    assert pre_truth["rf44_preupdate_state_materialized"] is True
    assert fw_truth["exact_backend_rebind_required"] is True
    assert fw_truth["manual_ExactCurrentI4NonlinearBackend_dataclass_construction_is_sufficient"] is False
    assert fw_truth["zero_or_synthetic_backend_fixture_is_scientific_correction_evidence"] is False
    assert pre_truth["final_normalized_momentum_gate"] == 1.0e-3
    assert fw_truth["final_normalized_momentum_gate"] == 1.0e-3
    assert pre_truth["final_normalized_divergence_gate"] == 1.0e-5
    assert fw_truth["final_normalized_divergence_gate"] == 1.0e-5
    assert pre_truth["residual_defined_free_forcing_allowed"] is False
    assert fw_truth["residual_defined_free_forcing_allowed"] is False


def test_exact_a2_1080_composite_remains_the_next_fail_closed_blocker() -> None:
    # This restack intentionally closes only the #1224 + #1226 topology gap.
    # Scientific RF30/RF39 promotion still requires the exact A2 #1080
    # composite runtime and its authenticated A1 dependency closure.
    assert importlib.util.find_spec(
        "openai_ns_reconstruction.kokuno_current_i4_leading_oscillatory_identity"
    ) is None
    assert mean.AGENT2_COMPOSITE_HEAD == "c40d8ddecd2971544a6e07dab093436b423cf326"
    assert mean.AGENT2_COMPOSITE_SOURCE_BLOB == "2a0a5aa5966b02da856bdcf51940f3c186042802"

    fw_truth = fw.truth_boundary()
    assert fw_truth["repository_candidate_rf34_rf39_correction_materializable_only_after_exact_rebind"] is True
    assert fw_truth["correction_applied_to_candidate"] is False
    assert fw_truth["rf44_rf49_nonlinear_remainder_recomputed"] is False
    assert fw_truth["cartesian_correction_velocity_materialized"] is False
    assert fw_truth["finite_correction_cycle_run"] is False
    assert fw_truth["heldout_ns_residual_assessed"] is False
    assert fw_truth["pde_validated"] is False
