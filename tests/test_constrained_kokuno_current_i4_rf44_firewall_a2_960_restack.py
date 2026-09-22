from __future__ import annotations

import hashlib
import inspect
from pathlib import Path

import numpy as np

import openai_ns_reconstruction.kokuno_current_i4_exact_backend_correction_firewall as firewall
import openai_ns_reconstruction.kokuno_current_i4_rf44_prestate_provider as prestate
from openai_ns_reconstruction import kokuno_oscillatory_batch_axis_safety as axis
from openai_ns_reconstruction import kokuno_oscillatory_batch_differentials as diff
from openai_ns_reconstruction import kokuno_public_candidate_velocity as public_candidate
from openai_ns_reconstruction import kokuno_public_carrier_resolved_velocity as carrier
from openai_ns_reconstruction import kokuno_public_oscillatory_vorticity_diagnostic as diagnostic
from openai_ns_reconstruction import kokuno_public_z_pullback_velocity as z_pullback
from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_DIFFERENTIAL_FUNCTION,
    AGENT2_DIFFERENTIAL_MODULE,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
)


PARENT_A2_1224_HEAD = "107cac08050fddeab778416c1caf019be0f5795b"
RF44_PRESTATE_BLOB = "0d8bb6471a4b7825b1b4b428d3b2c38d207729fd"
FIREWALL_BLOB = "0152a85f7fbaf9c32c24c9d5c7d0d11166e3b157"
A2_960_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
A3_1226_HEAD = "571abc47c98001886cdcc6042832635e7bd93f31"
EXPECTED_BLOBS = {
    "differential": "12df3bf6c949baeebf609b366f2973e7f515a157",
    "axis_safety": "599baa190ec742d0e532e5517078534124e68d6b",
    "vorticity_diagnostic": "4ae525bad1e4b83c9dcabc1a97d3931562d099a5",
    "z_pullback": "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653",
    "public_candidate": "a0e1b07981ba23a6bd547df14e8cdc89ad0d89ca",
    "carrier_resolved": "77edc85e4ac7aced997c90c2d93e7552e7d0c39b",
}


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    prefix = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(prefix + data).hexdigest()


def _module_blob(module: object) -> str:
    path = inspect.getsourcefile(module)
    assert path is not None
    return _git_blob_sha1(Path(path))


def test_rf44_firewall_and_exact_a2_960_runtime_coexist_unchanged() -> None:
    assert _module_blob(prestate) == RF44_PRESTATE_BLOB
    assert _module_blob(firewall) == FIREWALL_BLOB
    assert _module_blob(diff) == EXPECTED_BLOBS["differential"]
    assert _module_blob(axis) == EXPECTED_BLOBS["axis_safety"]
    assert _module_blob(diagnostic) == EXPECTED_BLOBS["vorticity_diagnostic"]
    assert _module_blob(z_pullback) == EXPECTED_BLOBS["z_pullback"]
    assert _module_blob(public_candidate) == EXPECTED_BLOBS["public_candidate"]
    assert _module_blob(carrier) == EXPECTED_BLOBS["carrier_resolved"]


def test_firewall_authentication_inputs_now_see_exact_a2_960_function() -> None:
    fn = diff.evaluate_oscillatory_batch_differentials
    assert fn.__module__ == AGENT2_DIFFERENTIAL_MODULE
    assert fn.__name__ == AGENT2_DIFFERENTIAL_FUNCTION
    source = inspect.getsourcefile(fn)
    assert source is not None
    assert _git_blob_sha1(Path(source)) == AGENT2_DIFFERENTIAL_SOURCE_BLOB
    assert AGENT2_DIFFERENTIAL_SOURCE_BLOB == EXPECTED_BLOBS["differential"]
    assert diff.FIXED_SPATIAL_STEP == 1.0e-3
    assert diff.public_contract()["caller_tunable_step"] is False
    assert diff.public_contract()["forbidden_inputs_present"] == []


def test_exact_a2_960_axis_and_support_exterior_stay_zero() -> None:
    t0, t1 = axis._time_interval()
    time = 0.5 * (t0 + t1)
    _, r1, _, z1 = axis._support()
    points = np.asarray(
        [
            (0.0, 0.0, 0.0),
            (2.0 * r1, 0.0, 0.0),
            (0.0, 0.0, 2.0 * max(abs(z1), 1.0)),
        ],
        dtype=float,
    )
    out = diff.evaluate_oscillatory_batch_differentials(points, time)
    assert not np.any(out.interior_mask)
    assert np.array_equal(out.velocity, np.zeros((3, 3)))
    assert np.array_equal(out.velocity_jacobian, np.zeros((3, 3, 3)))
    assert np.array_equal(out.divergence, np.zeros(3))
    assert np.array_equal(out.vorticity, np.zeros((3, 3)))


def test_restack_stops_before_a2_1080_and_real_exact_bind() -> None:
    a2 = prestate.truth_boundary()
    a3 = firewall.truth_boundary()
    assert a2["rf44_preupdate_state_materialized"] is True
    assert a2["rf44_postupdate_state_materialized"] is False
    assert a2["rf44_correction_increment_materialized"] is False
    assert a2["cartesian_correction_velocity_materialized"] is False
    assert a2["source_exact_rf44_prestate_claimed"] is False
    assert a3["exact_backend_rebind_required"] is True
    assert a3["correction_applied_to_candidate"] is False
    assert a3["cartesian_correction_velocity_materialized"] is False
    assert a3["heldout_ns_residual_assessed"] is False
    assert a3["pde_validated"] is False
    assert a2["final_normalized_momentum_gate"] == a3["final_normalized_momentum_gate"] == 1.0e-3
    assert a2["final_normalized_divergence_gate"] == a3["final_normalized_divergence_gate"] == 1.0e-5
    assert a2["residual_defined_free_forcing_allowed"] is False
    assert a3["residual_defined_free_forcing_allowed"] is False

    package_dir = Path(inspect.getsourcefile(prestate)).resolve().parent
    assert not (package_dir / "kokuno_current_i4_leading_oscillatory_identity.py").exists()

    receipt = {
        "task": "K2-OSC-115",
        "parent_a2_1224_head": PARENT_A2_1224_HEAD,
        "source_a2_960_head": A2_960_HEAD,
        "source_a3_1226_head": A3_1226_HEAD,
        "rf44_prestate_firewall_a2_960_same_lineage": True,
        "exact_a2_960_runtime_restored": True,
        "exact_a2_1080_runtime_restored": False,
        "real_exact_backend_bind_executed": False,
        "repository_candidate_scientific_correction_evidence": False,
        "rf44_postupdate_state_materialized": False,
        "cartesian_delta_u_materialized": False,
        "heldout_ns_residual_assessed": False,
        "pde_validated": False,
    }
    assert receipt["rf44_prestate_firewall_a2_960_same_lineage"] is True
    assert receipt["exact_a2_960_runtime_restored"] is True
    for key in (
        "exact_a2_1080_runtime_restored",
        "real_exact_backend_bind_executed",
        "repository_candidate_scientific_correction_evidence",
        "rf44_postupdate_state_materialized",
        "cartesian_delta_u_materialized",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ):
        assert receipt[key] is False
