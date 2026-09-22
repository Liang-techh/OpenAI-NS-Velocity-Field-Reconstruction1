from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction import st052_stable_identity_bound_streamline_trajectory as trajectory
from openai_ns_reconstruction import st052_streamline_inward_control_sensitivity as subject


def test_frozen_scope_and_truth_boundary() -> None:
    p = subject._protocol()
    assert p["source_streamline_pr"] == 1212
    assert p["source_control_pr"] == 1078
    assert p["channels"] == list(subject.CHANNELS)
    assert p["epsilons"] == [1.0e-3, 5.0e-4]
    assert p["primary_epsilon"] == 5.0e-4
    assert p["public_openai_numeric_target"] is None
    assert subject.TRUTH_BOUNDARY["candidate_changed"] is False
    assert subject.TRUTH_BOUNDARY["basis_dimension_changed"] is False
    assert subject.TRUTH_BOUNDARY["direct_visualization_fingerprint_improvement"] == 0.0
    assert subject.TRUTH_BOUNDARY["pde_validated"] is False
    assert subject.TRUTH_BOUNDARY["openai_field_identified"] is False


def test_baseline_transfer_gate_is_fail_closed() -> None:
    stable = np.zeros((3, 48, 3), dtype=float)
    same = stable.copy()
    passed = subject._compatibility(stable, same)
    assert passed["passes"] is True
    drifted = same.copy()
    drifted[1, 7, 2] = 2.0e-9
    failed = subject._compatibility(stable, drifted)
    assert failed["passes"] is False
    assert failed["max_absolute_velocity_mismatch"] == pytest.approx(2.0e-9)


def test_coherent_leverage_requires_sign_coherence_and_stability() -> None:
    stable = {"passes": True}
    yes, sign = subject._coherent_leverage(np.array([1.0, 2.0, 3.0]), stable)
    assert yes is True and sign == "negative"
    yes, sign = subject._coherent_leverage(np.array([-1.0, -2.0, -3.0]), stable)
    assert yes is True and sign == "positive"
    yes, sign = subject._coherent_leverage(np.array([1.0, -2.0, 3.0]), stable)
    assert yes is False and sign is None
    yes, sign = subject._coherent_leverage(np.array([1.0, 2.0, 3.0]), {"passes": False})
    assert yes is False and sign is None
    yes, sign = subject._coherent_leverage(np.array([1.0e-12, 2.0, 3.0]), stable)
    assert yes is False and sign is None


def test_local_linearity_gate_uses_frozen_drift_and_cosine() -> None:
    fine = np.array([1.0, 2.0, 3.0])
    close = 1.001 * fine
    row = subject._channel_stability(fine, close)
    assert row["passes"] is True
    far = np.array([1.0, -2.0, 3.0])
    row = subject._channel_stability(fine, far)
    assert row["passes"] is False


def _synthetic_grid(*, radial_sign: float) -> tuple[np.ndarray, np.ndarray]:
    axis = np.linspace(trajectory.BOX[0], trajectory.BOX[1], trajectory.GRID_RESOLUTION)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    rr = np.hypot(xx, yy)
    denom = np.maximum(rr, 0.2)
    erx = xx / denom
    ery = yy / denom
    etx = -yy / denom
    ety = xx / denom
    radial = 0.20 * float(radial_sign)
    swirl = 0.85
    axial = 0.25 * np.tanh(2.0 * zz)
    grid = np.stack(
        (radial * erx + swirl * etx, radial * ery + swirl * ety, axial), axis=-1
    )
    return axis, grid


def test_forward_metric_distinguishes_inward_from_outward_helices() -> None:
    axis, inward = _synthetic_grid(radial_sign=-1.0)
    m_in = subject._forward_metrics_from_grid(axis, inward)
    axis, outward = _synthetic_grid(radial_sign=1.0)
    m_out = subject._forward_metrics_from_grid(axis, outward)
    assert m_in["overall"]["mean_endpoint_radial_displacement"] < 0.0
    assert m_out["overall"]["mean_endpoint_radial_displacement"] > 0.0
    assert m_in["overall"]["mean_absolute_forward_turns"] > 0.0
    assert m_out["overall"]["mean_absolute_forward_turns"] > 0.0


def test_control_bundle_freezes_seed_grid_and_channel_contract(tmp_path: Path) -> None:
    seeds = trajectory.seed_points()
    axis = np.linspace(trajectory.BOX[0], trajectory.BOX[1], trajectory.GRID_RESOLUTION)
    shape = (trajectory.GRID_RESOLUTION,) * 3 + (3,)
    payload = {
        "seed_points": seeds,
        "baseline_seed_velocity": np.zeros((3, 48, 3)),
        "axis": axis,
    }
    payload.update({f"tangent_{name}": np.zeros(shape) for name in subject.CHANNELS})
    path = tmp_path / "controls.npz"
    np.savez_compressed(path, **payload)
    loaded = subject._load_control_bundle(path)
    assert set(loaded) == set(payload)

    bad = dict(payload)
    bad["seed_points"] = seeds.copy()
    bad["seed_points"][0, 0] += 1.0e-4
    bad_path = tmp_path / "bad.npz"
    np.savez_compressed(bad_path, **bad)
    with pytest.raises(ValueError, match="seed contract drift"):
        subject._load_control_bundle(bad_path)


def test_control_metadata_is_exact_head_bound(tmp_path: Path) -> None:
    meta = {
        "source_control_pr": 1078,
        "source_control_head": subject.SOURCE_CONTROL_HEAD,
        "channels": list(subject.CHANNELS),
        "times": list(subject.TIMES_FOR_TRANSFER),
        "time_for_tangent_grids": trajectory.TIME,
        "grid_resolution": trajectory.GRID_RESOLUTION,
        "note": "extra provenance is allowed",
    }
    path = tmp_path / "meta.json"
    path.write_text(json.dumps(meta))
    assert subject._load_control_metadata(path)["source_control_head"] == subject.SOURCE_CONTROL_HEAD
    meta["source_control_head"] = "0" * 40
    path.write_text(json.dumps(meta))
    with pytest.raises(ValueError, match="control metadata drift"):
        subject._load_control_metadata(path)


def test_routing_never_mutates_candidate_or_adds_basis() -> None:
    incompatible = {"passes": False}
    d = subject._decision(incompatible, None)
    assert d["local_five_control_trajectory_obstruction_established"] is False
    assert d["new_basis_added"] is False
    assert d["candidate_mutation_authorized"] is False

    compatible = {"passes": True}
    with_existing = {"existing_control_leverage_established": True}
    d = subject._decision(compatible, with_existing)
    assert d["existing_control_leverage_established"] is True
    assert d["local_five_control_trajectory_obstruction_established"] is False
    assert d["new_basis_added"] is False

    without_existing = {"existing_control_leverage_established": False}
    d = subject._decision(compatible, without_existing)
    assert d["local_five_control_trajectory_obstruction_established"] is True
    assert d["new_basis_authorized_in_this_increment"] is False
    assert d["actual_velocity_changed"] is False
    assert d["direct_visualization_fingerprint_improvement"] == 0.0
