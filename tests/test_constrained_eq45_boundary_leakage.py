import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE_PATH = ROOT / "artifacts/constrained/eq45_velocity_candidate_seed.json"
AUDIT_PATH = ROOT / "artifacts/constrained/eq45_boundary_leakage_audit_seed.json"


def _core_speed_rms(field, time, half_width=1.0, grid_size=31):
    dx = 2.0 * half_width / grid_size
    axis = -half_width + dx * (np.arange(grid_size, dtype=float) + 0.5)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((x, y, z), axis=-1).reshape(-1, 3)
    velocity = field.at_points(points, time)
    return float(np.sqrt(np.mean(np.sum(velocity * velocity, axis=1))))


def _cube_face_metrics(field, time, half_width, grid_size):
    axis = np.linspace(-half_width, half_width, grid_size)
    a, b = np.meshgrid(axis, axis, indexing="ij")
    faces = []
    normals = []
    for sign in (-1.0, 1.0):
        faces.append(
            np.stack((np.full_like(a, sign * half_width), a, b), axis=-1).reshape(-1, 3)
        )
        normals.append(np.array((sign, 0.0, 0.0)))
    for sign in (-1.0, 1.0):
        faces.append(
            np.stack((a, np.full_like(a, sign * half_width), b), axis=-1).reshape(-1, 3)
        )
        normals.append(np.array((0.0, sign, 0.0)))
    for sign in (-1.0, 1.0):
        faces.append(
            np.stack((a, b, np.full_like(a, sign * half_width)), axis=-1).reshape(-1, 3)
        )
        normals.append(np.array((0.0, 0.0, sign)))

    velocities = [field.at_points(points, time) for points in faces]
    speed_sq = np.concatenate([np.sum(value * value, axis=1) for value in velocities])
    normal_speed = np.concatenate(
        [velocities[index] @ normals[index] for index in range(len(velocities))]
    )
    return {
        "boundary_speed_rms": float(np.sqrt(np.mean(speed_sq))),
        "boundary_speed_max": float(np.sqrt(np.max(speed_sq))),
        "normal_speed_rms": float(np.sqrt(np.mean(normal_speed * normal_speed))),
    }


def test_eq45_boundary_leakage_checked_artifact_replays_public_velocity():
    report = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    field = Eq45VelocityCandidate.load_json(CANDIDATE_PATH)

    assert report["schema"] == "eq45_boundary_leakage_audit_v1"
    assert report["candidate_sha256"] == field.sha256
    assert report["claim_scope"] == "finite_window_boundary_tail_sampling_only"

    sampling = report["sampling"]
    assert sampling["face_grid_sizes"] == [33, 49, 65]
    assert sampling["cube_half_widths"] == [2.0, 2.5, 3.0]
    assert sampling["times"] == [0.25, 0.5, 0.75]

    computed = {}
    core_rms = {}
    for time in sampling["times"]:
        core_rms[time] = _core_speed_rms(
            field,
            time,
            sampling["core_reference_half_width"],
            sampling["core_reference_grid_size"],
        )
        for half_width in sampling["cube_half_widths"]:
            for grid_size in sampling["face_grid_sizes"]:
                computed[(time, half_width, grid_size)] = _cube_face_metrics(
                    field, time, half_width, grid_size
                )

    for expected in report["finest_grid_results"]:
        time = expected["time"]
        half_width = expected["half_width"]
        actual = computed[(time, half_width, 65)]
        assert np.isclose(core_rms[time], expected["core_speed_rms"], rtol=2e-10, atol=1e-13)
        for name in ("boundary_speed_rms", "boundary_speed_max", "normal_speed_rms"):
            assert np.isclose(actual[name], expected[name], rtol=2e-10, atol=1e-13)
        assert np.isclose(
            actual["boundary_speed_rms"] / core_rms[time],
            expected["boundary_rms_over_core"],
            rtol=2e-10,
            atol=1e-13,
        )
        assert np.isclose(
            actual["normal_speed_rms"] / core_rms[time],
            expected["normal_rms_over_core"],
            rtol=2e-10,
            atol=1e-13,
        )

    relative_changes = {name: [] for name in (
        "boundary_speed_rms", "boundary_speed_max", "normal_speed_rms"
    )}
    for time in sampling["times"]:
        finest_rms = []
        for half_width in sampling["cube_half_widths"]:
            medium = computed[(time, half_width, 49)]
            finest = computed[(time, half_width, 65)]
            finest_rms.append(finest["boundary_speed_rms"])
            for name in relative_changes:
                relative_changes[name].append(
                    abs(medium[name] - finest[name]) / max(abs(finest[name]), 1e-30)
                )
        assert finest_rms[0] > finest_rms[1] > finest_rms[2]

    summary = report["resolution_summary"]
    assert np.isclose(
        max(relative_changes["boundary_speed_rms"]),
        summary["max_relative_change_49_to_65_boundary_speed_rms"],
        rtol=2e-10,
        atol=1e-13,
    )
    assert np.isclose(
        max(relative_changes["normal_speed_rms"]),
        summary["max_relative_change_49_to_65_normal_speed_rms"],
        rtol=2e-10,
        atol=1e-13,
    )
    assert np.isclose(
        max(relative_changes["boundary_speed_max"]),
        summary["max_relative_change_49_to_65_boundary_speed_max"],
        rtol=2e-10,
        atol=1e-13,
    )
    assert max(relative_changes["boundary_speed_rms"]) < 0.01
    assert max(relative_changes["normal_speed_rms"]) < 0.01

    truth = report["truth_boundary"]
    assert truth["physical_support_validated"] is False
    assert truth["tail_beyond_largest_cube_resolved"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert truth["blowup_proved"] is False

    # This is a regression observation for the current frozen seed, not a universal threshold:
    # the registered L=2 plotting box has a visible early-time face tail, while an enlarged
    # L=3 window reduces the same RMS leakage by more than two orders of magnitude.
    early_l2 = computed[(0.25, 2.0, 65)]["boundary_speed_rms"] / core_rms[0.25]
    early_l3 = computed[(0.25, 3.0, 65)]["boundary_speed_rms"] / core_rms[0.25]
    assert early_l2 > 1e-2
    assert early_l3 < 1e-3
