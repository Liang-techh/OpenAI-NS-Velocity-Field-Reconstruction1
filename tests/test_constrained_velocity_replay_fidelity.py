import numpy as np
import pytest

from openai_ns_reconstruction.constrained_velocity_replay_fidelity import (
    ReplayLevel,
    audit_velocity_replay_fidelity,
)

SHA = "a" * 64
BOUNDS = np.array([[-0.8, 0.9], [-0.7, 0.6], [-0.5, 0.75]])
TIMES = np.array([0.25, 0.5, 0.75])


def reference(x, y, z, t):
    return np.column_stack((x + 2.0 * y + t, -0.5 * x + z - t, y - 2.0 * z + 0.25 * t))


def audit(levels, **kwargs):
    options = dict(
        bounds_xyz=BOUNDS,
        times=TIMES,
        points_per_time=64,
        seed=90210,
        candidate_sha256=SHA,
        provenance="frozen-grid replay audit fixture",
    )
    options.update(kwargs)
    return audit_velocity_replay_fidelity(reference, levels, **options)


def test_exact_replay_has_zero_error_and_read_only_common_probes():
    report = audit([ReplayLevel("32x36x40", (32, 36, 40), reference)])
    metric = report.metrics[0]
    assert metric.vector_rms_error == 0.0
    assert metric.vector_max_error == 0.0
    assert metric.normalized_vector_rms_error == 0.0
    assert np.array_equal(metric.component_rms_error, np.zeros(3))
    assert report.reference_rms_speed > 0.0
    assert len(report.probe_sha256) == 64
    for arr in (report.bounds_xyz, report.times, report.probe_points_xyz, metric.component_rms_error):
        assert not arr.flags.writeable
    assert report.metadata["derivatives_computed"] is False
    assert report.metadata["pde_acceptance_evidence"] is False
    assert report.metadata["automatic_resolution_selection"] is False
    assert report.metadata["visualization_ready"] is False
    assert report.metadata["pde_validated"] is False
    assert report.metadata["openai_field_identified"] is False


def test_known_offset_is_not_registered_or_scaled_away():
    offset = np.array([0.03, -0.04, 0.12])

    def shifted(x, y, z, t):
        return reference(x, y, z, t) + offset

    report = audit([ReplayLevel("coarse", (16, 18, 20), shifted)])
    metric = report.metrics[0]
    expected = float(np.linalg.norm(offset))
    assert metric.vector_rms_error == pytest.approx(expected, rel=1e-14)
    assert metric.vector_max_error == pytest.approx(expected, rel=1e-14)
    assert np.allclose(metric.component_rms_error, np.abs(offset), rtol=0, atol=1e-14)
    assert report.metadata["registration_or_camera_fit"] is False
    assert report.metadata["rescaling_or_component_fit"] is False


def test_refinement_trend_is_reported_without_automatic_selection():
    def level(error):
        def velocity(x, y, z, t):
            return reference(x, y, z, t) + np.column_stack(
                (error * np.ones_like(x), -0.5 * error * np.ones_like(x), 0.25 * error * np.ones_like(x))
            )
        return velocity

    report = audit(
        [
            ReplayLevel("coarse", (16, 18, 20), level(0.04)),
            ReplayLevel("medium", (32, 36, 40), level(0.01)),
            ReplayLevel("fine", (64, 72, 80), level(0.0025)),
        ]
    )
    rms = [m.vector_rms_error for m in report.metrics]
    assert rms[0] > rms[1] > rms[2] > 0.0
    assert report.metadata["automatic_resolution_selection"] is False


def test_all_callables_receive_identical_probes_and_seed_is_deterministic():
    calls = []

    def recorded(label):
        def velocity(x, y, z, t):
            calls.append((label, float(t), np.column_stack((x, y, z)).copy()))
            return reference(x, y, z, t)
        return velocity

    levels = [
        ReplayLevel("a", (16, 16, 16), recorded("a")),
        ReplayLevel("b", (32, 32, 32), recorded("b")),
    ]
    first = audit(levels)
    second = audit([ReplayLevel("a", (16, 16, 16), reference), ReplayLevel("b", (32, 32, 32), reference)])
    assert first.probe_sha256 == second.probe_sha256
    assert np.array_equal(first.probe_points_xyz, second.probe_points_xyz)
    by_label = {label: [] for label in ("a", "b")}
    for label, t, pts in calls:
        by_label[label].append((t, pts))
    assert [x[0] for x in by_label["a"]] == list(TIMES)
    assert [x[0] for x in by_label["b"]] == list(TIMES)
    for a, b in zip(by_label["a"], by_label["b"]):
        assert np.array_equal(a[1], b[1])
        assert np.array_equal(a[1], first.probe_points_xyz)


def test_fail_closed_on_zero_reference_and_bad_velocity_outputs():
    zero = lambda x, y, z, t: np.zeros((np.asarray(x).size, 3))
    with pytest.raises(ValueError, match="exact-zero"):
        audit_velocity_replay_fidelity(
            zero,
            [ReplayLevel("grid", (16, 16, 16), zero)],
            bounds_xyz=BOUNDS,
            times=TIMES,
            points_per_time=16,
            seed=1,
            candidate_sha256=SHA,
            provenance="zero reference",
        )

    def wrong_shape(x, y, z, t):
        return np.zeros((np.asarray(x).size, 2))

    with pytest.raises(ValueError, match="shape"):
        audit([ReplayLevel("bad", (16, 16, 16), wrong_shape)])

    def nonfinite(x, y, z, t):
        out = reference(x, y, z, t)
        out[0, 0] = np.nan
        return out

    with pytest.raises(ValueError, match="nonfinite"):
        audit([ReplayLevel("bad", (16, 16, 16), nonfinite)])


def test_fail_closed_on_invalid_contract_inputs():
    good = [ReplayLevel("grid", (16, 18, 20), reference)]
    bad_bounds = BOUNDS.copy()
    bad_bounds[0] = [1.0, -1.0]
    with pytest.raises(ValueError):
        audit(good, bounds_xyz=bad_bounds)
    with pytest.raises(ValueError):
        audit(good, times=np.array([0.5, 0.25]))
    with pytest.raises(ValueError):
        audit(good, points_per_time=24)
    with pytest.raises(ValueError):
        audit(good, seed=-1)
    with pytest.raises(ValueError):
        audit(good, candidate_sha256="bad")
    with pytest.raises(ValueError):
        audit(good, provenance="")
    with pytest.raises(ValueError):
        audit([ReplayLevel("a", (32, 32, 32), reference), ReplayLevel("b", (16, 64, 64), reference)])
    with pytest.raises(ValueError):
        audit([ReplayLevel("dup", (16, 16, 16), reference), ReplayLevel("dup", (32, 32, 32), reference)])
