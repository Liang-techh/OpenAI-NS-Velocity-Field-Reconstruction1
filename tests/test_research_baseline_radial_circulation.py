import json

import numpy as np
import pytest

import research_baseline.radial_circulation as rc


class FakeField:
    sha256 = rc.CANDIDATE_SHA256

    def at_points(self, points, time):
        p = np.asarray(points, dtype=float)
        x, y, _z = p.T
        r = np.hypot(x, y)
        c = np.divide(x, r, out=np.ones_like(x), where=r > 0)
        s = np.divide(y, r, out=np.zeros_like(y), where=r > 0)
        ur = -0.2 * r
        utheta = (1.0 + float(time)) * r
        uz = np.full_like(r, 0.3)
        return np.column_stack((ur * c - utheta * s, ur * s + utheta * c, uz))


def test_profile_contract_and_known_axisymmetric_field(monkeypatch):
    monkeypatch.setattr(rc, "load_best", lambda: FakeField())
    report = rc.measure_retained_st006_radial_circulation()
    assert report["candidate_sha256"] == rc.CANDIDATE_SHA256
    assert report["measurement_contract"]["times"] == [0.25, 0.5, 0.75]
    assert len(report["measurement_contract"]["radii"]) == 31
    assert report["interpretation_boundary"]["pde_validated"] is False
    assert report["interpretation_boundary"]["visual_correspondence_verified"] is False

    radii = np.asarray(report["measurement_contract"]["radii"])
    for entry, time in zip(report["per_time"], (0.25, 0.5, 0.75)):
        np.testing.assert_allclose(entry["mean_radial_velocity"], -0.2 * radii, atol=2e-15)
        np.testing.assert_allclose(
            entry["mean_circulating_velocity"], (1.0 + time) * radii, atol=2e-15
        )
        np.testing.assert_allclose(
            entry["circulation"], 2.0 * np.pi * (1.0 + time) * radii**2,
            rtol=2e-15,
            atol=2e-15,
        )
        assert entry["inward_mean_radial_ring_count"] == 31
        assert entry["outward_mean_radial_ring_count"] == 0
        assert entry["peak_abs_circulating_speed_radius"] == 1.5
        a, b = radii[0], radii[-1]
        expected_weighted_radius = ((b**3 - a**3) / 3.0) / ((b**2 - a**2) / 2.0)
        assert abs(entry["swirl_weighted_radius"] - expected_weighted_radius) < 1e-14
        assert entry["max_ring_component_deviation_over_velocity_rms"] < 1e-14

    assert report["endpoint_deltas"]["peak_abs_circulating_speed"] > 0.0
    payload = {k: v for k, v in report.items() if k != "report_sha256"}
    assert report["report_sha256"] == rc._report_digest(payload)


def test_bad_identity_fails_closed():
    field = FakeField()
    field.sha256 = "0" * 64
    with pytest.raises(ValueError, match="bound to retained ST006"):
        rc._sample_centerplane_profile(field)


def test_numerically_zero_swirl_fails_closed():
    class NoSwirl(FakeField):
        def at_points(self, points, time):
            p = np.asarray(points, dtype=float)
            x, y, _z = p.T
            r = np.hypot(x, y)
            c, s = x / r, y / r
            ur = -0.2 * r
            return np.column_stack((ur * c, ur * s, np.full_like(r, 0.3)))

    with pytest.raises(ValueError, match="swirl"):
        rc._sample_centerplane_profile(NoSwirl())


def test_nonfinite_and_bad_shape_fail_closed():
    class BadShape(FakeField):
        def at_points(self, points, time):
            return np.zeros((len(points), 2))

    with pytest.raises(ValueError, match="shape"):
        rc._sample_centerplane_profile(BadShape())

    class Nonfinite(FakeField):
        def at_points(self, points, time):
            out = np.zeros((len(points), 3))
            out[0, 0] = np.nan
            return out

    with pytest.raises(ValueError, match="finite"):
        rc._sample_centerplane_profile(Nonfinite())


def test_write_report_refuses_overwrite(monkeypatch, tmp_path):
    monkeypatch.setattr(rc, "load_best", lambda: FakeField())
    path = tmp_path / "radial.json"
    written = rc.write_report(path)
    assert json.loads(path.read_text())["report_sha256"] == written["report_sha256"]
    with pytest.raises(FileExistsError):
        rc.write_report(path)
