from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest

from research_baseline import CANDIDATE_SHA256
from research_baseline import vorticity_morphology as vm


def _rigid_rotation(points, time):
    points = np.asarray(points, dtype=float)
    return np.column_stack(
        (-2.0 * points[:, 1], 2.0 * points[:, 0], np.zeros(points.shape[0]))
    )


def _anisotropic_vorticity(points, time, scale=1.0):
    points = np.asarray(points, dtype=float)
    x = points[:, 0]
    return scale * np.column_stack(
        (np.zeros_like(x), x**3 / 3.0, np.zeros_like(x))
    )


def _row(velocity, *, time=0.5, resolution=17):
    return vm._diagnose_one_grid(
        velocity,
        time=time,
        resolution=resolution,
        candidate_sha256=CANDIDATE_SHA256,
        provenance="synthetic_test",
    )


def test_rigid_rotation_recovers_uniform_vorticity_and_isotropic_extent():
    row = _row(_rigid_rotation)
    assert row["vorticity_rms"] == pytest.approx(4.0, abs=1e-13)
    assert row["vorticity_max"] == pytest.approx(4.0, abs=1e-13)
    assert row["same_operator_divergence_rms"] == pytest.approx(0.0, abs=1e-14)
    np.testing.assert_allclose(row["enstrophy_weighted_centroid"], 0.0, atol=1e-14)
    assert row["principal_to_transverse_aspect_ratio"] == pytest.approx(1.0, abs=1e-13)
    eigenvalues = np.asarray(row["enstrophy_covariance_eigenvalues"])
    np.testing.assert_allclose(eigenvalues, eigenvalues[0], rtol=0.0, atol=1e-13)


def test_anisotropic_vorticity_has_intrinsic_principal_axis():
    row = _row(_anisotropic_vorticity)
    principal = np.asarray(row["principal_axis_registered_coordinates"])
    assert abs(principal[0]) > 1.0 - 1e-13
    assert abs(principal[1]) < 1e-13
    assert abs(principal[2]) < 1e-13
    assert row["principal_to_transverse_aspect_ratio"] > 1.4
    np.testing.assert_allclose(row["enstrophy_weighted_centroid"], 0.0, atol=1e-13)


def test_uniform_amplitude_scaling_preserves_geometry_and_scales_curl():
    base = _row(lambda points, time: _anisotropic_vorticity(points, time, 1.0))
    scaled = _row(lambda points, time: _anisotropic_vorticity(points, time, 3.25))
    for key in (
        "enstrophy_weighted_centroid",
        "enstrophy_covariance",
        "enstrophy_covariance_eigenvalues",
        "principal_axis_registered_coordinates",
    ):
        np.testing.assert_allclose(scaled[key], base[key], rtol=2e-13, atol=2e-13)
    assert scaled["principal_to_transverse_aspect_ratio"] == pytest.approx(
        base["principal_to_transverse_aspect_ratio"], rel=2e-13
    )
    assert scaled["vorticity_rms"] == pytest.approx(3.25 * base["vorticity_rms"], rel=2e-13)
    assert scaled["vorticity_max"] == pytest.approx(3.25 * base["vorticity_max"], rel=2e-13)


def test_retained_contract_report_is_deterministic_and_truth_bounded():
    class FakeField:
        sha256 = CANDIDATE_SHA256

        def at_points(self, points, time):
            return _anisotropic_vorticity(points, time, scale=1.0 + float(time))

    first = vm._measure_field(FakeField())
    second = vm._measure_field(FakeField())
    assert first == second
    assert first["candidate"] == "ST006"
    assert first["candidate_sha256"] == CANDIDATE_SHA256
    assert first["measurement_contract"]["times"] == [0.25, 0.5, 0.75]
    assert first["measurement_contract"]["resolutions"] == [17, 25, 33]
    assert first["interpretation_boundary"]["pde_validated"] is False
    assert first["interpretation_boundary"]["visualization_ready"] is False
    assert first["interpretation_boundary"]["visual_correspondence_verified"] is False
    assert first["interpretation_boundary"]["camera_fitted"] is False
    assert first["interpretation_boundary"]["used_for_pde_acceptance"] is False
    ratio = first["finest_resolution_endpoint_trend"]["vorticity_rms"]["ratio_final_over_initial"]
    assert ratio == pytest.approx(1.75 / 1.25, rel=2e-13)
    aspect_ratio = first["finest_resolution_endpoint_trend"][
        "principal_to_transverse_aspect_ratio"
    ]["ratio_final_over_initial"]
    assert aspect_ratio == pytest.approx(1.0, rel=2e-13)

    payload = dict(first)
    reported = payload.pop("report_sha256")
    expected = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    assert reported == expected


def test_fail_closed_inputs_and_report_overwrite(monkeypatch, tmp_path):
    with pytest.raises(ValueError, match="nonzero"):
        _row(lambda points, time: np.zeros((len(points), 3)))
    with pytest.raises(ValueError, match="shape"):
        _row(lambda points, time: np.zeros((len(points), 2)))
    with pytest.raises(ValueError, match="finite"):
        def bad(points, time):
            value = _rigid_rotation(points, time)
            value[0, 0] = np.nan
            return value
        _row(bad)
    with pytest.raises(ValueError, match="odd"):
        _row(_rigid_rotation, resolution=16)
    with pytest.raises(ValueError, match="registered"):
        _row(_rigid_rotation, time=0.2)
    with pytest.raises(ValueError, match="candidate_sha256"):
        vm._diagnose_one_grid(
            _rigid_rotation,
            time=0.5,
            resolution=17,
            candidate_sha256="bad",
            provenance="synthetic_test",
        )
    with pytest.raises(ValueError, match="provenance"):
        vm._diagnose_one_grid(
            _rigid_rotation,
            time=0.5,
            resolution=17,
            candidate_sha256=CANDIDATE_SHA256,
            provenance="",
        )

    class WrongField:
        sha256 = "0" * 64

    with pytest.raises(ValueError, match="retained ST006"):
        vm._measure_field(WrongField())

    report = {
        "schema": "synthetic",
        "report_sha256": "0" * 64,
        "interpretation_boundary": {"pde_validated": False},
    }
    monkeypatch.setattr(vm, "measure_retained_st006_vorticity_morphology", lambda: report)
    path = tmp_path / "report.json"
    written = vm.write_report(path)
    assert written == report
    assert json.loads(path.read_text()) == report
    with pytest.raises(FileExistsError, match="overwrite"):
        vm.write_report(path)
