from __future__ import annotations

import inspect
import runpy

import numpy as np

from openai_ns_reconstruction.kokuno_source_azimuthal_frame_oriented_multiharmonic_velocity import (
    KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity,
)
from openai_ns_reconstruction.kokuno_source_oriented_multiharmonic_diagnostics import (
    CHART_FD_STEPS,
    audit_oriented_multiharmonic_field,
    oriented_multiharmonic_diagnostics_contract,
    public_contract,
)


_ORIENTED_PARENT_TEST = (
    "tests/test_constrained_kokuno_source_azimuthal_frame_oriented_multiharmonic_velocity.py"
)


def _parent_field():
    ns = runpy.run_path(_ORIENTED_PARENT_TEST)
    return ns["_parent_field"]()


def _oriented(alpha=0.63):
    return KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        _parent_field(), frame_azimuth=alpha
    )


def test_alpha_zero_replay_has_exact_numerical_covariance_identity():
    report = audit_oriented_multiharmonic_field(_oriented(0.0))
    assert report["frame_azimuth"] == 0.0
    assert report["velocity_covariance_observation"]["base_probe_max_abs"] == 0.0
    assert report["offgrid_observation"]["velocity_covariance_max_abs"] == 0.0
    for row in report["resolutions"]:
        assert row["divergence_covariance_rms"] == 0.0
        assert row["gradient_covariance_relative_rms"] == 0.0
        assert row["vorticity_covariance_relative_rms"] == 0.0


def test_nonzero_orientation_three_resolution_covariance_is_finite_and_nontrivial():
    report = audit_oriented_multiharmonic_field(_oriented(0.63))
    assert [row["chart_step"] for row in report["resolutions"]] == list(CHART_FD_STEPS)
    assert report["velocity_covariance_observation"]["nonzero_oriented_probe_count"] > 0
    assert report["velocity_covariance_observation"]["base_probe_max_abs"] <= 5e-14
    assert report["offgrid_observation"]["velocity_covariance_max_abs"] <= 5e-14
    assert report["support_observation"]["common_axis_core_exact_zero"] is True
    assert report["support_observation"]["outer_spatial_support_exact_zero"] is True
    for row in report["resolutions"]:
        values = np.asarray(
            [
                row["reference_divergence_max_abs"],
                row["oriented_divergence_max_abs"],
                row["divergence_covariance_normalized_rms"],
                row["gradient_covariance_relative_rms"],
                row["vorticity_covariance_relative_rms"],
                row["oriented_vorticity_rms"],
                row["oriented_vorticity_max"],
                row["oriented_axial_vorticity_rms_fraction"],
            ]
        )
        assert np.all(np.isfinite(values))


def test_divergent_public_velocity_mutation_is_detected_by_independent_fd_path():
    field = _oriented(0.41)
    original = field.velocity

    def divergent_velocity(x, y, z, t):
        value = np.asarray(original(x, y, z, t), dtype=float).copy()
        x_b = np.broadcast_to(np.asarray(x, dtype=float), value.shape[:-1])
        value[..., 0] += 0.1 * x_b
        return value

    field.velocity = divergent_velocity
    report = audit_oriented_multiharmonic_field(field)
    finest = report["resolutions"][-1]
    assert finest["oriented_divergence_max_abs"] > 0.05
    assert finest["divergence_covariance_rms"] > 0.05


def test_report_is_deterministic_and_binds_exact_oriented_parent():
    field = _oriented(-0.37)
    a = audit_oriented_multiharmonic_field(field)
    b = audit_oriented_multiharmonic_field(field)
    assert a["diagnostic_sha256"] == b["diagnostic_sha256"]
    assert a == b
    assert a["parent"]["pr"] == 1170
    assert a["parent"]["head"] == "59d3e6d64fcdfeaf145cbee67e68b8a25e89489c"
    assert a["field_semantic_sha256"] == field.semantic_sha256


def test_truth_boundary_and_public_diagnostic_tuning_firewall():
    c = oriented_multiharmonic_diagnostics_contract()
    p = public_contract()
    assert c["oriented_public_velocity_three_resolution_diagnostic_materialized"] is True
    assert c["public_velocity_only_cartesian_derivatives_materialized"] is True
    assert c["rigid_rotation_divergence_covariance_numerically_assessed"] is True
    assert c["rigid_rotation_gradient_covariance_numerically_assessed"] is True
    assert c["rigid_rotation_vorticity_covariance_numerically_assessed"] is True
    assert c["provider_certified_axis_zero_support_checked"] is True
    assert c["provider_certified_outer_spatial_support_checked"] is True
    assert c["analytic_divergence_identity_proved_by_this_increment"] is False
    assert c["source_exact_orientation_recovered"] is False
    assert c["self_contained_velocity_xyzt_provider"] is False
    assert c["complete_ns_residual_assessed"] is False
    assert c["oscillation_before_after_ns_residual_compared"] is False
    assert c["residual_reduction_claimed"] is False
    assert c["paper_exact"] is False
    assert c["pde_validated"] is False
    assert p["forbidden_audit_inputs_present"] == []
    assert p["chart_fd_steps_frozen"] == list(CHART_FD_STEPS)
    assert set(inspect.signature(audit_oriented_multiharmonic_field).parameters) == {"field"}
