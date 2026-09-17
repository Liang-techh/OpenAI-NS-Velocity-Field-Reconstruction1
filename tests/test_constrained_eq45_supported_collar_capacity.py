import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_collar_capacity import (
    EXPECTED_PARENT_SHA256,
    PARAMETER_LABELS,
    TASK_ID,
    _with_phi_delta,
    audit_eq45_supported_collar_capacity,
)


def test_supported_collar_capacity_keeps_distinct_existing_controls(capsys):
    report = audit_eq45_supported_collar_capacity(
        coefficient_steps=(0.04, 0.02),
        times=(0.25, 0.50, 0.75),
        n_angles=8,
    )

    assert report["task_id"] == TASK_ID
    assert report["parent_sha256"] == EXPECTED_PARENT_SHA256
    assert report["supported_child_sha256"] != report["parent_sha256"]
    assert report["parameter_labels"] == list(PARAMETER_LABELS)
    assert report["coefficient_steps"] == pytest.approx([0.04, 0.02])
    assert report["times"] == pytest.approx([0.25, 0.50, 0.75])

    norms = np.asarray(report["finest_collar_response_norms"], dtype=float)
    assert norms.shape == (2,)
    assert np.all(np.isfinite(norms))
    assert np.all(norms > 1.0e-7)

    singular = np.asarray(report["singular_values"], dtype=float)
    assert singular.shape == (2,)
    assert np.all(np.isfinite(singular))
    assert report["numerical_rank"] == 2
    assert np.isfinite(report["condition_number"])
    assert report["condition_number"] < 1.0e6
    assert report["phi_1_2_novelty_outside_phi_1_0_span"] > 1.0e-3

    # The centered response is effectively linear in these bounded profile
    # coefficients.  A coarse->fine coefficient-step ladder must not reveal an
    # unstable local direction.
    refinement = np.asarray(report["step_refinement_relative_changes"], dtype=float)
    assert refinement.shape == (1,)
    assert np.all(np.isfinite(refinement))
    assert np.max(refinement) < 1.0e-8

    # The support transform is exactly the identity on the plateau, so the two
    # local basis-response directions must be inherited there to roundoff.
    plateau_change = np.asarray(
        report["plateau_supported_vs_parent_response_relative_change"], dtype=float
    )
    assert plateau_change.shape == (2,)
    assert np.max(plateau_change) < 1.0e-11

    # In the actual exterior collar the transform must genuinely change each
    # direction; otherwise the end/support connection has not reached that
    # basis channel at all.
    collar_change = np.asarray(
        report["collar_supported_vs_parent_response_relative_change"], dtype=float
    )
    assert collar_change.shape == (2,)
    assert np.all(np.isfinite(collar_change))
    assert np.all(collar_change > 1.0e-4)

    axial_share = np.asarray(report["axial_response_share"], dtype=float)
    assert axial_share.shape == (2,)
    assert np.all(np.isfinite(axial_share))
    assert np.all((axial_share > 0.0) & (axial_share < 1.0))

    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["production_coefficients_changed"] is False
    assert truth["new_basis_added"] is False
    assert truth["forcing_or_pressure_refit"] is False
    assert truth["physical_support_transform_changed"] is False
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False

    with capsys.disabled():
        print("AGENT7_SUPPORTED_COLLAR_CAPACITY_REPORT=" + json.dumps(report, sort_keys=True))


def test_supported_collar_capacity_fails_closed_on_identity_bounds_and_bad_inputs():
    parent = Eq45VelocityCandidate.seed()
    assert parent.sha256 == EXPECTED_PARENT_SHA256

    mutated = _with_phi_delta(parent, mode=(1, 0), delta=0.01)
    with pytest.raises(ValueError, match="frozen governed Eq45 parent identity"):
        audit_eq45_supported_collar_capacity(mutated)

    with pytest.raises(ValueError, match="strictly decreasing"):
        audit_eq45_supported_collar_capacity(coefficient_steps=(0.02, 0.04))
    with pytest.raises(ValueError, match="at least two"):
        audit_eq45_supported_collar_capacity(coefficient_steps=(0.02,))
    with pytest.raises(ValueError, match="inside the parent delivery interval"):
        audit_eq45_supported_collar_capacity(times=(0.25, 0.90))
    with pytest.raises(ValueError, match="n_angles"):
        audit_eq45_supported_collar_capacity(n_angles=3)
    with pytest.raises(ValueError, match="coefficient bound"):
        _with_phi_delta(parent, mode=(1, 0), delta=5.0)
    with pytest.raises(ValueError, match="mode must be"):
        _with_phi_delta(parent, mode=(0, 0), delta=0.01)
