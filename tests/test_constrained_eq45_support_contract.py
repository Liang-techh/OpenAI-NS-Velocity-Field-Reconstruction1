import math

import pytest

from openai_ns_reconstruction.constrained_eq45_support_contract import (
    audit_eq45_similarity_support,
)


def test_pr48_default_cutoffs_do_not_imply_cr001_physical_support():
    report = audit_eq45_similarity_support(
        X_cutoff=4.0,
        eta_cutoff=1.0,
        h=0.005,
        time_interval=(0.25, 0.75),
        physical_radial_support=2.0,
        physical_axial_support=2.0,
    )

    assert report.finite_physical_extent is False
    assert report.q_max is None
    assert report.radial_max is None
    assert report.axial_max is None
    assert report.eta_zero_radial_witness == pytest.approx(math.sqrt(6.0))
    assert report.max_X_cutoff_at_eta_zero == pytest.approx(8.0 / 3.0)
    assert report.similarity_support_alone_satisfies_physical_support is False
    assert "eta_zero_radial_extent_exceeds_registered_support" in report.reasons
    assert "eta_cutoff_reaches_or_crosses_unbounded_physical_branch_edge" in report.reasons
    assert report.visual_correspondence_verified is False
    assert report.pde_validated is False
    assert report.paper_exact is False

    with pytest.raises(ValueError, match="does not by itself satisfy"):
        audit_eq45_similarity_support(
            X_cutoff=4.0,
            eta_cutoff=1.0,
            h=0.005,
            time_interval=(0.25, 0.75),
            physical_radial_support=2.0,
            physical_axial_support=2.0,
            require_physical_support=True,
        )


def test_stricter_similarity_cutoffs_can_map_inside_registered_support():
    report = audit_eq45_similarity_support(
        X_cutoff=1.0,
        eta_cutoff=0.5,
        h=0.005,
        time_interval=(0.25, 0.75),
        physical_radial_support=2.0,
        physical_axial_support=2.0,
        require_physical_support=True,
    )

    assert report.q_max == pytest.approx(1.0)
    assert report.radial_max == pytest.approx(math.sqrt(2.0))
    assert report.axial_max == pytest.approx(0.5)
    assert report.similarity_support_alone_satisfies_physical_support is True


def test_finite_eta_cutoff_can_still_fail_radial_or_axial_support():
    radial = audit_eq45_similarity_support(
        X_cutoff=4.0,
        eta_cutoff=0.5,
        h=0.005,
        time_interval=(0.25, 0.75),
        physical_radial_support=2.0,
        physical_axial_support=2.0,
    )
    assert radial.finite_physical_extent is True
    assert radial.radial_max > 2.0
    assert radial.similarity_support_alone_satisfies_physical_support is False

    axial = audit_eq45_similarity_support(
        X_cutoff=0.1,
        eta_cutoff=0.95,
        h=0.005,
        time_interval=(0.25, 0.75),
        physical_radial_support=10.0,
        physical_axial_support=2.0,
    )
    assert axial.finite_physical_extent is True
    assert axial.axial_max > 2.0
    assert "axial_extent_exceeds_registered_support" in axial.reasons


def test_fail_closed_on_invalid_contract_inputs():
    common = dict(
        X_cutoff=1.0,
        eta_cutoff=0.5,
        h=0.005,
        time_interval=(0.25, 0.75),
        physical_radial_support=2.0,
        physical_axial_support=2.0,
    )

    for key, bad_value in (
        ("X_cutoff", 0.0),
        ("eta_cutoff", float("nan")),
        ("physical_radial_support", -1.0),
        ("physical_axial_support", float("inf")),
    ):
        values = dict(common)
        values[key] = bad_value
        with pytest.raises(ValueError):
            audit_eq45_similarity_support(**values)

    values = dict(common)
    values["h"] = 0.5
    with pytest.raises(ValueError, match="h must"):
        audit_eq45_similarity_support(**values)

    values = dict(common)
    values["time_interval"] = (0.75, 0.25)
    with pytest.raises(ValueError, match="time_interval"):
        audit_eq45_similarity_support(**values)
