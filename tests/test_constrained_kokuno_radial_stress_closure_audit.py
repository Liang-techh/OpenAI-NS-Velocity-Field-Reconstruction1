from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_actual_oscillatory_mean_stress import (
    _compact_radial_stress,
)
from openai_ns_reconstruction.kokuno_radial_stress_closure_audit import (
    ADMITTED_AGENT2_HEAD,
    INDEPENDENT_AGENT4_HEAD,
    RADIAL_COUNTS,
    _actual_ring_mean_sources,
    _audit_channel,
    audit_actual_radial_stress_closure,
)


def test_independent_radial_operator_closes_on_smooth_manufactured_source() -> None:
    radii = np.linspace(0.15, 1.35, 129)
    source = (radii - 0.15) ** 4 * (1.35 - radii) ** 4 * (
        1.0 + 0.2 * radii + 0.1 * radii**2
    )
    for exponent in (1, 2):
        audit = _audit_channel(
            radii,
            source,
            exponent=exponent,
            center=0.75,
            halfwidth=0.60,
        )
        assert np.isfinite(audit["operator_relative_rms"])
        assert np.isfinite(audit["operator_relative_max"])
        assert audit["operator_relative_rms"] < 5.0e-2
        assert audit["operator_relative_max"] < 1.0e-1
        assert audit["moment_complement_relative"] < 1.0e-12
        assert audit["edge_relative"] < 1.0e-12


def test_actual_radial_stress_audit_consumes_frozen_candidate_and_stays_fail_closed() -> None:
    receipt = audit_actual_radial_stress_closure()
    assert receipt["provenance"]["admitted_agent2_head"] == ADMITTED_AGENT2_HEAD
    assert receipt["provenance"]["independent_agent4_head"] == INDEPENDENT_AGENT4_HEAD
    assert [level["radial_count"] for level in receipt["levels"]] == list(RADIAL_COUNTS)
    assert isinstance(receipt["failed_guards"], list)

    for channel in ("theta_e2", "axial_e1"):
        convergence = receipt["convergence"][channel]
        assert len(convergence["relative_rms_by_radial_count"]) == len(RADIAL_COUNTS)
        assert len(convergence["refinement_ratios"]) == len(RADIAL_COUNTS) - 1
        assert np.all(np.isfinite(convergence["relative_rms_by_radial_count"]))
        assert np.all(np.asarray(convergence["refinement_ratios"]) > 0.0)
        assert np.isfinite(convergence["finest_relative_max"])
        assert convergence["finest_moment_complement_relative"] >= 0.0
        assert convergence["finest_edge_relative"] >= 0.0

    truth = receipt["truth_boundary"]
    assert truth["real_oscillatory_self_defect_component_consumed"] is True
    assert truth["compact_radial_stress_reconstruction_audited"] is True
    assert truth["full_same_cycle_composite_requested_stress_materialized"] is False
    assert truth["candidate_finite_head_mean_debt_materialized"] is False
    assert truth["signed_mean_inverse_input_ready"] is False
    assert truth["public_velocity_correction_materialized"] is False
    assert truth["finite_correction_cycle_rerun_allowed"] is False
    assert truth["finite_correction_cycle_run"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["pde_validated"] is False


def test_invalid_radial_resolution_fails_closed() -> None:
    with pytest.raises(ValueError):
        _actual_ring_mean_sources(15)


def test_existing_constructor_still_enforces_exact_zero_edges() -> None:
    radii = np.linspace(0.15, 1.35, 65)
    source = (radii - 0.15) ** 2 * (1.35 - radii) ** 2
    for exponent in (1, 2):
        built = _compact_radial_stress(
            radii,
            source,
            exponent=exponent,
            bump_center=0.75,
            bump_halfwidth=0.60,
        )
        assert built["stress_inner_edge"] == 0.0
        assert abs(built["stress_outer_edge"]) <= 2.0e-13
