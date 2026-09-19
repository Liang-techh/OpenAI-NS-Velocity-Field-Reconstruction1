from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_actual_oscillatory_mean_stress import (
    _compact_radial_stress,
)
from openai_ns_reconstruction.kokuno_radial_stress_closure_audit import (
    _actual_ring_mean_sources,
    _audit_channel,
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


def test_actual_candidate_source_can_feed_independent_radial_audit_without_surrogate() -> None:
    actual = _actual_ring_mean_sources(17)
    radii = np.asarray(actual["radii"], dtype=float)
    field = actual["field"]
    for source_name, exponent in (("theta_source", 2), ("axial_source", 1)):
        source = np.asarray(actual[source_name], dtype=float)
        assert source.shape == radii.shape
        assert np.all(np.isfinite(source))
        assert np.linalg.norm(source) > 0.0
        audit = _audit_channel(
            radii,
            source,
            exponent=exponent,
            center=field.radial_center,
            halfwidth=field.radial_halfwidth,
        )
        assert np.isfinite(audit["operator_relative_rms"])
        assert np.isfinite(audit["operator_relative_max"])
        assert audit["moment_complement_relative"] >= 0.0
        assert audit["edge_relative"] >= 0.0


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
