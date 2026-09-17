import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_optimized_energy import (
    DEFAULT_DOMAIN_LEVELS,
    DEFAULT_RESOLUTIONS,
    EXPECTED_CANONICAL_SHA256,
    audit_optimized_eq45_energy,
    finite_window_component_energy,
)


def test_optimized_eq45_energy_audit_replays_public_candidates_and_is_stable():
    report = audit_optimized_eq45_energy()

    assert report["schema"] == "eq45_optimized_energy_audit_v1"
    assert report["task_id"] == "CR007-EQ45-OPTIMIZED-ENERGY-AUDIT-015"
    assert report["canonical_candidate_sha256"] == EXPECTED_CANONICAL_SHA256
    assert report["optimized_candidate_sha256"] == (
        "b34902ee28b3965f245eceb5867df7409498063cb0ea80f1758931723b303c98"
    )
    assert report["source_optimization_task_id"] == "CR005-EQ45-PROFILE-FORCE-OPT-015"
    assert report["resolutions"] == list(DEFAULT_RESOLUTIONS)
    assert report["domain_levels"] == [list(level) for level in DEFAULT_DOMAIN_LEVELS]

    expected_totals = {
        0.25: (0.8942276, 1.1763595),
        0.50: (0.7335698, 0.9666628),
        0.75: (0.5225982, 0.6893733),
    }
    for row in report["finest_comparison"]:
        canonical_expected, optimized_expected = expected_totals[row["time"]]
        assert row["canonical_total_energy"] == pytest.approx(canonical_expected, rel=5e-5)
        assert row["optimized_total_energy"] == pytest.approx(optimized_expected, rel=5e-5)
        assert 0.30 < row["total_energy_relative_change"] < 0.34
        assert -0.15 < row["radial_energy_relative_change"] < -0.10
        assert 0.70 < row["swirl_energy_relative_change"] < 0.80
        assert 0.50 < row["axial_energy_relative_change"] < 0.60
        assert sum(row["canonical_fractions"].values()) == pytest.approx(1.0)
        assert sum(row["optimized_fractions"].values()) == pytest.approx(1.0)
        assert 0.41 < row["canonical_fractions"]["radial"] < 0.44
        assert 0.27 < row["optimized_fractions"]["radial"] < 0.30
        assert 0.24 < row["canonical_fractions"]["swirl"] < 0.26
        assert 0.32 < row["optimized_fractions"]["swirl"] < 0.34
        assert 0.32 < row["canonical_fractions"]["axial"] < 0.34
        assert 0.38 < row["optimized_fractions"]["axial"] < 0.40

    # The bulk integral is already stable on the last two midpoint grids; this
    # does not imply support validity or convergence of pointwise extrema.
    for row in report["resolution_stability"]:
        assert row["canonical_finest_relative_change"] < 1e-3
        assert row["optimized_finest_relative_change"] < 1e-3

    # At fixed dx=0.125 the L=2 box captures nearly all finite-window energy even
    # though the separate boundary audit still finds a nonzero early-time tail.
    for row in report["domain_sensitivity"]:
        assert 0.0 <= row["canonical_L2_vs_largest_relative_deficit"] < 2e-3
        assert 0.0 <= row["optimized_L2_vs_largest_relative_deficit"] < 2e-3

    assert report["velocity_changed_by_source_optimization"] is True
    assert report["velocity_changed_by_this_audit"] is False
    assert report["physical_support_validated"] is False
    assert report["cr001_energy_normalization_status"] == "pending_unknown"
    assert report["cr001_energy_gate_assessed"] is False
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False


def test_component_energy_matches_cartesian_energy_and_fails_closed():
    candidate = Eq45VelocityCandidate.seed()
    row = finite_window_component_energy(candidate, 0.5, half_width=1.0, resolution=8)
    assert row["total_energy"] > 0.0
    assert row["radial_energy"] + row["swirl_energy"] + row["axial_energy"] == pytest.approx(
        row["total_energy"], rel=1e-14, abs=1e-14
    )

    with pytest.raises(ValueError):
        finite_window_component_energy(candidate, 0.5, half_width=0.0, resolution=8)
    with pytest.raises(ValueError):
        finite_window_component_energy(candidate, 0.5, half_width=1.0, resolution=3)
    with pytest.raises(ValueError):
        audit_optimized_eq45_energy(resolutions=(16, 32))
    with pytest.raises(ValueError):
        audit_optimized_eq45_energy(domain_levels=((2.0, 32), (2.5, 32), (3.0, 48)))


def test_optimization_identity_drift_is_rejected(tmp_path):
    root = Path(__file__).resolve().parents[1]
    source = root / "artifacts/constrained/eq45_profile_force_optimization.json"
    payload = json.loads(source.read_text(encoding="utf-8"))
    payload["optimized_candidate"]["sha256"] = "0" * 64
    mutated = tmp_path / "optimization.json"
    mutated.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(RuntimeError, match="identity"):
        audit_optimized_eq45_energy(
            optimization_path=mutated,
            resolutions=(8, 10, 12),
            domain_levels=((1.0, 8), (1.25, 10), (1.5, 12)),
        )
