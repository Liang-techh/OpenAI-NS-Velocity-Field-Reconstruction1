import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_profile_separable_rank import (
    audit_selected_paper_core_separable_rank,
    diagnose_separable_profile_rank,
)


def _level(report, rank):
    return next(level for level in report.levels if level.rank == rank)


def test_exact_rank_two_value_and_derivative_capacity():
    X = np.linspace(0.0, 1.0, 19)
    eta = np.linspace(-1.0, 1.0, 23)
    f1 = 1.0 + X
    f2 = X * (1.0 - X)
    df1 = np.ones_like(X)
    df2 = 1.0 - 2.0 * X
    g1 = 1.0 + 0.3 * eta
    g2 = eta * eta - 0.2
    values = f1[:, None] * g1[None, :] + f2[:, None] * g2[None, :]
    derivatives = df1[:, None] * g1[None, :] + df2[:, None] * g2[None, :]

    report = diagnose_separable_profile_rank(
        values,
        derivatives,
        ranks=(1, 2, 3),
        rank_rtol=1e-12,
        diagnostic_tolerance=1e-10,
    )

    assert report.numerical_rank == 2
    assert _level(report, 1).worst_relative_error > 1e-3
    assert _level(report, 2).worst_relative_error < 2e-15
    assert report.smallest_rank_below_diagnostic_tolerance == 2
    assert report.pde_validated is False
    assert report.visual_correspondence_verified is False


def test_derivative_block_exposes_a_mode_hidden_by_value_only_fit():
    X = np.linspace(0.1, 1.0, 21)
    eta = np.linspace(-0.9, 0.9, 27)
    g1 = 1.0 + eta
    g2 = np.cos(np.pi * eta / 2.0)
    values = (1.0 + 0.2 * X)[:, None] * g1[None, :]
    derivatives = 0.2 * g1[None, :] + (0.8 * X)[:, None] * g2[None, :]

    report = diagnose_separable_profile_rank(
        values, derivatives, ranks=(1, 2), rank_rtol=1e-10
    )

    assert report.numerical_rank == 2
    assert _level(report, 1).radial_derivative_relative_error > 0.05
    assert _level(report, 2).worst_relative_error < 2e-15


def test_selected_paper_core_has_a_four_mode_knee_before_ill_conditioned_tail():
    report = audit_selected_paper_core_separable_rank(
        X_values=np.linspace(0.001, 0.409, 13),
        eta_values=np.linspace(-0.9, 0.9, 33),
        ranks=(1, 2, 3, 4, 6),
    )

    assert set(report.channels) == {"phi", "F", "U"}
    assert report.production_inner_seed_evaluated is True
    assert report.velocity_changed is False
    assert report.pde_validated is False
    assert report.visualization_ready is False
    for channel in report.channels.values():
        errors = [level.worst_relative_error for level in channel.levels]
        assert np.all(np.diff(errors) <= 1e-14)
        # On this deliberately coarser audit grid the singular value nearest
        # the 1e-3 cutoff can cross the threshold.  The full checked 25x65
        # calibration below carries the exact rank-4 statement.
        assert channel.numerical_rank in (4, 5)
        assert _level(channel, 4).condition_number < 500.0
        assert _level(channel, 6).condition_number > 1e4
        assert _level(channel, 6).worst_relative_error < 2e-5

    assert report.channels["F"].smallest_rank_below_diagnostic_tolerance == 4
    assert report.channels["phi"].smallest_rank_below_diagnostic_tolerance == 6
    assert report.channels["U"].smallest_rank_below_diagnostic_tolerance == 6


def test_checked_calibration_keeps_the_selected_seed_knee_visible():
    root = Path(__file__).resolve().parents[1]
    expected = json.loads(
        (root / "artifacts/constrained/profile_separable_rank_calibration.json").read_text()
    )
    report = audit_selected_paper_core_separable_rank()

    assert expected["inner_seed"] == {
        "sigma": 0.5,
        "maxdegree": 14,
        "eta_nodes": 257,
    }
    for name in ("phi", "F", "U"):
        current = report.channels[name]
        stored = expected["channels"][name]
        assert current.numerical_rank == stored["numerical_rank_at_rtol_1e-3"]
        assert current.smallest_rank_below_diagnostic_tolerance == stored[
            "smallest_rank_below_1e-3"
        ]
        for rank in (4, 6):
            current_level = _level(current, rank)
            stored_level = next(
                item for item in stored["levels"] if item["rank"] == rank
            )
            assert current_level.worst_relative_error == pytest.approx(
                stored_level["worst_relative_error"], rel=0.08, abs=2e-6
            )


def test_fail_closed_on_malformed_or_inactive_inputs():
    good = np.ones((4, 5))
    with pytest.raises(ValueError):
        diagnose_separable_profile_rank(good, np.ones((4, 4)))
    with pytest.raises(ValueError):
        diagnose_separable_profile_rank(np.full((4, 5), np.nan), good)
    with pytest.raises(ValueError):
        diagnose_separable_profile_rank(good, np.zeros((4, 5)))
    with pytest.raises(ValueError):
        diagnose_separable_profile_rank(good, good, ranks=(2, 2))
    with pytest.raises(ValueError):
        audit_selected_paper_core_separable_rank(
            X_values=(0.2, 0.1), eta_values=(-0.5, 0.5)
        )
