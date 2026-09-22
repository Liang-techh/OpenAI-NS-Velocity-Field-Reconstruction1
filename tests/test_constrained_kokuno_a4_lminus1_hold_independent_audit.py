from __future__ import annotations

from decimal import Decimal
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_lminus1_hold_independent_audit import (
    SEED,
    TIMES,
    _canonical_representability,
    _fd4,
    _heldout_points,
    _seam,
    _stage_axis,
    default_candidate,
    public_api_has_no_scientific_tuning_knobs,
)
from openai_ns_reconstruction.kokuno_pa16_current_cartesian_postswirl_lminus1_hold import (
    KokunoPA16CurrentCartesianPostSwirlLMinus1Hold,
)


class _SolenoidalManufactured:
    @staticmethod
    def velocity(x, y, z, t):
        vals = (Decimal.from_float(-float(y)), Decimal.from_float(float(x)), Decimal(0))
        return np.asarray(vals, dtype=object)


class _DivergentManufactured:
    @staticmethod
    def velocity(x, y, z, t):
        vals = (Decimal.from_float(1.0e-3 * float(x)), Decimal(0), Decimal(0))
        return np.asarray(vals, dtype=object)


def _manufactured_point():
    return {
        "x": 0.4,
        "y": 0.3,
        "z": 0.2,
        "t": 0.5,
        "r": 0.5,
        "q": 1.0,
        "rho_target": 0.5,
        "hold_s_target": 1.0,
        "role": "manufactured",
    }


def test_save_load_replays_exact_configuration_and_semantic_identity(tmp_path: Path):
    candidate = default_candidate()
    path = tmp_path / "candidate.json"
    before = candidate.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianPostSwirlLMinus1Hold.load_configuration(path)
    assert loaded.configuration() == before
    assert loaded.semantic_sha256 == candidate.semantic_sha256
    assert loaded.truth_boundary["source_l_minus1_hold_materialized"] is True
    assert loaded.truth_boundary["heldout_ns_residual_assessed"] is False
    assert loaded.truth_boundary["pde_validated"] is False


def test_heldout_protocol_is_deterministic_offgrid_and_strict_interior():
    candidate = default_candidate()
    a = _heldout_points(candidate)
    b = _heldout_points(candidate)
    assert a == b
    assert len(a) == 44
    assert {float(row["t"]) for row in a} == set(TIMES)
    assert all(0.0 < float(row["rho_target"]) < 1.0 for row in a)
    assert all(np.isfinite(float(row[k])) for row in a for k in ("x", "y", "z", "r"))
    assert all(float(row["r"]) > 0.0 for row in a)
    # The seeded azimuths must not collapse to the coordinate axes/cardinal grid.
    assert all(abs(np.sin(float(row["theta"]))) > 1.0e-6 for row in a)
    assert all(abs(np.cos(float(row["theta"]))) > 1.0e-6 for row in a)
    assert SEED == 9173951


def test_independent_fd4_calibrates_solenoidal_field_and_detects_mutation():
    point = _manufactured_point()
    J0 = _fd4(_SolenoidalManufactured(), point, 1.0e-6)
    div0 = J0[0][0] + J0[1][1] + J0[2][2]
    assert abs(div0) <= Decimal("1e-20")

    J1 = _fd4(_DivergentManufactured(), point, 1.0e-6)
    div1 = J1[0][0] + J1[1][1] + J1[2][2]
    assert abs(float(div1) - 1.0e-3) <= 2.0e-10


def test_canonical_absolute_fd_collapse_fails_closed_at_huge_coordinates():
    points = [{"x": 1.0e200, "y": -1.0e200, "z": 0.0}]
    receipt = _canonical_representability(points)
    assert receipt["all_representable"] is False
    assert any(not row["all_points_all_axes_representable"] for row in receipt["ladder"])


def test_hold_entry_seam_axis_and_stage_guards_are_public_field_checks():
    candidate = default_candidate()
    seam = _seam(candidate)
    firewall = _stage_axis(candidate)
    assert seam["all_exact"] is True
    assert firewall["exact_axis_transverse_velocity_zero"] is True
    assert firewall["before_hold_profile_fails_closed"] is True
    assert firewall["after_hold_profile_fails_closed"] is True


def test_public_audit_api_exposes_no_scientific_tuning_knobs():
    assert public_api_has_no_scientific_tuning_knobs() is True


def test_future_truth_states_remain_fail_closed():
    candidate = default_candidate()
    truth = candidate.truth_boundary
    for key in (
        "source_l_minus1_to_minus_h_transition_materialized",
        "source_terminal_multiplier_materialized",
        "source_exterior_heat_replacement_materialized",
        "outer_global_leading_velocity_materialized",
        "unified_global_cartesian_velocity_export_ready",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ):
        assert truth[key] is False
