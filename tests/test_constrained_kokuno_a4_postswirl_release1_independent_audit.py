from __future__ import annotations

import inspect
import json
from decimal import Decimal

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_postswirl_release1_independent_audit import (
    CANONICAL_ABSOLUTE_FD_LADDER,
    FINAL_PROJECT_DIVERGENCE_GATE,
    FINAL_PROJECT_MOMENTUM_GATE,
    NORMALIZED_DIVERGENCE_GATE,
    RELATIVE_FD_LADDER,
    SEED,
    TASK,
    TIMES,
    UPSTREAM_HEAD,
    _canonical_representability,
    _fd4,
    _heldout_points,
    _release_entry_seam_replay,
    audit_loaded_release1,
    default_candidate,
    public_api_has_no_scientific_tuning_knobs,
)
from openai_ns_reconstruction.kokuno_pa16_current_cartesian_postswirl_release1 import (
    KokunoPA16CurrentCartesianPostSwirlRelease1,
)


def test_frozen_protocol_and_project_gates_are_unchanged():
    assert TASK == "K4-VAL-123"
    assert UPSTREAM_HEAD == "f5b5b8b532e0802e181bfead5e748cd4008906bf"
    assert SEED == 9173941
    assert TIMES == (0.31, 0.47, 0.63, 0.71)
    assert RELATIVE_FD_LADDER == (2e-6, 1e-6, 5e-7)
    assert CANONICAL_ABSOLUTE_FD_LADDER == (0.02, 0.01, 0.005)
    assert NORMALIZED_DIVERGENCE_GATE == 1e-5
    assert FINAL_PROJECT_MOMENTUM_GATE == 1e-3
    assert FINAL_PROJECT_DIVERGENCE_GATE == 1e-5


class _LinearDecimalField:
    """Manufactured public Decimal field for the A4 FD4 operator only."""
    def __init__(self, divergent: bool = False):
        self.divergent = divergent

    def velocity(self, x, y, z, t):
        # Solenoidal: (x+2y, -y+3z, -4x), trace 1-1+0=0.
        # Mutation adds +1e-3*x to u_x so divergence becomes 1e-3.
        ux = x + 2.0*y + (1e-3*x if self.divergent else 0.0)
        uy = -y + 3.0*z
        uz = -4.0*x
        return np.asarray([Decimal.from_float(float(ux)), Decimal.from_float(float(uy)),
                           Decimal.from_float(float(uz))], dtype=object)


def _manufactured_point():
    return {"x": 1.2, "y": -0.7, "z": 0.4, "t": 0.5, "r": 1.4,
            "q": 1.0, "role": "manufactured", "s_target": 0.5}


def test_fd4_manufactured_solenoidal_and_divergent_mutation():
    p = _manufactured_point()
    J = _fd4(_LinearDecimalField(False), p, 1e-6)
    div = J[0][0] + J[1][1] + J[2][2]
    assert abs(float(div)) < 2e-10
    J_bad = _fd4(_LinearDecimalField(True), p, 1e-6)
    div_bad = J_bad[0][0] + J_bad[1][1] + J_bad[2][2]
    assert float(div_bad) == pytest.approx(1e-3, rel=2e-7, abs=2e-10)


def test_canonical_absolute_stencil_collapse_is_detected():
    p = {"x": 1e200, "y": -1e200, "z": 0.0, "t": 0.5, "r": 1e200,
         "q": 1.0, "role": "huge", "s_target": 0.5}
    receipt = _canonical_representability([p])
    assert receipt["all_representable"] is False
    assert any(not row["all_points_all_axes_representable"] for row in receipt["ladder"])


def test_heldout_points_are_deterministic_offgrid_and_strict_interior():
    c = default_candidate()
    a = _heldout_points(c)
    b = _heldout_points(c)
    assert a == b
    assert len(a) == 44
    roles = [p["role"] for p in a]
    assert sum(r.startswith("interior_band1") for r in roles) == 8
    assert sum(r.startswith("interior_band2") for r in roles) == 8
    assert sum(r.startswith("interior_band3") for r in roles) == 8
    assert sum(r.startswith("release_entry_seam") for r in roles) == 8
    assert sum(r.startswith("late_release") for r in roles) == 12
    for p in a:
        assert p["r"] > 0
        assert 0.0 < p["theta"] < 2.0*np.pi
        coords = c.similarity_coordinates_logX(p["x"], p["y"], p["z"], p["t"])
        got_s = float(np.asarray(coords["log_X"])) - c.log_X_rel
        assert got_s == pytest.approx(p["s_target"], abs=5e-11)
        assert float(np.asarray(coords["eta"])) == pytest.approx(0.0, abs=5e-13)


def test_release_entry_seam_is_exact_and_truth_remains_false():
    c = default_candidate()
    seam = _release_entry_seam_replay(c)
    assert seam["all_exact"] is True
    truth = c.truth_boundary
    assert truth["post_relative_swirl_release1_public_velocity_materialized"] is True
    for key in ("outer_global_leading_velocity_materialized", "matched_global_pressure_materialized",
                "restricted_forcing_materialized", "heldout_ns_residual_assessed", "pde_validated"):
        assert truth[key] is False


def test_save_load_rejects_truth_promotion(tmp_path):
    c = default_candidate()
    path = tmp_path / "candidate.json"
    payload = c.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianPostSwirlRelease1.load_configuration(path)
    assert loaded.configuration() == payload
    bad = json.loads(path.read_text())
    bad["truth_boundary"]["pde_validated"] = True
    path.write_text(json.dumps(bad))
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostSwirlRelease1.load_configuration(path)


def test_audit_public_api_has_no_scientific_tuning_knobs():
    assert public_api_has_no_scientific_tuning_knobs()
    sig = inspect.signature(audit_loaded_release1)
    assert set(sig.parameters) == {"loaded", "pre_serialization_reference"}
    forbidden = {"step","h","threshold","residual","forcing","force","pressure","viscosity",
                 "precision","gain","optimizer","tolerance"}
    assert not (forbidden & set(sig.parameters))
