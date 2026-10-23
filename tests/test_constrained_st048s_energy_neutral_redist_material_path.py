import copy

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_st048s_energy_neutral_redist_material_path import (
    AGENT7_REPORTED_ALPHA,
    REDISTRIBUTION_GAIN,
    TRUTH_BOUNDARY,
    RedistributionVelocity,
    _summarize_by_seed_radius,
    audit_truth_boundary,
    compact_interval_bump,
    redistribution_profile,
)


def _cyl(points, velocity):
    r = np.hypot(points[:, 0], points[:, 1])
    rx, ry = points[:, 0] / r, points[:, 1] / r
    return (
        rx*velocity[:, 0] + ry*velocity[:, 1],
        -ry*velocity[:, 0] + rx*velocity[:, 1],
        velocity[:, 2],
    )


def test_radial_profiles_put_capacity_inner_and_remove_outer():
    r = np.array([0.0, 0.3, 0.6, 0.9, 1.2, 1.85, 2.0])
    inner = compact_interval_bump(r, 0.30, 1.05)
    outer = compact_interval_bump(r, 0.95, 1.85)
    assert inner[2] > 0 and inner[3] > 0 and inner[4] == 0
    assert outer[2] == 0 and outer[3] == 0 and outer[4] > 0
    h = redistribution_profile(np.array([0.6, 0.9, 1.2]), AGENT7_REPORTED_ALPHA)
    assert h[0] > 0 and h[1] > 0 and h[2] < 0


def test_redistribution_changes_only_swirl_before_common_scale():
    points = np.array([[0.6, 0.0, 0.1], [0.0, 0.9, -0.2], [-1.2, 0.0, 0.3]])
    velocity = np.array([[0.2, 0.3, 0.4], [-0.5, 0.25, -0.1], [0.7, -0.2, 0.6]])
    field = RedistributionVelocity(lambda p, t: velocity, gain=REDISTRIBUTION_GAIN, alpha=AGENT7_REPORTED_ALPHA)
    before = _cyl(points, velocity)
    after = _cyl(points, field(points, 0.5))
    h = redistribution_profile(np.hypot(points[:,0], points[:,1]), AGENT7_REPORTED_ALPHA)
    assert np.allclose(after[0], before[0], rtol=0, atol=3e-16)
    assert np.allclose(after[1], before[1]*(1+REDISTRIBUTION_GAIN*h), rtol=3e-15, atol=3e-16)
    assert np.array_equal(after[2], before[2])
    with pytest.raises(ValueError):
        RedistributionVelocity(lambda p, t: velocity, gain=.051, alpha=AGENT7_REPORTED_ALPHA)


def test_radius_summary_uses_frozen_seed_metadata():
    per_path=[]
    pair_rows=[]
    for radius in (0.6,0.9,1.2):
        for i in range(16):
            per_path.append({"seed":{"radius":radius}, "absolute_turns":radius+i*1e-4, "radius_change":-0.01*radius})
        for i in range(8):
            pair_rows.append({"radius":radius, "axial_separation_change":0.02*radius})
    out=_summarize_by_seed_radius({"per_path":per_path,"pair_rows":pair_rows})
    assert set(out)=={"0.6","0.9","1.2"}
    assert all(row["path_count"]==16 for row in out.values())
    assert out["0.6"]["mean_pair_axial_separation_change"]==pytest.approx(.012)


def test_truth_boundary_rejects_promotion():
    report={"task_id":"CR-A9-050","candidate":{"pde_validated":False},
            "agent7_formula_crosscheck":{"crosscheck_passed":True},
            "truth_boundary":copy.deepcopy(TRUTH_BOUNDARY)}
    audit_truth_boundary(report)
    for key in ("visualization_ready","pde_validated","production_redistribution_gain_selected"):
        bad=copy.deepcopy(report); bad["truth_boundary"][key]=True
        with pytest.raises(ValueError):
            audit_truth_boundary(bad)
