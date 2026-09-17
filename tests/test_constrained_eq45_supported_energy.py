import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_energy import (
    audit_supported_eq45_energy,
    gauss_energy_and_components,
)


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"
CONSTRAINTS = ROOT / "configs" / "constraints.json"


class _ConstantField:
    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        out = np.empty_like(points)
        out[..., 0] = 1.0
        out[..., 1] = 2.0
        out[..., 2] = 2.0
        return out


def test_gauss_energy_integrator_matches_constant_field_exactly():
    report = gauss_energy_and_components(
        _ConstantField(), time=0.5, half_width=2.0, order=4
    )
    np.testing.assert_allclose(report["total"], 288.0, rtol=2e-14, atol=2e-14)
    assert report["component_closure_relative_error"] < 2e-15


def test_supported_eq45_energy_is_recomputed_after_serialized_support_transform(tmp_path):
    parent = Eq45VelocityCandidate.load_json(SEED)
    child_path = tmp_path / "supported_eq45.json"
    Eq45SupportedVelocityCandidate(parent=parent).save_json(child_path)
    child = Eq45SupportedVelocityCandidate.load_json(child_path)

    report = audit_supported_eq45_energy(
        child,
        constraints_path=CONSTRAINTS,
        times=(0.25, 0.5, 0.75),
    )

    assert report["candidate_sha256"] == child.sha256
    assert report["parent_sha256"] == parent.sha256
    assert report["quadrature_orders_per_axis"] == [24, 48, 96]
    assert report["reference_time"] == 0.25
    assert report["reference_energy"] == 1.0
    assert report["reference_energy_abs_tolerance"] == 0.001
    assert report["quadrature_relative_change_threshold"] == 0.001
    assert report["truth_boundary"]["physical_support_validated"] is False
    assert report["truth_boundary"]["visualization_ready"] is False
    assert report["truth_boundary"]["pde_validated"] is False

    assert len(report["rows"]) == 3
    for row in report["rows"]:
        totals = [entry["total"] for entry in row["quadrature"]]
        assert all(np.isfinite(value) and value > 0.0 for value in totals)
        assert np.isfinite(row["medium_to_finest_relative_change"])
        assert np.isfinite(row["support_transform_energy_change_fraction"])
        fractions = row["child_same_order_component_fractions"]
        np.testing.assert_allclose(
            fractions["radial"] + fractions["swirl"] + fractions["axial"],
            1.0,
            rtol=5e-12,
            atol=5e-12,
        )

    # One-shot calibration commit: force pytest to surface the deterministic
    # report in Actions. This assertion is removed immediately after harvest.
    raise AssertionError("CR007_SUPPORTED_ENERGY_CALIBRATION=" + json.dumps(report, sort_keys=True))
