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
AUDIT_RECORD = ROOT / "artifacts" / "constrained" / "eq45_supported_energy_audit.json"


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
    # |u|^2=9 and volume([-2,2]^3)=64, so E=0.5*9*64=288.
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
    expected = json.loads(AUDIT_RECORD.read_text(encoding="utf-8"))

    assert report["candidate_sha256"] == child.sha256 == expected["candidate_sha256"]
    assert report["parent_sha256"] == parent.sha256 == expected["parent_sha256"]
    assert report["quadrature_orders_per_axis"] == [24, 48, 96]
    assert report["reference_time"] == 0.25
    assert report["reference_energy"] == 1.0
    assert report["reference_energy_abs_tolerance"] == 0.001
    assert report["quadrature_relative_change_threshold"] == 0.001
    assert report["truth_boundary"]["physical_support_validated"] is False
    assert report["truth_boundary"]["visualization_ready"] is False
    assert report["truth_boundary"]["pde_validated"] is False

    # The support-connected child has a stable physical energy integral, but
    # the preregistered E(0.25)=1+-0.001 normalization is not met.
    assert report["all_quadrature_gates_pass"] is True
    assert report["all_energy_range_gates_pass"] is True
    assert report["reference_energy_gate_pass"] is False
    np.testing.assert_allclose(
        report["reference_observed_finest_energy"],
        expected["reference_observed_finest_energy"],
        rtol=2e-10,
        atol=2e-12,
    )
    np.testing.assert_allclose(
        report["reference_energy_abs_error"],
        expected["reference_energy_abs_error"],
        rtol=2e-10,
        atol=2e-12,
    )

    assert len(report["rows"]) == len(expected["rows"]) == 3
    for row, expected_row in zip(report["rows"], expected["rows"]):
        assert row["time"] == expected_row["time"]
        assert row["quadrature_gate_pass"] is True
        assert row["energy_range_gate_pass"] is True
        assert row["medium_to_finest_relative_change"] <= 0.001
        np.testing.assert_allclose(
            [entry["total"] for entry in row["quadrature"]],
            [
                expected_row["total_energy_by_order"][str(order)]
                for order in (24, 48, 96)
            ],
            rtol=2e-10,
            atol=2e-12,
        )
        np.testing.assert_allclose(
            row["support_transform_energy_change_fraction"],
            expected_row["support_transform_energy_change_fraction"],
            rtol=2e-10,
            atol=2e-12,
        )
        fractions = row["child_same_order_component_fractions"]
        np.testing.assert_allclose(
            fractions["radial"] + fractions["swirl"] + fractions["axial"],
            1.0,
            rtol=5e-12,
            atol=5e-12,
        )
        np.testing.assert_allclose(
            [fractions[name] for name in ("radial", "swirl", "axial")],
            [
                expected_row["child_same_order_component_fractions"][name]
                for name in ("radial", "swirl", "axial")
            ],
            rtol=2e-10,
            atol=2e-12,
        )
