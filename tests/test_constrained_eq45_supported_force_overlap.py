import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from openai_ns_reconstruction.constrained_axisymmetric_physical_taper import (
    AxisymmetricPhysicalTaper,
)
from openai_ns_reconstruction.constrained_force import RestrictedForce


RECEIPT = (
    Path(__file__).resolve().parents[1]
    / "artifacts"
    / "constrained"
    / "eq45_supported_force_overlap_audit.json"
)


def _cylindrical_l2_sq(force, order, radial_limit, axial_half_height):
    radial_nodes, radial_weights = leggauss(order)
    axial_nodes, axial_weights = leggauss(order)

    radius = 0.5 * radial_limit * (radial_nodes + 1.0)
    radial_weights = 0.5 * radial_limit * radial_weights
    z = axial_half_height * axial_nodes
    axial_weights = axial_half_height * axial_weights

    rr, zz = np.meshgrid(radius, z, indexing="ij")
    points = np.stack((rr, np.zeros_like(rr), zz), axis=-1)
    values = force(points, 0.5)
    density = np.sum(values * values, axis=-1)
    measure = (
        2.0
        * np.pi
        * rr
        * radial_weights[:, None]
        * axial_weights[None, :]
    )
    return float(np.sum(density * measure))


def _channel_metrics(force, order, taper):
    total = _cylindrical_l2_sq(
        force,
        order,
        taper.radial_support,
        taper.axial_half_height,
    )
    plateau_radius = taper.radial_support * np.sqrt(taper.radial_plateau_q)
    plateau_height = taper.axial_half_height * np.sqrt(taper.axial_plateau_q)
    plateau = _cylindrical_l2_sq(
        force,
        order,
        plateau_radius,
        plateau_height,
    )
    collar = total - plateau
    return {
        "total_l2_sq": total,
        "plateau_l2_sq": plateau,
        "collar_l2_sq": collar,
        "plateau_fraction": plateau / total,
        "collar_fraction": collar / total,
    }


def test_restricted_force_and_supported_velocity_share_outer_support():
    taper = AxisymmetricPhysicalTaper()
    assert taper.radial_support == 2.0
    assert taper.axial_half_height == 2.0
    assert taper.radial_plateau_q == 0.64
    assert taper.axial_plateau_q == 0.64

    boundary_and_exterior = np.array(
        [
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [0.0, 0.0, -2.0],
            [2.2, 0.0, 0.0],
            [0.0, 0.0, 2.2],
        ]
    )
    taper_value, _, _ = taper.factors(boundary_and_exterior)
    np.testing.assert_array_equal(taper_value, 0.0)

    for force in (RestrictedForce(a=1.0, c=0.0), RestrictedForce(a=0.0, c=1.0)):
        np.testing.assert_array_equal(force(boundary_and_exterior, 0.5), 0.0)

    interior = np.array([[1.0, 0.0, 0.5]])
    assert np.linalg.norm(RestrictedForce(a=1.0, c=0.0)(interior, 0.5)) > 0.0
    assert np.linalg.norm(RestrictedForce(a=0.0, c=1.0)(interior, 0.5)) > 0.0


def test_checked_overlap_receipt_replays_through_public_force():
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    assert receipt["schema"] == "eq45_supported_force_overlap_audit_v1"
    assert receipt["task_id"] == "CR004-EQ45-SUPPORT-FORCE-OVERLAP-018"

    taper = AxisymmetricPhysicalTaper()
    channels = {
        "a_poloidal": RestrictedForce(a=1.0, c=0.0),
        "c_swirl": RestrictedForce(a=0.0, c=1.0),
    }
    for order in receipt["quadrature"]["orders"]:
        for name, force in channels.items():
            measured = _channel_metrics(force, order, taper)
            expected = receipt["quadrature"]["basis_channels"][name][f"order_{order}"]
            for key, value in measured.items():
                np.testing.assert_allclose(value, expected[key], rtol=2e-12, atol=2e-12)

    a32 = receipt["quadrature"]["basis_channels"]["a_poloidal"]["order_32"]
    a48 = receipt["quadrature"]["basis_channels"]["a_poloidal"]["order_48"]
    c32 = receipt["quadrature"]["basis_channels"]["c_swirl"]["order_32"]
    c48 = receipt["quadrature"]["basis_channels"]["c_swirl"]["order_48"]

    assert abs(a48["collar_fraction"] - a32["collar_fraction"]) < 3e-6
    assert abs(c48["collar_fraction"] - c32["collar_fraction"]) < 3e-7
    assert 0.24 < a48["collar_fraction"] < 0.25
    assert 0.50 < c48["collar_fraction"] < 0.51


def test_overlap_audit_does_not_promote_force_or_pde_claims():
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    conclusion = receipt["conclusion"]
    truth = receipt["truth_boundary"]

    assert conclusion["outer_support_compatible"] is True
    assert conclusion["previous_force_coefficients_transferable_without_revalidation"] is False
    assert truth == {
        "velocity_changed": False,
        "force_family_changed": False,
        "force_refit_performed": False,
        "pde_validated": False,
        "visualization_ready_promoted": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
