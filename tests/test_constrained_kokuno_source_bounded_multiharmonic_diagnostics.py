from __future__ import annotations

import inspect
import math

import numpy as np

from openai_ns_reconstruction.kokuno_source_axis_safe_physical_velocity import (
    AxisSafeSourceHarmonicProvider,
)
from openai_ns_reconstruction.kokuno_source_bounded_multiharmonic_diagnostics import (
    CHART_FD_STEPS,
    audit_bounded_multiharmonic_field,
    bounded_multiharmonic_diagnostics_contract,
    public_contract,
)
from openai_ns_reconstruction.kokuno_source_bounded_multiharmonic_velocity import (
    BoundedHarmonicTerm,
    KokunoBoundedMultiHarmonicPhysicalVelocity,
)
from openai_ns_reconstruction.kokuno_source_physicalized_localized_harmonic import (
    source_chart_to_physical_rzt,
)
from openai_ns_reconstruction.kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizedZeroDataHarmonic,
)


def _chart_cartesian(theta: np.ndarray, cylindrical: np.ndarray) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    cylindrical = np.asarray(cylindrical, dtype=float)
    c = np.cos(theta)
    s = np.sin(theta)
    ur = cylindrical[..., 0]
    uth = cylindrical[..., 1]
    uz = cylindrical[..., 2]
    return np.stack((ur * c - uth * s, ur * s + uth * c, uz), axis=-1)


def _compact_toroidal_provider(
    *,
    provider_id: str,
    digest_char: str,
    core: float,
    width: float,
) -> AxisSafeSourceHarmonicProvider:
    """Manufactured compact curl A=(0,0,a(R)), u_theta=-d_R a.

    a=s^5(L-s)^5 on 0<s<L and zero outside.  Both a and the first four
    derivatives vanish at the support edges, so the certified axis core is
    compatible with the manufactured field rather than being a radius clip.
    """

    def evaluate(R, theta, Z, T, scaling):
        del Z, T, scaling
        R = np.asarray(R, dtype=float)
        theta = np.asarray(theta, dtype=float)
        s = R - core
        active = (s > 0.0) & (s < width)
        a_z = np.zeros_like(R)
        u_theta = np.zeros_like(R)
        sa = s[active]
        a_z[active] = sa**5 * (width - sa) ** 5
        u_theta[active] = -5.0 * sa**4 * (width - sa) ** 4 * (width - 2.0 * sa)

        potential = np.zeros(R.shape + (3,), dtype=np.complex128)
        velocity = np.zeros(R.shape + (3,), dtype=np.complex128)
        potential[..., 2] = a_z
        velocity[..., 1] = u_theta
        real_cyl = 2.0 * np.real(velocity)
        real_cart = _chart_cartesian(theta, real_cyl)
        zeros = np.zeros(R.shape + (3,), dtype=np.complex128)
        return SourceLocalizedZeroDataHarmonic(
            amplitude_sensitivity=None,
            phase_jet=None,
            coefficient_jet=None,
            complete_amplitude=None,
            localization_weight=np.where(active, 1.0, 0.0),
            dr_localization_weight=np.zeros_like(R),
            dz_localization_weight=np.zeros_like(R),
            localized_coefficient=zeros.copy(),
            localized_dr_coefficient=zeros.copy(),
            localized_dz_coefficient=zeros.copy(),
            complex_vector_potential=potential,
            complex_velocity_cylindrical=velocity,
            real_pair_velocity_cylindrical=real_cyl,
            real_pair_velocity_chart_cartesian=real_cart,
        )

    return AxisSafeSourceHarmonicProvider(
        provider_id=provider_id,
        provider_semantic_sha256=digest_char * 64,
        axis_zero_radius_chart=core,
        axis_zero_certified=True,
        evaluate_chart=evaluate,
    )


def _field() -> KokunoBoundedMultiHarmonicPhysicalVelocity:
    p1 = _compact_toroidal_provider(
        provider_id="manufactured-m1", digest_char="a", core=0.12, width=2.0
    )
    p2 = _compact_toroidal_provider(
        provider_id="manufactured-m2", digest_char="b", core=0.18, width=2.1
    )
    return KokunoBoundedMultiHarmonicPhysicalVelocity(
        (
            BoundedHarmonicTerm(p1, amplitude=0.38, phase_offset=0.23),
            BoundedHarmonicTerm(p2, amplitude=0.47, phase_offset=-0.41),
        ),
        ell=8,
        h=0.004,
    )


def test_three_resolution_diagnostic_sees_solenoidal_manufactured_curls_and_vorticity():
    report = audit_bounded_multiharmonic_field(_field())
    assert report["schema"].endswith("diagnostics-v1")
    assert report["protocol"]["chart_fd_steps"] == list(CHART_FD_STEPS)
    assert report["protocol"]["probe_count"] == 8
    assert report["protocol"]["derivatives_consume_public_velocity_only"] is True
    assert report["protocol"]["scientific_acceptance_gate"] is None
    assert len(report["resolutions"]) == 3
    assert report["velocity_observation"]["nonzero_probe_count"] >= 6
    assert report["axis_support_observation"]["common_axis_core_exact_zero"] is True
    assert report["axis_support_observation"]["outer_support_assessed"] is False
    assert report["offgrid_observation"]["finite"] is True
    assert report["offgrid_observation"]["speed_rms"] > 0.0

    normalized = [row["normalized_divergence_rms"] for row in report["resolutions"]]
    vorticity = [row["vorticity_rms"] for row in report["resolutions"]]
    axial_fraction = [row["axial_vorticity_rms_fraction"] for row in report["resolutions"]]
    assert all(np.isfinite(normalized))
    assert all(value < 2.0e-2 for value in normalized)
    assert normalized[-1] <= max(normalized[0], 1.0e-12)
    assert all(value > 0.0 and np.isfinite(value) for value in vorticity)
    assert all(0.0 <= value <= 1.0 + 1.0e-12 for value in axial_fraction)
    assert report["resolution_stability_observation"]["medium_to_fine_vorticity_relative_rms_drift"] < 5.0e-2
    assert len(report["diagnostic_sha256"]) == 64


def test_manufactured_provider_has_real_outer_support_but_diagnostic_does_not_launder_it_into_contract():
    field = _field()
    # Both manufactured terms are zero beyond their private test-only support.
    R = np.asarray((2.6, 2.8, 3.0))
    theta = np.asarray((0.2, 1.0, 2.4))
    r, z, t = source_chart_to_physical_rzt(
        R,
        np.zeros(3),
        np.ones(3),
        ell=field.scaling.ell,
        h=field.scaling.h,
    )
    velocity = field.velocity(r * np.cos(theta), r * np.sin(theta), z, t)
    assert np.array_equal(velocity, np.zeros((3, 3)))

    report = audit_bounded_multiharmonic_field(field)
    assert report["axis_support_observation"]["outer_support_assessed"] is False
    assert "provider contract exposes axis-zero core" in report["axis_support_observation"]["reason_outer_support_not_assessed"]


def test_public_velocity_fd_detects_an_explicit_divergent_mutation():
    class DivergentField(KokunoBoundedMultiHarmonicPhysicalVelocity):
        def velocity(self, x, y, z, t):
            base = np.array(super().velocity(x, y, z, t), dtype=float, copy=True)
            x_b = np.broadcast_arrays(
                np.asarray(x, dtype=float),
                np.asarray(y, dtype=float),
                np.asarray(z, dtype=float),
                np.asarray(t, dtype=float),
            )[0]
            base[..., 0] += 100.0 * x_b
            return base

    base = _field()
    mutated = DivergentField(base.terms, ell=base.scaling.ell, h=base.scaling.h)
    clean = audit_bounded_multiharmonic_field(base)
    bad = audit_bounded_multiharmonic_field(mutated)
    clean_fine = clean["resolutions"][-1]["divergence_max_abs"]
    bad_fine = bad["resolutions"][-1]["divergence_max_abs"]
    assert bad_fine > clean_fine + 50.0
    assert bad["resolutions"][-1]["normalized_divergence_rms"] > 5.0e-2


def test_axis_common_core_is_exact_zero_at_public_velocity_surface():
    field = _field()
    common_core = min(term.provider.axis_zero_radius_chart for term in field.terms)
    R = np.asarray((0.0, 0.25 * common_core, 0.75 * common_core))
    r, z, t = source_chart_to_physical_rzt(
        R,
        np.asarray((0.0, 0.1, -0.1)),
        np.asarray((0.8, 0.9, 1.0)),
        ell=field.scaling.ell,
        h=field.scaling.h,
    )
    got = field.velocity(r, np.zeros_like(r), z, t)
    assert np.array_equal(got, np.zeros((3, 3)))


def test_truth_boundary_keeps_residual_source_and_outer_support_claims_false():
    contract = bounded_multiharmonic_diagnostics_contract()
    assert contract["public_velocity_only_cartesian_derivatives_materialized"] is True
    assert contract["three_resolution_divergence_observation_materialized"] is True
    assert contract["three_resolution_vorticity_morphology_observation_materialized"] is True
    assert contract["certified_axis_zero_support_checked"] is True
    assert contract["outer_support_boundary_assessed"] is False
    assert contract["analytic_divergence_identity_proved_by_this_increment"] is False
    assert contract["source_provider_self_contained"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["oscillation_before_after_ns_residual_compared"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_public_surface_has_no_scientific_tuning_inputs():
    assert list(inspect.signature(audit_bounded_multiharmonic_field).parameters) == ["field"]
    public = public_contract()
    assert public["audit_inputs"] == ["field"]
    assert public["forbidden_audit_inputs_present"] == []
    assert public["chart_fd_steps"] == list(CHART_FD_STEPS)


def test_diagnostic_is_deterministic_for_same_field_identity():
    field = _field()
    first = audit_bounded_multiharmonic_field(field)
    second = audit_bounded_multiharmonic_field(field)
    assert first["field_semantic_sha256"] == second["field_semantic_sha256"]
    assert first["diagnostic_sha256"] == second["diagnostic_sha256"]
    assert first["resolutions"] == second["resolutions"]
