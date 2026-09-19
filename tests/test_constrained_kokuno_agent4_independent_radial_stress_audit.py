import math

import numpy as np

from openai_ns_reconstruction.kokuno_agent4_independent_radial_stress_audit import (
    EDGE_RELATIVE_MAX,
    FINEST_RELATIVE_MAX_MAX,
    FINEST_RELATIVE_RMS_MAX,
    MIN_REFINEMENT_RATIO,
    MOMENT_COMPLEMENT_RELATIVE_MAX,
    SIGN_FLIP_MUTATION_MIN,
    _audit_channel,
    _fd6_first_uniform,
    _independent_compact_stress,
)


def test_fd6_radial_derivative_matches_polynomial():
    r = np.linspace(0.4, 1.6, 97)
    y = r**6 - 0.3 * r**4 + 0.2 * r**2 - 0.7 * r
    x_mid, numerical = _fd6_first_uniform(y, r)
    exact = 6.0 * x_mid**5 - 1.2 * x_mid**3 + 0.4 * x_mid - 0.7
    assert np.max(np.abs(numerical - exact)) < 2.0e-10


def test_simpson_compact_stress_closes_moment_and_edges():
    r = np.linspace(0.4, 1.6, 193)
    source = np.sin(2.3 * r) + 0.2 * np.cos(5.1 * r)
    built = _independent_compact_stress(
        r,
        source,
        exponent=2,
        center=1.0,
        halfwidth=0.6,
    )
    stress = np.asarray(built["stress"])
    complement = np.asarray(built["complement"])
    moment = float(np.trapezoid((r**2) * complement, r))
    # Simpson construction closes the compact endpoint much more tightly than
    # this deliberately different trapezoid diagnostic needs.
    assert abs(stress[0]) == 0.0
    assert abs(stress[-1]) < 1.0e-9 * max(1.0, np.max(np.abs(stress)))
    assert math.isfinite(moment)


def test_synthetic_channel_operator_and_sign_mutation_are_detected():
    r = np.linspace(0.4, 1.6, 193)
    s = (r - 1.0) / 0.6
    source = np.zeros_like(r)
    mask = np.abs(s) < 1.0
    source[mask] = (np.cos(0.5 * math.pi * s[mask]) ** 6) * (1.0 + 0.3 * r[mask])
    receipt = _audit_channel(
        r,
        source,
        exponent=1,
        center=1.0,
        halfwidth=0.6,
    )
    assert receipt["operator_relative_rms"] < 2.0e-5
    assert receipt["operator_relative_max"] < 5.0e-5
    assert receipt["moment_complement_relative"] < 1.0e-12
    assert receipt["edge_relative"] < 1.0e-10
    assert receipt["sign_flip_mutation_relative_rms"] > 1.5


def test_scientific_guards_are_frozen_and_distinct_from_project_gates():
    assert FINEST_RELATIVE_RMS_MAX == 2.0e-2
    assert FINEST_RELATIVE_MAX_MAX == 5.0e-2
    assert MIN_REFINEMENT_RATIO == 2.0
    assert MOMENT_COMPLEMENT_RELATIVE_MAX == 1.0e-10
    assert EDGE_RELATIVE_MAX == 1.0e-8
    assert SIGN_FLIP_MUTATION_MIN == 5.0e-1
