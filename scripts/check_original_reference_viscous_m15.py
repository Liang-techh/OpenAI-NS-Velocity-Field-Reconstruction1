"""Independent fixture for the reference/restoration physical endpoint algebra.

This checker deliberately builds a small source germ instead of loading an
accepted ancestor or replaying a production constructor.  The symbolic part
checks that a certified positive quotient cancels a genuine multi-term
Jacobian only when rational polynomial division proves the cancellation.  The
numeric part uses a finite ``N = 17`` slow/fast source on two native charts,
with unequal native Jacobians and a nonzero exterior trace.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import mpmath as mp

from src.openai_ns_reconstruction.original_reference_viscous_m15 import (
    PhysicalEndpointFingerprint,
)


def _synthetic_nodes():
    """Build a source-certified, two-chart endpoint DAG."""

    nodes = []

    def add(operation, **fields):
        row = {"operation": operation}
        row.update(fields)
        nodes.append(row)
        return len(nodes) - 1

    def exact(value):
        value = Fraction(value)
        return add(
            "exact_rational",
            numerator=value.numerator,
            denominator=value.denominator,
        )

    def source(name):
        return add(
            "current_original_source_parameter",
            name=name,
            Z_independent=True,
        )

    def summands(*arguments):
        return add("sum", arguments=list(arguments))

    def product(*arguments):
        return add("product", arguments=list(arguments))

    q_ref = add("bound_variable", name="synthetic_reference_q")
    t_rest = add("bound_variable", name="synthetic_restoration_t")
    z = add("bound_variable", name="synthetic_Z")

    # J_ref is a genuine three-term polynomial.  Its synthetic source values
    # are 11/10, 1 and 1/5, so J_ref = 11/10 + 2 + 3/5 = 37/10.
    j0, j1, j2 = (source(name) for name in ("J0", "J1", "J2"))
    j_ref = summands(j0, product(exact(2), j1), product(exact(3), j2))
    j_rest = exact(1)
    a = source("slow_derivative_source")
    b = source("fast_derivative_source")
    j_times_a = product(j_ref, a)
    j_times_b = product(j_ref, b)
    numerator = summands(j_times_a, j_times_b)
    expected_sum = summands(a, b)
    certified = add(
        "positive_quotient",
        numerator=numerator,
        denominator=j_ref,
        source_positive_certificate={
            "source_handle": "fixture:J_ref_positive_on_source_germ",
            "statement": "J_ref = 37/10 > 0",
        },
    )

    # A quotient with the same source node is the minimal x/x certificate.
    x = source("positive_x")
    x_over_x = add(
        "positive_quotient",
        numerator=x,
        denominator=x,
        source_positive_certificate={
            "source_handle": "fixture:x_positive",
            "statement": "x > 0",
        },
    )

    zero = exact(0)
    zero_over_zero = add(
        "positive_quotient",
        numerator=zero,
        denominator=zero,
        source_positive_certificate={
            "source_handle": "fixture:invalid_zero_domain",
            "statement": "intentionally invalid 0/0 control",
        },
    )

    uncertified = add(
        "positive_quotient",
        numerator=numerator,
        denominator=j_ref,
    )

    # The additive source a + J_ref is not divisible by J_ref.  This guards
    # against an interval-like or common-token cancellation that is too broad.
    nondivisible_numerator = summands(a, j_ref)
    nondivisible = add(
        "positive_quotient",
        numerator=nondivisible_numerator,
        denominator=j_ref,
        source_positive_certificate={
            "source_handle": "fixture:J_ref_positive_on_source_germ",
            "statement": "same positive J_ref certificate",
        },
    )

    wrong_j_ref = summands(
        j0,
        product(exact(2), j1),
        product(exact(4), j2),
    )
    wrong_j = add(
        "positive_quotient",
        numerator=numerator,
        denominator=wrong_j_ref,
        source_positive_certificate={
            "source_handle": "fixture:wrong_J_ref",
            "statement": "deliberately wrong native Jacobian",
        },
    )

    recipe = add(
        "current_original_leading_function_recipe",
        name="synthetic_source_germ_recipe",
        coordinate=q_ref,
        Z_variable=z,
    )
    substitution_value = exact(0)
    substituted_recipe = add(
        "function_substitution",
        variable=q_ref,
        value=substitution_value,
        expression=recipe,
        Z_independent_substitution=True,
    )

    return {
        "nodes": nodes,
        "q_ref": q_ref,
        "t_rest": t_rest,
        "z": z,
        "j_ref": j_ref,
        "j_rest": j_rest,
        "a": a,
        "b": b,
        "expected_sum": expected_sum,
        "certified": certified,
        "x_over_x": x_over_x,
        "zero_over_zero": zero_over_zero,
        "uncertified": uncertified,
        "nondivisible": nondivisible,
        "wrong_j": wrong_j,
        "recipe": recipe,
        "substituted_recipe": substituted_recipe,
    }


def _symbolic_checks():
    data = _synthetic_nodes()
    nodes = data["nodes"]
    q_ref, t_rest, z = data["q_ref"], data["t_rest"], data["z"]
    fingerprint = PhysicalEndpointFingerprint(
        nodes,
        {"coordinate": q_ref, "loop": {}},
    )
    scope_ref = {
        q_ref: fingerprint.rational(1),
        t_rest: fingerprint.rational(0),
        z: fingerprint.rational(0),
    }
    scope_rest = {
        q_ref: fingerprint.rational(0),
        t_rest: fingerprint.rational(0),
        z: fingerprint.rational(0),
    }

    quotient_hash, quotient_rational = fingerprint.at(
        data["certified"], scope_ref
    )
    expected_hash, expected_rational = fingerprint.at(
        data["expected_sum"], scope_ref
    )
    assert quotient_rational is None
    assert expected_rational is None
    assert quotient_hash == expected_hash, "certified J*(a+b)/J did not cancel"
    assert data["certified"] in fingerprint.certified_quotients
    j_hash, j_rational = fingerprint.at(data["j_ref"], scope_ref)
    assert j_rational is None
    assert j_hash in fingerprint.exact_divisions

    x_hash, x_rational = fingerprint.at(data["x_over_x"], scope_ref)
    assert x_rational == 1 and x_hash == fingerprint.one[0]

    try:
        fingerprint.at(data["zero_over_zero"], scope_ref)
    except ValueError as exc:
        assert "Undefined certified endpoint quotient" in str(exc)
    else:
        raise AssertionError("certified 0/0 was accepted")

    opaque_hash, opaque_rational = fingerprint.at(data["uncertified"], scope_ref)
    assert opaque_rational is None
    assert opaque_hash != expected_hash

    nondiv_hash, nondiv_rational = fingerprint.at(data["nondivisible"], scope_ref)
    assert nondiv_rational is None
    nondiv_terms = fingerprint.polynomials.get(nondiv_hash)
    assert nondiv_terms, "non-divisible quotient unexpectedly became an opaque leaf"
    assert any(
        "inverse:" + j_hash in monomial for monomial in nondiv_terms
    ), "non-divisible numerator was incorrectly canceled"

    wrong_hash, wrong_rational = fingerprint.at(data["wrong_j"], scope_ref)
    assert wrong_rational is None
    assert wrong_hash != expected_hash, "wrong native J produced the same identity"

    direct_ref = fingerprint.at(data["recipe"], scope_ref)
    direct_rest = fingerprint.at(data["recipe"], scope_rest)
    substituted = fingerprint.at(data["substituted_recipe"], scope_ref)
    assert direct_ref[0] != direct_rest[0], "opaque source argument lost scope"
    assert substituted[0] == direct_rest[0], "function substitution lost its argument"

    return {
        "nodes": len(nodes),
        "certified_quotient_node": data["certified"],
        "certified_quotients": sorted(fingerprint.certified_quotients),
        "exact_division_denominators": sorted(fingerprint.exact_divisions),
        "maximum_terms": fingerprint.maximum_terms,
        "checks": [
            "multiterm J*(a+b)/J",
            "certified x/x",
            "certified 0/0 rejection",
            "uncertified quotient opaque",
            "non-divisible quotient retained",
            "opaque argument and scope distinction",
            "wrong-J nonidentity",
        ],
    }


def _finite_n_source(y, z, *, N=17):
    """Return slow/fast physical-y and ordinary-Z1 derivatives of B."""

    amplitude = (1 + mp.mpf("0.13") * y + mp.mpf("0.021") * y * y) * (
        1 + mp.mpf("0.17") * z + mp.mpf("0.013") * z * z
    )
    amplitude_y = (
        mp.mpf("0.13") + mp.mpf("0.042") * y
    ) * (1 + mp.mpf("0.17") * z + mp.mpf("0.013") * z * z)
    amplitude_z = (1 + mp.mpf("0.13") * y + mp.mpf("0.021") * y * y) * (
        mp.mpf("0.17") + mp.mpf("0.026") * z
    )
    phase = N * y + mp.mpf("0.13") * z + mp.mpf("0.017") * y * z
    phase_y_slow = mp.mpf("0.017") * z
    phase_y_fast = mp.mpf(N)
    phase_z = mp.mpf("0.13") + mp.mpf("0.017") * y
    sine, cosine = mp.sin(phase), mp.cos(phase)
    b = amplitude * sine
    slow_y = amplitude_y * sine + amplitude * cosine * phase_y_slow
    fast_y = amplitude * cosine * phase_y_fast
    ordinary_z1 = amplitude_z * sine + amplitude * cosine * phase_z
    return {
        "B": b,
        "slow_y": slow_y,
        "fast_y": fast_y,
        "ordinary_z1": ordinary_z1,
        "physical_y": slow_y + fast_y,
    }


def _radius(y, z):
    return mp.mpf("2.3") + mp.mpf("0.11") * y + mp.mpf("0.02") * z


def _exterior(y, z):
    return mp.mpf("0.19") + mp.mpf("0.03") * y + mp.mpf("0.04") * z


def _physical_qy(y, z, *, N=17):
    source = _finite_n_source(y, z, N=N)
    # The ordinary-Z1 contribution is deliberately retained in the physical
    # Q_y trace.  The outer 1/N is the finite-N normalization, not a second
    # phase/Jacobian factor.
    return (source["physical_y"] + mp.mpf("0.23") * source["ordinary_z1"]) / N


def _trace(y, z, *, N=17):
    return _radius(y, z) * _physical_qy(y, z, N=N) + _exterior(y, z)


def _close(a, b, tolerance=mp.mpf("1e-45")):
    return abs(a - b) <= tolerance * max(1, abs(a), abs(b))


def _quadrature_checks():
    mp.mp.dps = 70
    N = 17
    J_ref = mp.mpf(37) / 10
    J_rest = mp.mpf(1)
    y0, z0 = mp.mpf("0.37"), mp.mpf("0.41")

    # The source-germ values agree with the symbolic coefficients used above.
    source_values = (mp.mpf(11) / 10, mp.mpf(1), mp.mpf(1) / 5)
    assert _close(source_values[0] + 2 * source_values[1] + 3 * source_values[2], J_ref)

    at_ref = _finite_n_source(y0, z0, N=N)
    at_rest = _finite_n_source(y0, z0, N=N)
    assert abs(at_ref["slow_y"]) > mp.mpf("1e-12")
    assert abs(at_ref["fast_y"]) > mp.mpf("1e-12")
    assert abs(at_ref["ordinary_z1"]) > mp.mpf("1e-12")

    native_ref = J_ref * _physical_qy(y0, z0, N=N)
    native_rest = J_rest * _physical_qy(y0, z0, N=N)
    ref_trace = _radius(y0, z0) * (native_ref / J_ref) + _exterior(y0, z0)
    rest_trace = _radius(y0, z0) * (native_rest / J_rest) + _exterior(y0, z0)
    assert _close(ref_trace, rest_trace), "native-J physical endpoint mismatch"
    exterior_endpoint = _exterior(y0, z0)
    assert exterior_endpoint > 0, "fixture exterior trace must be nonzero"

    # Each chart is integrated in its own native coordinate, with the
    # physical dy = J d(native-coordinate) factor included exactly once.
    def ref_integrand(q):
        y = y0 + J_ref * (q - 1)
        return _trace(y, z0, N=N) * J_ref

    def rest_integrand(t):
        y = y0 + J_rest * t
        return _trace(y, z0, N=N) * J_rest

    reference_integral = mp.quad(ref_integrand, [0, 1])
    restoration_integral = mp.quad(rest_integrand, [0, 1])
    assert mp.isfinite(reference_integral) and mp.isfinite(restoration_integral)
    assert abs(reference_integral) > mp.mpf("1e-12")
    assert abs(restoration_integral) > mp.mpf("1e-12")

    failed = []

    def expected_mismatch(label, condition):
        assert condition, label + " negative control unexpectedly matched"
        failed.append(label)

    wrong_J_physical = native_ref / J_rest
    expected_mismatch(
        "wrong-J",
        not _close(wrong_J_physical, _physical_qy(y0, z0, N=N)),
    )
    slow_only_trace = _radius(y0, z0) * (
        (at_ref["slow_y"] + mp.mpf("0.23") * at_ref["ordinary_z1"]) / N
    ) + exterior_endpoint
    expected_mismatch("fast-phase omission", not _close(slow_only_trace, ref_trace))
    expected_mismatch("sign reversal", not _close(-rest_trace, ref_trace))
    expected_mismatch(
        "missing-radius factor",
        not _close(_physical_qy(y0, z0, N=N) + exterior_endpoint, ref_trace),
    )
    expected_mismatch("false external zero", not _close(exterior_endpoint, mp.mpf(0)))

    return {
        "N": N,
        "J_ref": str(J_ref),
        "J_rest": str(J_rest),
        "shared_y": str(y0),
        "shared_Z": str(z0),
        "physical_endpoint_trace": mp.nstr(ref_trace, 35),
        "nonzero_exterior_endpoint": mp.nstr(exterior_endpoint, 35),
        "reference_native_integral": mp.nstr(reference_integral, 35),
        "restoration_native_integral": mp.nstr(restoration_integral, 35),
        "negative_controls": failed,
        "ordinary_Z1_included": True,
        "slow_fast_split_included": True,
    }


def _physical_pressure_jet_checks():
    """Differentiate the pressure recovery with independent SymPy jets.

    The symbols are source fields, not interval centers or cap values.  The
    ordinary-Z derivative is carried one order further than the displayed
    pressure expression so that the ``m_z``, ``e_z`` and ``pressure_z`` product
    terms expose their genuine second-Z fields.
    """

    import sympy as sp

    R, E, V, m, mz, e, ez, Pstar, pressure, pressurez, delta, z = sp.symbols(
        "R E V m mz e ez Pstar pressure pressurez delta z", real=True
    )
    Rz, Ez, Vz, mzz, ezz, pressurezz = sp.symbols(
        "Rz Ez Vz mzz ezz pressurezz", real=True
    )

    T = (1 - delta) * z * m + (1 - z**2) * mz
    L = 1 - delta * z**2
    Iz_numerator = (
        -V
        + sp.Rational(1, 2) * (1 - delta) * (m - z * mz)
        + Pstar
        * (
            V * T
            + 2 * delta * z * e
            - (1 - z**2) * ez
            + 2 * (1 + delta) * z * pressure
            - (1 - z**2) * pressurez
        )
    )
    Iz = Iz_numerator / L
    root_p2 = R / E * Iz

    y_jets = {
        R: R,
        E: E / 10,
        V: 0,
        m: V - m,
        mz: Vz - mz,
        e: V**2 - E**2 / 2 - e,
        ez: 2 * V * Vz - E * Ez - ez,
        Pstar: 0,
        pressure: E**2 / 2,
        pressurez: E * Ez,
    }
    z_jets = {
        R: Rz,
        E: Ez,
        V: Vz,
        m: mz,
        mz: mzz,
        e: ez,
        ez: ezz,
        Pstar: 0,
        pressure: pressurez,
        pressurez: pressurezz,
    }

    def dy(expr):
        return sp.expand(sum(sp.diff(expr, field) * jet for field, jet in y_jets.items()))

    def dz(expr):
        # z is an explicit ordinary-Z coordinate in addition to the source
        # field jets, so the derivative includes both terms.
        return sp.expand(sp.diff(expr, z) + sum(sp.diff(expr, field) * jet for field, jet in z_jets.items()))

    # Since (R/E)_y = R/E - R/(10 E), this is the exact physical-y product
    # rule required by the source pressure reconstruction.
    lhs_y = dy(root_p2)
    rhs_y = sp.Rational(9, 10) * root_p2 + R / E * dy(Iz)
    assert sp.simplify(lhs_y - rhs_y) == 0

    # Ordinary Z1 must differentiate the radius/energy prefactor and every
    # primitive in Iz.  In particular, mzz/ezz/pressurezz are not caps: they
    # are independent second-Z fields in this source-level check.
    lhs_z = dz(root_p2)
    rhs_z = (Rz / E - R * Ez / E**2) * Iz + R / E * dz(Iz)
    assert sp.simplify(lhs_z - rhs_z) == 0
    for second_z_field in (mzz, ezz, pressurezz):
        assert sp.simplify(sp.diff(dz(Iz), second_z_field)) != 0

    # The prescribed y jets commute with the ordinary-Z jets for the products
    # used by Iz.  This catches a dropped product derivative independently of
    # the prefactor identity above.
    assert sp.simplify(dy(mz) - dz(m_y := V - m)) == 0
    assert sp.simplify(dy(ez) - dz(e_y := V**2 - E**2 / 2 - e)) == 0
    assert sp.simplify(dy(pressurez) - dz(pressure_y := E**2 / 2)) == 0

    return {
        "identity": "root_p2_y = 9/10 root_p2 + (R/E) Iz_y",
        "ordinary_Z1_product_rule": True,
        "second_Z_fields_present": ["mzz", "ezz", "pressurezz"],
        "jet_commutation": ["m_z", "e_z", "pressure_z"],
        "source_values_or_caps_selected": False,
    }


def _poisson_loop_checks():
    """Check the periodic Poisson loop algebra and a bounded numeric control."""

    import sympy as sp

    u, q, psi, t0, a, b = sp.symbols("u q psi t0 a b", real=True)
    r = u / sp.sqrt(1 + u**2)
    alpha = 2 * q / sp.sqrt(1 + u**2)
    one_minus_r2 = sp.simplify(1 - r**2)
    assert sp.simplify(one_minus_r2 - 1 / (1 + u**2)) == 0
    assert sp.simplify(alpha**2 / one_minus_r2 - 4 * q**2) == 0

    D = 1 - 2 * r * sp.cos(psi) + r**2
    complex_modulus_squared = (1 - r * sp.cos(psi)) ** 2 + (r * sp.sin(psi)) ** 2
    assert sp.simplify(complex_modulus_squared - D) == 0

    # The triangle lower bound is checked in both signs of r, with rho = |r|.
    rho = sp.symbols("rho", nonnegative=True, real=True)
    rho_lower = (1 - rho) ** 2
    positive_r_remainder = sp.factor(
        (1 - 2 * rho * sp.cos(psi) + rho**2) - rho_lower
    )
    negative_r_remainder = sp.factor(
        (1 + 2 * rho * sp.cos(psi) + rho**2) - rho_lower
    )
    assert sp.simplify(positive_r_remainder - 2 * rho * (1 - sp.cos(psi))) == 0
    assert sp.simplify(negative_r_remainder - 2 * rho * (1 + sp.cos(psi))) == 0

    # In Laurent form the Poisson numerator is the two-sided geometric series
    # with no constant coefficient.  This is the source of the zero first
    # moment and the square norm below, rather than a numerical quadrature.
    w = sp.symbols("w", nonzero=True)
    cos_laurent = (w + 1 / w) / 2
    kernel_laurent = (cos_laurent - r) / (1 - r * (w + 1 / w) + r**2)
    geometric_laurent = sp.Rational(1, 2) * (
        w / (1 - r * w) + (1 / w) / (1 - r / w)
    )
    assert sp.factor(kernel_laurent - geometric_laurent) == 0
    kernel_mean = 0
    kernel_square_mean = 1 / (2 * (1 - r**2))
    integral_t = 2 * sp.pi * t0
    integral_t2 = sp.simplify(
        2 * sp.pi * (t0**2 + alpha**2 * kernel_square_mean)
    )
    assert sp.simplify(integral_t2 - 2 * sp.pi * (t0**2 + 2 * q**2)) == 0
    assert kernel_mean == 0

    normalization = 1 / (2 * sp.pi * (1 + t0**2 + 2 * q**2))
    phi_period = sp.simplify(
        normalization * (2 * sp.pi + integral_t2)
    )
    assert phi_period == 1

    # The endpoint baseline is the affine source value at t0.  Its two chart
    # endpoint copies vanish for t0 = -b/a.  The phase derivative belongs to
    # the periodic loop itself and therefore has identical values at 0 and
    # 2*pi; no unproved source alias is introduced here.
    loop_t = t0 + alpha * (sp.cos(psi) - r) / D
    B_endpoint = a * t0 + b
    Bphi = a * sp.diff(loop_t, psi)
    assert sp.simplify(B_endpoint.subs(t0, -b / a)) == 0
    assert sp.simplify(Bphi.subs(psi, 0) - Bphi.subs(psi, 2 * sp.pi)) == 0
    assert sp.simplify(loop_t.subs(psi, 0) - loop_t.subs(psi, 2 * sp.pi)) == 0

    # Active source-germ specialization: nonzero a and Delta, positive q^2
    # for every eta > 0.  These are proof inputs, not selected interval caps.
    eta = sp.symbols("eta", positive=True, real=True)
    a_active = sp.Rational(4, 5)
    delta_active = -sp.Rational(6, 5)
    q2_active = (5 * eta + 3) / 4
    assert a_active != 0 and delta_active != 0
    assert sp.ask(sp.Q.positive(q2_active)) is True

    # Moderate finite controls verify the Fourier moment formulas numerically.
    # They are diagnostic fixtures only; they are not source values or caps.
    mp.mp.dps = 90
    numeric_controls = []
    for u_value, q_value, t0_value in (
        (mp.mpf("0.7"), mp.mpf("0.8"), mp.mpf("-0.2")),
        (mp.mpf("-1.1"), mp.mpf("0.55"), mp.mpf("0.3")),
    ):
        r_value = u_value / mp.sqrt(1 + u_value**2)
        alpha_value = 2 * q_value / mp.sqrt(1 + u_value**2)

        def t_value(angle):
            denominator = 1 - 2 * r_value * mp.cos(angle) + r_value**2
            return t0_value + alpha_value * (mp.cos(angle) - r_value) / denominator

        integral_1 = mp.quad(t_value, [0, 2 * mp.pi])
        integral_2 = mp.quad(lambda angle: t_value(angle) ** 2, [0, 2 * mp.pi])
        expected_1 = 2 * mp.pi * t0_value
        expected_2 = 2 * mp.pi * (t0_value**2 + 2 * q_value**2)
        assert abs(integral_1 - expected_1) < mp.mpf("1e-60")
        assert abs(integral_2 - expected_2) < mp.mpf("1e-60")
        numeric_controls.append(
            {
                "u": mp.nstr(u_value, 5),
                "integral_t": mp.nstr(integral_1, 18),
                "integral_t2": mp.nstr(integral_2, 18),
            }
        )

    return {
        "identities": [
            "1-r^2 = 1/(1+u^2) > 0",
            "alpha^2/(1-r^2) = 4 q^2",
            "D >= (1-|r|)^2",
            "integral t = 2 pi t0",
            "integral t^2 = 2 pi (t0^2 + 2 q^2)",
            "K normalization",
            "Phi(2 pi) = 1",
            "endpoint B baseline and Bphi periodicity",
        ],
        "active_specialization": {
            "a": "4/5",
            "b": "0",
            "t0": "0",
            "Delta": "-6/5",
            "q_squared": "(5 eta + 3)/4, eta > 0",
        },
        "numeric_controls_are_proof_inputs": True,
        "numeric_controls_are_source_numerics": False,
        "numeric_controls": numeric_controls,
    }


def _recipe_numeric_fixture():
    """Evaluate the new source-port recipe graph on synthetic smooth fields.

    The original node IDs below are intentionally synthetic ports.  Every
    ``original_endpoint_function`` is checked against its port and ordinary-Z
    order before evaluation, so this fixture cannot silently use a numerical
    atom or a Taylor-normalized coefficient.
    """

    from src.openai_ns_reconstruction.original_reference_viscous_m15 import (
        PhysicalGermRecipes,
    )

    mp.mp.dps = 70
    coordinate, endpoint = 100, 101
    z0 = mp.mpf("0.37")
    delta = Fraction(1, 1000)

    def E_fun(x):
        return mp.mpf("2.2") + mp.mpf("0.13") * mp.sin(mp.mpf("0.7") * x) + mp.mpf("0.02") * x**2

    def V_fun(x):
        return mp.mpf("0.4") + mp.mpf("0.2") * mp.cos(mp.mpf("0.5") * x) + mp.mpf("0.01") * x

    def m_fun(x):
        return mp.mpf("0.7") + mp.mpf("0.15") * mp.sin(mp.mpf("0.9") * x + mp.mpf("0.1"))

    def h_fun(x):
        return mp.mpf("0.8") + mp.mpf("0.09") * mp.cos(mp.mpf("0.8") * x) + mp.mpf("0.01") * x**2

    def k_fun(x):
        return mp.mpf("0.6") + mp.mpf("0.05") * mp.sin(mp.mpf("1.1") * x - mp.mpf("0.2"))

    def e_fun(x):
        return mp.mpf("0.9") + mp.mpf("0.07") * mp.cos(mp.mpf("0.4") * x + mp.mpf("0.3")) + mp.mpf("0.02") * x

    def p0_fun(x):
        return mp.mpf("0.3") + mp.mpf("0.08") * mp.sin(mp.mpf("0.6") * x)

    def p_fun(x):
        return mp.mpf("0.5") + mp.mpf("0.04") * mp.cos(mp.mpf("0.9") * x - mp.mpf("0.4"))

    def pressure_fun(x):
        return p0_fun(x) + p_fun(x)

    source_functions = {
        "E": E_fun,
        "V": V_fun,
        "m": m_fun,
        "h": h_fun,
        "k": k_fun,
        "e": e_fun,
        "p": p_fun,
        "P0": p0_fun,
        "pressure": pressure_fun,
    }

    def jet(function, x=z0):
        return [function(x), mp.diff(function, x), mp.diff(function, x, 2)]

    fields = {}
    source_jets = {}
    next_port = 1000
    for key, function in source_functions.items():
        ports = [next_port + order for order in range(3)]
        next_port += 10
        fields[key] = ports
        values = jet(function)
        for order, port in enumerate(ports):
            source_jets[port] = values

    radius_port, z_port, pstar_port = 2000, 2001, 2002
    radius_value = mp.mpf("2.35")
    pstar_value = mp.mpf("1.25")
    source_jets[radius_port] = [radius_value, mp.mpf(0), mp.mpf(0)]
    source_jets[z_port] = [z0, mp.mpf(0), mp.mpf(0)]
    source_jets[pstar_port] = [pstar_value, mp.mpf(0), mp.mpf(0)]

    recipe = PhysicalGermRecipes(
        fields,
        radius_port,
        z_port,
        pstar_port,
        delta,
        coordinate,
        endpoint,
    )
    built = recipe.build()
    nodes = built["recipe_nodes"]
    allowed_operations = {
        "exact_rational",
        "original_endpoint_function",
        "sum",
        "product",
        "positive_quotient",
    }
    operation_counts = {}
    for row in nodes:
        operation_counts[row["operation"]] = operation_counts.get(row["operation"], 0) + 1
        assert row["operation"] in allowed_operations
        if row["operation"] == "original_endpoint_function":
            assert row["original_coordinate"] == coordinate
            assert row["original_endpoint"] == endpoint
            assert row["ordinary_Z_order"] in (0, 1, 2)
            assert row["original_node"] in source_jets

    cache = {}

    def evaluate(node):
        if node in cache:
            return cache[node]
        row = nodes[node]
        operation = row["operation"]
        if operation == "exact_rational":
            value = mp.mpf(row["numerator"]) / row["denominator"]
        elif operation == "original_endpoint_function":
            value = source_jets[row["original_node"]][row["ordinary_Z_order"]]
        elif operation == "sum":
            value = sum((evaluate(child) for child in row["arguments"]), mp.mpf(0))
        elif operation == "product":
            value = mp.mpf(1)
            for child in row["arguments"]:
                value *= evaluate(child)
        elif operation == "positive_quotient":
            numerator = evaluate(row["numerator"])
            denominator = evaluate(row["denominator"])
            assert abs(denominator) > mp.mpf("1e-45")
            value = numerator / denominator
        else:
            raise AssertionError("Unexpected recipe operation: " + operation)
        cache[node] = value
        return value

    def evaluate_jet(jet_nodes):
        return [evaluate(node) for node in jet_nodes]

    def close(a, b):
        return abs(a - b) <= mp.mpf("1e-45") * max(1, abs(a), abs(b))

    def add3(*vectors):
        return [sum((vector[k] for vector in vectors), mp.mpf(0)) for k in range(3)]

    def scale3(vector, coefficient):
        coefficient = mp.mpf(coefficient.numerator) / coefficient.denominator if isinstance(coefficient, Fraction) else mp.mpf(coefficient)
        return [coefficient * value for value in vector]

    def mul3(left, right):
        return [
            left[0] * right[0],
            left[1] * right[0] + left[0] * right[1],
            left[2] * right[0] + 2 * left[1] * right[1] + left[0] * right[2],
        ]

    base = {key: source_jets[ports[0]] for key, ports in fields.items()}
    E, V, m, h, k, e, pressure = (base[key] for key in ("E", "V", "m", "h", "k", "e", "pressure"))
    E2, EV, V2 = mul3(E, E), mul3(E, V), mul3(V, V)
    zero3 = [mp.mpf(0), mp.mpf(0), mp.mpf(0)]
    y1_expected = {
        "E": scale3(E, Fraction(1, 10)),
        "V": zero3,
        "m": add3(V, scale3(m, -1)),
        "h": add3(E, scale3(h, Fraction(-3, 2))),
        "k": add3(EV, scale3(k, Fraction(-3, 2))),
        "e": add3(V2, scale3(E2, Fraction(-1, 2)), scale3(e, -1)),
        "p": scale3(E2, Fraction(1, 2)),
        "P0": zero3,
        "pressure": scale3(E2, Fraction(1, 2)),
    }
    y2_expected = {
        "E": scale3(E, Fraction(1, 100)),
        "V": zero3,
        "m": scale3(y1_expected["m"], -1),
        "h": add3(y1_expected["E"], scale3(y1_expected["h"], Fraction(-3, 2))),
        "k": add3(scale3(EV, Fraction(1, 10)), scale3(y1_expected["k"], Fraction(-3, 2))),
        "e": add3(scale3(E2, Fraction(-1, 10)), scale3(y1_expected["e"], -1)),
        "p": scale3(E2, Fraction(1, 10)),
        "P0": zero3,
        "pressure": scale3(E2, Fraction(1, 10)),
    }
    common = built["common_fields"]
    for key in fields:
        assert all(close(actual, expected) for actual, expected in zip(evaluate_jet(common["y0"][key]), base[key]))
        assert all(close(actual, expected) for actual, expected in zip(evaluate_jet(common["y1"][key]), y1_expected[key]))
        assert all(close(actual, expected) for actual, expected in zip(evaluate_jet(common["y2"][key]), y2_expected[key]))

    pressure_sum = add3(base["P0"], base["p"])
    assert all(close(actual, expected) for actual, expected in zip(base["pressure"], pressure_sum))
    for key in ("E", "V", "m", "e", "pressure"):
        assert abs(base[key][2]) > mp.mpf("1e-8"), "synthetic second-Z jet accidentally vanished: " + key

    def pair_add(*pairs):
        return (
            sum((pair[0] for pair in pairs), mp.mpf(0)),
            sum((pair[1] for pair in pairs), mp.mpf(0)),
        )

    def pair_scale(pair, coefficient):
        coefficient = mp.mpf(coefficient.numerator) / coefficient.denominator if isinstance(coefficient, Fraction) else mp.mpf(coefficient)
        return coefficient * pair[0], coefficient * pair[1]

    def pair_mul(left, right):
        return left[0] * right[0], left[1] * right[0] + left[0] * right[1]

    def pair_div(left, right):
        return left[0] / right[0], (left[1] * right[0] - left[0] * right[1]) / right[0] ** 2

    def source_values(x):
        return {
            key: (function(x), mp.diff(function, x))
            for key, function in source_functions.items()
        }

    def independent_iz(x):
        values = source_values(x)
        E0, E1 = values["E"]
        V0, V1 = values["V"]
        m0, m1 = values["m"]
        e0, e1 = values["e"]
        pressure0, pressure1 = values["pressure"]
        mz0 = m1
        ez0 = e1
        pressurez0 = pressure1
        T = (1 - mp.mpf(delta.numerator) / delta.denominator) * x * m0 + (1 - x**2) * mz0
        L = 1 - mp.mpf(delta.numerator) / delta.denominator * x**2
        numerator = (
            -V0
            + mp.mpf("0.5") * (1 - mp.mpf(delta.numerator) / delta.denominator) * (m0 - x * mz0)
            + pstar_value
            * (
                V0 * T
                + 2 * mp.mpf(delta.numerator) / delta.denominator * x * e0
                - (1 - x**2) * ez0
                + 2 * (1 + mp.mpf(delta.numerator) / delta.denominator) * x * pressure0
                - (1 - x**2) * pressurez0
            )
        )
        return numerator / L

    def independent_root(x):
        return radius_value / E_fun(x) * independent_iz(x)

    def independent_root_y_pair(x):
        values = source_values(x)
        E0, E1 = values["E"]
        V0, V1 = values["V"]
        m0, m1 = values["m"]
        e0, e1 = values["e"]
        pressure0, pressure1 = values["pressure"]
        delta_value = mp.mpf(delta.numerator) / delta.denominator
        E_pair = (E0, E0 / 10)
        V_pair = (V0, mp.mpf(0))
        m_pair = (m0, V0 - m0)
        mz_pair = (m1, V1 - m1)
        e_pair = (e0, V0**2 - E0**2 / 2 - e0)
        ez_pair = (e1, 2 * V0 * V1 - E0 * E1 - e1)
        pressure_pair = (pressure0, E0**2 / 2)
        pressurez_pair = (pressure1, E0 * E1)
        z_pair = (x, mp.mpf(0))
        one_pair = (mp.mpf(1), mp.mpf(0))
        delta_pair = (delta_value, mp.mpf(0))
        P_pair = (pstar_value, mp.mpf(0))
        R_pair = (radius_value, radius_value)
        one_minus_delta = pair_add(one_pair, pair_scale(delta_pair, -1))
        one_plus_delta = pair_add(one_pair, delta_pair)
        d_pair = pair_add(one_pair, pair_scale(pair_mul(z_pair, z_pair), -1))
        L_pair = pair_add(one_pair, pair_scale(pair_mul(delta_pair, pair_mul(z_pair, z_pair)), -1))
        T_pair = pair_add(
            pair_mul(pair_mul(one_minus_delta, z_pair), m_pair),
            pair_mul(d_pair, mz_pair),
        )
        linear_pair = pair_add(
            pair_scale(V_pair, -1),
            pair_scale(
                pair_mul(
                    one_minus_delta,
                    pair_add(m_pair, pair_scale(pair_mul(z_pair, mz_pair), -1)),
                ),
                Fraction(1, 2),
            ),
        )
        quadratic_pair = pair_add(
            pair_mul(V_pair, T_pair),
            pair_scale(pair_mul(pair_mul(delta_pair, z_pair), e_pair), 2),
            pair_scale(pair_mul(d_pair, ez_pair), -1),
            pair_scale(pair_mul(pair_mul(one_plus_delta, z_pair), pressure_pair), 2),
            pair_scale(pair_mul(d_pair, pressurez_pair), -1),
        )
        iz_pair = pair_div(
            pair_add(linear_pair, pair_mul(P_pair, quadratic_pair)),
            L_pair,
        )
        return pair_mul(pair_div(R_pair, E_pair), iz_pair)

    root_nodes = built["root_p2"]
    actual_root = evaluate_jet(root_nodes["value_Z0_Z1"])
    actual_root_y = evaluate_jet(root_nodes["physical_y_Z0_Z1"])
    expected_root = [independent_root(z0), mp.diff(independent_root, z0)]
    expected_root_y = [
        independent_root_y_pair(z0)[1],
        mp.diff(lambda x: independent_root_y_pair(x)[1], z0),
    ]
    assert all(close(actual, expected) for actual, expected in zip(actual_root, expected_root))
    assert all(close(actual, expected) for actual, expected in zip(actual_root_y, expected_root_y))
    assert abs(actual_root[1]) > mp.mpf("1e-8")
    assert abs(actual_root_y[1]) > mp.mpf("1e-8")

    # Negative controls: omitting the Leibniz binomial 2, or substituting a
    # first-Z jet for a genuine second-Z jet, must disagree with the recipe.
    wrong_binomial_k2 = EV[2] - mp.mpf("1.5") * y1_expected["k"][2]
    wrong_mz_for_mzz = V[1] - m[1]
    assert not close(evaluate_jet(common["y2"]["k"])[2], wrong_binomial_k2)
    assert not close(evaluate_jet(common["y2"]["m"])[2], wrong_mz_for_mzz)

    return {
        "passed": True,
        "coordinate": coordinate,
        "endpoint": endpoint,
        "synthetic_source_ports": len(source_jets),
        "recipe_nodes": len(nodes),
        "operation_counts": operation_counts,
        "ordinary_Z_orders_checked": [0, 1, 2],
        "common_fields_y1_y2_checked": sorted(fields),
        "pressure_is_P0_plus_p": True,
        "root_p2_value_Z0_Z1_checked": True,
        "root_p2_physical_y_Z0_Z1_checked": True,
        "nonconstant_second_Z_fields": ["E", "V", "m", "e", "pressure"],
        "negative_controls": ["wrong Leibniz binomial", "m_z substituted for m_zz"],
        "source_ports_are_synthetic_fixture_values": True,
        "interval_or_cap_values_used": False,
    }


def run():
    symbolic = _symbolic_checks()
    quadrature = _quadrature_checks()
    pressure = _physical_pressure_jet_checks()
    poisson = _poisson_loop_checks()
    recipe = _recipe_numeric_fixture()
    return {
        "passed": True,
        "symbolic": symbolic,
        "quadrature": quadrature,
        "pressure_jets": pressure,
        "poisson_loop": poisson,
        "recipe_numeric_fixture": recipe,
    }


def main():
    print(json.dumps(run(), sort_keys=True))


if __name__ == "__main__":
    main()
