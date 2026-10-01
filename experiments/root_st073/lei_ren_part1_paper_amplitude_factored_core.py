"""Exact radial core recurrence in the gauge ``F = F0(Z) A``.

The unfactored recurrence in :mod:`lei_ren_part1_paper_core_recursion`
propagates the radial coefficients of ``F`` directly.  Here the axis value
is factored out.  The input ``ell_Z_taylor`` is the Taylor jet of
``ell = F0'/F0`` and ``S_Z_taylor`` is the jet of ``S = F0**2``.  The
returned ``A`` rows therefore have ``A[0] = 1`` and the physical swirl and
pressure terms retain their factors of ``S``.

All rows use ordinary Taylor coefficients in ``Z``.  The available axial
derivative depth is consumed in the same way as ``core_coefficients``:
the initial depth is the minimum of the four input jet lengths, and one
coefficient is removed after each radial step.  This is an exact algebraic
rewrite of the finite recurrence; it does not add a remainder or a global
core claim.
"""

from __future__ import annotations

import mpmath as mp


def factored_core_coefficients(
    Z,
    delta,
    *,
    ell_Z_taylor,
    S_Z_taylor,
    U0_Z_taylor,
    P0_Z_taylor,
    radial_degree=3,
    precision=160,
    scalar_converter=None,
):
    """Return radial rows for ``A``, ``Uz``, and physical ``P``.

    ``ell_Z_taylor[k]`` and ``S_Z_taylor[k]`` are the coefficients of
    ``F0'/F0`` and ``F0**2``.  ``A`` is normalized so its axis row is one.
    For every radial index ``n >= 1``, the returned pressure row is physical:

        P_n = S * sum(A_i A_(n-1-i)) / n.

    ``P[0]`` is copied from ``P0_Z_taylor``.  The ``scalar_converter`` hook
    follows the component arithmetic hook in ``core_coefficients``; all
    coefficients are converted before recurrence arithmetic.
    """

    with mp.workdps(precision):
        degree = int(radial_degree)
        if degree < 1 or degree != radial_degree:
            raise ValueError("Positive integer radial_degree required")
        K = min(
            len(ell_Z_taylor),
            len(S_Z_taylor),
            len(U0_Z_taylor),
            len(P0_Z_taylor),
        )
        initial_count = K
        if K < degree + 2:
            raise ValueError("Insufficient axis Z Taylor coefficients")

        convert = mp.mpf if scalar_converter is None else scalar_converter
        zero = convert(0)

        def pad(values):
            return list(map(convert, values)) + [convert(0)] * (K - len(values))

        def const(value):
            return pad([value])

        def add(*args):
            return [sum(arg[k] for arg in args) for k in range(K)]

        def scale(values, scalar):
            return [value * scalar for value in values]

        def mul(left, right):
            return [
                sum(left[i] * right[k - i] for i in range(k + 1))
                for k in range(K)
            ]

        def diff(values):
            return [
                (k + 1) * values[k + 1] if k + 1 < len(values) else zero
                for k in range(K)
            ]

        def inverse(values):
            quotient = [convert(1) / values[0]]
            for k in range(1, K):
                quotient.append(
                    -sum(
                        values[i] * quotient[k - i]
                        for i in range(1, k + 1)
                    )
                    / values[0]
                )
            return quotient

        z = pad([
            Z if hasattr(Z, "_mpi_") else mp.mpf(str(Z)),
            1,
        ])
        delta_value = mp.mpf(str(delta))
        dt = convert(delta_value)
        one = const(1)
        d = add(one, scale(mul(z, z), -1))
        L = add(one, scale(mul(z, z), -dt))
        invL = inverse(L)
        ell = pad(ell_Z_taylor)
        S = pad(S_Z_taylor)
        a = [const(1)]
        u = [pad(U0_Z_taylor)]
        pressure = [pad(P0_Z_taylor)]

        for n in range(degree):
            # One Z derivative is consumed by each radial recurrence, exactly
            # as in core_coefficients.
            K = initial_count - n - 1
            W = []
            H = []
            for i in range(n + 1):
                wi = scale(
                    add(
                        scale(mul(z, u[i]), 1 - dt),
                        mul(d, diff(u[i])),
                    ),
                    -convert(1) / (i + 1),
                )
                if i == 0:
                    wi = add(one, wi)
                W.append(wi)
                hi = mul(d, u[i])
                if i == 0:
                    hi = add(hi, scale(z, (1 - dt) / 2))
                H.append(hi)

            aa_u = const(0)
            uu = const(0)
            rhs_a = const(0)
            rhs_u = const(0)
            for i in range(n + 1):
                # The former F*Uz convolution is now A*Uz after removing F0.
                aa_u = add(aa_u, mul(u[i], a[n - i]))
                uu = add(uu, mul(u[i], u[n - i]))

                # F_Z becomes F0*(A_Z + ell*A), so F0 cancels from this
                # equation before the new A row is solved.
                a_derivative = add(
                    diff(a[n - i]),
                    mul(ell, a[n - i]),
                )
                rhs_a = add(
                    rhs_a,
                    scale(mul(W[i], a[n - i]), n - i + 1),
                    mul(H[i], a_derivative),
                )
                rhs_u = add(
                    rhs_u,
                    scale(mul(W[i], u[n - i]), n - i),
                    mul(H[i], diff(u[n - i])),
                )

            rhs_a = add(
                rhs_a,
                scale(add(a[n], scale(mul(z, aa_u), -2)), dt / 2),
            )
            rhs_u = add(
                rhs_u,
                scale(add(u[n], scale(mul(z, uu), -2)), (1 + dt) / 2),
                mul(d, diff(pressure[n])),
                scale(mul(z, pressure[n]), -2 * (1 + dt)),
            )

            # The axial swirl source is F^2 = S*A^2.  It enters only after
            # the first radial step, as in the unfactored recurrence.
            if n:
                aa_previous = const(0)
                for i in range(n):
                    aa_previous = add(
                        aa_previous,
                        mul(a[i], a[n - 1 - i]),
                    )
                rhs_u = add(
                    rhs_u,
                    scale(mul(z, mul(S, aa_previous)), -2),
                )

            a.append(
                scale(
                    mul(invL, rhs_a),
                    convert(1) / (2 * (n + 1) * (n + 2)),
                )
            )
            u.append(
                scale(
                    mul(invL, rhs_u),
                    convert(1) / (2 * (n + 1) ** 2),
                )
            )

            # P[0] is the supplied physical axis pressure.  Every generated
            # P[n+1] is the physical S times the A*A convolution.
            aa_current = const(0)
            for i in range(n + 1):
                aa_current = add(
                    aa_current,
                    mul(a[i], a[n - i]),
                )
            pressure.append(
                scale(
                    mul(S, aa_current),
                    convert(1) / (n + 1),
                )
            )

        return {
            "A": a,
            "Uz": u,
            "P": pressure,
            "Z": Z if hasattr(Z, "_mpi_") else mp.mpf(str(Z)),
            "delta": delta_value,
            "radial_degree": degree,
            "initial_Z_degree": initial_count - 1,
            "precision": precision,
            "component_scalar_arithmetic": scalar_converter is not None,
            "gauge_factored": True,
            "axis_amplitude_normalized": True,
            "physical_pressure_rows": True,
            "scope": "Local exact-equation radial jets in F=F0*A gauge; no temporal recursion or global matching.",
        }


# A descriptive alias for callers that prefer the module-level name.
amplitude_factored_core_coefficients = factored_core_coefficients


__all__ = [
    "factored_core_coefficients",
    "amplitude_factored_core_coefficients",
]
