"""Directed interval integral for the finite angular flattening stage.

The stage is evaluated pointwise at one declared ``Z``.  It encloses the
normalized pressure integral and its first two Z derivatives while retaining
the stored schedule amplitude and the finite ``y_v`` to ``y_f`` interval.
This is a local radial datum, not a uniform-Z or whole-source certificate.
"""

import mpmath as mp

from lei_ren_part1_paper_high_order_preheat_integrals import switch_taylor
from lei_ren_part1_paper_interval_taylor import (
    IntervalTaylor,
    constant,
    integrate_symmetric,
)
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def _as_interval(ctx, value):
    return value if hasattr(value, "_mpi_") else ctx.mpf(str(value))


def _width(value):
    lo, hi = endpoints(value)
    return mp.mpf(hi) - mp.mpf(lo)


def _density_jets(ctx, value, cell, *, base, length, slope, z, q, log_q, order):
    """Return density, first-Z, and second-Z Taylor jets in local u."""

    center_u = IntervalTaylor.variable(ctx, value, order)
    cell_u = IntervalTaylor.variable(ctx, cell, order)
    sigma_center = switch_taylor(ctx, value, order)
    sigma_cell = switch_taylor(ctx, cell, order)
    beta_center = 2 * (1 - sigma_center)
    beta_cell = 2 * (1 - sigma_cell)
    ell_center = constant(ctx, base, order) + center_u * (slope * length)
    ell_cell = constant(ctx, base, order) + cell_u * (slope * length)
    log_two = ctx.ln(2)

    def make_density(ell, beta):
        logarithm = 2 * ell + (beta - 3) * log_two - beta * log_q
        density = logarithm.exp()
        first_factor = -2 * beta * z / q
        second_factor = -2 * beta / q + 4 * beta * (beta + 1) * z * z / (q * q)
        return density, density * first_factor, density * second_factor

    return make_density(ell_center, beta_center), make_density(ell_cell, beta_cell), sigma_cell


def _edge_integral(ctx, left, right, *, base, length, slope, z, q, log_q, beta_base, beta_cell):
    """Integrate the affine exponential edge and bound its beta correction."""

    width = right - left
    rate = 2 * slope * length
    if endpoints(rate)[1] >= 0:
        raise ValueError("flatten affine edge requires a strictly negative rate")
    # Integral in u of exp(2*base + rate*u), evaluated without subtraction
    # of nearly equal endpoint exponentials.
    baseline = (
        ctx.exp(2 * base + rate * left)
        * ctx.expm1(rate * width)
        / rate
    )
    baseline *= ctx.exp((beta_base - 3) * ctx.ln(2) - beta_base * log_q)
    beta_delta = beta_cell - beta_base
    correction = ctx.exp(beta_delta * (ctx.ln(2) - log_q))
    value_u = baseline * correction
    value = length * value_u
    first_factor = -2 * beta_cell * z / q
    second_factor = -2 * beta_cell / q + 4 * beta_cell * (beta_cell + 1) * z * z / (q * q)
    return value, value * first_factor, value * second_factor


def integrate_flatten_core(
    enclosures,
    *,
    base,
    length,
    mu,
    Z=".3",
    panels=128,
    order=12,
):
    """Integrate a flatten interval from explicit ``base``, ``length`` and ``mu``.

    ``base`` is the interval enclosure for ``ell(y_v) = log(A(y_v)/Pstar)``;
    ``length`` is the stored y-length and may be any positive resolved value
    for fixture or local use.  The public schedule wrapper below supplies the
    actual stored ``y_f-y_v`` length.
    """

    if isinstance(panels, bool) or not isinstance(panels, int) or panels < 4:
        raise ValueError("Require integer panels>=4")
    if isinstance(order, bool) or not isinstance(order, int) or order < 2:
        raise ValueError("Require integer order>=2")
    iv = enclosures.iv
    base = _as_interval(iv, base)
    length = _as_interval(iv, length)
    mu = _as_interval(iv, mu)
    z = _as_interval(iv, Z)
    if endpoints(length)[0] <= 0:
        raise ValueError("flatten length must be positive")
    if endpoints(mu)[0] <= 0:
        raise ValueError("flatten mu must be positive")
    zlo, zhi = endpoints(z)
    if zlo < -1 or zhi > 1:
        raise ValueError("Require pointwise Z in [-1,1]")

    slope = -iv.mpf(".5") - mu
    q = iv.mpf(1) + z * z
    log_q = iv.ln(q)
    h = iv.mpf(1) / panels
    radius = h / 2
    total = iv.mpf(0)
    first = iv.mpf(0)
    second = iv.mpf(0)
    edge_records = []

    for index in range(panels):
        left = h * index
        right = h * (index + 1)
        cell = iv.mpf([endpoints(left)[0], endpoints(right)[1]])
        if index == 0 or index == panels - 1:
            # Do not differentiate sigma at a flat endpoint.  The affine
            # amplitude is integrated exactly and only the tiny beta change
            # is enclosed by a monotone interval correction.
            beta_cell = 2 * (1 - enclosures.sigma_interval(cell))
            beta_base = iv.mpf(2 if index == 0 else 0)
            value, value_first, value_second = _edge_integral(
                iv,
                left,
                right,
                base=base,
                length=length,
                slope=slope,
                z=z,
                q=q,
                log_q=log_q,
                beta_base=beta_base,
                beta_cell=beta_cell,
            )
            edge_records.append({"index": index, "beta_interval": beta_cell, "method": "analytic_affine_edge"})
        else:
            center = (left + right) / 2
            (density_center, first_center, second_center), (density_cell, first_cell, second_cell), _ = _density_jets(
                iv,
                center,
                cell,
                base=base,
                length=length,
                slope=slope,
                z=z,
                q=q,
                log_q=log_q,
                order=order,
            )
            value = integrate_symmetric(density_center, density_cell, radius)
            value_first = integrate_symmetric(first_center, first_cell, radius)
            value_second = integrate_symmetric(second_center, second_cell, radius)
            value *= length
            value_first *= length
            value_second *= length
        total += value
        first += value_first
        second += value_second

    return {
        "integrals": {"value": total, "first_Z": first, "second_Z": second},
        "normalized_mass_interval": total,
        "normalized_mass_at_Z_interval": total,
        "normalized_mass_Z_interval": first,
        "normalized_mass_ZZ_interval": second,
        "normalized_pressure_integral_interval": total,
        "normalized_pressure_integral_Z_interval": first,
        "normalized_pressure_integral_ZZ_interval": second,
        "value_interval": total,
        "first_Z_interval": first,
        "second_Z_interval": second,
        "panels": panels,
        "Taylor_order": order,
        "Z": z,
        "stored_flatten_length": length,
        "ell_yv_interval": base,
        "mu_interval": mu,
        "edge_panels": edge_records,
        "edge_method": "analytic affine amplitude with interval beta correction",
        "edge_precision_limitation": True,
        "pointwise_Z_only": True,
        "uniform_Z_bound": False,
        "directed_interval_arithmetic": True,
        "source_scope": "stored schedule parameters",
        "original_parameter_errors_enclosed": False,
        "quadrature_error_enclosed": True,
        "full_pressure_error_enclosed": False,
        "five_defect_interval_closure": False,
        "interval_widths": {
            "value": _width(total),
            "first_Z": _width(first),
            "second_Z": _width(second),
        },
    }


def integrate_flatten(enclosures, *, Z=".3", panels=128, order=12):
    """Integrate the stored schedule's ``z_flatten`` stage pointwise in Z."""

    schedule = enclosures.schedule
    left, right = schedule._stage_bounds["z_flatten"]
    enclosures.stage_pressure_upper("z_flatten")
    if right - left != 100:
        raise ValueError("stored z_flatten length must equal 100")
    base = enclosures.log_amplitude_ratio(left)["interval"]
    length = enclosures.scalar(right) - enclosures.scalar(left)
    result = integrate_flatten_core(
        enclosures,
        base=base,
        length=length,
        mu=enclosures.scalar(schedule.mu),
        Z=Z,
        panels=panels,
        order=order,
    )
    result.update(
        stage="z_flatten",
        stored_left=left,
        stored_right=right,
        stored_length_expected=100,
        stored_length_is_100=(right - left == 100),
    )
    return result


__all__ = ["integrate_flatten", "integrate_flatten_core"]
