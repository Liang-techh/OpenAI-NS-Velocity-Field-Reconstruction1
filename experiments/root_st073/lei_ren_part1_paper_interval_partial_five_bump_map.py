"""Directed partial five-bump map on the physical interval ``[1, 2]``.

This companion uses the frozen directed partial-integral backend and keeps
all five Section 10.2 rows in the coefficient order
``(c1, c2, xi1, xi2, xi3)``.  The map accepts first-order
``IntervalTaylor`` coefficients and returns the corresponding C1 rows
``L(x) h + Q(x, h)``.  Full-support mass rows use the exact algebraic
identity ``int gamma = 1``; every other full-support coefficient is loaded
from the immutable directed weight receipt.

The module is a fixed bump-map backend only.  It does not install a field,
solve controls, or certify a cone or source error.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

import mpmath as mp

from lei_ren_part1_paper_bump_integral_enclosures import (
    CENTERS,
    RADIUS,
    _fraction,
)
from lei_ren_part1_paper_bump_partial_enclosures import (
    FULL_MODULE,
    FULL_RECEIPT,
    PaperBumpPartialEnclosures,
)
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_interval_five_bump_inverse import weights as parent_weights
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


HERE = Path(__file__).resolve().parent
COEFFICIENT_ORDER = ("c1", "c2", "xi1", "xi2", "xi3")
CENTER_FRACTIONS = tuple(_fraction(value, "center") for value in CENTERS)
RADIUS_FRACTION = _fraction(RADIUS, "radius")


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _pack(value: Any) -> Any:
    if isinstance(value, IntervalTaylor):
        return [_pack(item) for item in value.coefficients]
    if isinstance(value, dict):
        return {key: _pack(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_pack(item) for item in value]
    return value


def _overlap(left: Any, right: Any) -> bool:
    ll, lh = endpoints(left)
    rl, rh = endpoints(right)
    return ll <= rh and rl <= lh


def _as_vector(value: Any, name: str, length: int = 5) -> tuple[Any, ...]:
    if hasattr(value, "rows") and hasattr(value, "cols"):
        if int(value.rows) * int(value.cols) != length:
            raise ValueError(f"{name} must contain exactly {length} entries")
        return tuple(value[index] for index in range(length))
    try:
        result = tuple(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an iterable of {length} entries") from exc
    if len(result) != length:
        raise ValueError(f"{name} must contain exactly {length} entries")
    return result


def _restore_in_context(ctx: Any, value: Any) -> Any:
    """Restore scalar/encoded intervals into this map's interval context."""

    if isinstance(value, IntervalTaylor):
        if value.ctx is ctx:
            return value
        return IntervalTaylor(ctx, [_restore_in_context(ctx, item) for item in value.coefficients])
    if isinstance(value, dict):
        return restore_value(ctx, value)
    if hasattr(value, "_mpi_") or hasattr(value, "_mpf_"):
        return restore_value(ctx, value)
    return value


def _coerce_jet(ctx: Any, value: Any, name: str) -> IntervalTaylor:
    if isinstance(value, IntervalTaylor):
        jet = _restore_in_context(ctx, value)
    elif isinstance(value, (list, tuple)):
        jet = IntervalTaylor(ctx, [_restore_in_context(ctx, item) for item in value])
    else:
        jet = IntervalTaylor(ctx, [_restore_in_context(ctx, value)])
    if jet.order < 1:
        raise ValueError(f"{name} must contain value and first C1 coefficient")
    return jet


class IntervalPartialFiveBumpMap:
    """Directed partial five-bump map with C1 coefficient arithmetic."""

    coefficient_order = COEFFICIENT_ORDER
    centers = CENTER_FRACTIONS
    radius = RADIUS_FRACTION

    def __init__(
        self,
        ctx: Any | None = None,
        *,
        precision: int = 80,
        order: int = 12,
        tolerance: Any = "1e-20",
        receipt_path: str | Path | None = None,
    ) -> None:
        self.backend = PaperBumpPartialEnclosures(
            precision=precision,
            order=order,
            tolerance=tolerance,
            receipt_path=receipt_path,
        )
        # The frozen partial backend remains at the receipt precision, while
        # callers such as the repaired reference field may carry a higher
        # precision target context.  Restore every backend interval into that
        # target context before it enters IntervalTaylor arithmetic.
        if ctx is not None and not hasattr(ctx, "mpf"):
            raise TypeError("ctx must provide an mpmath interval mpf method")
        self.ctx = ctx if ctx is not None else self.backend.ctx
        self.precision = self.backend.precision
        self.order = self.backend.order
        self.tolerance = self.backend.tolerance
        self.normalization = _restore_in_context(self.ctx, self.backend.normalization)
        self.full_module_sha256 = _hash(HERE / FULL_MODULE)
        self.full_receipt_sha256 = _hash(self.backend.full_receipt_path)
        self.metadata = dict(self.backend.metadata)
        self.metadata.update(
            {
                "coefficient_order": list(COEFFICIENT_ORDER),
                "partial_map": True,
                "full_support_mass_exact_one": True,
                "support_disjoint_mixed_products_only_0_2": True,
                "no_midpoint_weights": True,
                "map_enclosure": True,
                "field_installed": False,
                "controls_solved": False,
                "cone_certified": False,
            }
        )

    def _validate_x(self, x: Any) -> Fraction:
        value = _fraction(x, "x")
        if value < 1 or value > 2:
            raise ValueError("x must be an exact rational in [1, 2]")
        return value

    def _weight(self, center: Fraction, power: Any, x: Fraction, multiplicity: int) -> Any:
        value = _restore_in_context(
            self.ctx, self.backend.weight(center, power, x, multiplicity)
        )
        # The first row is the normalized mass identity.  Use it exactly on
        # completed supports, while retaining the saved directed intervals
        # for all weighted and quadratic coefficients.
        if (
            multiplicity == 1
            and _fraction(power, "power") == 0
            and x >= center + self.radius
        ):
            return self.ctx.mpf(1)
        return value

    def weights_at(self, x: Any) -> dict[str, Any]:
        """Return the directed L, fg, gg, ff, and ff_over_x weights at x."""

        x_fraction = self._validate_x(x)
        c = self.ctx
        zero = c.mpf(0)

        def w(index: int, power: Any, multiplicity: int = 1) -> Any:
            return self._weight(
                CENTER_FRACTIONS[index], power, x_fraction, multiplicity
            )

        L = [[zero for _ in range(5)] for _ in range(5)]
        # b1=gamma_1 and b2=gamma_3.
        L[0][0] = w(0, 0)
        L[0][1] = w(2, 0)
        L[1][0] = w(0, "3/5")
        L[1][1] = w(2, "3/5")
        for index in range(3):
            column = 2 + index
            L[2][column] = w(index, "1/2")
            L[3][column] = -w(index, "1/10")
            L[4][column] = w(index, "-9/10")
        return dict(
            L=L,
            fg=[w(0, "1/2", 2), w(2, "1/2", 2)],
            gg=[w(0, 0, 2), w(2, 0, 2)],
            ff=[w(index, 0, 2) for index in range(3)],
            ff_over_x=[w(index, -1, 2) for index in range(3)],
        )

    def bump_values(self, x: Any) -> dict[str, Any]:
        """Return directed beta and d beta/dx values at exact rational x.

        Support endpoints are handled before forming ``1-t**2``; therefore
        the singular-looking formula is never evaluated at ``|t| = 1``.
        """

        x_fraction = self._validate_x(x)
        c = self.ctx
        values = []
        derivatives = []
        support_status = []
        radius = c.mpf(RADIUS_FRACTION.numerator) / c.mpf(RADIUS_FRACTION.denominator)
        for center in CENTER_FRACTIONS:
            left, right = center - RADIUS_FRACTION, center + RADIUS_FRACTION
            if x_fraction <= left or x_fraction >= right:
                values.append(c.mpf(0))
                derivatives.append(c.mpf(0))
                support_status.append("outside_or_endpoint_exact_zero")
                continue
            x_point = c.mpf(x_fraction.numerator) / c.mpf(x_fraction.denominator)
            center_point = c.mpf(center.numerator) / c.mpf(center.denominator)
            t = (x_point - center_point) / radius
            q = 1 - t * t
            if endpoints(q)[0] <= 0:
                raise ArithmeticError("interior bump coordinate lost q positivity")
            raw = c.exp(-1 / q)
            gamma = raw / (radius * self.normalization)
            derivative = gamma * (-2 * t / (radius * q * q))
            values.append(gamma)
            derivatives.append(derivative)
            support_status.append("interior_directed")
        return dict(
            x=str(x_fraction),
            values=values,
            derivatives=derivatives,
            centers=[str(center) for center in CENTER_FRACTIONS],
            support_status=support_status,
        )

    def apply(self, h: Any, invAm2: Any, x: Any) -> list[Any]:
        """Return the five directed C1 rows ``L(x)h + Q(x,h)``."""

        weights = self.weights_at(x)
        coefficients = tuple(
            _coerce_jet(self.ctx, item, f"h[{index}]")
            for index, item in enumerate(_as_vector(h, "h"))
        )
        inv = _coerce_jet(self.ctx, invAm2, "invAm2")
        c1, c2, x1, x2, x3 = coefficients
        linear = [
            sum((coefficients[column] * weights["L"][row][column] for column in range(5)), coefficients[0] * 0)
            for row in range(5)
        ]
        quadratic = [
            coefficients[0] * 0,
            (c1 * x1) * weights["fg"][0] + (c2 * x3) * weights["fg"][1],
            coefficients[0] * 0,
            inv
            * sum(
                ((value * value) * weights["gg"][index] for index, value in enumerate((c1, c2))),
                coefficients[0] * 0,
            )
            - sum(
                ((value * value) * weights["ff"][index] for index, value in enumerate((x1, x2, x3))),
                coefficients[0] * 0,
            )
            / 2,
            sum(
                ((value * value) * weights["ff_over_x"][index] for index, value in enumerate((x1, x2, x3))),
                coefficients[0] * 0,
            )
            / 2,
        ]
        return [linear[row] + quadratic[row] for row in range(5)]


def bump_values(backend: IntervalPartialFiveBumpMap, x: Any) -> dict[str, Any]:
    """Convenience wrapper for the class method."""

    return backend.bump_values(x)


def _independent_bump(t: mp.mpf) -> mp.mpf:
    if abs(t) >= 1:
        return mp.mpf(0)
    return mp.exp(-1 / (1 - t * t))


def _independent_weights(x: Fraction) -> dict[str, Any]:
    """Independent high-precision quadrature fixture, separate from backend."""

    with mp.workdps(140):
        xx = mp.mpf(x.numerator) / x.denominator
        radius = mp.mpf(1) / 40
        centers = [mp.mpf("1.25"), mp.mpf("1.5"), mp.mpf("1.75")]
        norm = mp.quad(_independent_bump, [-1, 0, 1])

        def gamma(center, y):
            return _independent_bump((y - center) / radius) / (radius * norm)

        def integrate(center, power, multiplicity=1):
            left = max(mp.mpf(1), center - radius)
            right = min(xx, center + radius)
            if right <= left:
                return mp.mpf(0)
            return mp.quad(
                lambda y: y**power * gamma(center, y) ** multiplicity,
                [left, (left + right) / 2, right],
            )

        def w(index, power, multiplicity=1):
            return integrate(centers[index], mp.mpf(power), multiplicity)

        L = [[mp.mpf(0) for _ in range(5)] for _ in range(5)]
        L[0][0] = w(0, 0)
        L[0][1] = w(2, 0)
        L[1][0] = w(0, mp.mpf(".6"))
        L[1][1] = w(2, mp.mpf(".6"))
        for index in range(3):
            L[2][2 + index] = w(index, mp.mpf(".5"))
            L[3][2 + index] = -w(index, mp.mpf(".1"))
            L[4][2 + index] = w(index, mp.mpf("-.9"))
        return dict(
            L=L,
            fg=[w(0, mp.mpf(".5"), 2), w(2, mp.mpf(".5"), 2)],
            gg=[w(0, 0, 2), w(2, 0, 2)],
            ff=[w(index, 0, 2) for index in range(3)],
            ff_over_x=[w(index, -1, 2) for index in range(3)],
        )


def _independent_apply(W, h, inv):
    c1, c2, x1, x2, x3 = h
    linear = [
        sum(W["L"][row][column] * h[column] for column in range(5))
        for row in range(5)
    ]
    quadratic = [
        0,
        W["fg"][0] * c1 * x1 + W["fg"][1] * c2 * x3,
        0,
        inv * (W["gg"][0] * c1 * c1 + W["gg"][1] * c2 * c2)
        - sum(W["ff"][i] * h[2 + i] ** 2 for i in range(3)) / 2,
        sum(W["ff_over_x"][i] * h[2 + i] ** 2 for i in range(3)) / 2,
    ]
    return [linear[i] + quadratic[i] for i in range(5)]


def _independent_apply_derivative(W, h, dh, inv, dinv):
    """Differentiate the independent quadratic map at ``(h, inv)``.

    This is intentionally separate from ``_independent_apply``: evaluating
    the latter at ``dh`` would square the derivative vector and is not the
    derivative of the quadratic terms.
    """

    c1, c2, x1, x2, x3 = h
    dc1, dc2, dx1, dx2, dx3 = dh
    linear = [
        sum(W["L"][row][column] * dh[column] for column in range(5))
        for row in range(5)
    ]
    quadratic_derivative = [
        0,
        W["fg"][0] * (dc1 * x1 + c1 * dx1)
        + W["fg"][1] * (dc2 * x3 + c2 * dx3),
        0,
        dinv * (W["gg"][0] * c1 * c1 + W["gg"][1] * c2 * c2)
        + inv
        * (
            W["gg"][0] * 2 * c1 * dc1
            + W["gg"][1] * 2 * c2 * dc2
        )
        - sum(
            W["ff"][i] * 2 * h[2 + i] * dh[2 + i]
            for i in range(3)
        )
        / 2,
        sum(
            W["ff_over_x"][i] * 2 * h[2 + i] * dh[2 + i]
            for i in range(3)
        )
        / 2,
    ]
    return [linear[i] + quadratic_derivative[i] for i in range(5)]


def run() -> dict[str, Any]:
    """Run endpoint, bump, full-map, and independent C1 fixture checks."""

    backend = IntervalPartialFiveBumpMap()
    c = backend.ctx
    parent_raw = json.loads(
        (HERE / "lei_ren_part1_paper_bump_integral_enclosures_check.json").read_text(
            encoding="utf-8"
        )
    )
    parent_W = parent_weights(c, parent_raw)
    full_W = backend.weights_at(Fraction(2))
    overlap_count = 0
    for name in ("L", "fg", "gg", "ff", "ff_over_x"):
        left = full_W[name]
        right = parent_W[name]
        if isinstance(left, list):
            left_values = [item for row in left for item in (row if isinstance(row, list) else [row])]
            right_values = [item for row in right for item in (row if isinstance(row, list) else [row])]
        else:
            left_values, right_values = left, right
        for a, b in zip(left_values, right_values):
            if not _overlap(a, b):
                raise AssertionError(f"full x=2 coefficient does not overlap parent W: {name}")
            overlap_count += 1
    if endpoints(full_W["L"][0][0]) != endpoints(c.mpf(1)):
        raise AssertionError("completed first axial mass is not exact one")

    edge_checks = 0
    for center in CENTER_FRACTIONS:
        for edge in (center - RADIUS_FRACTION, center + RADIUS_FRACTION):
            packet = backend.bump_values(edge)
            index = CENTER_FRACTIONS.index(center)
            if endpoints(packet["values"][index]) != endpoints(c.mpf(0)):
                raise AssertionError("bump value did not vanish at support endpoint")
            if endpoints(packet["derivatives"][index]) != endpoints(c.mpf(0)):
                raise AssertionError("bump derivative did not vanish at support endpoint")
            edge_checks += 2

    with mp.workdps(140):
        x = Fraction(3, 2)
        W = _independent_weights(x)
        h_value = [mp.mpf("1e-14"), mp.mpf("-2e-14"), mp.mpf("3e-32"), mp.mpf("-4e-32"), mp.mpf("5e-32")]
        h_derivative = [mp.mpf("2e-15"), mp.mpf("-3e-15"), mp.mpf("4e-33"), mp.mpf("-5e-33"), mp.mpf("6e-33")]
        inv_value, inv_derivative = mp.mpf(".7"), mp.mpf(".03")
        expected_value = _independent_apply(W, h_value, inv_value)
        expected_derivative = _independent_apply_derivative(
            W, h_value, h_derivative, inv_value, inv_derivative
        )
        h = [IntervalTaylor(c, [c.mpf(value), c.mpf(derivative)]) for value, derivative in zip(h_value, h_derivative)]
        inv = IntervalTaylor(c, [c.mpf(inv_value), c.mpf(inv_derivative)])
        actual = backend.apply(h, inv, x)
        c1_checks = 0
        for row in range(5):
            for order, reference in enumerate((expected_value[row], expected_derivative[row])):
                lo, hi = endpoints(actual[row][order])
                if not lo <= reference <= hi:
                    raise AssertionError(("independent C1 map coefficient excluded", row, order))
                c1_checks += 1

    report = dict(
        precision=backend.precision,
        order=backend.order,
        tolerance="1e-20",
        coefficient_order=list(COEFFICIENT_ORDER),
        normalization=backend.normalization,
        full_module_sha256=backend.full_module_sha256,
        full_receipt_sha256=backend.full_receipt_sha256,
        full_x2_parent_coefficient_overlap_count=overlap_count,
        full_x2_first_mass_exact_one=True,
        support_endpoint_value_derivative_zero_checks=edge_checks,
        independent_clipped_center="3/2",
        independent_c1_map_coefficients_contained=c1_checks,
        directed_partial_backend=True,
        independent_quadrature_used=True,
        no_midpoint_weights=True,
        source_error_enclosed=False,
        parameter_error_enclosed=False,
        cone_certified=False,
        field_installed=False,
        controls_solved=False,
    )
    report["input_hashes"] = {
        name: _hash(HERE / name)
        for name in (
            Path(__file__).name,
            FULL_MODULE,
            FULL_RECEIPT,
            "lei_ren_part1_paper_bump_partial_enclosures.py",
            "lei_ren_part1_paper_interval_five_bump_inverse.py",
        )
    }
    output = HERE / "lei_ren_part1_paper_interval_partial_five_bump_map.json"
    output.write_text(json.dumps(encode(_pack(report)), indent=2) + "\n", encoding="utf-8")
    print(
        "Directed partial five-bump map checks:",
        overlap_count,
        "full coefficients overlapped;",
        c1_checks,
        "C1 coefficients contained",
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
