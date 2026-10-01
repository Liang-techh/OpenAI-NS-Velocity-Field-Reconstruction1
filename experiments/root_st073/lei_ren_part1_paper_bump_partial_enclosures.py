"""Directed partial Section 10.2 bump integrals.

This adapter keeps the frozen full-support enclosure backend unchanged.  It
loads the exact full-weight and normalizer intervals from its receipt and
returns those intervals verbatim whenever the requested upper endpoint has
passed a bump support.  A clipped support is integrated in the transformed
coordinate ``t=(x-center)/(1/40)`` using the frozen backend's directed Taylor
cell primitive.  Endpoint strips use a positive analytic bound, so the raw
bump formula is never evaluated at ``|t|=1``.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import operator
from pathlib import Path
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_bump_integral_enclosures import (
    CENTERS,
    RADIUS,
    PaperBumpIntegralEnclosures,
    _fraction,
    _interval_endpoints,
    _record_interval,
)


FULL_MODULE = "lei_ren_part1_paper_bump_integral_enclosures.py"
FULL_RECEIPT = "lei_ren_part1_paper_bump_integral_enclosures_check.json"


def _index(value: Any, name: str) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return int(result)


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _positive_interval(ctx: MPIntervalContext, value: Any) -> Any:
    """Return a directed ``[0, upper(value)]`` interval."""

    upper = _interval_endpoints(value)[1]
    return ctx.mpf([ctx.mpf(0), upper])


def _decode_interval(ctx: MPIntervalContext, record: Any) -> Any:
    if isinstance(record, dict) and "lower_exact_mpf_tuple" in record:
        lower = mp.make_mpf(tuple(record["lower_exact_mpf_tuple"]))
        upper = mp.make_mpf(tuple(record["upper_exact_mpf_tuple"]))
        return ctx.mpf([lower, upper])
    if isinstance(record, dict) and "lower" in record and "upper" in record:
        return ctx.mpf([ctx.mpf(str(record["lower"])), ctx.mpf(str(record["upper"]))])
    raise ValueError("interval record is missing directed endpoints")


def _overlap(left: Any, right: Any) -> bool:
    left_lo, left_hi = _interval_endpoints(left)
    right_lo, right_hi = _interval_endpoints(right)
    return left_lo <= right_hi and right_lo <= left_hi


class PaperBumpPartialEnclosures:
    """Directed ``int_1^x`` weights for the fixed paper bumps."""

    def __init__(
        self,
        *,
        precision: int = 80,
        order: int = 12,
        tolerance: Any = "1e-20",
        receipt_path: str | Path | None = None,
    ) -> None:
        self.precision = _index(precision, "precision")
        self.order = _index(order, "order")
        self.tolerance = _fraction(tolerance, "tolerance")
        if self.tolerance <= 0:
            raise ValueError("tolerance must be positive")
        base = Path(__file__).parent
        self.full_module_path = base / FULL_MODULE
        self.full_receipt_path = Path(receipt_path) if receipt_path else base / FULL_RECEIPT
        if not self.full_module_path.exists() or not self.full_receipt_path.exists():
            raise FileNotFoundError("frozen full-weight module/receipt is unavailable")
        self.receipt = json.loads(self.full_receipt_path.read_text(encoding="utf-8"))
        module_hash = hashlib.sha256(self.full_module_path.read_bytes()).hexdigest()
        if self.receipt.get("module_sha256") != module_hash:
            raise ValueError("full-weight receipt does not match frozen module hash")
        settings = self.receipt.get("settings", {})
        if int(settings.get("precision", self.precision)) != self.precision:
            raise ValueError("partial backend precision must match frozen receipt precision")
        if int(settings.get("Taylor_order", self.order)) != self.order:
            raise ValueError("partial backend order must match frozen receipt Taylor order")
        self.ctx = MPIntervalContext()
        self.ctx.dps = self.precision
        self.radius = self._point(RADIUS)
        self.centers = tuple(_fraction(value, "center") for value in CENTERS)
        self.supports = {
            center: (center - RADIUS, center + RADIUS) for center in self.centers
        }
        self.full_weights = self._load_full_weights()
        self.normalization = _decode_interval(self.ctx, self.receipt["normalization"])
        # Reuse the frozen Taylor and derivative-cell implementation for
        # clipped cells.  Its normalizer is independently recomputed only as
        # a consistency check; full-support calls below never use it.
        self._backend = PaperBumpIntegralEnclosures(
            precision=self.precision,
            order=self.order,
            tolerance=_fraction_text(self.tolerance),
            initial_cells=int(settings.get("initial_cells", 16)),
            max_depth=int(settings.get("max_depth", 24)),
        )
        if not _overlap(self.normalization, self._backend.normalization):
            raise ValueError("saved and recomputed normalizer intervals are disjoint")
        self.tail_width_fraction = self._backend.tail_width_fraction
        self.tail_cut_fraction = self._backend.tail_cut_fraction
        self.tail_width = self._backend.tail_width
        self._cache: dict[tuple[Fraction, Fraction, Fraction, int], Any] = {}
        self._records: dict[tuple[Fraction, Fraction, Fraction, int], dict[str, Any]] = {}
        self.metadata = {
            "source_section": "Part I Section 10.2, equations (10.7)--(10.10)",
            "full_module": FULL_MODULE,
            "full_receipt": str(self.full_receipt_path.name),
            "full_module_sha256": module_hash,
            "normalization_loaded_from_receipt": True,
            "full_support_loaded_from_receipt": True,
            "clipped_cells_reuse_frozen_backend": True,
            "endpoint_tail_enclosed": True,
            "directed_interval_arithmetic": True,
            "source_error_enclosed": False,
            "parameter_error_enclosed": False,
            "cone_certified": False,
        }

    def _point(self, value: Any) -> Any:
        rational = value if isinstance(value, Fraction) else _fraction(value)
        return self.ctx.mpf(rational.numerator) / self.ctx.mpf(rational.denominator)

    def _load_full_weights(self) -> dict[tuple[Fraction, Fraction, int], Any]:
        weights: dict[tuple[Fraction, Fraction, int], Any] = {}
        records = self.receipt.get("weight_records", {})
        for key, item in records.items():
            try:
                center_text, power_text, multiplicity_text = key.split("|", 2)
                multiplicity = int(multiplicity_text.removeprefix("k"))
                center = _fraction(center_text, "center")
                power = _fraction(power_text, "power")
            except (ValueError, TypeError) as exc:
                raise ValueError(f"invalid full-weight receipt key: {key!r}") from exc
            weights[(center, power, multiplicity)] = _decode_interval(
                self.ctx, item["weight_interval"]
            )
        if not weights:
            raise ValueError("full-weight receipt contains no weight records")
        return weights

    @staticmethod
    def _key(center: Fraction, power: Fraction, x: Fraction, multiplicity: int):
        return center, power, x, multiplicity

    def _tail_segment_bound(
        self,
        center: Fraction,
        power: Fraction,
        multiplicity: int,
        left: Fraction,
        right: Fraction,
    ) -> Any:
        if right <= left:
            return self.ctx.mpf(0)
        x_left = self._point(center) + self.radius * self._point(left)
        x_right = self._point(center) + self.radius * self._point(right)
        x_lo = min(_interval_endpoints(x_left)[0], _interval_endpoints(x_right)[0])
        x_hi = max(_interval_endpoints(x_left)[1], _interval_endpoints(x_right)[1])
        x_interval = self.ctx.mpf([x_lo, x_hi])
        weight = self._backend._power(x_interval, power)
        bump_bound = self.ctx.exp(
            -self.ctx.mpf(multiplicity) / (2 * self.tail_width)
        )
        width = self._point(right - left)
        return _positive_interval(self.ctx, width * weight * bump_bound)

    def _integrate_clipped_raw(
        self,
        center: Fraction,
        power: Fraction,
        multiplicity: int,
        endpoint: Fraction,
    ) -> tuple[Any, dict[str, Any]]:
        cut = self.tail_cut_fraction
        total = self.ctx.mpf(0)
        tail_parts = []
        stats = {
            "accepted_cells": 0,
            "splits": 0,
            "max_depth": 0,
            "remainder_interval": self.ctx.mpf(0),
        }
        if endpoint <= -cut:
            tail = self._tail_segment_bound(
                center, power, multiplicity, Fraction(-1), endpoint
            )
            total += tail
            tail_parts.append((Fraction(-1), endpoint, tail))
            return total, {
                "tail_parts": tail_parts,
                "interior_cells": stats,
            }
        # The full left endpoint strip is included once the clipped endpoint
        # reaches the interior.  The right strip is included only when clipped.
        left_tail = self._tail_segment_bound(
            center, power, multiplicity, Fraction(-1), -cut
        )
        total += left_tail
        tail_parts.append((Fraction(-1), -cut, left_tail))
        interior_end = min(endpoint, cut)
        interior_length = interior_end + cut
        if interior_length > 0:
            full_interior = 2 * cut
            target = self._backend.tolerance / 2
            for index in range(self._backend.initial_cells):
                left = -cut + full_interior * index / self._backend.initial_cells
                right = -cut + full_interior * (index + 1) / self._backend.initial_cells
                clipped_left = max(left, -cut)
                clipped_right = min(right, interior_end)
                if clipped_right <= clipped_left:
                    continue
                budget = target * (
                    self._point(clipped_right - clipped_left)
                    / self._point(interior_length)
                )
                total += self._backend._adaptive_cell(
                    center,
                    power,
                    multiplicity,
                    clipped_left,
                    clipped_right,
                    budget,
                    0,
                    stats,
                )
        if endpoint > cut:
            right_tail = self._tail_segment_bound(
                center, power, multiplicity, cut, endpoint
            )
            total += right_tail
            tail_parts.append((cut, endpoint, right_tail))
        return total, {
            "tail_parts": tail_parts,
            "interior_cells": stats,
        }

    def _normalize(self, raw: Any, multiplicity: int) -> Any:
        if multiplicity == 1:
            return raw / self.normalization
        return raw / (self.radius * self.normalization * self.normalization)

    def weight(
        self,
        center: Any,
        power: Any,
        x: Any,
        multiplicity: int = 1,
    ) -> Any:
        """Return a directed enclosure of ``int_1^x y**power*gamma**k dy``."""

        center_fraction = _fraction(center, "center")
        power_fraction = _fraction(power, "power")
        x_fraction = _fraction(x, "x")
        multiplicity = _index(multiplicity, "multiplicity")
        if center_fraction not in self.supports:
            raise ValueError("center is not one of the three paper bumps")
        if multiplicity not in (1, 2):
            raise ValueError("multiplicity must be 1 or 2")
        left_support, right_support = self.supports[center_fraction]
        # The caller may pass a radius beyond [1,2]; the paper partial map is
        # naturally clamped at its declared outer endpoint.
        clipped_x = max(Fraction(1), min(Fraction(2), x_fraction))
        key = self._key(center_fraction, power_fraction, clipped_x, multiplicity)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        if clipped_x <= left_support:
            value = self.ctx.mpf(0)
            record = {
                "interval": _record_interval(value),
                "status": "lower_empty_exact_zero",
                "x_input": _fraction_text(x_fraction),
                "x_clipped": _fraction_text(clipped_x),
                "full_support_saved": False,
            }
        elif clipped_x >= right_support:
            full_key = (center_fraction, power_fraction, multiplicity)
            if full_key not in self.full_weights:
                raise KeyError(
                    "full-support weight is absent from the frozen receipt; "
                    "cannot safely return a recomputed full interval"
                )
            value = self.full_weights[full_key]
            record = {
                "interval": _record_interval(value),
                "status": "full_support_saved_interval",
                "x_input": _fraction_text(x_fraction),
                "x_clipped": _fraction_text(clipped_x),
                "full_support_saved": True,
                "full_receipt_key": f"{_fraction_text(center_fraction)}|{_fraction_text(power_fraction)}|k{multiplicity}",
            }
        else:
            endpoint = (clipped_x - center_fraction) / RADIUS
            raw, raw_record = self._integrate_clipped_raw(
                center_fraction, power_fraction, multiplicity, endpoint
            )
            value = self._normalize(raw, multiplicity)
            # The integrand is nonnegative; safely remove only any tiny
            # negative lower endpoint introduced by interval Taylor padding.
            lower, upper = _interval_endpoints(value)
            value = self.ctx.mpf([max(mp.mpf(0), lower), upper])
            record = {
                "interval": _record_interval(value),
                "status": "clipped_taylor_interval",
                "x_input": _fraction_text(x_fraction),
                "x_clipped": _fraction_text(clipped_x),
                "t_endpoint": _fraction_text(endpoint),
                "full_support_saved": False,
                "raw_interval": _record_interval(raw),
                "raw_details": raw_record,
            }
        self._cache[key] = value
        self._records[key] = record
        return value

    def integral_record(
        self,
        center: Any,
        power: Any,
        x: Any,
        multiplicity: int = 1,
    ) -> dict[str, Any]:
        center_fraction = _fraction(center, "center")
        power_fraction = _fraction(power, "power")
        x_fraction = max(Fraction(1), min(Fraction(2), _fraction(x, "x")))
        key = self._key(center_fraction, power_fraction, x_fraction, _index(multiplicity, "multiplicity"))
        self.weight(center, power, x, multiplicity)
        return dict(self._records[key])


__all__ = ["PaperBumpPartialEnclosures"]

