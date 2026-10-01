"""Directed finite radial-core regeneration with a coherent pressure error.

The accepted exact inlet JSON supplies the finite axis rows at ``Z=.3``.  This
module reruns only the finite radial recursion with a rectangular interval
pressure ring.  The supplied error receipt is the combined normalized
preheat pressure Taylor error; conservatively, its physical bound is applied
independently to pressure powers 0 and 1 of every initial pressure row.

This is conditional on the stored schedule parameters and finite datum.  It
does not enclose source-parameter, axis-data generation, pressure-order,
radial-series, or later collar/ODE errors.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp

try:  # package imports
    from .lei_ren_part1_paper_component_pressure_core import PressurePolynomial
    from .lei_ren_part1_paper_core_recursion import core_coefficients
    from .lei_ren_part1_paper_interval_collar_core_inlet import (
        DEFAULT_CACHE as DEFAULT_INPUT,
        _endpoints,
        _read_cache,
    )
    from .lei_ren_part1_paper_interval_difference import IntervalDifference
    from .lei_ren_part1_paper_interval_pressure_width_jet import (
        IntervalPressureWidthJet,
    )
except (ImportError, ValueError):  # direct experiment execution
    from lei_ren_part1_paper_component_pressure_core import PressurePolynomial
    from lei_ren_part1_paper_core_recursion import core_coefficients
    from lei_ren_part1_paper_interval_collar_core_inlet import (
        DEFAULT_CACHE as DEFAULT_INPUT,
        _endpoints,
        _read_cache,
    )
    from lei_ren_part1_paper_interval_difference import IntervalDifference
    from lei_ren_part1_paper_interval_pressure_width_jet import (
        IntervalPressureWidthJet,
    )


MINIMUM_PRECISION = 473
PRESSURE_ORDER = 9
WIDTH_ORDER = 2
RADIAL_DEGREE = 18
DEFAULT_ERROR_RECEIPT = Path(__file__).with_name(
    "lei_ren_part1_paper_coherent_uniform_pressure_high_derivatives.json"
)


def _decode_mpf(value: Any) -> mp.mpf:
    """Decode one exact MP value from the endpoint receipt."""

    if isinstance(value, Mapping):
        if "exact_mpf_tuple" in value:
            return mp.make_mpf(tuple(value["exact_mpf_tuple"]))
        if "value" in value and isinstance(value["value"], str):
            return mp.mpf(value["value"])
    if isinstance(value, (list, tuple)) and len(value) == 4:
        return mp.make_mpf(tuple(value))
    return mp.mpf(value)


def _read_error_receipt(path: Any) -> tuple[dict[str, Any], list[mp.mpf], Path, str]:
    receipt_path = Path(path).resolve()
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    derivatives = receipt.get("derivatives")
    if not isinstance(derivatives, list):
        raise ValueError("pressure error receipt lacks derivative envelopes")
    bounds: list[mp.mpf] = []
    for order, row in enumerate(derivatives):
        if not isinstance(row, Mapping):
            raise ValueError(f"invalid pressure error row {order}")
        item = row.get("normalized_Taylor_coefficient_error_upper")
        if item is None:
            item = row.get("normalized_derivative_error_upper")
        if item is None:
            raise ValueError(f"pressure error row {order} lacks normalized bound")
        bound = _decode_mpf(item)
        if bound < 0:
            raise ValueError(f"pressure error row {order} is negative")
        bounds.append(bound)
    digest = hashlib.sha256(receipt_path.read_bytes()).hexdigest()
    return dict(receipt), bounds, receipt_path, digest


def _pair(ctx: Any, nominal: Any, difference: Any = 0) -> IntervalDifference:
    return IntervalDifference(ctx, ctx.mpf(nominal), ctx.mpf(difference))


def _ring(
    ctx: Any,
    atoms: Mapping[tuple[int, int], IntervalDifference] | Any = 0,
) -> IntervalPressureWidthJet:
    if isinstance(atoms, Mapping):
        return IntervalPressureWidthJet._from_atoms(
            ctx, atoms, PRESSURE_ORDER, WIDTH_ORDER
        )
    return IntervalPressureWidthJet(
        ctx,
        atoms,
        pressure_order=PRESSURE_ORDER,
        width_order=WIDTH_ORDER,
    )


def _positive_bound(ctx: Any, normalized_bound: mp.mpf) -> Any:
    """Return the directed upper physical error ``exp(28) * B[k]``."""

    physical = ctx.exp(ctx.mpf(28)) * ctx.mpf(normalized_bound)
    _, upper = _endpoints(physical)
    return upper


def _polynomial_ring(
    ctx: Any,
    polynomial: PressurePolynomial,
    *,
    error_bound: Any = 0,
) -> IntervalPressureWidthJet:
    """Convert one stored pressure polynomial, retaining both P powers."""

    if not hasattr(polynomial, "atoms"):
        polynomial = PressurePolynomial(polynomial)
    atoms: dict[tuple[int, int], IntervalDifference] = {}
    error = ctx.mpf(error_bound)
    if hasattr(error_bound, "_mpi_"):
        _, error = _endpoints(error_bound)
        error = ctx.mpf(error)
    powers = {int(power) for power in polynomial.atoms}
    if error:
        # The combined bound applies independently to prefix and post-Rv
        # subsets.  A missing stored atom still receives its conservative
        # subset error around zero.
        powers.update((0, 1))
    for power in sorted(powers):
        coefficient = polynomial.atoms.get(power, mp.mpf(0))
        difference = (
            ctx.mpf([-error, error]) if error and power in (0, 1) else ctx.mpf(0)
        )
        atoms[(power, 0)] = _pair(ctx, coefficient, difference)
    return _ring(ctx, atoms)


def _converter(ctx: Any):
    def convert(value: Any) -> IntervalPressureWidthJet:
        if isinstance(value, IntervalPressureWidthJet):
            if value.ctx is not ctx:
                raise ValueError("interval pressure contexts must match")
            if value.orders != (PRESSURE_ORDER, WIDTH_ORDER):
                raise ValueError("interval pressure orders must match")
            return value
        if isinstance(value, PressurePolynomial) or hasattr(value, "atoms"):
            return _polynomial_ring(ctx, value)
        return _ring(ctx, value)

    return convert


def _source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def regenerate_core(
    cache_path: Any = DEFAULT_INPUT,
    *,
    error_receipt: Any = DEFAULT_ERROR_RECEIPT,
    precision: int = MINIMUM_PRECISION,
    radial_degree: int = RADIAL_DEGREE,
) -> dict[str, Any]:
    """Regenerate the finite radial core with a coherent interval pressure tail.

    ``core`` is directly consumable by
    ``build_inlet_from_interval_core(core, ctx, Lambda, delta)``.  The returned
    ``raw_source`` retains the exact stored ``PressurePolynomial`` rows for
    independent comparison and receipt generation.
    """

    precision = int(precision)
    radial_degree = int(radial_degree)
    if precision < MINIMUM_PRECISION:
        raise ValueError(f"precision must be at least {MINIMUM_PRECISION}")
    if radial_degree < 1:
        raise ValueError("radial_degree must be positive")

    raw_source, source_path, source_digest = _read_cache(cache_path)
    receipt, normalized_bounds, receipt_path, receipt_digest = _read_error_receipt(
        error_receipt
    )
    coefficients = raw_source["coefficients"]
    initial_rows = min(
        len(coefficients["F"][0]),
        len(coefficients["Uz"][0]),
        len(coefficients["P"][0]),
    )
    required = radial_degree + 2
    if initial_rows < required:
        raise ValueError(
            f"stored axis rows provide {initial_rows} terms; {required} required"
        )
    if len(normalized_bounds) < initial_rows:
        raise ValueError(
            f"pressure receipt provides {len(normalized_bounds)} bounds; "
            f"{initial_rows} required"
        )

    ctx = mp.ctx_iv.MPIntervalContext()
    ctx.dps = precision
    convert = _converter(ctx)

    # The error receipt is normalized to Pstar**2.  Apply its physical bound
    # independently to the stored prefix and post-Rv pressure powers.  This
    # is conservative because each subset is bounded by the combined total.
    pressure_axis = []
    physical_bounds = []
    with mp.workdps(precision + 40):
        for k, polynomial in enumerate(coefficients["P"][0]):
            bound = _positive_bound(ctx, normalized_bounds[k])
            physical_bounds.append(bound)
            pressure_axis.append(_polynomial_ring(ctx, polynomial, error_bound=bound))

        f_axis = [convert(item) for item in coefficients["F"][0]]
        uz_axis = [convert(item) for item in coefficients["Uz"][0]]
        z = ctx.mpf(str(raw_source["Z"]))
        delta = raw_source["delta"]
        core = core_coefficients(
            z,
            delta,
            F0_Z_taylor=f_axis,
            U0_Z_taylor=uz_axis,
            P0_Z_taylor=pressure_axis,
            radial_degree=radial_degree,
            precision=precision,
            scalar_converter=convert,
        )

    # Keep the original finite metadata while making the conditional error
    # provenance explicit.  This is metadata only; no source or ODE claim is
    # inferred from the interval recursion.
    core.update(
        interval_pressure_error_propagated=True,
        pressure_error_receipt=receipt_path.name,
        pressure_error_receipt_sha256=receipt_digest,
        pressure_error_normalization="combined normalized Pstar^2 bound applied separately to P powers 0 and 1",
        source_parameter_errors_enclosed=False,
        axis_generation_roundoff_enclosed=False,
        pressure_order_remainder_enclosed=False,
        radial_series_remainder_enclosed=False,
        source_remainder_enclosed=False,
        collar_RK_remainder_enclosed=False,
        temporal_recursion=False,
        accepted_source_path="accepted coherent-waiting schedule; stored exact finite inlet input",
    )
    return {
        "ctx": ctx,
        "raw_source": raw_source,
        "core": core,
        "regenerated_core": core,
        "receipt": receipt,
        "receipt_metadata": {
            "source_path": str(source_path),
            "source_sha256": source_digest,
            "error_receipt_path": str(receipt_path),
            "error_receipt_sha256": receipt_digest,
            "accepted_source_path": "accepted coherent-waiting schedule; stored exact finite inlet input",
            "precision": precision,
            "radial_degree": radial_degree,
            "pressure_order": PRESSURE_ORDER,
            "width_order": WIDTH_ORDER,
            "initial_axis_terms": initial_rows,
            "normalized_error_bound_count": len(normalized_bounds),
            "physical_error_bounds": physical_bounds,
            "pressure_powers_perturbed": [0, 1],
            "axis_perturbation": False,
            "source_parameter_errors_enclosed": False,
            "axis_generation_roundoff_enclosed": False,
            "pressure_order_remainder_enclosed": False,
            "radial_series_remainder_enclosed": False,
            "collar_RK_remainder_enclosed": False,
            "temporal_recursion": False,
        },
        # Keep these as the exact stored MP scalars: the parent consumer checks
        # identity against the accepted input before using the interval rows.
        "Lambda": raw_source["Lambda"],
        "delta": raw_source["delta"],
        "Z": raw_source["Z"],
    }


def _encode_interval(value: Any) -> dict[str, Any]:
    lower, upper = _endpoints(value)
    return {
        "lower": mp.nstr(lower, 60),
        "upper": mp.nstr(upper, 60),
        "lower_exact_mpf_tuple": list(value._mpi_[0]),
        "upper_exact_mpf_tuple": list(value._mpi_[1]),
    }


def encode_receipt_metadata(metadata: Mapping[str, Any]) -> dict[str, Any]:
    """Encode metadata while preserving exact interval bound endpoints."""

    result: dict[str, Any] = {}
    for key, value in metadata.items():
        if hasattr(value, "_mpi_"):
            result[key] = _encode_interval(value)
        elif isinstance(value, Mapping):
            result[key] = encode_receipt_metadata(value)
        elif isinstance(value, (list, tuple)):
            result[key] = [
                encode_receipt_metadata({'item':item})['item']
                for item in value
            ]
        elif isinstance(value, mp.mpf):
            result[key] = {
                "value": mp.nstr(value, 60),
                "exact_mpf_tuple": list(value._mpf_),
            }
        else:
            result[key] = value
    return result


__all__ = [
    "MINIMUM_PRECISION",
    "PRESSURE_ORDER",
    "WIDTH_ORDER",
    "RADIAL_DEGREE",
    "DEFAULT_INPUT",
    "DEFAULT_ERROR_RECEIPT",
    "regenerate_core",
    "encode_receipt_metadata",
]
