"""First-collar ``dx`` measure and signed source-cell adapter.

The first collar uses the physical coordinate

``x(s) = 2 exp(hb*s/2)``.

Consequently ``dx/ds = (hb/2)*x``.  For the direct coordinate ``q=s/sc``,
``dx/dq = (hb*sc/2)*x``.  The source code also uses a normalized coordinate
``q_source = 4*s/sc``; that chart has ``s=sc*q_source/4`` and
``dx/dq_source = (hb*sc/8)*x``.  The ``q_normalization`` argument keeps the
two conventions explicit.  This leaf keeps those identities in caller-owned
directed intervals and reweights signed source sectors cell by cell.  It is
deliberately a measure adapter: it does not evaluate a source producer, prove
that two source names are identical, or claim the legacy first-collar
identity.

There are two representations of the same physical measure in the result.
The default ``formal_weight`` is the analytic cell enclosure
``sc*x_hull*(q_right-q_left)/(2*q_normalization)`` and contains no division
by ``hb``.  ``physical_weight = hb*formal_weight`` contains the physical
factor, while the emitted formal sector increments its ``hb`` power by one.
Thus a closed physical ``hb`` interval (including zero as a geometric bound)
does not need to be materialized or inverted, and it is not also inserted as
a coefficient beside the formal atom.  The ownership record makes the
distinction visible to downstream callers.  A caller may opt out of formal
sectors with ``formal_hb_measure=False``.

The existing source giving the radial chart is
``collar_pressure_regular_singular_proof.py`` (``raw_rho = 4*exp(hb*s)``
near line 462 and the ``d rho/ds = hb*rho`` contract near line 626).  The
square-root chart and its ``dx`` identity are also recorded by
``original_first_collar_transport.py:3,137,275``.  Those references are
provenance, not imported runtime dependencies.
"""

from __future__ import annotations

import inspect
from collections.abc import Mapping, Sequence
from typing import Any, Callable


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
SCOPE = "first-collar signed source dx reweighting only"

__all__ = [
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "SCOPE",
    "x",
    "dx_ds",
    "dx_dq",
    "map_first_collar",
    "source_cell_weights",
    "ordinary_z1_chain",
    "build",
    "build_collar_green_measure",
    "build_g4",
    "integrate_source",
    "integrate",
    "green_measure",
    "adapt",
    "reweight",
]


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _interval(c: Any, value: Any, label: str = "interval") -> Any:
    """Coerce a caller value without midpoint or endpoint selection."""

    if _is_interval(value):
        return c.make_mpf((_raw_endpoint(value, False), _raw_endpoint(value, True)))
    if isinstance(value, Mapping):
        if "value" in value:
            return _interval(c, value["value"], label)
        if "interval" in value:
            return _interval(c, value["interval"], label)
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            return c.make_mpf(
                (
                    tuple(value["lower_exact_mpf_tuple"]),
                    tuple(value["upper_exact_mpf_tuple"]),
                )
            )
        if "lower" in value and "upper" in value:
            return c.mpf([value["lower"], value["upper"]])
        raise TypeError(f"{label} needs lower/upper endpoints")
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return c.mpf([value[0], value[1]])
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"unsupported {label}") from exc


def _lower(value: Any) -> Any:
    if _is_interval(value):
        return value.a
    return value


def _upper(value: Any) -> Any:
    if _is_interval(value):
        return value.b
    return value


def _hull(c: Any, *values: Any) -> Any:
    if not values:
        raise ValueError("at least one value is required for an interval hull")
    lowers = [_raw_endpoint(value, False) for value in values]
    uppers = [_raw_endpoint(value, True) for value in values]
    # Endpoint objects from MPIntervalContext are comparable scalar mpf
    # values.  Keeping the orientation separately is important; this hull is
    # only the range used to evaluate a monotone exponential.
    lo = min(lowers)
    hi = max(uppers)
    return c.make_mpf((lo, hi))


def _zero(c: Any) -> Any:
    return c.mpf(0)


def _contains_zero(value: Any) -> bool:
    return _lower(value) <= 0 <= _upper(value)


def _same_interval(c: Any, left: Any, right: Any) -> bool:
    left, right = _interval(c, left), _interval(c, right)
    return _lower(left) == _lower(right) and _upper(left) == _upper(right)


def _require_nonzero(value: Any, label: str) -> None:
    if _contains_zero(value):
        raise ValueError(f"{label} interval contains zero")


def _owner(value: Any, default: Any = None) -> Any:
    if isinstance(value, Mapping):
        return value.get("owner", default)
    return default


def _source_atoms(source_atoms: Any) -> None:
    if source_atoms is None:
        return
    if isinstance(source_atoms, Mapping):
        names = set(source_atoms)
    elif isinstance(source_atoms, Sequence) and not isinstance(source_atoms, (str, bytes, bytearray)):
        names = set()
        for item in source_atoms:
            if isinstance(item, Mapping):
                item_name = item.get("name", item.get("atom", item.get("power")))
                if item_name is not None:
                    names.add(item_name)
            else:
                names.add(item)
    else:
        raise TypeError("source_atoms must name the formal source atoms")
    if names != set(SOURCE_POWER_BASIS):
        raise ValueError("source_atoms must contain exactly hb, F0_anchor, Pstar_squared, Lambda")


def _z_independence_certificate(certificate: Any) -> dict[str, Any]:
    if not isinstance(certificate, Mapping):
        raise ValueError("caller ordinary-Z independence certificate is required")
    hb_sc = certificate.get("hb_sc_Z_independent") is True
    hb_sc = hb_sc or (
        certificate.get("hb_Z_independent") is True
        and certificate.get("sc_Z_independent") is True
    )
    if not hb_sc:
        raise ValueError("caller certificate must establish hb/sc Z independence")
    return {
        "hb_sc_Z_independent": True,
        "rho_radius_Z_independent": certificate.get("rho_radius_Z_independent") is True,
        "certificate": dict(certificate),
    }


def _check_geometry_owners(
    *,
    geometry_owner: Any,
    rho: Any,
    radius: Any,
    rho_owner: Any,
    radius_owner: Any,
) -> dict[str, Any]:
    resolved_rho = _owner(rho, rho_owner)
    resolved_radius = _owner(radius, radius_owner)
    if geometry_owner is not None:
        if resolved_rho is not None and resolved_rho != geometry_owner:
            raise ValueError("rho owner disagrees with geometry_owner")
        if resolved_radius is not None and resolved_radius != geometry_owner:
            raise ValueError("radius owner disagrees with geometry_owner")
        resolved_rho = resolved_rho if resolved_rho is not None else geometry_owner
        resolved_radius = resolved_radius if resolved_radius is not None else geometry_owner
    if rho is not None and radius is not None:
        if resolved_rho is None or resolved_radius is None:
            raise ValueError("rho and radius require one shared owner")
        if resolved_rho != resolved_radius:
            raise ValueError("rho and radius must have the same owner")
    return {
        "rho_owner": resolved_rho,
        "radius_owner": resolved_radius,
        "same_owner": bool(
            rho is not None
            and radius is not None
            and resolved_rho is not None
            and resolved_rho == resolved_radius
        ),
    }


def x(c: Any, hb: Any, s: Any) -> Any:
    """Return the directed interval enclosure of ``2*exp(hb*s/2)``."""

    hb, s = _interval(c, hb, "hb"), _interval(c, s, "s")
    return c.mpf(2) * c.exp(hb * s / 2)


def dx_ds(c: Any, hb: Any, s: Any, *, x_value: Any = None) -> Any:
    """Return ``(hb/2)*x(s)`` without choosing a point value."""

    hb = _interval(c, hb, "hb")
    value = x(c, hb, s) if x_value is None else _interval(c, x_value, "x(s)")
    return hb * value / 2


def _q_normalization(c: Any, value: Any) -> Any:
    result = _interval(c, value, "q_normalization")
    _require_nonzero(result, "q_normalization")
    return result


def dx_dq(
    c: Any,
    hb: Any,
    sc: Any,
    q: Any,
    *,
    q_normalization: Any = 1,
    x_value: Any = None,
) -> Any:
    """Return ``(hb*sc/(2*q_normalization))*x``.

    ``q_normalization=1`` is the direct ``q=s/sc`` chart.  The source
    normalized chart uses ``q_normalization=4``.
    """

    hb, sc, q = (
        _interval(c, value, label)
        for value, label in ((hb, "hb"), (sc, "sc"), (q, "q"))
    )
    q_normalization = _q_normalization(c, q_normalization)
    value = (
        x(c, hb, sc * q / q_normalization)
        if x_value is None
        else _interval(c, x_value, "x(s)")
    )
    return hb * sc * value / (2 * q_normalization)


def map_first_collar(
    c: Any,
    *,
    hb: Any,
    sc: Any,
    q: Any,
    q_normalization: Any = 1,
) -> dict[str, Any]:
    """Return the point/range map and both Jacobians for one ``q`` range."""

    hb, sc, q = (
        _interval(c, value, label)
        for value, label in ((hb, "hb"), (sc, "sc"), (q, "q"))
    )
    q_normalization = _q_normalization(c, q_normalization)
    s = sc * q / q_normalization
    xv = x(c, hb, s)
    return {
        "q": q,
        "s": s,
        "x": xv,
        "dx_ds": dx_ds(c, hb, s, x_value=xv),
        "dx_dq": dx_dq(c, hb, sc, q, q_normalization=q_normalization, x_value=xv),
        "q_normalization": q_normalization,
        "mapping": "s=sc*q/q_normalization; x=2*exp(hb*s/2)",
        "jacobian_identities": {
            "dx_ds": "(hb/2)*x",
            "dx_dq": "(hb*sc/(2*q_normalization))*x",
        },
    }


def _partition(c: Any, count: Any) -> list[tuple[Any, Any]]:
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        raise ValueError("cell_count must be a positive integer")
    return [(c.mpf(index) / count, c.mpf(index + 1) / count) for index in range(count)]


def _cell_bounds(c: Any, raw: Any, index: int) -> tuple[Any, Any, Any]:
    """Normalize one cell while retaining directed left/right endpoints."""

    if isinstance(raw, Mapping):
        nested = raw.get("q_interval", raw.get("interval"))
        left = raw.get("q_left", raw.get("left"))
        right = raw.get("q_right", raw.get("right"))
        if nested is not None and (left is None or right is None):
            if not isinstance(nested, Sequence) or len(nested) != 2:
                raise TypeError(f"cells[{index}].q_interval must have two endpoints")
            left, right = nested
        if left is None or right is None:
            raise ValueError(f"cells[{index}] needs q_left and q_right")
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)) and len(raw) == 2:
        left, right = raw
    else:
        raise TypeError(f"cells[{index}] must be an endpoint pair or mapping")
    left_iv = _interval(c, left, f"cells[{index}].q_left")
    right_iv = _interval(c, right, f"cells[{index}].q_right")
    q_hull = _hull(c, left_iv, right_iv)
    return left_iv, right_iv, q_hull


def _normalize_cells(c: Any, cells: Any, cell_count: int) -> list[dict[str, Any]]:
    if cells is None:
        raw_cells: list[Any] = _partition(c, cell_count)
    elif isinstance(cells, int):
        raw_cells = _partition(c, cells)
    elif isinstance(cells, Mapping) and "cells" in cells:
        raw_cells = list(cells["cells"])
    elif isinstance(cells, Sequence) and not isinstance(cells, (str, bytes, bytearray)):
        raw_cells = list(cells)
    else:
        raise TypeError("cells must be a positive count or a sequence of cells")
    if not raw_cells:
        raise ValueError("at least one source cell is required")

    result: list[dict[str, Any]] = []
    for index, raw in enumerate(raw_cells):
        left, right, q_hull = _cell_bounds(c, raw, index)
        result.append(
            {
                "index": index,
                "q_left": left,
                "q_right": right,
                "q_interval": q_hull,
                "raw": raw,
            }
        )
    return result


def _cell_geometry(
    c: Any,
    hb: Any,
    sc: Any,
    cell: Mapping[str, Any],
    q_normalization: Any = 1,
) -> dict[str, Any]:
    q_normalization = _q_normalization(c, q_normalization)
    q_left, q_right, q_hull = cell["q_left"], cell["q_right"], cell["q_interval"]
    s_left, s_right = sc * q_left / q_normalization, sc * q_right / q_normalization
    s_hull = _hull(c, s_left, s_right)
    x_left, x_right = x(c, hb, s_left), x(c, hb, s_right)
    x_hull = x(c, hb, s_hull)
    q_span = q_right - q_left
    # The hull encloses the positive/negative x values on this directed cell.
    # Multiplying by the directed q span preserves orientation and avoids
    # dividing by a physical hb enclosure that may contain zero.
    formal_weight = sc * x_hull * q_span / (2 * q_normalization)
    physical_weight = hb * formal_weight
    endpoint_weight = x_right - x_left
    return {
        "index": cell["index"],
        "q_left": q_left,
        "q_right": q_right,
        "q_interval": q_hull,
        "s_left": s_left,
        "s_right": s_right,
        "s_interval": s_hull,
        "x_left": x_left,
        "x_right": x_right,
        "x_interval": x_hull,
        "dx_ds": hb * x_hull / 2,
        "dx_dq": hb * sc * x_hull / (2 * q_normalization),
        "q_normalization": q_normalization,
        "q_span": q_span,
        "physical_weight": physical_weight,
        "formal_weight": formal_weight,
        "endpoint_weight": endpoint_weight,
        "weight": physical_weight,
        "measure_formula": "physical_weight=hb*sc*x_hull*(q_right-q_left)/(2*q_normalization)",
        "formal_measure_formula": "sc*x_hull*(q_right-q_left)/(2*q_normalization)",
        "weight_is_enclosure": True,
        "endpoint_weight_is_diagnostic_only": True,
        "exact_endpoint_weight_claimed": False,
        "orientation_preserved": True,
    }


def source_cell_weights(
    c: Any,
    *,
    sc: Any,
    hb: Any,
    cells: Any = None,
    cell_count: int = 1,
    q_normalization: Any = 1,
) -> list[dict[str, Any]]:
    """Return directed ``dx`` weights for every source cell.

    ``sc`` may be negative.  Cell left/right orientation is kept in the
    directed q span; no absolute width or sorted endpoint is substituted.
    ``endpoint_weight`` is retained only as a diagnostic; the certified
    output weight is the analytic enclosure built from ``x_interval``.
    """

    hb, sc = _interval(c, hb, "hb"), _interval(c, sc, "sc")
    _require_nonzero(sc, "sc")
    normalized = _normalize_cells(c, cells, cell_count)
    return [_cell_geometry(c, hb, sc, cell, q_normalization) for cell in normalized]


def _powers(value: Any, label: str) -> tuple[int, int, int, int]:
    powers = value.get("powers") if isinstance(value, Mapping) else None
    if powers is None and isinstance(value, Mapping):
        powers = (
            value.get("hb_power", 0),
            value.get("F0_power", 0),
            value.get("Pstar_squared_power", value.get("Pstar_power", 0)),
            value.get("Lambda_power", 0),
        )
    if not isinstance(powers, (list, tuple)) or len(powers) != 4:
        raise ValueError(f"{label} needs four source powers")
    if any(not isinstance(power, int) or isinstance(power, bool) for power in powers):
        raise TypeError(f"{label} source powers must be integers")
    return tuple(int(power) for power in powers)


def _term(c: Any, value: Any, label: str, *, default_origin: Any) -> dict[str, Any]:
    if isinstance(value, Mapping) and (
        "powers" in value
        or "coefficient" in value
        or "coefficient_interval" in value
    ):
        powers = _powers(value, label)
        coefficient_raw = value.get("coefficient", value.get("coefficient_interval"))
        if coefficient_raw is None:
            raise ValueError(f"{label} coefficient is missing")
        coefficient = _interval(c, coefficient_raw, f"{label}.coefficient")
        origin = value.get("origin", default_origin)
        result = dict(value)
        result.update({"powers": powers, "coefficient": coefficient, "origin": origin})
        return result
    return {
        "powers": (0, 0, 0, 0),
        "coefficient": _interval(c, value, f"{label}.coefficient"),
        "origin": default_origin,
    }


def _terms(c: Any, value: Any, label: str) -> list[dict[str, Any]]:
    if value is None:
        return []
    if isinstance(value, Mapping) and (
        "powers" in value
        or "coefficient" in value
        or "coefficient_interval" in value
    ):
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        values = [value]
    return [
        _term(c, item, f"{label}[{index}]", default_origin=("caller_source", label, index))
        for index, item in enumerate(values)
    ]


def _rows(c: Any, value: Any, label: str) -> dict[int, list[dict[str, Any]]]:
    """Read source rows ``Z0`` and ``Z1`` without manufacturing a derivative."""

    if value is None:
        return {0: [], 1: []}
    if isinstance(value, Mapping) and not (
        "powers" in value
        or "coefficient" in value
        or "coefficient_interval" in value
    ):
        row_values: dict[int, Any] = {}
        for key, item in value.items():
            if key in (0, "0", "Z0", "z0", "value"):
                order = 0
            elif key in (1, "1", "Z1", "z1", "derivative", "ordinary_Z1"):
                order = 1
            else:
                continue
            if order in row_values:
                raise ValueError(f"duplicate ordinary-Z row in {label}")
            row_values[order] = item
        if row_values:
            return {order: _terms(c, row_values.get(order), f"{label}.Z{order}") for order in (0, 1)}
    # A plain value is a Z0 source.  Z1 must be supplied separately by the
    # caller; returning an empty row prevents a fabricated derivative.
    return {0: _terms(c, value, f"{label}.Z0"), 1: []}


def _source_for_cell(c: Any, cell: Mapping[str, Any], source_value: Any = None) -> dict[int, list[dict[str, Any]]]:
    raw = source_value
    if raw is None:
        candidate = cell.get("raw")
        if isinstance(candidate, Mapping):
            if "source" in candidate:
                raw = candidate["source"]
            elif "source_Z" in candidate:
                raw = candidate["source_Z"]
            elif any(key in candidate for key in (0, 1, "0", "1", "Z0", "Z1", "z0", "z1")):
                raw = candidate
            elif "terms" in candidate:
                raw = candidate["terms"]
            # Explicit aliases are accepted when a cell has a separate Z1.
            if raw is not candidate and "Z1" in candidate:
                raw = {"Z0": raw, "Z1": candidate["Z1"]}
            if raw is None and "z0" in candidate:
                raw = {"Z0": candidate["z0"], "Z1": candidate.get("z1")}
    return _rows(c, raw, f"source_cell[{cell['index']}]")


def _call_source(callback: Callable[..., Any], cell: Mapping[str, Any]) -> Any:
    """Call a source callback using the documented one/two argument forms."""

    try:
        parameters = list(inspect.signature(callback).parameters.values())
    except (TypeError, ValueError):
        parameters = []
    positional = [
        parameter
        for parameter in parameters
        if parameter.kind
        in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD)
    ]
    if len(positional) >= 2:
        return callback(cell["q_interval"], cell)
    if len(positional) == 1:
        return callback(cell)
    return callback()


def _weighted_terms(
    c: Any,
    terms: Sequence[Mapping[str, Any]],
    cell: Mapping[str, Any],
    *,
    formal_hb_measure: bool,
    source_row: int,
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    factor = cell["formal_weight"] if formal_hb_measure else cell["physical_weight"]
    for index, raw in enumerate(terms):
        term = dict(raw)
        origin_value = term.get("origin")
        already_measured = (
            term.get("formal_hb_measure_attached") is True
            or term.get("measure_hb_attached") is True
            or term.get("hb_measure_attached") is True
            or (
                isinstance(origin_value, Mapping)
                and origin_value.get("formal_hb_measure_attached") is True
            )
        )
        if formal_hb_measure and already_measured:
            raise ValueError("source term already carries a formal hb measure")
        if formal_hb_measure and (
            term.get("physical_hb_coefficient") is True
            or term.get("hb_is_physical_coefficient") is True
            or term.get("coefficient_contains_physical_hb") is True
        ):
            raise ValueError(
                "source term declares a physical hb coefficient while formal hb measure is enabled"
            )
        powers = list(_powers(term, f"source term {cell['index']}:{source_row}:{index}"))
        if formal_hb_measure:
            powers[0] += 1
        origin = term.get("origin", ("caller_source", cell["index"], source_row, index))
        origin_metadata = {
            "operation": "directed_dx_source_cell_reweight",
            "cell_index": cell["index"],
            "ordinary_Z_order": source_row,
            "source_term_index": index,
            "formal_hb_measure_attached": bool(formal_hb_measure),
            "formal_hb_power_added": 1 if formal_hb_measure else 0,
            "physical_hb_used_only_in_geometry": True,
            "physical_hb_in_coefficient": False,
            "physical_weight_is_signed": True,
            "weight_is_enclosure": True,
            "exact_endpoint_weight_claimed": False,
            "formal_weight_formula": "sc*x_hull*(q_right-q_left)/(2*q_normalization)",
            "origin": origin,
        }
        result.append(
            {
                "powers": tuple(powers),
                "coefficient": term["coefficient"] * factor,
                "origin": origin_metadata,
                "source_cell": cell["index"],
                "ordinary_Z_order": source_row,
                "input_origin": origin,
            }
        )
    return result


def ordinary_z1_chain(
    c: Any,
    *,
    cell: Mapping[str, Any],
    rows: Mapping[int, Sequence[Mapping[str, Any]]],
    formal_hb_measure: bool,
) -> dict[str, Any]:
    """Build the ordinary-Z1 row while keeping the measure derivative zero.

    This is a chain record, not a proof of source identity.  The caller's
    ``hb/sc`` independence certificate is retained by :func:`build`.
    """

    zero_rows = rows.get(0, rows.get("Z0", rows.get("z0", ())))
    one_rows = rows.get(1, rows.get("Z1", rows.get("z1", ())))
    weighted_zero = _weighted_terms(
        c, zero_rows, cell, formal_hb_measure=formal_hb_measure, source_row=0
    )
    weighted_one = _weighted_terms(
        c, one_rows, cell, formal_hb_measure=formal_hb_measure, source_row=1
    )
    return {
        "Z0": weighted_zero,
        "Z1": weighted_one,
        "weight_Z1": _zero(c),
        "formula": "d_Z(source*weight)=source_Z1*weight+source_Z0*weight_Z1",
        "weight_Z1_is_zero": True,
        "geometry_not_differentiated": True,
        "ordinary_Z_semantics": "fixed directed q/sc/hb geometry; caller-supplied source Z1",
    }


def build(
    c: Any,
    *,
    sc: Any,
    hb: Any = None,
    physical_hb: Any = None,
    epsilon: Any = None,
    cells: Any = None,
    cell_count: int = 1,
    source_cells: Any = None,
    source: Callable[..., Any] | None = None,
    source_Z: Any = None,
    g4: Any = None,
    forcing: Any = None,
    source_family: Any = None,
    source_atoms: Any = None,
    geometry_certificate: Any = None,
    certificate: Any = None,
    formal_hb_measure: bool = True,
    formal_hb_owner: Any = "green_measure_adapter",
    physical_hb_owner: Any = "caller_geometry",
    epsilon_owner: Any = None,
    geometry_owner: Any = None,
    rho: Any = None,
    radius: Any = None,
    rho_owner: Any = None,
    radius_owner: Any = None,
    include_z1: bool = True,
    q_normalization: Any = None,
    scope: Any = SCOPE,
) -> dict[str, Any]:
    """Build directed first-collar map weights and optional signed source rows.

    ``source`` is called as ``source(q_interval, cell_record)`` when it has
    two positional parameters, or as ``source(cell_record)`` when it has one.
    It returns either a Z0 value or a mapping with ``Z0`` and ``Z1`` rows.
    A cell mapping may instead contain ``source``/``source_Z``/``Z0``/``Z1``
    directly.  Terms use the usual four formal powers
    ``(hb,F0_anchor,Pstar_squared,Lambda)``.

    The default formal mode attaches one formal ``hb`` power and multiplies
    the coefficient by the analytic ``formal_weight`` enclosure.  This is the
    single explicit measure factor and remains valid when the physical ``hb``
    interval contains zero.  In physical mode (``formal_hb_measure=False``),
    the full signed ``physical_weight`` enclosure is put in the coefficient
    and no formal power is added.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf") or not hasattr(c, "exp"):
        raise TypeError("c must be the caller's interval context")
    if hb is not None and physical_hb is not None and not _same_interval(c, hb, physical_hb):
        raise ValueError("hb and physical_hb disagree")
    if hb is None:
        hb = physical_hb
    if hb is None:
        raise TypeError("physical hb interval is required for the chart")
    hb, sc = _interval(c, hb, "hb"), _interval(c, sc, "sc")
    _require_nonzero(sc, "sc")
    if not isinstance(formal_hb_measure, bool):
        raise TypeError("formal_hb_measure must be boolean")
    if formal_hb_measure and formal_hb_owner is None:
        raise ValueError("formal_hb_owner is required when formal hb measure is enabled")
    if epsilon is not None:
        epsilon = _interval(c, epsilon, "epsilon")
        if epsilon_owner is None:
            raise ValueError("epsilon_owner is required for a physical epsilon")
    _source_atoms(source_atoms)
    if geometry_certificate is not None and certificate is not None and geometry_certificate != certificate:
        raise ValueError("geometry_certificate and certificate disagree")
    if geometry_certificate is None:
        geometry_certificate = certificate

    aliases = [value is not None for value in (g4, forcing)]
    if sum(aliases) > 1:
        raise ValueError("g4 and forcing aliases are mutually exclusive")
    if g4 is not None:
        forcing = g4
    if forcing is not None:
        if any(value is not None for value in (source_cells, source, source_Z)):
            raise ValueError("forcing/g4 cannot be combined with source inputs")
        # A row map is one source for every cell; a callback or cell-index map
        # remains a per-cell source provider.
        if isinstance(forcing, Mapping) and any(
            key in forcing for key in (0, 1, "0", "1", "Z0", "Z1", "z0", "z1")
        ):
            source_Z = forcing
        elif (
            g4 is not None
            and isinstance(forcing, Sequence)
            and not isinstance(forcing, (str, bytes, bytearray))
            and len(forcing) == 2
            and cells is None
            and cell_count == 1
        ):
            source_Z = {"Z0": forcing[0], "Z1": forcing[1]}
        else:
            source_cells = forcing
    if source is not None and source_cells is not None:
        raise ValueError("source and source_cells are mutually exclusive")
    if source is not None and not callable(source):
        raise TypeError("source must be callable")
    if source_family is not None and isinstance(source_family, str) and not source_family.strip():
        raise ValueError("source_family must not be empty")
    if g4 is not None and source_family is None:
        raise ValueError("same-source g4 forcing requires source_family")
    if g4 is not None and source_atoms is None:
        raise ValueError("same-source g4 forcing requires source_atoms")

    if q_normalization is None:
        # The native g4/source chart uses q_source=4*s/sc.  Generic callers
        # retain the direct q=s/sc convention unless they opt into q=4s/sc.
        q_normalization = 4 if g4 is not None or forcing is not None else 1
    q_normalization = _q_normalization(c, q_normalization)
    owner_info = _check_geometry_owners(
        geometry_owner=geometry_owner,
        rho=rho,
        radius=radius,
        rho_owner=rho_owner,
        radius_owner=radius_owner,
    )
    rho_iv = None if rho is None else _interval(c, rho, "rho")
    radius_iv = None if radius is None else _interval(c, radius, "radius")

    raw_cells = _normalize_cells(c, cells, cell_count)
    weights = [_cell_geometry(c, hb, sc, cell, q_normalization) for cell in raw_cells]
    has_sources = source is not None or source_cells is not None or source_Z is not None or any(
        isinstance(cell.get("raw"), Mapping)
        and any(
            key in cell["raw"]
            for key in ("source", "source_Z", "terms", "Z0", "Z1", "z0", "z1")
        )
        for cell in raw_cells
    )
    if has_sources and source_family is None:
        raise ValueError("source integration requires a same-source family label")
    z_certificate = None
    if include_z1 and has_sources:
        z_certificate = _z_independence_certificate(geometry_certificate)

    callback = source
    if source_cells is not None and callable(source_cells):
        callback = source_cells
        source_cells = None
    if callback is None and callable(source_Z):
        callback = source_Z
        source_Z = None

    weighted_cells: list[dict[str, Any]] = []
    integral_rows: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    z1_rows: list[dict[str, Any]] = []
    for index, (cell, geometry) in enumerate(zip(raw_cells, weights)):
        cell_with_geometry = dict(cell)
        cell_with_geometry.update(geometry)
        supplied = None
        if callback is not None:
            supplied = _call_source(callback, cell_with_geometry)
        elif source_cells is not None:
            if isinstance(source_cells, Mapping):
                supplied = source_cells.get(index, source_cells.get(str(index)))
            elif isinstance(source_cells, Sequence) and not isinstance(source_cells, (str, bytes, bytearray)):
                if index >= len(source_cells):
                    raise ValueError("source_cells does not cover every geometry cell")
                supplied = source_cells[index]
            else:
                raise TypeError("source_cells must be a callback, mapping, or sequence")
        elif source_Z is not None:
            supplied = source_Z

        rows = _source_for_cell(c, cell_with_geometry, supplied)
        if supplied is None and source_Z is not None:
            rows = _rows(c, source_Z, f"source_Z[{index}]")
        if include_z1 and has_sources and not rows[1]:
            raise ValueError(f"source cell {index} must supply an ordinary Z1 row")
        chain = ordinary_z1_chain(
            c,
            cell=geometry,
            rows=rows,
            formal_hb_measure=formal_hb_measure,
        )
        integral_rows[0].extend(chain["Z0"])
        integral_rows[1].extend(chain["Z1"])
        z1_rows.append(chain)
        weighted_cells.append(
            {
                "index": index,
                "geometry": geometry,
                "source_rows": {"Z0": chain["Z0"], "Z1": chain["Z1"]},
                "ordinary_Z1_chain": chain,
            }
        )

    if not include_z1:
        z_certificate = None
        integral_rows[1] = []
    return {
        "source_family": source_family,
        "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "sc": sc,
        "hb": hb,
        "q_normalization": q_normalization,
        "q_definition": "q=s/sc"
        if _lower(q_normalization) == 1 and _upper(q_normalization) == 1
        else "q=4*s/sc",
        "s_from_q": "s=sc*q"
        if _lower(q_normalization) == 1 and _upper(q_normalization) == 1
        else "s=sc*q/4",
        "dx_dq_formula": "(hb*sc/2)*x"
        if _lower(q_normalization) == 1 and _upper(q_normalization) == 1
        else "(hb*sc/8)*x",
        "epsilon": epsilon,
        "rho": rho_iv,
        "radius": radius_iv,
        "formal_hb_measure": bool(formal_hb_measure),
        "source_cell_weights": weights,
        "cell_weights": weights,
        "weights": weights,
        "cells": weighted_cells,
        "source_integral": {"Z0": integral_rows[0], "Z1": integral_rows[1]},
        "forcing": {"Z0": integral_rows[0], "Z1": integral_rows[1]},
        "g4": {"Z0": integral_rows[0], "Z1": integral_rows[1]} if g4 is not None else None,
        "ordinary_Z1": z1_rows,
        "ordinary_Z1_chain": z1_rows,
        "Z1_chain": z1_rows,
        "physical_epsilon_owner": epsilon_owner,
        "physical_hb_owner": physical_hb_owner,
        "formal_hb_owner": formal_hb_owner,
        "ownership": {
            "physical_epsilon_owner": epsilon_owner,
            "physical_hb_owner": physical_hb_owner,
            "formal_hb_owner": formal_hb_owner,
            "physical_hb_used_for_chart": True,
            "physical_hb_inserted_as_coefficient": not formal_hb_measure,
            "formal_hb_attached_once": bool(formal_hb_measure),
            "epsilon_is_not_automatically_multiplied": True,
            "double_count_guard": "physical hb is chart geometry; formal hb is one measure power",
            "geometry_owner": geometry_owner,
            **owner_info,
        },
        "geometry_certificate": geometry_certificate,
        "ordinary_Z1_certificate": z_certificate,
        "geometry_certificate_independently_verified": False,
        "source_identity_proved": False,
        "same_source_identity_proved": False,
        "m15_computed": False,
        "complete_m15_computed": False,
        "scope": scope,
        "legacy_identity": False,
        "legacy_identity_claimed": False,
        "scope_exclusions": (
            "source identity",
            "legacy first-collar equality",
            "tail/history closure",
            "epsilon multiplication or physical source selection",
        ),
        "mapping_reference": (
            "src/openai_ns_reconstruction/collar_pressure_regular_singular_proof.py:462",
            "src/openai_ns_reconstruction/collar_pressure_regular_singular_proof.py:626",
            "src/openai_ns_reconstruction/original_first_collar_transport.py:3",
            "src/openai_ns_reconstruction/original_first_collar_transport.py:137",
            "src/openai_ns_reconstruction/original_first_collar_transport.py:275",
        ),
        "source_terms_merged": False,
        "single_measure_attachment": True,
        "weight_semantics": "directed analytic x-hull enclosure, not exact endpoint subtraction",
        "signed_cell_orientation_preserved": True,
        "point_value_selected": False,
    }


build_collar_green_measure = build
build_g4 = build
integrate_source = build
integrate = build
green_measure = build
adapt = build
reweight = build
