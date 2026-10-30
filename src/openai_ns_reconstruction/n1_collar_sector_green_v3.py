"""Standalone first-collar first-potential correction by signed sectors.

This leaf is intentionally *not* the previously saved bare ``K0`` Green
response.  It computes only the new finite correction

    ``Dg = (-2*chi) K0 o K0 g5``

at the physical endpoint ``q*=1`` of the first collar.  The source is supplied
as whole-cell, signed ``g5`` sectors through ordinary ``Z0`` and ``Z1``.  No
ancestor source producer, receipt, baseline Green response, all-word tail, or
physical m15 assembly is consulted here.

The collar chart is

    ``x(q) = 2 exp(hb*sc*q/8)`` and ``dx = hb*sc*x(q)/8 dq``.

Writing ``c=sc/8`` and ``phi(z)=(1-exp(-z))/z``, the two positive kernel
factors used by the nested integral are

    ``K0(x(q),x(r)) dx/dr = hb**2*c**2*x(r)**2*(q-r)*phi(2*hb*c*(q-r))``

and

    ``d_x K0(x(q),x(r)) dx/dr = hb*c*x(r)*exp(-3*hb*c*(q-r))``.

The returned mass intervals factor the displayed ``hb`` powers out of their
coefficients.  Thus every physical ``dx`` is owned exactly once, while each
vanishing K0 factor remains a formal ``hb`` atom: the F1 correction adds four
``hb`` powers and the radial derivative correction adds three.  Geometry uses
the caller's directed ``hb`` interval; no midpoint or numerical ``hb`` value
is selected.

The result is a conditional finite source-sector enclosure, not a claim that
the complete physical field has been evaluated.  The caller must provide an
explicit ``geometry_certificate`` proving both
``hb_sc_Z_independent=True`` and ``rho_radius_Z_independent=True`` (or the
equivalent paired flags) before fixed-geometry ordinary-Z1 transport is
accepted.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
FIELDS = ("F1", "d_x_F1")
SOURCE_ROW = "5"
CELL_COUNT = 4

__all__ = [
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "FIELDS",
    "SOURCE_ROW",
    "CELL_COUNT",
    "x",
    "dx_dq",
    "phi_nonnegative",
    "endpoint_kernel_masses",
    "kernel_masses",
    "physical_kernel_masses",
    "build",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _raw_endpoint(value: Any, upper: bool) -> Any:
    """Return the exact low-level endpoint for interval or scalar mpf input."""

    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    # ``_bounds`` returns ordinary mpmath mpf scalars.  ``make_mpf`` expects
    # their low-level ``_mpf_`` tuples, not the scalar object itself.
    raw = getattr(value, "_mpf_", None)
    if raw is not None:
        return raw
    return value


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _interval(c: Any, value: Any, label: str = "interval") -> Any:
    """Coerce a caller value without selecting a midpoint."""

    if _is_interval(value):
        return c.make_mpf((_raw_endpoint(value, False), _raw_endpoint(value, True)))
    if isinstance(value, Mapping):
        for key in ("value", "interval", "box", "enclosure"):
            if key in value:
                return _interval(c, value[key], label)
        if "lower" in value and "upper" in value:
            return c.mpf([value["lower"], value["upper"]])
        raise TypeError(f"{label} requires an interval-compatible value")
    if _is_sequence(value) and len(value) == 2:
        return c.mpf([value[0], value[1]])
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{label} is not interval-compatible") from exc


def _bounds(value: Any) -> tuple[Any, Any]:
    raw = getattr(value, "_mpi_", None)
    if raw is None:
        return value, value
    # The conversion is used only for ordering/validation.  All arithmetic
    # stays in the caller's directed interval context.
    from mpmath import mp

    return mp.make_mpf(raw[0]), mp.make_mpf(raw[1])


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _point(c: Any, value: Any) -> Any:
    value = _interval(c, value)
    raw = value._mpi_[0]
    return c.make_mpf((raw, raw))


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _interval(c, value, label)
    lo, hi = _bounds(value)
    if hi < 0:
        raise ValueError(f"{label} is negative")
    if lo < 0:
        lo = type(hi)(0)
    return _make_interval(c, lo, hi)


def _positive(c: Any, value: Any, label: str) -> None:
    lo, _ = _bounds(_interval(c, value, label))
    if lo <= 0:
        raise ValueError(f"{label} must be strictly positive")


def _geometry_evidence(value: Any) -> dict[str, Any]:
    """Require explicit ordinary-Z independence for the physical chart."""

    if not isinstance(value, Mapping):
        raise ValueError(
            "geometry_certificate with hb_sc_Z_independent and "
            "rho_radius_Z_independent is required"
        )
    evidence = value.get("geometry", value)
    if not isinstance(evidence, Mapping):
        raise ValueError("geometry_certificate.geometry must be a mapping")
    hb_sc = evidence.get("hb_sc_Z_independent") is True
    hb_sc = hb_sc or (
        evidence.get("hb_Z_independent") is True
        and evidence.get("sc_Z_independent") is True
    )
    rho_radius = evidence.get("rho_radius_Z_independent") is True
    rho_radius = rho_radius or (
        evidence.get("rho_Z_independent") is True
        and evidence.get("radius_Z_independent") is True
    )
    if not hb_sc or not rho_radius:
        missing = []
        if not hb_sc:
            missing.append("hb_sc_Z_independent")
        if not rho_radius:
            missing.append("rho_radius_Z_independent")
        raise ValueError("geometry_certificate lacks explicit true flags: " + ", ".join(missing))
    return {
        "hb_sc_Z_independent": True,
        "rho_radius_Z_independent": True,
        "raw": dict(value),
    }


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    for value in values:
        if not any(value == previous for previous in result):
            result.append(value)
    return result


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    if _is_sequence(value):
        return list(value)
    return [value]


def _origin_info(value: Any) -> dict[str, Any]:
    """Collect source provenance without changing the caller's metadata."""

    if not isinstance(value, Mapping):
        return {"source_ids": [], "error_classes": [], "uncertainty": False}
    ids: list[Any] = []
    for key in ("source_ids", "sources", "source_id", "source", "input_id"):
        if key in value:
            ids.extend(_items(value[key]))
    classes: list[Any] = []
    for key in ("error_classes", "error_class", "uncertainty_class"):
        if key in value:
            classes.extend(_items(value[key]))
    uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
    operation = str(value.get("operation", "")).lower()
    uncertain = uncertain or any(token in operation for token in ("cauchy", "symmetric", "uncertainty"))
    for key in ("origin", "metadata", "parents"):
        child = value.get(key)
        if child is None or child is value:
            continue
        children = child if _is_sequence(child) else [child]
        for nested in children:
            info = _origin_info(nested)
            ids.extend(info["source_ids"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or info["uncertainty"]
    return {
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": bool(uncertain),
    }


def _powers(raw: Mapping[str, Any], label: str) -> tuple[int, int, int, int]:
    value = raw.get("powers")
    if value is None:
        value = (
            raw.get("hb_power", 0),
            raw.get("F0_power", 0),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)),
            raw.get("Lambda_power", 0),
        )
    if not _is_sequence(value) or len(value) != 4:
        raise ValueError(f"{label}.powers must contain four integer formal powers")
    if any(not isinstance(item, int) or isinstance(item, bool) for item in value):
        raise ValueError(f"{label}.powers must contain four integer formal powers")
    return tuple(int(item) for item in value)


def _term(c: Any, raw: Any, *, label: str, default_source_id: Any,
          require_source_id: bool = True) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        if require_source_id:
            raise TypeError(f"{label} must be a sector mapping with source IDs")
        raw = {"powers": (0, 0, 0, 0), "coefficient": raw,
               "origin": {"source_ids": [default_source_id]}}
    coefficient_raw = raw.get("coefficient", raw.get("coefficient_interval"))
    if coefficient_raw is None:
        raise ValueError(f"{label}.coefficient is required")
    powers = _powers(raw, label)
    origin = dict(raw.get("origin", {})) if isinstance(raw.get("origin"), Mapping) else {}
    info = _origin_info(raw)
    ids = list(info["source_ids"])
    if not ids and default_source_id is not None and not require_source_id:
        ids = [default_source_id]
    if not ids:
        raise ValueError(f"{label} must retain at least one source ID")
    origin["source_ids"] = _unique(ids)
    origin.setdefault("error_classes", list(info["error_classes"]))
    origin.setdefault("uncertainty", bool(info["uncertainty"]))
    result = dict(raw)
    result.update({
        "powers": powers,
        "coefficient": _interval(c, coefficient_raw, f"{label}.coefficient"),
        "origin": origin,
        "source_ids": list(origin["source_ids"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(info["uncertainty"]),
    })
    return result


def _term_list(c: Any, raw: Any, *, label: str, default_source_id: Any,
               require_source_id: bool = True) -> list[dict[str, Any]]:
    if raw is None:
        return []
    if isinstance(raw, Mapping) and (
        "coefficient" in raw or "coefficient_interval" in raw or "powers" in raw
    ):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        values = [raw]
    return [
        _term(c, value, label=f"{label}[{index}]", default_source_id=default_source_id,
              require_source_id=require_source_id)
        for index, value in enumerate(values)
    ]


_ORDER_ALIASES = {
    0: (0, "0", "Z0", "z0", "value"),
    1: (1, "1", "Z1", "z1", "derivative", "ordinary_Z1"),
}


def _lookup_order(mapping: Mapping[Any, Any], order: int) -> tuple[bool, Any]:
    for key in _ORDER_ALIASES[order]:
        if key in mapping:
            return True, mapping[key]
    return False, None


def _rows(c: Any, raw: Any, *, label: str, default_source_id: Any,
          require_source_id: bool = True) -> dict[int, list[dict[str, Any]]]:
    if isinstance(raw, Mapping):
        for wrapper in ("ordinary_Z", "source_Z", "source_jets", "jets", "rows", "g5"):
            if wrapper in raw:
                raw = raw[wrapper]
                break
    if isinstance(raw, Mapping):
        found0, z0 = _lookup_order(raw, 0)
        found1, z1 = _lookup_order(raw, 1)
    elif _is_sequence(raw) and len(raw) == 2 and not all(
        isinstance(item, Mapping) and ("powers" in item or "coefficient" in item)
        for item in raw
    ):
        found0, z0 = True, raw[0]
        found1, z1 = True, raw[1]
    else:
        found0, z0 = True, raw
        found1, z1 = False, None
    if not found0 or not found1:
        raise ValueError(f"{label} must provide explicit ordinary Z0 and Z1 rows")
    return {
        0: _term_list(c, z0, label=f"{label}.Z0", default_source_id=f"{default_source_id}:Z0",
                      require_source_id=require_source_id),
        1: _term_list(c, z1, label=f"{label}.Z1", default_source_id=f"{default_source_id}:Z1",
                      require_source_id=require_source_id),
    }


def _looks_like_cell(value: Any) -> bool:
    if not isinstance(value, Mapping):
        return False
    if any(key in value for key in ("q_left", "q_right", "q_interval", "phase_interval", "cell")):
        return True
    return any(key in value for key in ("ordinary_Z", "source_Z", "source_jets", "g5", "Z0", "Z1", "z0", "z1"))


def _default_cells(c: Any) -> list[dict[str, Any]]:
    return [
        {"index": index, "q_left": c.mpf(index) / CELL_COUNT,
         "q_right": c.mpf(index + 1) / CELL_COUNT}
        for index in range(CELL_COUNT)
    ]


def _normalise_cells(c: Any, cells: Any) -> list[dict[str, Any]]:
    if cells is None:
        raw_cells = _default_cells(c)
    elif isinstance(cells, int):
        if isinstance(cells, bool) or cells < 1:
            raise ValueError("cells must be a positive integer")
        raw_cells = [
            {"index": index, "q_left": c.mpf(index) / cells,
             "q_right": c.mpf(index + 1) / cells}
            for index in range(cells)
        ]
    elif isinstance(cells, Mapping) and "cells" in cells:
        raw_cells = list(cells["cells"])
    elif _is_sequence(cells):
        raw_cells = list(cells)
    else:
        raise TypeError("cells must be a positive count or cell sequence")
    if not raw_cells:
        raise ValueError("at least one phase cell is required")
    result: list[dict[str, Any]] = []
    previous = None
    for index, raw in enumerate(raw_cells):
        if isinstance(raw, Mapping):
            nested = raw.get("q_interval", raw.get("phase_interval", raw.get("interval")))
            left = raw.get("q_left", raw.get("left"))
            right = raw.get("q_right", raw.get("right"))
            if nested is not None and (left is None or right is None):
                if not _is_sequence(nested) or len(nested) != 2:
                    raise TypeError(f"cells[{index}].q_interval must have two endpoints")
                left, right = nested
        elif _is_sequence(raw) and len(raw) == 2:
            left, right = raw
        else:
            raise TypeError(f"cells[{index}] must specify q_left and q_right")
        if left is None or right is None:
            raise ValueError(f"cells[{index}] needs q_left and q_right")
        left, right = _interval(c, left, f"cells[{index}].q_left"), _interval(c, right, f"cells[{index}].q_right")
        lo, _ = _bounds(left)
        _, hi = _bounds(right)
        if lo < 0 or hi > 1 or hi <= lo:
            raise ValueError(f"cells[{index}] must lie in [0,1] with positive width")
        if previous is not None:
            previous_hi = _bounds(previous)[1]
            if lo != previous_hi:
                raise ValueError("phase cells must be contiguous and ordered")
        previous = right
        result.append({"index": index, "q_left": left, "q_right": right,
                       "q_interval": _make_interval(c, left, right),
                       "q_width": right - left})
    if _bounds(result[0]["q_left"])[0] != 0 or _bounds(result[-1]["q_right"])[1] != 1:
        raise ValueError("phase cells must cover exactly q in [0,1]")
    return result


def _normalise_source_cells(c: Any, raw: Any, cells: Sequence[Mapping[str, Any]]) -> list[dict[int, list[dict[str, Any]]]]:
    count = len(cells)
    if isinstance(raw, Mapping) and "cells" in raw:
        raw = raw["cells"]
    per_cell: list[Any] | None = None
    if isinstance(raw, Mapping) and not any(_lookup_order(raw, order)[0] for order in ORDINARY_Z_ORDERS):
        if all(index in raw or str(index) in raw for index in range(count)):
            per_cell = [raw.get(index, raw.get(str(index))) for index in range(count)]
    elif _is_sequence(raw) and len(raw) == count and any(_looks_like_cell(item) for item in raw):
        per_cell = list(raw)
    if per_cell is None:
        per_cell = [raw for _ in range(count)]
    if len(per_cell) != count:
        raise ValueError("source cells do not match phase-cell count")
    return [
        _rows(c, item, label=f"g5_cell[{index}]", default_source_id=f"g5_cell:{index}",
              require_source_id=True)
        for index, item in enumerate(per_cell)
    ]


def _uniform_packet(raw: Any, count: int) -> bool:
    """Return whether *raw* is one whole-collar packet, not cellwise data."""

    if isinstance(raw, Mapping):
        if "cells" in raw:
            return False
        if any(_lookup_order(raw, order)[0] for order in ORDINARY_Z_ORDERS):
            return True
        if all(index in raw or str(index) in raw for index in range(count)):
            return False
        return True
    if _is_sequence(raw):
        return not (len(raw) == count and any(_looks_like_cell(item) for item in raw))
    return True


def _normalise_chi(c: Any, chi: Any, chi_Z: Any, cells: Sequence[Mapping[str, Any]]) -> list[dict[int, list[dict[str, Any]]]]:
    count = len(cells)
    if chi is None:
        raise ValueError("chi is required for the first-potential correction")
    # A separately supplied chi_Z augments a scalar or Z0-only chi packet.
    if chi_Z is not None and isinstance(chi, Mapping) and _lookup_order(chi, 1)[0]:
        raise ValueError("chi_Z was supplied twice")
    if _is_sequence(chi) and len(chi) == count and any(_looks_like_cell(item) for item in chi):
        values = list(chi)
    elif isinstance(chi, Mapping) and "cells" in chi:
        values = list(chi["cells"])
    elif isinstance(chi, Mapping) and not any(_lookup_order(chi, order)[0] for order in ORDINARY_Z_ORDERS):
        if all(index in chi or str(index) in chi for index in range(count)):
            values = [chi.get(index, chi.get(str(index))) for index in range(count)]
        else:
            values = [chi] * count
    else:
        values = [chi] * count
    if len(values) != count:
        raise ValueError("chi cells do not match phase-cell count")
    result: list[dict[int, list[dict[str, Any]]]] = []
    for index, value in enumerate(values):
        if chi_Z is not None:
            if isinstance(value, Mapping):
                z0_found, z0 = _lookup_order(value, 0)
            else:
                z0_found, z0 = True, value
            if not z0_found:
                raise ValueError(f"chi cell {index} has no Z0 value")
            rows = {
                0: _term_list(c, z0, label=f"chi_cell[{index}].Z0", default_source_id=f"chi:{index}:Z0", require_source_id=False),
                1: _term_list(c, chi_Z if not (_is_sequence(chi_Z) and len(chi_Z) == count) else chi_Z[index],
                              label=f"chi_cell[{index}].Z1", default_source_id=f"chi:{index}:Z1", require_source_id=False),
            }
        else:
            if isinstance(value, Mapping):
                candidate = value
                for wrapper in ("ordinary_Z", "chi_Z", "chi_jets", "jets", "rows"):
                    if wrapper in candidate:
                        candidate = candidate[wrapper]
                        break
                found0, _ = _lookup_order(candidate, 0)
                found1, _ = _lookup_order(candidate, 1)
                if found0 or found1:
                    if not found0:
                        raise ValueError(f"chi cell {index} has no Z0 value")
                    if not found1:
                        raise ValueError(
                            f"chi cell {index} requires an explicit chi_Z/Z1 row"
                        )
                    rows = _rows(c, value, label=f"chi_cell[{index}]", default_source_id=f"chi:{index}", require_source_id=False)
                else:
                    raise ValueError(
                        f"chi cell {index} requires explicit chi_Z or an ordinary Z1 row"
                    )
            else:
                raise ValueError(
                    f"chi cell {index} requires an explicit chi_Z row"
                )
        result.append(rows)
    return result


def x(c: Any, hb: Any, sc: Any, q: Any) -> Any:
    """Directed geometry ``x(q)=2*exp(hb*sc*q/8)``."""

    hb, sc, q = _interval(c, hb, "hb"), _interval(c, sc, "sc"), _interval(c, q, "q")
    return c.mpf(2) * c.exp(hb * sc * q / 8)


def dx_dq(c: Any, hb: Any, sc: Any, q: Any) -> Any:
    """Return ``hb*sc*x(q)/8`` as a directed interval."""

    return _interval(c, hb, "hb") * _interval(c, sc, "sc") * x(c, hb, sc, q) / 8


def phi_nonnegative(c: Any, z: Any) -> Any:
    """Safe interval bound for ``phi(z)=(1-exp(-z))/z`` on ``z>=0``.

    The value at zero is handled by the theorem ``exp(-z) <= phi(z) <= 1``;
    no interval division at a removable zero is attempted.
    """

    z = _nonnegative(c, z, "phi argument")
    lower = c.exp(-z)
    return _make_interval(c, lower.a, c.mpf(1))


def _geometry_mass(c: Any, hb: Any, sc: Any, q_star: Any,
                   outer_left: Any, outer_right: Any,
                   source_left: Any, source_right: Any) -> dict[str, Any]:
    hb, sc, q_star = (_interval(c, value, name) for value, name in (
        (hb, "hb"), (sc, "sc"), (q_star, "q_star")))
    outer_left, outer_right = _interval(c, outer_left, "outer_left"), _interval(c, outer_right, "outer_right")
    source_left, source_right = _interval(c, source_left, "source_left"), _interval(c, source_right, "source_right")
    outer = _make_interval(c, outer_left, outer_right)
    source = _make_interval(c, source_left, source_right)
    outer_width = _nonnegative(c, outer_right - outer_left, "outer width")
    source_width = _nonnegative(c, source_right - source_left, "source width")
    d_outer = _nonnegative(c, q_star - outer, "outer target distance")
    # The source cell is either strictly before the outer cell or its diagonal
    # triangular cell.  In the latter case the rectangle bound uses the full
    # nonnegative difference range [0, width], which safely contains p<=r.
    d_inner_raw = outer - source
    _, d_inner_hi = _bounds(d_inner_raw)
    d_inner = _make_interval(c, max(_bounds(d_inner_raw)[0], type(d_inner_hi)(0)), d_inner_hi)
    cphase = sc / 8
    hphase = hb * cphase
    x_outer = c.mpf(2) * c.exp(hphase * outer)
    x_source = c.mpf(2) * c.exp(hphase * source)
    phi_outer = phi_nonnegative(c, 2 * hphase * d_outer)
    phi_inner = phi_nonnegative(c, 2 * hphase * d_inner)
    decay_outer = c.exp(-3 * hphase * d_outer)
    # Formal hb factors are deliberately absent from these coefficients.
    f1_mass = (
        cphase**4 * x_outer**2 * d_outer * phi_outer
        * x_source**2 * d_inner * phi_inner
        * outer_width * source_width
    )
    dx_mass = (
        cphase**3 * x_outer * decay_outer
        * x_source**2 * d_inner * phi_inner
        * outer_width * source_width
    )
    return {
        "F1": _nonnegative(c, f1_mass, "F1 physical kernel mass"),
        "d_x_F1": _nonnegative(c, dx_mass, "d_x_F1 physical kernel mass"),
        "outer_q_interval": outer,
        "source_q_interval": source,
        "q_star": q_star,
        "outer_width": outer_width,
        "source_width": source_width,
        "outer_distance": d_outer,
        "inner_distance": d_inner,
        "x_outer_interval": x_outer,
        "x_source_interval": x_source,
        "phi_outer_interval": phi_outer,
        "phi_inner_interval": phi_inner,
        "decay_outer_interval": decay_outer,
        "formal_hb_power": {"F1": 4, "d_x_F1": 3},
        "physical_dx_count": {"F1": 2, "d_x_F1": 2},
        "vanishing_K0_count": {"F1": 2, "d_x_F1": 1},
        "measure_formula": "two physical dx factors are included once in K0 o K0",
        "coefficient_excludes_formal_hb": True,
        "safe_rectangle_geometry": True,
    }


def _whole_domain_mass(c: Any, hb: Any, sc: Any, q_star: Any = 1) -> dict[str, Any]:
    """Safe full ``0<=p<=r<=q*=1`` masses for a whole-collar source box.

    The bounds factor out the formal ``hb`` powers.  With ``c=sc/8`` and
    ``A=sc/16`` the positive F1 mass is bounded by
    ``(2/3)c**4*[exp(-4A),exp(4A)]``.  For the target derivative we use the
    safe ``sinhc(A)<=exp(A)`` consequence, giving
    ``(4/3)c**3*[1,exp(3A)]``.  These are coefficient bounds, not selected
    physical values.
    """

    hb, sc, q_star = (_interval(c, value, name) for value, name in (
        (hb, "hb"), (sc, "sc"), (q_star, "q_star")))
    q_lo, q_hi = _bounds(q_star)
    if q_lo != 1 or q_hi != 1:
        raise ValueError("whole-domain mass is defined only at q*=1")
    cphase = sc / 8
    A = sc / 16
    exp_symmetric = _make_interval(c, c.exp(-4 * A).a, c.exp(4 * A).b)
    f1 = c.mpf(2) / 3 * cphase**4 * exp_symmetric
    dx = c.mpf(4) / 3 * cphase**3 * _make_interval(c, c.mpf(1), c.exp(3 * A).b)
    zero = c.mpf(0)
    one = c.mpf(1)
    return {
        "F1": _nonnegative(c, f1, "whole-collar F1 mass"),
        "d_x_F1": _nonnegative(c, dx, "whole-collar d_x_F1 mass"),
        "outer_q_interval": _make_interval(c, zero, one),
        "source_q_interval": _make_interval(c, zero, one),
        "q_star": q_star,
        "outer_width": one,
        "source_width": one,
        "outer_distance": _make_interval(c, zero, one),
        "inner_distance": _make_interval(c, zero, one),
        "x_outer_interval": x(c, hb, sc, _make_interval(c, zero, one)),
        "x_source_interval": x(c, hb, sc, _make_interval(c, zero, one)),
        "phi_outer_interval": _make_interval(c, c.mpf(0), c.mpf(1)),
        "phi_inner_interval": _make_interval(c, c.mpf(0), c.mpf(1)),
        "decay_outer_interval": _make_interval(c, c.exp(-3 * hb * cphase).a, c.mpf(1)),
        "formal_hb_power": {"F1": 4, "d_x_F1": 3},
        "physical_dx_count": {"F1": 2, "d_x_F1": 2},
        "vanishing_K0_count": {"F1": 2, "d_x_F1": 1},
        "measure_formula": "whole-collar positive mass bounds; two physical dx factors once",
        "coefficient_excludes_formal_hb": True,
        "safe_rectangle_geometry": False,
        "whole_domain_geometry": True,
        "full_domain_mass_bounds": {
            "F1": "(2/3)*(sc/8)^4*[exp(-4A),exp(4A)]",
            "d_x_F1": "(4/3)*(sc/8)^3*[1,exp(3A)]",
            "A": "sc/16",
        },
    }


def endpoint_kernel_masses(
    c: Any,
    *,
    hb: Any,
    sc: Any,
    outer_q: Any,
    source_q: Any,
    q_star: Any = 1,
) -> dict[str, Any]:
    """Return safe positive coefficient masses for one nested cell pair.

    ``outer_q`` and ``source_q`` are endpoint pairs.  The formal powers of
    ``hb`` are returned separately and are not materialized in either mass.
    """

    return _geometry_mass(c, hb, sc, q_star, outer_q[0], outer_q[1], source_q[0], source_q[1])


kernel_masses = endpoint_kernel_masses
physical_kernel_masses = endpoint_kernel_masses


def _combine_origin(*terms: Mapping[str, Any], **metadata: Any) -> dict[str, Any]:
    ids: list[Any] = []
    classes: list[Any] = []
    uncertain = bool(metadata.pop("uncertainty", False))
    for term in terms:
        info = _origin_info(term)
        ids.extend(info["source_ids"])
        classes.extend(info["error_classes"])
        uncertain = uncertain or info["uncertainty"]
    ids.extend(_items(metadata.pop("source_ids", [])))
    classes.extend(_items(metadata.pop("error_classes", [])))
    result = {
        "operation": metadata.pop("operation", "Dg_sector_product"),
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": bool(uncertain),
        "symbolic_cancellation_claimed": False,
    }
    result.update(metadata)
    return result


def _output_term(c: Any, chi_term: Mapping[str, Any], source_term: Mapping[str, Any],
                 mass: Any, *, target_cell: int, source_cell: int,
                 chi_order: int, source_order: int, field: str,
                 leibniz_multiplier: int = 1) -> dict[str, Any]:
    powers = tuple(a + b for a, b in zip(chi_term["powers"], source_term["powers"]))
    hb_added = 4 if field == "F1" else 3
    powers = (powers[0] + hb_added, powers[1], powers[2], powers[3])
    coefficient = -2 * leibniz_multiplier * chi_term["coefficient"] * source_term["coefficient"] * mass
    origin = _combine_origin(
        chi_term.get("origin", chi_term), source_term.get("origin", source_term),
        source_ids=[f"outer_cell:{target_cell}", f"source_cell:{source_cell}"],
        operation="first_collar_Dg_sector",
        target_cell=target_cell,
        source_cell=source_cell,
        source_row=SOURCE_ROW,
        chi_order=chi_order,
        source_order=source_order,
        field=field,
        leibniz_multiplier=leibniz_multiplier,
        physical_dx_count=2,
        formal_hb_power_added=hb_added,
        coefficient_contains_formal_hb=False,
        baseline_K0_reentered=False,
    )
    return {
        "powers": powers,
        "coefficient": coefficient,
        "origin": origin,
        "source_ids": list(origin["source_ids"]),
        "error_classes": list(origin["error_classes"]),
        "uncertainty": bool(origin["uncertainty"]),
        "target_cell": target_cell,
        "source_cell": source_cell,
        "source_row": SOURCE_ROW,
        "chi_order": chi_order,
        "source_order": source_order,
        "leibniz_multiplier": leibniz_multiplier,
        "formal_hb_power_added": hb_added,
        "physical_dx_count": 2,
        "vanishing_K0_count": 2 if field == "F1" else 1,
        "kernel": "(-2*chi) K0 o K0",
        "kernel_mass": mass,
        "kernel_mass_is_coefficient_only": True,
        "coefficient_contains_formal_hb": False,
        "formal_hb_measure_attached": True,
        "measure_attached_once": True,
        "baseline_K0_included": False,
    }


def build(
    c: Any,
    *,
    hb: Any,
    sc: Any,
    g5: Any = None,
    g5_source: Any = None,
    source_Z: Any = None,
    source: Any = None,
    source_cells: Any = None,
    g: Any = None,
    chi: Any = None,
    chi_Z: Any = None,
    chi_jets: Any = None,
    cells: Any = None,
    q_star: Any = 1,
    source_family: Any = None,
    source_atoms: Any = None,
    certificate: Any = None,
    geometry_certificate: Any = None,
) -> dict[str, Any]:
    """Build the new first-potential correction at endpoint ``q*=1``.

    ``g5``/``g5_source``/``source_Z``/``source``/``source_cells``/``g`` are
    aliases.  A source may be a single ordinary ``{Z0, Z1}`` packet reused as a certified whole-cell box,
    or a sequence/mapping with one such packet per phase cell.  Every source
    term must carry four integer powers and at least one source ID.

    ``chi`` may be scalar with an explicit ``chi_Z``, ordinary ``{Z0,Z1}``,
    or cellwise packets carrying explicit Z1 rows.  An omitted ``chi_Z`` is
    rejected unless every supplied chi packet explicitly contains its Z1 row.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf") or not hasattr(c, "exp"):
        raise TypeError("c must be the caller's directed interval context")
    hb, sc = _interval(c, hb, "hb"), _interval(c, sc, "sc")
    hb_lo, hb_hi = _bounds(hb)
    if hb_lo < 0 or hb_hi > type(hb_hi)(1) / 2:
        raise ValueError("hb must be certified in [0, 1/2]")
    _positive(c, sc, "sc")
    q_star = _interval(c, q_star, "q_star")
    q_lo, q_hi = _bounds(q_star)
    if q_lo != 1 or q_hi != 1:
        raise ValueError("this leaf computes only the exact endpoint q*=1")
    cert = geometry_certificate if geometry_certificate is not None else certificate
    if geometry_certificate is not None and certificate is not None and geometry_certificate != certificate:
        raise ValueError("geometry_certificate and certificate disagree")
    geometry_evidence = _geometry_evidence(cert)
    cert = geometry_evidence["raw"]
    if source_atoms is not None:
        names = set(source_atoms) if isinstance(source_atoms, Mapping) else set(source_atoms)
        if names != set(SOURCE_POWER_BASIS):
            raise ValueError("source_atoms must contain exactly the formal source basis")
    atoms = source_atoms if source_atoms is not None else list(SOURCE_POWER_BASIS)
    cells_norm = _normalise_cells(c, cells)
    aliases = [value is not None for value in (g5, g5_source, source_Z, source, source_cells, g)]
    if sum(aliases) > 1:
        raise ValueError("g5/g5_source/source_Z/source/source_cells/g are aliases and only one may be supplied")
    raw_source = next((value for value in (g5, g5_source, source_Z, source, source_cells, g) if value is not None), None)
    if raw_source is None:
        raise ValueError("whole-cell row-5 source sectors are required")
    source_rows = _normalise_source_cells(c, raw_source, cells_norm)
    if chi_jets is not None:
        if chi is not None or chi_Z is not None:
            raise ValueError("chi_jets cannot be combined with chi or chi_Z")
        chi = chi_jets
    chi_rows = _normalise_chi(c, chi, chi_Z, cells_norm)

    all_terms = {field: {order: [] for order in ORDINARY_Z_ORDERS} for field in FIELDS}
    cell_records: list[dict[str, Any]] = []
    mass_records: list[dict[str, Any]] = []
    source_uniform = _uniform_packet(raw_source, len(cells_norm))
    chi_uniform = _uniform_packet(chi, len(cells_norm))
    if source_uniform and chi_uniform:
        pair_specs = [
            ("whole", "whole", None, None, _whole_domain_mass(c, hb, sc, q_star),
             chi_rows[0], source_rows[0], False)
        ]
    else:
        pair_specs = []
        for outer_index, outer in enumerate(cells_norm):
            for source_index in range(outer_index + 1):
                source_cell = cells_norm[source_index]
                pair_specs.append((
                    outer_index, source_index, outer, source_cell,
                    _geometry_mass(
                        c, hb, sc, q_star,
                        outer["q_left"], outer["q_right"],
                        source_cell["q_left"], source_cell["q_right"],
                    ),
                    chi_rows[outer_index], source_rows[source_index],
                    outer_index == source_index,
                ))
    for outer_index, source_index, outer, source_cell, masses, chi_row, source_row, triangular in pair_specs:
        record = {
            "target_q": q_star,
            "outer_cell": outer_index,
            "source_cell": source_index,
            "outer_q_interval": masses["outer_q_interval"],
            "source_q_interval": masses["source_q_interval"],
            "triangular_inner_domain": triangular,
            "whole_domain_source_box": outer_index == "whole",
            "source_row": SOURCE_ROW,
            "kernel_masses": masses,
            "signed_terms": {field: {order: [] for order in ORDINARY_Z_ORDERS} for field in FIELDS},
        }
        # Z0: chi*g.  Z1: chi_Z*g + chi*g_Z at fixed certified geometry.
        for field in FIELDS:
            for chi_order, source_order, multiplier in ((0, 0, 1), (1, 0, 1), (0, 1, 1)):
                output_order = 0 if chi_order == 0 and source_order == 0 else 1
                mass = masses[field]
                for chi_term in chi_row[chi_order]:
                    for source_term in source_row[source_order]:
                        term = _output_term(
                            c, chi_term, source_term, mass,
                            target_cell=outer_index, source_cell=source_index,
                            chi_order=chi_order, source_order=source_order,
                            field=field, leibniz_multiplier=multiplier,
                        )
                        record["signed_terms"][field][output_order].append(term)
                        all_terms[field][output_order].append(term)
        cell_records.append(record)
        mass_records.append({
            "outer_cell": outer_index,
            "source_cell": source_index,
            "target_q": q_star,
            "triangular_inner_domain": triangular,
            "whole_domain_source_box": outer_index == "whole",
            **masses,
        })

    geometry = {
        "q_interval": "[0,1]",
        "q_star": q_star,
        "x_formula": "2*exp(hb*sc*q/8)",
        "dx_dq_formula": "hb*sc*x(q)/8",
        "x_q_star": x(c, hb, sc, q_star),
        "dx_dq_q_star": dx_dq(c, hb, sc, q_star),
        "hb_interval": hb,
        "sc_interval": sc,
        "hb_certified_in_closed_half": True,
        "sc_strictly_positive": True,
        "ordinary_Z_geometry_fixed": geometry_evidence["hb_sc_Z_independent"] and geometry_evidence["rho_radius_Z_independent"],
        "formal_hb_not_materialized": True,
        "formal_hb_measure_count": 2,
        "physical_dx_each_attached_once": True,
    }
    response = {
        field: {f"Z{order}": list(all_terms[field][order]) for order in ORDINARY_Z_ORDERS}
        for field in FIELDS
    }
    return {
        "schema": "n1_collar_sector_green_v3",
        "status": "conditional_partial_first_potential_correction",
        "source_family": source_family,
        "source_atoms": atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "source_row": SOURCE_ROW,
        "fields": list(FIELDS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "geometry": geometry,
        "cells": cells_norm,
        "target_q": q_star,
        "endpoint_q_star": q_star,
        "kernel": "Dg=(-2*chi) K0 o K0 g5",
        "kernel_mass_definition": "safe positive whole-cell interval geometry at q*=1",
        "physical_kernel_masses": mass_records,
        "cell_terms": cell_records,
        "whole_domain_source_box": bool(source_uniform and chi_uniform),
        "source_partition_required": False,
        "whole_domain_mass_formula": {
            "F1": "(2/3)*(sc/8)^4*[exp(-4*sc/16),exp(4*sc/16)]",
            "d_x_F1": "(4/3)*(sc/8)^3*[1,exp(3*sc/16)]",
        },
        "response": response,
        "Dg": response,
        "F1": response["F1"],
        "d_x_F1": response["d_x_F1"],
        "signed_sectors": True,
        "source_ids_preserved": True,
        "formal_powers_preserved": True,
        "formal_hb_powers_added": {"F1": 4, "d_x_F1": 3},
        "ordinary_Z1_product_rule": "chi_Z*g5_Z0 + chi*g5_Z1 at fixed q geometry",
        "ordinary_Z_geometry_differentiated": False,
        "baseline_K0_reentered": False,
        "baseline_K0_reference_only": True,
        "allword_tail_added": False,
        "Kchi_resummation_added": False,
        "m15_computed": False,
        "complete_physical_field": False,
        "finite_physical_field_claimed": False,
        "ancestor_imports": False,
        "receipt_reads": False,
        "certificate_passthrough": cert,
        "geometry_certificate": geometry_evidence,
        "conditional_on_geometry_certificate": True,
    }


build_collar_sector_green = build
build_Dg = build
