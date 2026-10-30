"""Local row-5 ``g5`` tail for the collar Kchi potential.

This is a source-sector-only leaf for the omitted words

    ``(Kchi - K0 - D) g5``
    ``= sum(j>=2) (-2*chi*K0)**j K0 g5``.

The caller supplies the already saved, whole-collar local ``g5`` sectors.  The
leaf does not evaluate the source, re-enter the first word ``D``, or consume a
previous receipt.  It returns only a symmetric enclosure for the omitted
words.  A signed ``g5`` sector is therefore retained as provenance and as an
absolute source radius, but is never turned into a signed tail value.

Put ``c=sc/8``, ``A=sc/16`` and ``E=exp(A)``.  On
``x(q)=2*exp(hb*c*q)`` with ``dx=hb*c*x(q)*dq``, one initial physical ``K0``
map has the whole-collar coefficient mass

    ``m0 = 2*c**2*E**2``

after its formal ``hb**2`` is factored out.  With
``C0=sup|chi|`` and ``C1=sup|chi_Z|``, use

    ``k0=8*c**2*E**2*C0``       ``kz=8*c**2*E**2*C1``
    ``kx=4*c*E*C0``              ``kxz=4*c*E*C1``
    ``T=exp(k0/4)``.

For source radii ``S0`` and ``S1`` the directed radii are

    ``F_Z0 = S0*hb**6*m0*k0**2*T/2``
    ``F_Z1 = S1*hb**6*m0*k0**2*T/2
             + S0*hb**6*m0*kz*k0*T``
    ``Fx_Z0 = S0*hb**5*m0*kx*k0*T/2``
    ``Fx_Z1 = S1*hb**5*m0*kx*k0*T/2
              + S0*hb**5*m0*(kxz*k0+kx*kz)*T/2``.

The ``+6`` and ``+5`` formal ``hb`` increments count the initial ``K0 g5``
and the two Kchi words at the start of the ``j>=2`` tail.  Every physical
``dx`` is attached once: one to the initial source map and one per ordered
Volterra word.  The interval ``T`` uses only the certified geometry bound
``0<hb<1/2``; ``hb`` remains a formal source atom and no midpoint is chosen.
Ordinary Z1 is the fixed-geometry product rule
``chi_Z*g5_Z0 + chi*g5_Z1``.  No ``P1`` field, m15 term, feedback term, or
terminal closure is represented here.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_collar_local_g5_j2_tail_v1"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
FIELDS = ("F", "Fx")
LEGACY_FIELDS = ("F1", "d_x_F1")
BRANCHES = ("F_from_g5", "Fx_from_g5")
FORMAL_HB_POWER_ADDED = {"F_from_g5": 6, "Fx_from_g5": 5}

__all__ = [
    "SCHEMA", "SOURCE_POWER_BASIS", "ORDINARY_Z_ORDERS", "FIELDS",
    "LEGACY_FIELDS", "BRANCHES", "FORMAL_HB_POWER_ADDED",
    "kernel_bounds", "tail_kernel_bounds", "build", "build_tail",
    "local_g5_tail", "evaluate",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    return list(value) if _is_sequence(value) else [value]


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    for value in values:
        if not any(value == old for old in result):
            result.append(value)
    return result


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    raw = getattr(value, "_mpf_", None)
    return raw if raw is not None else value


def _make(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw(lower, False), _raw(upper, True)))


def _iv(c: Any, value: Any, label: str = "interval") -> Any:
    """Coerce an interval-compatible value without choosing a midpoint."""

    if _is_interval(value):
        return _make(c, value, value)
    if isinstance(value, Mapping):
        for key in ("value", "interval", "box", "enclosure", "data"):
            if key in value:
                return _iv(c, value[key], label)
        for left, right in (("lower", "upper"), ("lo", "hi"), ("left", "right")):
            if left in value and right in value:
                return _make(c, value[left], value[right])
        raise TypeError(f"{label} requires an interval-compatible value")
    if _is_sequence(value) and len(value) == 2:
        return _make(c, value[0], value[1])
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{label} is not interval-compatible") from exc


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    value = _iv(c, value)
    lo, hi = value._mpi_
    return c.make_mpf((lo, lo)), c.make_mpf((hi, hi))


def _lo(c: Any, value: Any) -> Any:
    return _bounds(c, value)[0]


def _hi(c: Any, value: Any) -> Any:
    return _bounds(c, value)[1]


def _zero(c: Any) -> Any:
    return c.mpf(0)


def _abs_upper(c: Any, value: Any) -> Any:
    lo, hi = _bounds(c, value)
    if lo >= 0:
        return hi
    if hi <= 0:
        return -lo
    return max(-lo, hi)


def _make_nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lo, hi = _bounds(c, value)
    if hi < 0:
        raise ValueError(f"{label} is negative")
    if lo < 0:
        lo = _zero(c)
    return _make(c, lo, hi)


def _symmetric(c: Any, radius: Any) -> Any:
    radius = _hi(c, radius)
    if radius < 0:
        raise ValueError("symmetric radius must be nonnegative")
    return _make(c, -radius, radius)


def _positive(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    if _lo(c, value) <= 0:
        raise ValueError(f"{label} must be strictly positive")
    return value


def _strict_hb(c: Any, value: Any, *, allow_closed: bool) -> Any:
    value = _iv(c, value, "hb")
    lower, upper = _bounds(c, value)
    half = c.mpf(1) / 2
    if lower < 0 or upper > half:
        raise ValueError("hb must be enclosed in [0,1/2]")
    if not allow_closed and (lower <= 0 or upper >= half):
        raise ValueError("hb must enclose strict 0<hb<1/2")
    return value


def _meta(value: Any, *, seen: set[int] | None = None) -> dict[str, Any]:
    """Collect provenance while leaving source metadata structurally intact."""

    result = {
        "source_ids": [], "error_classes": [], "uncertainty": False,
        "error_only": False,
    }
    if not isinstance(value, Mapping):
        return result
    if seen is None:
        seen = set()
    if id(value) in seen:
        return result
    seen.add(id(value))
    for key in ("source_ids", "sources", "source_id", "source", "input_id"):
        if key in value:
            result["source_ids"].extend(_items(value[key]))
    for key in ("error_classes", "error_class", "uncertainty_class"):
        if key in value:
            result["error_classes"].extend(_items(value[key]))
    operation = str(value.get("operation", "")).lower()
    kind = str(value.get("kind", value.get("error_kind", ""))).lower()
    result["uncertainty"] = bool(
        value.get("uncertainty", value.get("is_uncertainty", False))
    ) or any(token in operation or token in kind for token in ("error", "uncertainty", "symmetric"))
    result["error_only"] = bool(value.get("error_only", value.get("symmetric_error", False)))
    for key in ("origin", "metadata", "parents"):
        child = value.get(key)
        if child is None or child is value:
            continue
        for nested in (_items(child) if _is_sequence(child) else [child]):
            info = _meta(nested, seen=seen)
            result["source_ids"].extend(info["source_ids"])
            result["error_classes"].extend(info["error_classes"])
            result["uncertainty"] = result["uncertainty"] or info["uncertainty"]
            result["error_only"] = result["error_only"] or info["error_only"]
    result["source_ids"] = _unique(result["source_ids"])
    result["error_classes"] = _unique(result["error_classes"])
    return result


def _origin_descriptor(value: Any) -> Any:
    if isinstance(value, str):
        return value
    if isinstance(value, Mapping):
        return dict(value)
    return None


def _powers(raw: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a sector mapping")
    powers = raw.get("powers")
    if powers is None:
        powers = (
            raw.get("hb_power", 0),
            raw.get("F0_power", raw.get("F0_anchor_power", 0)),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)),
            raw.get("Lambda_power", raw.get("lambda_power", 0)),
        )
    if not _is_sequence(powers) or len(powers) != 4 or any(
        isinstance(item, bool) or not isinstance(item, int) for item in powers
    ):
        raise ValueError(f"{label}.powers must contain four integer formal powers")
    return tuple(int(item) for item in powers)


def _term(c: Any, raw: Any, label: str, default_id: str) -> dict[str, Any]:
    # Scalar chi rows are accepted only as caller-owned values with the
    # deterministic fallback ID below.  Saved g5 sectors are mappings and
    # carry their own source IDs; those IDs are never replaced.
    item = dict(raw) if isinstance(raw, Mapping) else {"coefficient": raw}
    coefficient = next(
        (item[key] for key in ("coefficient", "coefficient_interval", "value", "interval") if key in item),
        None,
    )
    if coefficient is None:
        raise ValueError(f"{label} is missing coefficient")
    info = _meta(item)
    ids = list(info["source_ids"])
    if not ids:
        ids = [default_id]
    classes = list(info["error_classes"])
    error_only = bool(info["error_only"])
    if error_only and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    origin_raw = item.get("origin")
    origin_metadata = {
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": bool(info["uncertainty"] or error_only),
    }
    if isinstance(origin_raw, str):
        origin: Any = origin_raw
        origin_text = origin_raw
    elif isinstance(origin_raw, Mapping):
        origin = dict(origin_raw)
        origin.setdefault("source_ids", list(origin_metadata["source_ids"]))
        origin.setdefault("error_classes", list(origin_metadata["error_classes"]))
        origin.setdefault("uncertainty", origin_metadata["uncertainty"])
        origin_text = None
    else:
        origin = dict(origin_metadata)
        origin_text = None
    result = dict(item)
    result.update({
        "powers": _powers(item, label),
        "coefficient": _iv(c, coefficient, label + ".coefficient"),
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": origin_metadata["uncertainty"],
        "origin": origin,
    })
    if origin_text is not None:
        result["origin_text"] = origin_text
    if error_only:
        result.update({"error_only": True, "symmetric_error": True})
    return result


def _term_list(c: Any, raw: Any, label: str, default_id: str) -> list[dict[str, Any]]:
    if raw is None:
        raise ValueError(f"{label} must not be empty")
    if isinstance(raw, Mapping) and any(
        key in raw for key in ("coefficient", "coefficient_interval", "value", "interval", "powers")
    ):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        values = [raw]
    if not values:
        raise ValueError(f"{label} must not be empty")
    return [_term(c, value, f"{label}[{index}]", f"{default_id}:{index}") for index, value in enumerate(values)]


def _order(key: Any) -> int | None:
    if isinstance(key, int) and not isinstance(key, bool):
        return key if key in ORDINARY_Z_ORDERS else None
    if isinstance(key, str):
        text = key.lower()
        if text in {"0", "z0", "value", "ordinary_z0"}:
            return 0
        if text in {"1", "z1", "derivative", "ordinary_z1", "ordinary_z"}:
            return 1
    return None


def _rows(c: Any, raw: Any, label: str, default_id: str) -> dict[int, list[dict[str, Any]]]:
    if isinstance(raw, Mapping):
        for wrapper in ("ordinary_Z", "source_Z", "source_jets", "jets", "rows", "g5"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)):
                raw = raw[wrapper]
                break
    if not isinstance(raw, Mapping):
        if _is_sequence(raw) and len(raw) == 2 and not all(
            isinstance(item, Mapping) and any(key in item for key in ("coefficient", "powers", "value", "interval"))
            for item in raw
        ):
            return {
                index: _term_list(c, raw[index], f"{label}.Z{index}", f"{default_id}:Z{index}")
                for index in ORDINARY_Z_ORDERS
            }
        raise ValueError(f"{label} requires explicit ordinary Z0 and Z1 rows")
    found: dict[int, Any] = {}
    for key, value in raw.items():
        index = _order(key)
        if index is not None:
            if index in found:
                raise ValueError(f"{label} supplies Z{index} twice")
            found[index] = value
    if set(found) != {0, 1}:
        raise ValueError(f"{label} requires explicit ordinary Z0 and Z1 rows")
    return {
        index: _term_list(c, found[index], f"{label}.Z{index}", f"{default_id}:Z{index}")
        for index in ORDINARY_Z_ORDERS
    }


def _chi(c: Any, chi: Any, chi_Z: Any) -> dict[int, list[dict[str, Any]]]:
    if chi is None:
        raise ValueError("chi is required")
    raw = chi
    if isinstance(raw, Mapping):
        for wrapper in ("ordinary_Z", "chi_jets", "jets", "rows"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)):
                raw = raw[wrapper]
                break
    if chi_Z is not None:
        if isinstance(raw, Mapping) and any(_order(key) == 1 for key in raw):
            raise ValueError("chi_Z supplied twice")
        if isinstance(raw, Mapping) and not any(
            key in raw for key in ("coefficient", "coefficient_interval", "value", "interval", "powers")
        ):
            z0 = next((raw[key] for key in ("Z0", "z0", 0, "0", "value") if key in raw), None)
            if z0 is None:
                raise ValueError("chi Z0 is required")
        else:
            z0 = raw
        return {0: _term_list(c, z0, "chi.Z0", "chi:Z0"), 1: _term_list(c, chi_Z, "chi.Z1", "chi:Z1")}
    return _rows(c, raw, "chi", "chi")


def _validate_chi_rows(chi: Mapping[int, Sequence[Mapping[str, Any]]]) -> None:
    for order in ORDINARY_Z_ORDERS:
        terms = chi.get(order)
        if terms is None or len(terms) != 1:
            raise ValueError(f"chi Z{order} requires exactly one source-bound sector")
        if tuple(terms[0].get("powers", ())) != (0, 0, 0, 0):
            raise ValueError("chi formal powers must be zero when folded into the tail majorant")
        if not terms[0].get("source_ids"):
            raise ValueError(f"chi Z{order} must retain a source ID")


def _geometry(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError("explicit geometry certificate is required")
    evidence = value.get("geometry", value)
    if not isinstance(evidence, Mapping):
        raise ValueError("geometry certificate must be a mapping")
    hb_sc = evidence.get("hb_sc_Z_independent") is True or (
        evidence.get("hb_Z_independent") is True and evidence.get("sc_Z_independent") is True
    )
    rho = evidence.get("rho_radius_Z_independent") is True or (
        evidence.get("rho_Z_independent") is True and evidence.get("radius_Z_independent") is True
    )
    strict = any(
        evidence.get(key) is True or value.get(key) is True
        for key in ("strict_hb_admitted", "strict_hb_admissibility", "hb_strictly_between_zero_and_half")
    )
    proof = evidence.get(
        "strict_hb_admission_proof",
        evidence.get("strict_hb_proof", value.get("strict_hb_admission_proof")),
    )
    if proof is True:
        strict = True
    elif isinstance(proof, Mapping):
        strict = strict or proof.get("strict_hb_admitted") is True or (
            proof.get("lower_strict") is True and proof.get("upper_strict") is True
        )
    elif isinstance(proof, str):
        strict = strict or proof.replace(" ", "").lower() in {
            "0<hb<1/2", "0<hb<.5", "strict:0<hb<1/2"
        }
    if not hb_sc or not rho or not strict:
        missing = []
        if not hb_sc:
            missing.append("hb_sc_Z_independent")
        if not rho:
            missing.append("rho_radius_Z_independent")
        if not strict:
            missing.append("strict_hb_admitted")
        raise ValueError("geometry certificate lacks: " + ", ".join(missing))
    return {
        "hb_sc_Z_independent": True,
        "rho_radius_Z_independent": True,
        "strict_hb_admitted": True,
        "raw": dict(value),
    }


def _atoms(value: Any) -> list[Any]:
    if value is None:
        return list(SOURCE_POWER_BASIS)
    result = list(value.keys()) if isinstance(value, Mapping) else list(value)
    if set(result) != set(SOURCE_POWER_BASIS):
        raise ValueError("source_atoms must contain exactly four formal atoms")
    return result


def _sup(c: Any, terms: Sequence[Mapping[str, Any]]) -> Any:
    if len(terms) != 1:
        raise ValueError("the folded chi majorant requires exactly one sector")
    return _make(c, _zero(c), _abs_upper(c, terms[0]["coefficient"]))


def kernel_bounds(
    c: Any,
    *,
    sc: Any,
    chi_Z0: Sequence[Mapping[str, Any]],
    chi_Z1: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Return directed whole-collar coefficient majorants."""

    _validate_chi_rows({0: chi_Z0, 1: chi_Z1})
    sc = _positive(c, sc, "sc")
    phase, A = sc / 8, sc / 16
    E = _make_nonnegative(c, c.exp(A), "E")
    C0, C1 = _sup(c, chi_Z0), _sup(c, chi_Z1)
    m0 = _make_nonnegative(c, 2 * phase**2 * E**2, "m0")
    k0 = _make_nonnegative(c, 8 * phase**2 * E**2 * C0, "k0")
    kz = _make_nonnegative(c, 8 * phase**2 * E**2 * C1, "kz")
    kx = _make_nonnegative(c, 4 * phase * E * C0, "kx")
    kxz = _make_nonnegative(c, 4 * phase * E * C1, "kxz")
    T = _make_nonnegative(c, c.exp(k0 / 4), "T")
    return {
        "c": phase, "A": A, "E": E, "C0": C0, "C1": C1,
        "m0": m0, "k0": k0, "kz": kz, "kx": kx, "kxz": kxz, "T": T,
        "formulas": {
            "c": "sc/8", "A": "sc/16", "E": "exp(A)",
            "m0": "2*c**2*E**2", "C0": "sup|chi|", "C1": "sup|chi_Z|",
            "k0": "8*c**2*E**2*C0", "kz": "8*c**2*E**2*C1",
            "kx": "4*c*E*C0", "kxz": "4*c*E*C1", "T": "exp(k0/4)",
        },
        "m0_formal_hb_power": 2,
        "m0_physical_dx_count": 1,
        "tail_starts_at_j": 2,
    }


tail_kernel_bounds = kernel_bounds


def _source_radius(c: Any, term: Mapping[str, Any]) -> Any:
    return _make_nonnegative(c, _abs_upper(c, term["coefficient"]), "g5 source radius")


def _combine_ids(*values: Any) -> list[Any]:
    ids: list[Any] = []
    for value in values:
        if isinstance(value, Mapping):
            ids.extend(_items(value.get("source_ids")))
        elif _is_sequence(value):
            ids.extend(value)
        elif value is not None:
            ids.append(value)
    return _unique(ids)


def _combine_classes(*values: Any) -> list[Any]:
    classes: list[Any] = []
    for value in values:
        if isinstance(value, Mapping):
            classes.extend(_items(value.get("error_classes")))
        elif _is_sequence(value):
            classes.extend(value)
        elif value is not None:
            classes.append(value)
    return _unique(classes)


def _output_term(
    c: Any,
    *,
    source: Mapping[str, Any],
    chi: Mapping[str, Any],
    chi_provenance: Sequence[Mapping[str, Any]] | None = None,
    branch: str,
    source_order: int,
    chi_order: int,
    radius_factor: Any,
) -> dict[str, Any]:
    source_radius = _source_radius(c, source)
    radius = _make_nonnegative(c, source_radius * radius_factor, "local g5 j2 tail radius")
    coefficient = _symmetric(c, radius)
    chi_terms = list(chi_provenance) if chi_provenance is not None else [chi]
    if not chi_terms:
        raise ValueError("chi provenance must not be empty")
    # Only the emitted chi factor contributes to formal powers.  The extra
    # Z0 term below is provenance for the kz*k0 majorant, not a second
    # multiplicative chi factor.
    powers0 = tuple(a + b for a, b in zip(source["powers"], chi["powers"]))
    powers = list(powers0)
    powers[0] += FORMAL_HB_POWER_ADDED[branch]
    source_ids = _combine_ids(source)
    chi_ids = _combine_ids(*chi_terms)
    chi_classes = _combine_classes(*chi_terms)
    source_ids = _unique(source_ids + chi_ids)
    classes = _combine_classes(source, *chi_terms)
    classes = _unique(classes + ["j2_kernel_remainder", "symmetric_error"])
    piece = (
        "chi_Z*g5_Z0" if (chi_order, source_order) == (1, 0)
        else "chi*g5_Z1" if (chi_order, source_order) == (0, 1)
        else "chi*g5_Z0"
    )
    field = "F" if branch.startswith("F_") else "Fx"
    origin_texts = [
        text for text in (
            source.get("origin_text"),
            *(term.get("origin_text") for term in chi_terms),
        )
        if isinstance(text, str)
    ]
    chi_origins = [_origin_descriptor(term.get("origin")) for term in chi_terms]
    origin = {
        "operation": "local_g5_j>=2_Kchi_tail_v1",
        "branch": branch, "field": field,
        "chi_order": chi_order, "source_order": source_order,
        "ordinary_Z_piece": piece,
        "source_ids": source_ids, "chi_source_ids": chi_ids,
        "chi_source_classes": chi_classes, "chi_origins": chi_origins,
        "chi_provenance_count": len(chi_terms),
        "chi_factor_count": 1,
        "chi_formal_powers_multiplied": [list(chi["powers"])],
        "chi_provenance_not_multiplied": len(chi_terms) > 1,
        "g5_source_ids": list(source.get("source_ids", [])),
        "error_classes": classes, "uncertainty": True,
        "error_only": True, "symmetric_error": True,
        "tail_starts_at_j": 2, "operator_word_sign": "(-1)^j",
        "operator_word_coefficient": "(-2)^j",
        "absolute_majorant_drops_alternating_sign": True,
        "tail_word_sign_sequence": "j=2:+1, j=3:-1, alternating",
        "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch],
        "coefficient_contains_formal_hb": False,
        "initial_K0_attached_once": True,
        "physical_dx_count_minimum": 3,
        "signed_source_preserved": True,
        "signed_source_coefficient_preserved": True,
        "signed_tail_value": False,
        "baseline_K0_included": False, "j1_included": False,
        "source_origin": _origin_descriptor(source.get("origin")),
        "chi_origin": _origin_descriptor(chi.get("origin")),
        "origin_texts": origin_texts,
    }
    return {
        "powers": tuple(powers), "coefficient": coefficient,
        "signed_source_coefficient": source["coefficient"],
        "signed_source_coefficient_preserved": True,
        "source_ids": source_ids, "chi_source_ids": chi_ids,
        "chi_source_classes": chi_classes, "chi_origins": chi_origins,
        "chi_provenance_count": len(chi_terms), "chi_factor_count": 1,
        "chi_formal_powers_multiplied": [list(chi["powers"])],
        "chi_provenance_not_multiplied": len(chi_terms) > 1,
        "g5_source_ids": list(source.get("source_ids", [])),
        "error_classes": classes, "uncertainty": True, "origin": origin,
        "origin_texts": origin_texts,
        "branch": branch, "field": field,
        "chi_order": chi_order, "source_order": source_order,
        "ordinary_Z_piece": piece,
        "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch],
        "tail_starts_at_j": 2, "kernel_tail_order": "j>=2",
        "radius": radius, "radius_factor": radius_factor,
        "coefficient_is_symmetric_radius": True,
        "coefficient_contains_formal_hb": False,
        "initial_K0_formal_hb_power": 2,
        "initial_K0_physical_dx_count": 1,
        "tail_physical_dx_count_minimum": 2,
        "physical_dx_count_minimum": 3,
        "physical_dx_count": ">=3 per omitted Volterra word",
        "physical_dx_count_per_word": ">=3",
        "physical_measure_count": ">=3 per omitted Volterra word",
        "physical_measure_attached_once": True,
        "signed_source_term": False, "signed_tail_value": False,
        "error_only": True, "symmetric_error": True,
        "baseline_K0_included": False, "j1_included": False,
        "operator_word_sign": "(-1)^j",
        "operator_word_coefficient": "(-2)^j",
        "absolute_majorant_drops_alternating_sign": True,
        "tail_word_sign_sequence": "j=2:+1, j=3:-1, alternating",
        "double_integration": False,
    }


def _branch_factors(c: Any, constants: Mapping[str, Any], branch: str) -> dict[str, Any]:
    m0 = constants["m0"]
    k0, kz, kx, kxz, T = (
        constants[name] for name in ("k0", "kz", "kx", "kxz", "T")
    )
    if branch == "F_from_g5":
        z0 = m0 * k0 * k0 * T / 2
        z1_source = z0
        z1_chi = m0 * kz * k0 * T
    elif branch == "Fx_from_g5":
        z0 = m0 * kx * k0 * T / 2
        z1_source = z0
        z1_chi = m0 * (kxz * k0 + kx * kz) * T / 2
    else:
        raise ValueError(f"unknown branch {branch}")
    return {
        "branch": branch, "field": "F" if branch.startswith("F_") else "Fx",
        "Z0": _make_nonnegative(c, z0, branch + ".Z0_factor"),
        "Z1_source": _make_nonnegative(c, z1_source, branch + ".Z1_source_factor"),
        "Z1_chi": _make_nonnegative(c, z1_chi, branch + ".Z1_chi_factor"),
        "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch],
        "formulas": {
            "Z0": "m0*k0**2*T/2" if branch == "F_from_g5" else "m0*kx*k0*T/2",
            "Z1_source": "same as Z0",
            "Z1_chi": "m0*kz*k0*T" if branch == "F_from_g5" else "m0*(kxz*k0+kx*kz)*T/2",
        },
    }


def _branch_set(
    c: Any,
    source_rows: Mapping[int, Sequence[Mapping[str, Any]]],
    chi_rows: Mapping[int, Sequence[Mapping[str, Any]]],
    constants: Mapping[str, Any],
    branch: str,
) -> dict[str, Any]:
    factors = _branch_factors(c, constants, branch)
    rows = {"Z0": [], "Z1": []}
    for source in source_rows[0]:
        rows["Z0"].append(_output_term(
            c, source=source, chi=chi_rows[0][0], branch=branch,
            source_order=0, chi_order=0, radius_factor=factors["Z0"],
        ))
    for source in source_rows[1]:
        rows["Z1"].append(_output_term(
            c, source=source, chi=chi_rows[0][0], branch=branch,
            source_order=1, chi_order=0, radius_factor=factors["Z1_source"],
        ))
    for source in source_rows[0]:
        rows["Z1"].append(_output_term(
            c, source=source, chi=chi_rows[1][0],
            chi_provenance=(chi_rows[0][0], chi_rows[1][0]), branch=branch,
            source_order=0, chi_order=1, radius_factor=factors["Z1_chi"],
        ))
    return {
        "branch": branch, "field": "F" if branch.startswith("F_") else "Fx",
        "operator": "(Kchi-K0-D)g5",
        "formula": "sum(j>=2)(-2*chi*K0)^j K0*g5",
        "response": rows, "Z0": rows["Z0"], "Z1": rows["Z1"],
        "radius_factors": factors,
        "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch],
        "tail_starts_at_j": 2, "signed_sectors": False,
        "physical_dx_count_per_word": ">=3",
        "physical_measure_count": ">=3 per omitted Volterra word",
        "error_only_sectors_separate": True, "all_terms_error_only": True,
        "baseline_K0_included": False, "j1_included": False,
        "initial_K0_attached_once": True, "double_integration": False,
    }


def build(
    c: Any,
    *,
    hb: Any,
    sc: Any,
    g5: Any = None,
    g5_source: Any = None,
    source: Any = None,
    chi: Any = None,
    chi_Z: Any = None,
    chi_jets: Any = None,
    source_family: Any = None,
    source_atoms: Any = None,
    geometry_certificate: Any = None,
    certificate: Any = None,
    q_star: Any = 1,
) -> dict[str, Any]:
    """Build the local whole-collar ``g5`` homogeneous ``j>=2`` tail."""

    if not all(hasattr(c, key) for key in ("mpf", "make_mpf", "exp")):
        raise TypeError("c must be a directed interval context")
    if chi is not None and chi_jets is not None:
        raise ValueError("chi and chi_jets were supplied twice")
    chi = chi if chi is not None else chi_jets
    if geometry_certificate is not None and certificate is not None and geometry_certificate != certificate:
        raise ValueError("geometry_certificate and certificate disagree")
    cert = _geometry(geometry_certificate if geometry_certificate is not None else certificate)
    atoms = _atoms(source_atoms)
    q = _iv(c, q_star, "q_star")
    if _lo(c, q) != 1 or _hi(c, q) != 1:
        raise ValueError("q_star must be exact endpoint 1")
    hb = _strict_hb(c, hb, allow_closed=cert["strict_hb_admitted"])
    sc = _positive(c, sc, "sc")
    chi_rows = _chi(c, chi, chi_Z)
    _validate_chi_rows(chi_rows)
    aliases = [value is not None for value in (g5, g5_source, source)]
    if sum(aliases) > 1:
        raise ValueError("g5, g5_source and source are aliases and only one may be supplied")
    raw_source = next((value for value in (g5, g5_source, source) if value is not None), None)
    if raw_source is None:
        raise ValueError("whole-collar g5 sectors are required")
    source_rows = _rows(c, raw_source, "g5", "g5")
    constants = kernel_bounds(c, sc=sc, chi_Z0=chi_rows[0], chi_Z1=chi_rows[1])
    grouped = {
        branch: _branch_set(c, source_rows, chi_rows, constants, branch)
        for branch in BRANCHES
    }
    response = {
        "F": {"Z0": list(grouped["F_from_g5"]["Z0"]), "Z1": list(grouped["F_from_g5"]["Z1"])},
        "Fx": {"Z0": list(grouped["Fx_from_g5"]["Z0"]), "Z1": list(grouped["Fx_from_g5"]["Z1"])},
    }
    endpoint = c.mpf(2) * c.exp(hb * sc / 8)
    geometry = {
        "q_star": q, "q_interval": "[0,1]", "x_endpoint": endpoint,
        "x_formula": "2*exp(hb*sc*q/8)", "dx_dq_formula": "hb*sc*x(q)/8",
        "c": constants["c"], "A": constants["A"], "E": constants["E"],
        "hb_interval": hb, "sc_interval": sc,
        "admissible_hb": "0<hb<1/2", "strict_hb_admitted": True,
        "closed_hb_enclosure_allowed": True,
        "hb_sc_Z_independent": True, "rho_radius_Z_independent": True,
        "ordinary_Z_geometry_fixed": True, "ordinary_Z_geometry_differentiated": False,
        "formal_hb_not_materialized": True, "midpoint_values_used": False,
        "initial_K0_formal_hb_power": 2, "initial_K0_physical_dx_count": 1,
        "tail_physical_dx_count_minimum": 2, "physical_dx_count_minimum": 3,
        "physical_dx_double_counted": False,
    }
    component_radii = {
        "F": {
            "S0_formula": "|g5_Z0|", "S1_formula": "|g5_Z1|",
            "R_F_Z0": "S0*hb**6*m0*k0**2*T/2",
            "R_F_Z1": "S1*hb**6*m0*k0**2*T/2 + S0*hb**6*m0*kz*k0*T",
        },
        "Fx": {
            "S0_formula": "|g5_Z0|", "S1_formula": "|g5_Z1|",
            "R_Fx_Z0": "S0*hb**5*m0*kx*k0*T/2",
            "R_Fx_Z1": "S1*hb**5*m0*kx*k0*T/2 + S0*hb**5*m0*(kxz*k0+kx*kz)*T/2",
        },
    }
    return {
        "schema": SCHEMA,
        "status": "conditional_local_g5_j2_kchi_tail",
        "source_family": source_family, "source_atoms": atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "fields": list(FIELDS), "legacy_fields": list(LEGACY_FIELDS),
        "ordinary_Z_orders": [0, 1], "branches": grouped,
        "branch_names": list(BRANCHES), "response": response,
        "tail": response, "remainder": response,
        "F": response["F"], "Fx": response["Fx"],
        "F1": response["F"], "d_x_F1": response["Fx"],
        "g5_jets": {f"Z{order}": list(source_rows[order]) for order in ORDINARY_Z_ORDERS},
        "chi_jets": {f"Z{order}": list(chi_rows[order]) for order in ORDINARY_Z_ORDERS},
        "geometry": geometry, "geometry_certificate": cert,
        "certificate_passthrough": cert["raw"],
        "kernel": {
            "operator": "(Kchi-K0-D)g5",
            "series": "sum_{j>=2}(-2*chi*K0)^j K0*g5",
            "tail_order": "j>=2", "homogeneous": True,
            "constants": constants, "component_error_radii": component_radii,
            "causal_positive_majorant": True,
        },
        "kernel_bounds": constants, "constants": constants,
        "component_error_radii": component_radii,
        "chi_Z_required": True, "chi_source_sector_count_required": 1,
        "chi_formal_powers_required": (0, 0, 0, 0),
        "ordinary_Z1_product_rule": "chi_Z*g5_Z0 + chi*g5_Z1",
        "ordinary_Z1_chi_majorant_uses_both_chi_rows": True,
        "ordinary_Z1_chi_factor_count": 1,
        "ordinary_Z1_chi_formal_powers_multiplied": "chi_Z only; both chi rows retained as provenance",
        "ordinary_Z1_geometry_differentiated": False,
        "source_ids_preserved": True, "source_signs_preserved_as_provenance": True,
        "chi_ids_preserved": True, "Z_provenance_preserved": True,
        "error_classes_preserved": True, "formal_powers_preserved": True,
        "formal_hb_powers_added": dict(FORMAL_HB_POWER_ADDED),
        "symmetric_error_outputs": True, "error_only_outputs_explicit": True,
        "no_signed_tail_value": True, "operator_word_sign": "(-1)^j",
        "operator_word_coefficient": "(-2)^j",
        "absolute_majorant_drops_alternating_sign": True,
        "tail_word_sign_sequence": "j=2:+1, j=3:-1, alternating",
        "initial_K0_attached_once": True, "physical_dx_each_attached_once": True,
        "physical_dx_double_counted": False,
        "baseline_K0_included": False, "j1_included": False,
        "D_recomputed": False, "first_potential_reference_only": True,
        "Kchi_higher_words_included": True, "j2_kernel_remainder_included": True,
        "allword_tail_added": False, "Q_feedback_included": False,
        "P_feedback_included": False, "P1_included": False,
        "m15_computed": False, "appended_to_ledger": False,
        "complete_collar_field": False, "complete_collar": False,
        "finite_physical_field_claimed": False, "row5_local_source_only": True,
        "old_builder_called": False, "ancestor_imports": False,
        "receipt_reads": False, "automatic_code_execution": False,
    }


build_tail = build
local_g5_tail = build
evaluate = build

