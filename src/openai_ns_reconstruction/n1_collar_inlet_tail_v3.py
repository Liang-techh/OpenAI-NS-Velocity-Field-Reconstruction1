"""Formal-sector ``j >= 2`` homogeneous collar Kchi tail.

This module is a deliberately small, caller-owned adapter.  It accepts the
partial inlet packets ``G0g``, ``core_Dg`` and ``core_symmetric_tail`` and
returns only symmetric error sectors for the omitted homogeneous Kchi words.
The finite K0 inlet and the ``j = 1`` potential word are not recomputed here.

The collar chart is ``x(q) = 2 exp(hb*sc*q/8)``, with
``c = sc/8``, ``A = sc/16`` and ``E = exp(A)``.  For
``C0 = sup |chi|`` and ``C1 = sup |chi_Z|`` the majorants are
``k0 = 8*c**2*E**2*C0``, ``kz = 8*c**2*E**2*C1``,
``kx = 4*c*E*C0`` and ``kxz = 4*c*E*C1`` with ``T = exp(k0/4)``.

The resulting sectors retain the caller's formal powers and add the physical
increments ``(+4,+5,+3,+4)`` for
``(F<-Fc,F<-Fxc,Fx<-Fc,Fx<-Fxc)``.  Every returned coefficient is an
explicit symmetric ``[-R,R]`` enclosure; no signed tail value is inferred.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_collar_inlet_tail_v3"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
FIELDS = ("F", "Fx")
INLET_COMPONENTS = ("Fc", "Fxc")
BRANCHES = ("F_from_Fc", "F_from_Fxc", "Fx_from_Fc", "Fx_from_Fxc")
FORMAL_HB_POWER_ADDED = {
    "F_from_Fc": 4,
    "F_from_Fxc": 5,
    "Fx_from_Fc": 3,
    "Fx_from_Fxc": 4,
}
GROUP_ALIASES = {
    "G0g": ("G0g", "g0g"),
    "core_Dg": ("core_Dg", "core_dg", "Dg", "d_g", "core"),
    "core_symmetric_tail": (
        "core_symmetric_tail", "core_symmetric", "symmetric_tail", "tail",
    ),
}

__all__ = [
    "SCHEMA", "SOURCE_POWER_BASIS", "ORDINARY_Z_ORDERS", "FIELDS",
    "INLET_COMPONENTS", "BRANCHES", "FORMAL_HB_POWER_ADDED",
    "kernel_bounds", "tail_kernel_bounds", "build", "build_inlet_tail",
    "build_tail", "inlet_tail", "evaluate",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


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
    """Coerce a value into the caller's directed interval context."""

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


def _symmetric(c: Any, radius: Any) -> Any:
    radius = _hi(c, radius)
    if radius < 0:
        raise ValueError("symmetric radius must be nonnegative")
    return _make(c, -radius, radius)


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lo, hi = value._mpi_
    if c.make_mpf((hi, hi)) < 0:
        raise ValueError(f"{label} is negative")
    if c.make_mpf((lo, lo)) < 0:
        lo = _zero(c)._mpi_[0]
    return c.make_mpf((lo, hi))


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
    """Collect nested provenance without evaluating any source expression."""

    result = {"source_ids": [], "error_classes": [], "uncertainty": False, "error_only": False}
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
            nested_info = _meta(nested, seen=seen)
            result["source_ids"].extend(nested_info["source_ids"])
            result["error_classes"].extend(nested_info["error_classes"])
            result["uncertainty"] = result["uncertainty"] or nested_info["uncertainty"]
            result["error_only"] = result["error_only"] or nested_info["error_only"]
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
        return (0, 0, 0, 0)
    value = raw.get("powers")
    if value is None:
        value = (
            raw.get("hb_power", 0),
            raw.get("F0_power", raw.get("F0_anchor_power", 0)),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)),
            raw.get("Lambda_power", raw.get("lambda_power", 0)),
        )
    if not _is_sequence(value) or len(value) != 4 or any(
        isinstance(item, bool) or not isinstance(item, int) for item in value
    ):
        raise ValueError(f"{label} must carry four integer formal powers")
    return tuple(int(item) for item in value)


def _term(c: Any, raw: Any, label: str, default_id: str) -> dict[str, Any]:
    item = dict(raw) if isinstance(raw, Mapping) else {"coefficient": raw}
    coefficient = next(
        (item[key] for key in ("coefficient", "coefficient_interval", "value", "interval") if key in item),
        None,
    )
    if coefficient is None:
        raise ValueError(f"{label} is missing coefficient")
    info = _meta(item)
    ids = list(info["source_ids"]) or [default_id]
    classes = list(info["error_classes"])
    error_only = bool(info["error_only"])
    if error_only and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    origin_raw = item.get("origin")
    origin_metadata = {
        "source_ids": _unique(ids), "error_classes": _unique(classes),
        "uncertainty": bool(info["uncertainty"] or error_only),
    }
    result = dict(item)
    result.update({
        "powers": _powers(item, label),
        "coefficient": _iv(c, coefficient, label + ".coefficient"),
        "source_ids": _unique(ids), "error_classes": _unique(classes),
        "uncertainty": origin_metadata["uncertainty"],
        "origin_metadata": origin_metadata,
    })
    if isinstance(origin_raw, str):
        result["origin"] = origin_raw
        result["origin_text"] = origin_raw
    elif isinstance(origin_raw, Mapping):
        mapped = dict(origin_raw)
        mapped.setdefault("source_ids", list(origin_metadata["source_ids"]))
        mapped.setdefault("error_classes", list(origin_metadata["error_classes"]))
        mapped.setdefault("uncertainty", origin_metadata["uncertainty"])
        result["origin"] = mapped
    else:
        result["origin"] = dict(origin_metadata)
    if error_only:
        result.update({"error_only": True, "symmetric_error": True})
    return result


def _term_like(value: Any) -> bool:
    return isinstance(value, Mapping) and any(
        key in value for key in ("coefficient", "coefficient_interval", "value", "powers")
    )


def _terms(c: Any, raw: Any, label: str, default_id: str) -> list[dict[str, Any]]:
    if raw is None:
        raise TypeError(f"{label} is required")
    if _term_like(raw):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        values = [raw]
    if not values:
        raise ValueError(f"{label} must not be empty")
    return [_term(c, value, f"{label}[{i}]", f"{default_id}:{i}") for i, value in enumerate(values)]


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
        for wrapper in ("ordinary_Z", "source_Z", "source_jets", "jets", "rows", "incoming"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)):
                raw = raw[wrapper]
                break
    if isinstance(raw, Mapping):
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
            index: _terms(c, found[index], f"{label}.Z{index}", f"{default_id}:Z{index}")
            for index in ORDINARY_Z_ORDERS
        }
    if _is_sequence(raw) and len(raw) == 2 and not all(_term_like(item) for item in raw):
        return {
            index: _terms(c, raw[index], f"{label}.Z{index}", f"{default_id}:Z{index}")
            for index in ORDINARY_Z_ORDERS
        }
    raise ValueError(f"{label} requires explicit ordinary Z0 and Z1 rows")


def _lookup(packet: Mapping[Any, Any], names: Sequence[str]) -> Any:
    for name in names:
        if name in packet:
            return packet[name]
    lowered = {str(key).lower(): value for key, value in packet.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return None


def _incoming(c: Any, packet: Any, label: str) -> dict[str, dict[int, list[dict[str, Any]]]]:
    if not isinstance(packet, Mapping):
        raise TypeError(f"{label} must map Fc and Fxc")
    fc = _lookup(packet, ("Fc", "F_c", "F", "incoming_Fc"))
    fxc = _lookup(packet, ("Fxc", "Fx_c", "F_xc", "Fx", "incoming_Fxc"))
    if fc is None or fxc is None:
        raise ValueError(f"{label} requires both Fc and Fxc")
    return {
        "Fc": _rows(c, fc, label + ".Fc", label + ":Fc"),
        "Fxc": _rows(c, fxc, label + ".Fxc", label + ":Fxc"),
    }


def _group_name(value: Any) -> str | None:
    text = str(value).lower().replace("-", "_").replace(" ", "_")
    for canonical, aliases in GROUP_ALIASES.items():
        normalized = {str(alias).lower().replace("-", "_").replace(" ", "_") for alias in aliases}
        if text in normalized:
            return canonical
    # Preserve arbitrary, uniquely named caller groups.  The only names that
    # are rejected are empty names; aliases above still canonicalize legacy
    # spellings for the three established packets.
    return str(value) if str(value) else None


def _groups(
    c: Any, incoming: Any, *, Fc: Any, Fxc: Any, incoming_Fc: Any,
    incoming_Fxc: Any, inlet_groups: Any, G0g: Any, core_Dg: Any,
    core_symmetric_tail: Any,
) -> dict[str, dict[str, dict[int, list[dict[str, Any]]]]]:
    packets: dict[str, Any] = {}
    for name, packet in (("G0g", G0g), ("core_Dg", core_Dg), ("core_symmetric_tail", core_symmetric_tail)):
        if packet is not None:
            packets[name] = packet
    if inlet_groups is not None:
        if not isinstance(inlet_groups, Mapping):
            raise TypeError("inlet_groups must be a mapping")
        for name, packet in inlet_groups.items():
            canonical = _group_name(name)
            if canonical is None or canonical in packets:
                raise ValueError(f"invalid or duplicate inlet group {name!r}")
            packets[canonical] = packet
    direct = any(value is not None for value in (Fc, Fxc, incoming_Fc, incoming_Fxc))
    if direct:
        if incoming is not None or packets:
            raise ValueError("direct Fc/Fxc cannot be combined with grouped incoming")
        packets["core_Dg"] = {"Fc": Fc if Fc is not None else incoming_Fc, "Fxc": Fxc if Fxc is not None else incoming_Fxc}
    elif incoming is not None:
        if not isinstance(incoming, Mapping):
            raise TypeError("incoming must be a mapping")
        nested = incoming.get("groups", incoming.get("inlet_groups"))
        direct_packet = _lookup(incoming, ("Fc", "Fxc", "F", "Fx", "F_c", "F_xc"))
        if nested is not None:
            if direct_packet is not None or not isinstance(nested, Mapping):
                raise ValueError("incoming groups packet is malformed")
            for name, packet in nested.items():
                canonical = _group_name(name)
                if canonical is None or canonical in packets:
                    raise ValueError(f"invalid or duplicate inlet group {name!r}")
                packets[canonical] = packet
        elif direct_packet is None and any(_group_name(key) is not None for key in incoming):
            for name, packet in incoming.items():
                canonical = _group_name(name)
                if canonical is not None:
                    if canonical in packets:
                        raise ValueError(f"duplicate inlet group {canonical}")
                    packets[canonical] = packet
        elif packets:
            raise ValueError("direct incoming cannot be combined with grouped incoming")
        else:
            packets["core_Dg"] = incoming
    if not packets:
        raise ValueError("at least one caller-owned inlet packet is required")
    return {name: _incoming(c, packet, name) for name, packet in packets.items()}


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
        if isinstance(raw, Mapping) and not _term_like(raw):
            z0 = _lookup(raw, ("Z0", "z0", "0", "value"))
            if z0 is None:
                raise ValueError("chi Z0 is required")
        else:
            z0 = raw
        return {0: _terms(c, z0, "chi.Z0", "chi:Z0"), 1: _terms(c, chi_Z, "chi.Z1", "chi:Z1")}
    return _rows(c, raw, "chi", "chi")


def _validate_chi_rows(chi: Mapping[int, Sequence[Mapping[str, Any]]]) -> None:
    """Validate the source-bound axis contract used by the nonlinear majorant.

    The supplied ``C0``/``C1`` coefficients are folded into ``k0``/``kz``
    and then into ``T``.  Consequently they cannot themselves carry formal
    atoms or represent several additive sectors.  Keeping one explicit,
    source-bound term per ordinary-Z row prevents a formal power from being
    silently counted once in the coefficient bound and again in the output
    sector.
    """

    for order in ORDINARY_Z_ORDERS:
        terms = chi.get(order)
        if terms is None or len(terms) != 1:
            raise ValueError(f"chi Z{order} requires exactly one source-bound sector")
        term = terms[0]
        if tuple(term.get("powers", ())) != (0, 0, 0, 0):
            raise ValueError("chi formal powers must be zero when folded into k0/kz")
        if not term.get("source_ids"):
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
    proof = evidence.get("strict_hb_admission_proof", evidence.get("strict_hb_proof", value.get("strict_hb_admission_proof")))
    if proof is True:
        strict = True
    elif isinstance(proof, Mapping):
        strict = strict or proof.get("strict_hb_admitted") is True or (
            proof.get("lower_strict") is True and proof.get("upper_strict") is True
        )
    elif isinstance(proof, str):
        strict = strict or proof.replace(" ", "").lower() in {"0<hb<1/2", "0<hb<.5", "strict:0<hb<1/2"}
    if not hb_sc or not rho or not strict:
        missing = []
        if not hb_sc:
            missing.append("hb_sc_Z_independent")
        if not rho:
            missing.append("rho_radius_Z_independent")
        if not strict:
            missing.append("strict_hb_admitted")
        raise ValueError("geometry certificate lacks: " + ", ".join(missing))
    return {"hb_sc_Z_independent": True, "rho_radius_Z_independent": True, "strict_hb_admitted": True, "raw": dict(value)}


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


def kernel_bounds(c: Any, *, sc: Any, chi_Z0: Sequence[Mapping[str, Any]], chi_Z1: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Return the directed majorants used by the ``j >= 2`` tail."""
    _validate_chi_rows({0: chi_Z0, 1: chi_Z1})
    sc = _positive(c, sc, "sc")
    C0, C1 = _sup(c, chi_Z0), _sup(c, chi_Z1)
    phase, A = sc / 8, sc / 16
    E = _nonnegative(c, c.exp(A), "E")
    k0 = _nonnegative(c, 8 * phase**2 * E**2 * C0, "k0")
    kz = _nonnegative(c, 8 * phase**2 * E**2 * C1, "kz")
    kx = _nonnegative(c, 4 * phase * E * C0, "kx")
    kxz = _nonnegative(c, 4 * phase * E * C1, "kxz")
    T = _nonnegative(c, c.exp(k0 / 4), "T")
    return {
        "c": phase, "A": A, "E": E, "C0": C0, "C1": C1,
        "k0": k0, "kz": kz, "kx": kx, "kxz": kxz, "T": T,
        "formulas": {
            "c": "sc/8", "A": "sc/16", "E": "exp(A)", "C0": "sup|chi|", "C1": "sup|chi_Z|",
            "k0": "8*c**2*E**2*C0", "kz": "8*c**2*E**2*C1", "kx": "4*c*E*C0", "kxz": "4*c*E*C1", "T": "exp(k0/4)",
        },
    }


tail_kernel_bounds = kernel_bounds


def _branch_factors(c: Any, bounds: Mapping[str, Any], sc: Any, branch: str) -> dict[str, Any]:
    """Return route factors with all formal hb powers removed."""
    k0, kz, kx, kxz, T = bounds["k0"], bounds["kz"], bounds["kx"], bounds["kxz"], bounds["T"]
    component = "Fxc" if branch.endswith("Fxc") else "Fc"
    source_scale = sc / 4 if component == "Fxc" else c.mpf(1)
    if branch.startswith("F_"):
        z0_factor = source_scale * k0 * k0 * T / 2
        z1_chi_factor = source_scale * kz * k0 * T
        field = "F"
    else:
        z0_factor = source_scale * kx * k0 * T / 2
        z1_chi_factor = source_scale * (kxz * k0 + kx * kz) * T / 2
        field = "Fx"
    return {
        "field": field, "incoming_component": component, "source_scale": source_scale,
        "Z0": _nonnegative(c, z0_factor, branch + ".Z0_factor"),
        "Z1_source": _nonnegative(c, z0_factor, branch + ".Z1_source_factor"),
        "Z1_chi": _nonnegative(c, z1_chi_factor, branch + ".Z1_chi_factor"),
    }


def _output_term(c: Any, *, source: Mapping[str, Any], chi: Mapping[str, Any], chi_provenance: Sequence[Mapping[str, Any]] | None = None, branch: str, group: str, source_order: int, chi_order: int, radius_factor: Any) -> dict[str, Any]:
    chi_terms = list(chi_provenance) if chi_provenance is not None else [chi]
    chi_ids = _unique(
        source_id
        for chi_term in chi_terms
        for source_id in chi_term.get("source_ids", [])
    )
    chi_classes = _unique(
        error_class
        for chi_term in chi_terms
        for error_class in chi_term.get("error_classes", [])
    )
    source_ids = _unique(list(source.get("source_ids", [])) + chi_ids)
    error_classes = _unique(list(source.get("error_classes", [])) + chi_classes)
    if "j2_kernel_remainder" not in error_classes:
        error_classes.append("j2_kernel_remainder")
    if "symmetric_error" not in error_classes:
        error_classes.append("symmetric_error")
    source_radius = _abs_upper(c, source["coefficient"])
    radius = _nonnegative(c, source_radius * radius_factor, "tail radius")
    coefficient = _symmetric(c, radius)
    powers0 = tuple(a + b for a, b in zip(source["powers"], chi["powers"]))
    powers = list(powers0)
    powers[0] += FORMAL_HB_POWER_ADDED[branch]
    piece = "chi_Z*inlet_Z0" if (chi_order, source_order) == (1, 0) else "chi*inlet_Z1" if (chi_order, source_order) == (0, 1) else "chi*inlet_Z0"
    origin_texts = [text for text in (chi.get("origin_text"), source.get("origin_text")) if isinstance(text, str)]
    origin = {
        "operation": "j>=2_Kchi_collar_homogeneous_tail_v3", "group": group, "branch": branch,
        "field": "F" if branch.startswith("F_") else "Fx", "incoming_component": "Fxc" if branch.endswith("Fxc") else "Fc",
        "chi_order": chi_order, "source_order": source_order, "incoming_order": source_order, "leibniz_piece": piece,
        "source_ids": source_ids, "chi_source_ids": chi_ids, "incoming_source_ids": list(source.get("source_ids", [])),
        "chi_source_classes": chi_classes,
        "error_classes": error_classes, "uncertainty": True, "error_only": True, "symmetric_error": True,
        "tail_starts_at_j": 2, "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch], "coefficient_contains_formal_hb": False,
        "chi_origin": _origin_descriptor(chi.get("origin")),
        "chi_origins": [_origin_descriptor(item.get("origin")) for item in chi_terms],
        "incoming_origin": _origin_descriptor(source.get("origin")),
    }
    return {
        "powers": tuple(powers), "coefficient": coefficient, "source_ids": source_ids,
        "chi_source_ids": chi_ids, "chi_source_classes": chi_classes,
        "incoming_source_ids": list(source.get("source_ids", [])),
        "error_classes": error_classes, "uncertainty": True, "origin": origin, "origin_texts": origin_texts,
        "group": group, "branch": branch, "field": "F" if branch.startswith("F_") else "Fx",
        "incoming_component": "Fxc" if branch.endswith("Fxc") else "Fc", "chi_order": chi_order,
        "source_order": source_order, "incoming_order": source_order, "leibniz_piece": piece,
        "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch], "tail_starts_at_j": 2, "kernel_tail_order": "j>=2",
        "radius": radius, "radius_factor": radius_factor, "coefficient_is_symmetric_radius": True,
        "coefficient_contains_formal_hb": False, "physical_measure_count": ">=2 per Volterra word", "physical_dx_count": ">=2 per Volterra word",
        "physical_dx_count_per_word": ">=2",
        "signed_source_term": False, "error_only": True, "symmetric_error": True,
        "baseline_K0_included": False, "j1_included": False,
    }


def _branch_set(c: Any, incoming: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]], chi: Mapping[int, Sequence[Mapping[str, Any]]], bounds: Mapping[str, Any], sc: Any, group: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    components = {"F_from_Fc": ("F", "Fc"), "F_from_Fxc": ("F", "Fxc"), "Fx_from_Fc": ("Fx", "Fc"), "Fx_from_Fxc": ("Fx", "Fxc")}
    for branch in BRANCHES:
        field, component = components[branch]
        factors = _branch_factors(c, bounds, sc, branch)
        rows = {"Z0": [], "Z1": []}
        for chi_term in chi[0]:
            for source_term in incoming[component][0]:
                rows["Z0"].append(_output_term(c, source=source_term, chi=chi_term, branch=branch, group=group, source_order=0, chi_order=0, radius_factor=factors["Z0"]))
        for source_term in incoming[component][1]:
            for chi_term in chi[0]:
                rows["Z1"].append(_output_term(c, source=source_term, chi=chi_term, branch=branch, group=group, source_order=1, chi_order=0, radius_factor=factors["Z1_source"]))
        for source_term in incoming[component][0]:
            for chi_term in chi[1]:
                rows["Z1"].append(_output_term(c, source=source_term, chi=chi_term, chi_provenance=(chi[0][0], chi[1][0]), branch=branch, group=group, source_order=0, chi_order=1, radius_factor=factors["Z1_chi"]))
        result[branch] = {
            "branch": branch, "field": field, "incoming_component": component, "response": rows,
            "remainder": rows, "tail": rows, "Z0": rows["Z0"], "Z1": rows["Z1"], "radius_factors": factors,
            "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch], "tail_starts_at_j": 2, "kernel_tail_order": "j>=2",
            "signed_sectors": False, "error_only_sectors_separate": True, "all_terms_error_only": True,
            "baseline_K0_included": False, "j1_included": False,
        }
    return result


def _aggregate(branches: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, list[Any]]]:
    """Concatenate sectors only; never interval-add error sectors."""
    return {
        "F": {"Z0": list(branches["F_from_Fc"]["Z0"]) + list(branches["F_from_Fxc"]["Z0"]), "Z1": list(branches["F_from_Fc"]["Z1"]) + list(branches["F_from_Fxc"]["Z1"])},
        "Fx": {"Z0": list(branches["Fx_from_Fc"]["Z0"]) + list(branches["Fx_from_Fxc"]["Z0"]), "Z1": list(branches["Fx_from_Fc"]["Z1"]) + list(branches["Fx_from_Fxc"]["Z1"])},
    }


def build(
    c: Any, incoming: Any = None, chi: Any = None, *, hb: Any, sc: Any,
    Fc: Any = None, Fxc: Any = None, incoming_Fc: Any = None, incoming_Fxc: Any = None,
    inlet_groups: Any = None, G0g: Any = None, core_Dg: Any = None, core_symmetric_tail: Any = None,
    chi_Z: Any = None, chi_jets: Any = None, source_family: Any = None, source_atoms: Any = None,
    geometry_certificate: Any = None, certificate: Any = None, q_star: Any = 1,
) -> dict[str, Any]:
    """Build the formal-sector symmetric ``j >= 2`` collar tail."""
    if not all(hasattr(c, name) for name in ("mpf", "make_mpf", "exp")):
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
    groups = _groups(c, incoming, Fc=Fc, Fxc=Fxc, incoming_Fc=incoming_Fc, incoming_Fxc=incoming_Fxc, inlet_groups=inlet_groups, G0g=G0g, core_Dg=core_Dg, core_symmetric_tail=core_symmetric_tail)
    constants = kernel_bounds(c, sc=sc, chi_Z0=chi_rows[0], chi_Z1=chi_rows[1])
    grouped = {name: _branch_set(c, packet, chi_rows, constants, sc, name) for name, packet in groups.items()}
    primary = "core_Dg" if "core_Dg" in grouped else next(iter(grouped))
    branches, aggregate = grouped[primary], _aggregate(grouped[primary])
    endpoint = c.mpf(2) * c.exp(hb * sc / 8)
    component_bounds = {
        "F": {"B_formula": "|Fc| + hb*(sc/4)*|Fxc|", "BZ_formula": "|Fc_Z| + hb*(sc/4)*|Fxc_Z|", "R_F": "B*hb**4*k0**2*T/2", "R_FZ": "BZ*hb**4*k0**2*T/2 + B*hb**4*kz*k0*T"},
        "Fx": {"B_formula": "|Fc| + hb*(sc/4)*|Fxc|", "BZ_formula": "|Fc_Z| + hb*(sc/4)*|Fxc_Z|", "R_Fx": "B*hb**3*kx*k0*T/2", "R_FxZ": "BZ*hb**3*kx*k0*T/2 + B*hb**3*(kxz*k0+kx*kz)*T/2"},
    }
    geometry = {
        "q_star": q, "q_interval": "[0,1]", "x_endpoint": endpoint, "x_formula": "2*exp(hb*sc*q/8)",
        "c": constants["c"], "A": constants["A"], "E": constants["E"], "hb_interval": hb, "sc_interval": sc,
        "admissible_hb": "0<hb<1/2", "strict_hb_admitted": True, "closed_hb_enclosure_allowed": True,
        "hb_sc_Z_independent": True, "rho_radius_Z_independent": True, "ordinary_Z_geometry_fixed": True,
        "ordinary_Z_geometry_differentiated": False, "formal_hb_not_materialized": True, "midpoint_values_used": False,
        "physical_dt_count": ">=2 per Volterra word", "physical_dx_count_per_word": ">=2",
    }
    return {
        "schema": SCHEMA, "status": "conditional_partial_j2_kchi_collar_inlet_tail", "source_family": source_family,
        "source_atoms": atoms, "source_power_basis": list(SOURCE_POWER_BASIS), "fields": list(FIELDS),
        "inlet_components": list(INLET_COMPONENTS), "ordinary_Z_orders": [0, 1], "branches": branches,
        "branch_names": list(BRANCHES), "response": aggregate, "tail": aggregate, "remainder": aggregate,
        "homogeneous_tail": aggregate, "D": aggregate, "D_x": aggregate["Fx"], "F": aggregate["F"], "Fx": aggregate["Fx"],
        "inlet_groups": list(grouped), "primary_inlet_group": primary, "group_results": grouped,
        "grouped_inlet_packets_separate": True, "grouped_inlet_packets_merged": False,
        "all_branches": {f"{name}:{branch}": item for name, packet in grouped.items() for branch, item in packet.items()},
        "chi_jets": {f"Z{order}": list(chi_rows[order]) for order in ORDINARY_Z_ORDERS},
        "incoming_jets": {component: {f"Z{order}": list(groups[primary][component][order]) for order in ORDINARY_Z_ORDERS} for component in INLET_COMPONENTS},
        "geometry": geometry, "geometry_certificate": cert, "certificate_passthrough": cert["raw"],
        "kernel": {"Kchi_tail": "sum_{j>=2}(-2*chi)^j*T^j*K0", "tail_order": "j>=2", "homogeneous": True, "constants": constants, "component_error_radii": component_bounds, "causal_positive_majorant": True},
        "kernel_bounds": constants, "constants": constants, "component_error_radii": component_bounds, "chi_Z_required": True,
        "chi_source_sector_count_required": 1, "chi_formal_powers_required": (0, 0, 0, 0),
        "ordinary_Z1_product_rule": "chi_Z*inlet_Z0 + chi*inlet_Z1", "ordinary_Z1_geometry_differentiated": False,
        "source_ids_preserved": True, "chi_ids_preserved": True, "Z_provenance_preserved": True, "error_classes_preserved": True,
        "formal_powers_preserved": True, "formal_hb_powers_added": dict(FORMAL_HB_POWER_ADDED), "symmetric_error_outputs": True,
        "error_only_outputs_explicit": True, "no_signed_cancellation": True, "finite_partial_inlet_only": True,
        "j2_kernel_remainder_included": True, "j2_remainder_included": True, "Kchi_higher_words_included": True,
        "baseline_K0_reentered": False, "baseline_K0_included": False, "j1_recomputed": False, "j1_included": False,
        "allword_tail_added": False, "Q_feedback_included": False, "P_feedback_included": False, "m15_computed": False,
        "complete_collar_field": False, "complete_collar": False, "complete_inlet": False, "complete_Kchi_inlet": False,
        "finite_physical_field_claimed": False, "old_builder_called": False, "ancestor_imports": False,
        "receipt_reads": False, "automatic_code_execution": False,
    }


build_inlet_tail = build
build_tail = build
inlet_tail = build
evaluate = build
