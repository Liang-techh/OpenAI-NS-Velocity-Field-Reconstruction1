"""Formal-sector homogeneous Q/P collar transport at the first endpoint.

This module is deliberately a small, independent consumer.  It accepts caller
owned six-component inlet sectors and returns three disjoint pieces:

``K0``
    the no-potential homogeneous transport of ``(F,U,K,P,Fx,Ux)``;
``first_potential``
    the four finite ``-2 chi K0 H0`` branches; and
``j2_tail``
    a separate, symmetric enclosure for the omitted potential words ``j>=2``.

The four source powers ``(hb,F0_anchor,Pstar_squared,Lambda)`` stay on the
sector records.  They are never evaluated.  In particular, factors which
vanish at the seam are represented by an added formal ``hb`` power and a
coefficient interval.  No old builder, checker, receipt, or source artifact
is imported here.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_collar_qp_transport_v2"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
FIELDS = ("F", "U", "K", "P", "Fx", "Ux")
K0_FIELDS = FIELDS
POTENTIAL_FIELDS = ("F", "Fx")
INLET_COMPONENTS = ("F", "U", "K", "P", "Fx", "Ux")
BRANCHES = ("F_from_Fc", "F_from_Fxc", "Fx_from_Fc", "Fx_from_Fxc")

# The first potential word adds the indicated formal hb power.  The same
# increments are used by the j>=2 majorant, with two more vanishing K0 factors.
FIRST_HB_POWER_ADDED = {
    "F_from_Fc": 2,
    "F_from_Fxc": 3,
    "Fx_from_Fc": 1,
    "Fx_from_Fxc": 2,
}
J2_HB_POWER_ADDED = {
    "F_from_Fc": 4,
    "F_from_Fxc": 5,
    "Fx_from_Fc": 3,
    "Fx_from_Fxc": 4,
}

K0_BRANCHES = (
    "F_from_Fc",
    "F_from_Fxc",
    "U_from_Uc",
    "U_from_Uxc",
    "K_from_Kc",
    "K_from_Uxc",
    "P_from_Pc",
    "Fx_from_Fxc",
    "Ux_from_Uxc",
)

MASS_FORMULAS = {
    "F_from_Fc": "2*c**2*[exp(-2*A),exp(2*A)]",
    "F_from_Fxc": "(4/3)c**3*[exp(-4*A),exp(2*A)]",
    "Fx_from_Fc": "2*c*[exp(-3*A),exp(A)]",
    "Fx_from_Fxc": "2*c**2*[exp(-5*A),exp(A)]",
}

_COMPONENT_ALIASES = {
    "F": ("F", "Fc", "F_c", "F1"),
    "U": ("U", "Uc", "U_c", "Uz", "Uz1"),
    "K": ("K", "Kc", "K_c", "K1"),
    "P": ("P", "Pc", "P_c", "P1"),
    "Fx": ("Fx", "Fxc", "F_xc", "d_x_F", "d_x_F1"),
    "Ux": ("Ux", "Uxc", "U_xc", "d_x_U", "d_x_Uz1"),
}


__all__ = [
    "SCHEMA",
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "FIELDS",
    "K0_FIELDS",
    "POTENTIAL_FIELDS",
    "INLET_COMPONENTS",
    "BRANCHES",
    "K0_BRANCHES",
    "FIRST_HB_POWER_ADDED",
    "J2_HB_POWER_ADDED",
    "MASS_FORMULAS",
    "inlet_kernel_masses",
    "kernel_masses",
    "kernel_bounds",
    "build",
    "build_qp_transport",
    "transport",
    "build_transport",
    "evaluate",
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
    if _is_interval(value):
        return _make(c, value, value)
    if isinstance(value, Mapping):
        for key in ("value", "interval", "box", "enclosure", "data"):
            if key in value:
                return _iv(c, value[key], label)
        if "lower" in value and "upper" in value:
            return _make(c, value["lower"], value["upper"])
        if "lo" in value and "hi" in value:
            return _make(c, value["lo"], value["hi"])
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
    return max(-lo, hi)


def _symmetric(c: Any, value: Any) -> Any:
    radius = _abs_upper(c, value)
    return _make(c, -radius, radius)


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lo, hi = _bounds(c, value)
    if hi < 0:
        raise ValueError(f"{label} is negative")
    if lo < 0:
        lo = _zero(c)
    return _make(c, lo, hi)


def _positive(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lo, _ = _bounds(c, value)
    if lo <= 0:
        raise ValueError(f"{label} must be strictly positive")
    return value


def _strict_hb(c: Any, value: Any, *, allow_closed: bool) -> Any:
    value = _iv(c, value, "hb")
    lo, hi = _bounds(c, value)
    half = c.mpf(1) / 2
    if lo < 0 or hi > half:
        raise ValueError("hb must be contained in [0,1/2]")
    if not allow_closed and (lo <= 0 or hi >= half):
        raise ValueError("hb must enclose 0<hb<1/2")
    return value


def _expbox(c: Any, left: Any, right: Any) -> Any:
    left = _lo(c, left)
    right = _hi(c, right)
    lower = c.exp(left)
    upper = c.exp(right)
    return _make(c, lower, upper)


def _decay_bound(c: Any, sc: Any, multiplier: Any) -> Any:
    """Conservative ``[exp(-multiplier*sc/8),1]`` coefficient box."""

    sc_hi = _hi(c, sc)
    lower = c.exp(-multiplier * sc_hi / 8)
    return _make(c, lower, _zero(c) + 1)


def _meta(value: Any, *, seen: set[int] | None = None) -> dict[str, Any]:
    empty = {"source_ids": [], "error_classes": [], "uncertainty": False, "error_only": False}
    if not isinstance(value, Mapping):
        return empty
    if seen is None:
        seen = set()
    if id(value) in seen:
        return empty
    seen.add(id(value))
    ids: list[Any] = []
    for key in ("source_ids", "sources", "source_id", "source", "input_id"):
        if key in value:
            ids.extend(_items(value[key]))
    classes: list[Any] = []
    for key in ("error_classes", "error_class", "uncertainty_class"):
        if key in value:
            classes.extend(_items(value[key]))
    operation = str(value.get("operation", "")).lower()
    kind = str(value.get("kind", value.get("error_kind", ""))).lower()
    uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
    uncertain = uncertain or any(token in operation or token in kind for token in ("error", "uncertainty", "symmetric"))
    error_only = bool(value.get("error_only", value.get("symmetric_error", value.get("is_error_only", False))))
    for key in ("origin", "metadata", "parents"):
        child = value.get(key)
        if child is None or child is value:
            continue
        children = _items(child) if _is_sequence(child) else [child]
        for item in children:
            if not isinstance(item, Mapping):
                continue
            nested = _meta(item, seen=seen)
            ids.extend(nested["source_ids"])
            classes.extend(nested["error_classes"])
            uncertain = uncertain or nested["uncertainty"]
            error_only = error_only or nested["error_only"]
    return {
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": uncertain,
        "error_only": error_only,
    }


def _origin_descriptor(value: Any) -> Any:
    if isinstance(value, str):
        return value
    if isinstance(value, Mapping):
        return dict(value)
    return None


def _powers(raw: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(raw, Mapping):
        return (0, 0, 0, 0)
    powers = raw.get("powers", raw.get("formal_powers", raw.get("source_powers")))
    if powers is None:
        powers = (
            raw.get("hb_power", 0),
            raw.get("F0_power", raw.get("F0_anchor_power", raw.get("F0", 0))),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", raw.get("Pstar2_power", 0))),
            raw.get("Lambda_power", raw.get("lambda_power", 0)),
        )
    if not _is_sequence(powers) or len(powers) != 4 or any(isinstance(x, bool) or not isinstance(x, int) for x in powers):
        raise ValueError(f"{label} must carry four integer formal powers")
    return tuple(int(x) for x in powers)


def _term(c: Any, raw: Any, label: str, default_id: str) -> dict[str, Any]:
    item = dict(raw) if isinstance(raw, Mapping) else {"coefficient": raw}
    coefficient = next((item[key] for key in ("coefficient", "coefficient_interval", "value", "interval") if key in item), None)
    if coefficient is None:
        raise ValueError(f"{label} is missing coefficient")
    info = _meta(item)
    ids = list(info["source_ids"]) or [default_id]
    classes = list(info["error_classes"])
    error_only = bool(info["error_only"])
    if error_only and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    origin_raw = item.get("origin")
    origin_meta = {
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": bool(info["uncertainty"] or error_only),
    }
    result = dict(item)
    result.update(
        {
            "powers": _powers(item, label),
            "coefficient": _iv(c, coefficient, label + ".coefficient"),
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": origin_meta["uncertainty"],
            "origin_metadata": origin_meta,
        }
    )
    if isinstance(origin_raw, str):
        result["origin"] = origin_raw
        result["origin_text"] = origin_raw
    elif isinstance(origin_raw, Mapping):
        mapped = dict(origin_raw)
        mapped.update(origin_meta)
        result["origin"] = mapped
    else:
        result["origin"] = origin_meta
    if error_only:
        result.update({"error_only": True, "symmetric_error": True})
    return result


def _term_like(value: Any) -> bool:
    return isinstance(value, Mapping) and any(key in value for key in ("coefficient", "coefficient_interval", "value", "interval", "powers"))


def _one_term(c: Any, raw: Any, label: str, default_id: str) -> dict[str, Any] | None:
    if raw is None:
        raise ValueError(f"{label} is required")
    if _term_like(raw):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        values = [raw]
    # A caller may explicitly carry a structurally zero state row.  It is
    # represented by an empty list, never by a fabricated zero sector.
    if not values:
        return None
    if len(values) != 1:
        raise ValueError(f"{label} must contain exactly one sector")
    return _term(c, values[0], label, default_id)


def _term_list(c: Any, raw: Any, label: str, default_id: str) -> list[dict[str, Any]]:
    """Normalize one finite row while retaining each caller-owned sector."""

    if raw is None:
        raise ValueError(f"{label} is required")
    if _term_like(raw):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        values = [raw]
    return [
        _term(c, value, f"{label}[{index}]", f"{default_id}:{index}")
        for index, value in enumerate(values)
    ]


def _order(value: Any) -> int | None:
    if isinstance(value, int) and not isinstance(value, bool):
        return value if value in ORDINARY_Z_ORDERS else None
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"0", "z0", "value"}:
            return 0
        if text in {"1", "z1", "derivative", "ordinary_z1"}:
            return 1
    return None


def _rows(c: Any, raw: Any, label: str, default_id: str) -> dict[int, list[dict[str, Any]]]:
    if isinstance(raw, Mapping):
        for wrapper in ("ordinary_Z", "source_Z", "source_jets", "jets", "rows", "incoming", "state"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)):
                raw = raw[wrapper]
                break
    if not isinstance(raw, Mapping):
        raise ValueError(f"{label} requires explicit ordinary Z0 and Z1 rows")
    found: dict[int, Any] = {}
    for key, value in raw.items():
        index = _order(key)
        if index is not None:
            if index in found:
                raise ValueError(f"{label} supplies Z{index} twice")
            found[index] = value
    if set(found) != {0, 1}:
        # A two-entry sequence is accepted only after it has been wrapped in a
        # mapping; this guard makes missing ordinary-Z rows explicit.
        raise ValueError(f"{label} requires explicit ordinary Z0 and Z1 rows")
    return {
        index: _term_list(c, found[index], f"{label}.Z{index}", f"{default_id}:Z{index}")
        for index in ORDINARY_Z_ORDERS
    }


def _lookup(mapping: Mapping[Any, Any], names: Sequence[str]) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    lowered = {str(key).lower(): value for key, value in mapping.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return None


def _packet_wrapper(raw: Any) -> Any:
    if not isinstance(raw, Mapping):
        return raw
    for key in ("state", "six_state", "components", "fields", "incoming", "inlet"):
        if key in raw and isinstance(raw[key], Mapping):
            return raw[key]
    return raw


def _packet(c: Any, raw: Any, label: str) -> tuple[dict[str, dict[int, list[dict[str, Any]]]], bool, dict[str, Any]]:
    outer = raw
    raw = _packet_wrapper(raw)
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a six-component mapping")
    result: dict[str, dict[int, list[dict[str, Any]]]] = {}
    for component in FIELDS:
        value = _lookup(raw, _COMPONENT_ALIASES[component])
        if value is None:
            raise ValueError(f"{label} is missing six-state component {component}")
        result[component] = _rows(c, value, f"{label}.{component}", f"{label}:{component}")
    metadata: dict[str, Any] = {}
    for container in (outer, raw):
        if isinstance(container, Mapping) and isinstance(container.get("metadata"), Mapping):
            metadata.update(dict(container["metadata"]))
    skip_values = []
    for container in (outer, raw):
        if isinstance(container, Mapping):
            skip_values.extend(
                value
                for value in (
                    _lookup(container, ("skip_K0", "skip_k0", "skip_k0_transport")),
                    container.get("metadata", {}).get("skip_K0") if isinstance(container.get("metadata"), Mapping) else None,
                    container.get("metadata", {}).get("skip_k0") if isinstance(container.get("metadata"), Mapping) else None,
                )
                if value is not None
            )
    skip = bool(
        any(skip_values)
    )
    metadata["skip_K0"] = skip
    return result, skip, metadata


def _looks_like_packet(raw: Any) -> bool:
    if not isinstance(raw, Mapping):
        return False
    keys = {str(key).lower() for key in raw}
    aliases = {alias.lower() for values in _COMPONENT_ALIASES.values() for alias in values}
    if any(key in keys for key in ("state", "six_state", "components", "fields")):
        return True
    # A top-level group is allowed to be named ``P`` or ``F``.  Distinguish
    # that from a direct six-state packet by checking that the alias value is
    # itself a row mapping (rather than another packet mapping).
    for key, value in raw.items():
        if str(key).lower() not in aliases:
            continue
        if isinstance(value, Mapping):
            row_keys = {_order(item) for item in value}
            if {0, 1}.issubset(row_keys):
                return True
        elif _is_sequence(value):
            return True
    return False


def _groups(
    c: Any,
    incoming: Any,
    *,
    inlet_groups: Any,
    groups: Any,
    incoming_groups: Any,
    direct: Mapping[str, Any],
) -> dict[str, tuple[dict[str, dict[int, list[dict[str, Any]]]], bool, dict[str, Any]]]:
    if any(value is not None for value in direct.values()) and any(
        value is not None for value in (inlet_groups, groups, incoming_groups, incoming)
    ):
        raise ValueError("direct six-component aliases cannot be mixed with grouped inlet input")
    packets: list[tuple[Any, Any]] = []

    def add_mapping(raw: Any, label: str) -> None:
        if raw is None:
            return
        if not isinstance(raw, Mapping):
            raise TypeError(f"{label} must be a mapping")
        nested = None
        for wrapper in ("groups", "inlet_groups", "incoming_groups", "packets"):
            candidate = raw.get(wrapper)
            if not isinstance(candidate, Mapping) or _looks_like_packet(candidate):
                continue
            if any(_looks_like_packet(packet) for packet in candidate.values()):
                nested = candidate
                break
        if nested is not None:
            if not isinstance(nested, Mapping):
                raise TypeError(f"{label}.groups must be a mapping")
            for name, packet in nested.items():
                packets.append((name, packet))
        elif _looks_like_packet(raw):
            packets.append(("inlet", raw))
        else:
            for name, packet in raw.items():
                packets.append((name, packet))

    add_mapping(inlet_groups, "inlet_groups")
    add_mapping(groups, "groups")
    add_mapping(incoming_groups, "incoming_groups")
    add_mapping(incoming, "incoming")
    if not packets and any(value is not None for value in direct.values()):
        packet = {key: value for key, value in direct.items() if value is not None}
        packets.append(("inlet", packet))
    if not packets:
        raise ValueError("at least one caller-owned six-component inlet group is required")

    result: dict[str, tuple[dict[str, dict[int, list[dict[str, Any]]]], bool, dict[str, Any]]] = {}
    for name, packet in packets:
        if not isinstance(name, str) or not name.strip():
            raise ValueError("inlet group names must be nonempty strings")
        if name in result:
            raise ValueError(f"duplicate inlet group {name!r}")
        result[name] = _packet(c, packet, f"group[{name}]")
    return result


def _geometry(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError("explicit geometry certificate is required")
    evidence = value.get("geometry", value)
    if not isinstance(evidence, Mapping):
        raise ValueError("geometry certificate must be a mapping")
    hb_sc = (
        evidence.get("hb_sc_Z_independent") is True
        or (evidence.get("hb_Z_independent") is True and evidence.get("sc_Z_independent") is True)
        or (evidence.get("sc_phase_Z_independent") is True and evidence.get("hb_epsilon_Z_independent") is True)
    )
    rho = (
        evidence.get("rho_radius_Z_independent") is True
        or (evidence.get("rho_Z_independent") is True and evidence.get("radius_Z_independent") is True)
        or evidence.get("rho_radius_fixed") is True
    )
    strict = any(
        evidence.get(key) is True or value.get(key) is True
        for key in ("strict_hb_admitted", "strict_hb_admissibility", "hb_strictly_between_zero_and_half", "strict_hb")
    )
    proof = evidence.get("strict_hb_admission_proof", evidence.get("strict_hb_proof", value.get("strict_hb_admission_proof")))
    if proof is True:
        strict = True
    elif isinstance(proof, Mapping):
        strict = strict or proof.get("strict_hb_admitted") is True or (proof.get("lower_strict") is True and proof.get("upper_strict") is True)
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


def _chi(c: Any, chi: Any, chi_Z: Any) -> dict[int, dict[str, Any]]:
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
        rows = {
            0: _one_term(c, z0, "chi.Z0", "chi:Z0"),
            1: _one_term(c, chi_Z, "chi.Z1", "chi:Z1"),
        }
    else:
        if not isinstance(raw, Mapping):
            raise ValueError("chi requires explicit ordinary Z0 and Z1 rows")
        found: dict[int, Any] = {}
        for key, value in raw.items():
            index = _order(key)
            if index is not None:
                if index in found:
                    raise ValueError(f"chi supplies Z{index} twice")
                found[index] = value
        if set(found) != {0, 1}:
            raise ValueError("chi requires explicit ordinary Z0 and Z1 rows")
        rows = {
            index: _one_term(c, found[index], f"chi.Z{index}", f"chi:Z{index}")
            for index in ORDINARY_Z_ORDERS
        }
    for order, term in rows.items():
        if term is None:
            raise ValueError(f"chi Z{order} cannot be structurally zero")
        if tuple(term["powers"]) != (0, 0, 0, 0):
            raise ValueError(f"chi Z{order} must have zero formal powers")
        if not term.get("source_ids"):
            raise ValueError(f"chi Z{order} must retain a source ID")
    return rows


def inlet_kernel_masses(c: Any, *, hb: Any, sc: Any, strict_hb_admitted: bool = False) -> dict[str, Any]:
    """Return the four first-potential coefficient masses.

    ``hb`` is validated only; the returned masses contain no evaluated hb.
    """

    _strict_hb(c, hb, allow_closed=bool(strict_hb_admitted))
    sc = _positive(c, sc, "sc")
    phase, A = sc / 8, sc / 16
    result = {
        "F_from_Fc": _nonnegative(c, 2 * phase**2 * _expbox(c, -2 * A, 2 * A), "F<-Fc mass"),
        "F_from_Fxc": _nonnegative(c, c.mpf(4) / 3 * phase**3 * _expbox(c, -4 * A, 2 * A), "F<-Fxc mass"),
        "Fx_from_Fc": _nonnegative(c, 2 * phase * _expbox(c, -3 * A, A), "Fx<-Fc mass"),
        "Fx_from_Fxc": _nonnegative(c, 2 * phase**2 * _expbox(c, -5 * A, A), "Fx<-Fxc mass"),
    }
    result.update(
        {
            "Fc_to_F": result["F_from_Fc"],
            "Fxc_to_F": result["F_from_Fxc"],
            "Fc_to_Fx": result["Fx_from_Fc"],
            "Fxc_to_Fx": result["Fx_from_Fxc"],
        }
    )
    return result


kernel_masses = inlet_kernel_masses


def kernel_bounds(c: Any, *, sc: Any, chi: Any, chi_Z: Any) -> dict[str, Any]:
    """Return the source-free nonnegative majorants used by ``j>=2``."""

    sc = _positive(c, sc, "sc")
    C0 = _abs_upper(c, chi["coefficient"])
    C1 = _abs_upper(c, chi_Z["coefficient"])
    phase, A = sc / 8, sc / 16
    E = _nonnegative(c, _expbox(c, A, A), "E")
    k0 = _nonnegative(c, 8 * phase**2 * E**2 * C0, "k0")
    kz = _nonnegative(c, 8 * phase**2 * E**2 * C1, "kz")
    kx = _nonnegative(c, 4 * phase * E * C0, "kx")
    kxz = _nonnegative(c, 4 * phase * E * C1, "kxz")
    T = _nonnegative(c, _expbox(c, k0 / 4, k0 / 4), "T")
    return {
        "c": phase,
        "A": A,
        "E": E,
        "C0": C0,
        "C1": C1,
        "k0": k0,
        "kz": kz,
        "kx": kx,
        "kxz": kxz,
        "T": T,
        "formulas": {
            "c": "sc/8",
            "A": "sc/16",
            "E": "exp(A)",
            "C0": "sup|chi|",
            "C1": "sup|chi_Z|",
            "k0": "8*c**2*exp(2*A)*sup|chi|",
            "kz": "8*c**2*exp(2*A)*sup|chi_Z|",
            "kx": "4*c*exp(A)*sup|chi|",
            "kxz": "4*c*exp(A)*sup|chi_Z|",
            "T": "exp(k0/4)",
        },
    }


def _powers_with_hb(term: Mapping[str, Any], added: int) -> tuple[int, int, int, int]:
    base = tuple(term["powers"])
    return (base[0] + added, base[1], base[2], base[3])


def _derived_term(
    c: Any,
    source: Mapping[str, Any],
    factor: Any,
    *,
    added: int,
    group: str,
    operation: str,
    branch: str,
    field: str,
    source_component: str,
    error_only_override: bool | None = None,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    source_meta = _meta(source)
    error_only = bool(source_meta["error_only"] or source.get("symmetric_error"))
    if error_only_override is not None:
        error_only = bool(error_only_override)
    coefficient = source["coefficient"] * factor
    if error_only:
        coefficient = _symmetric(c, coefficient)
    classes = list(source_meta["error_classes"])
    if error_only and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    result = {
        "powers": _powers_with_hb(source, added),
        "coefficient": coefficient,
        "source_ids": list(source_meta["source_ids"]),
        "incoming_source_ids": list(source_meta["source_ids"]),
        "input_ref": {"group": group, "source_ids": list(source_meta["source_ids"])},
        "error_classes": _unique(classes),
        "uncertainty": bool(source_meta["uncertainty"] or error_only),
        "origin": {
            "operation": operation,
            "group": group,
            "branch": branch,
            "field": field,
            "source_component": source_component,
            "source_ids": list(source_meta["source_ids"]),
            "error_classes": _unique(classes),
            "uncertainty": bool(source_meta["uncertainty"] or error_only),
            "error_only": error_only,
            "symmetric_error": error_only,
            "formal_hb_power_added": added,
            "coefficient_contains_formal_hb": False,
            "source_origin": _origin_descriptor(source.get("origin")),
        },
        "group": group,
        "branch": branch,
        "field": field,
        "source_component": source_component,
        "incoming_component": source_component,
        "formal_hb_power_added": added,
        "coefficient_contains_formal_hb": False,
        "signed_source_term": not error_only,
        "error_only": error_only,
        "symmetric_error": error_only,
    }
    if isinstance(source.get("origin"), str):
        result["origin_texts"] = [source["origin"]]
    else:
        result["origin_texts"] = []
    if extra:
        result.update(dict(extra))
        origin = result.get("origin")
        if isinstance(origin, Mapping):
            merged = dict(origin)
            merged.update(dict(extra.get("origin", {})) if isinstance(extra.get("origin"), Mapping) else {})
            result["origin"] = merged
    return result


def _k0_term(
    c: Any,
    source: Mapping[str, Any],
    factor: Any,
    *,
    added: int,
    group: str,
    branch: str,
    field: str,
    source_component: str,
) -> dict[str, Any]:
    return _derived_term(
        c,
        source,
        factor,
        added=added,
        group=group,
        operation="K0_homogeneous_inlet_transport_v2",
        branch=branch,
        field=field,
        source_component=source_component,
        extra={"stage": "K0", "baseline_K0": True, "baseline_K0_included": True},
    )


def _empty_fields() -> dict[str, dict[str, list[Any]]]:
    return {field: {"Z0": [], "Z1": []} for field in FIELDS}


def _empty_branches(names: Sequence[str]) -> dict[str, dict[str, Any]]:
    return {
        name: {"branch": name, "Z0": [], "Z1": [], "response": {"Z0": [], "Z1": []}}
        for name in names
    }


def _aggregate_fields(branches: Mapping[str, Mapping[str, Any]], fields: Sequence[str]) -> dict[str, dict[str, list[Any]]]:
    result = {field: {"Z0": [], "Z1": []} for field in FIELDS}
    for branch in branches.values():
        field = branch.get("field")
        if field not in fields:
            continue
        for order in ORDINARY_Z_ORDERS:
            result[field][f"Z{order}"].extend(list(branch.get(f"Z{order}", [])))
    return result


def _k0_stage(c: Any, packet: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]], *, group: str, sc: Any, skip: bool) -> dict[str, Any]:
    branches = _empty_branches(K0_BRANCHES)
    field_for_branch = {
        "F_from_Fc": ("F", "F", "Fc", 0, "1"),
        "F_from_Fxc": ("F", "Fx", "Fxc", 1, "hb*(sc/4)*phi(d)"),
        "U_from_Uc": ("U", "U", "Uc", 0, "1"),
        "U_from_Uxc": ("U", "Ux", "Uxc", 1, "d=hb*(sc/4)"),
        "K_from_Kc": ("K", "K", "Kc", 0, "exp(-d)"),
        "K_from_Uxc": ("K", "Ux", "Uxc", 1, "-hb*(sc/4)*phi(d)"),
        "P_from_Pc": ("P", "P", "Pc", 0, "1"),
        "Fx_from_Fxc": ("Fx", "Fx", "Fxc", 0, "exp(-3*d/2)"),
        "Ux_from_Uxc": ("Ux", "Ux", "Uxc", 0, "exp(-d/2)"),
    }
    decay = _decay_bound(c, sc, c.mpf(1))
    decay_3 = _decay_bound(c, sc, c.mpf(3) / 2)
    decay_half = _decay_bound(c, sc, c.mpf(1) / 2)
    phi = decay
    b_factor = sc / 4 * phi
    minus_b_factor = -b_factor
    d_factor = sc / 4
    factors = {
        "F_from_Fc": (_zero(c) + 1, 0),
        "F_from_Fxc": (b_factor, 1),
        "U_from_Uc": (_zero(c) + 1, 0),
        "U_from_Uxc": (d_factor, 1),
        "K_from_Kc": (decay, 0),
        "K_from_Uxc": (minus_b_factor, 1),
        "P_from_Pc": (_zero(c) + 1, 0),
        "Fx_from_Fxc": (decay_3, 0),
        "Ux_from_Uxc": (decay_half, 0),
    }
    for branch, (field, source_field, source_component, delta, formula) in field_for_branch.items():
        factor, added = factors[branch]
        branches[branch].update(
            {
                "field": field,
                "source_component": source_component,
                "formal_hb_power_added": added,
                "factor": factor,
                "factor_formula": formula,
                "signed_sectors": not skip,
                "error_only_sectors_separate": True,
                "baseline_K0_included": not skip,
            }
        )
        if skip:
            continue
        for order in ORDINARY_Z_ORDERS:
            for source in packet[source_field][order]:
                output = _k0_term(
                    c,
                    source,
                    factor,
                    added=added,
                    group=group,
                    branch=branch,
                    field=field,
                    source_component=source_component,
                )
                output["source_order"] = order
                output["incoming_order"] = order
                output["stage"] = "K0"
                if isinstance(output.get("origin"), Mapping):
                    output["origin"].update({"source_order": order, "incoming_order": order, "stage": "K0"})
                branches[branch][f"Z{order}"].append(output)
                branches[branch]["response"][f"Z{order}"].append(output)
    fields = _aggregate_fields(branches, FIELDS)
    return {
        "stage": "K0",
        "fields": fields,
        "response": fields,
        "branches": branches,
        "skip_K0": bool(skip),
        "baseline_K0_included": not skip,
        "baseline_K0_recomputed": not skip,
        "structural_zero_fields": list(FIELDS) if skip else [],
        "homogeneous": True,
        "formal_hb_not_materialized": True,
        "K0_formula": {
            "F": "Fc + (1-exp(-d))*Fxc",
            "Fx": "exp(-3*d/2)*Fxc",
            "U": "Uc + d*Uxc",
            "K": "exp(-d)*Kc - (1-exp(-d))*Uxc",
            "P": "Pc",
            "Ux": "exp(-d/2)*Uxc",
            "d": "hb*sc/4",
        },
        "factor_intervals": {
            "phi": phi,
            "B_over_hb": sc / 4 * phi,
            "exp_minus_d": decay,
            "exp_minus_3d_over_2": decay_3,
            "exp_minus_d_over_2": decay_half,
        },
    }


def _first_term(
    c: Any,
    chi: Mapping[str, Any],
    source: Mapping[str, Any],
    *,
    branch: str,
    field: str,
    source_component: str,
    chi_order: int,
    source_order: int,
    mass: Any,
    group: str,
) -> dict[str, Any]:
    coefficient = -2 * chi["coefficient"] * source["coefficient"] * mass
    error_only = bool(_meta(chi)["error_only"] or _meta(source)["error_only"] or chi.get("symmetric_error") or source.get("symmetric_error"))
    if error_only:
        coefficient = _symmetric(c, coefficient)
    ids = _unique(list(_meta(chi)["source_ids"]) + list(_meta(source)["source_ids"]))
    classes = _unique(list(_meta(chi)["error_classes"]) + list(_meta(source)["error_classes"]))
    if error_only and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    piece = "chi_Z*inlet_Z0" if (chi_order, source_order) == (1, 0) else "chi*inlet_Z1" if (chi_order, source_order) == (0, 1) else "chi*inlet_Z0"
    added = FIRST_HB_POWER_ADDED[branch]
    origin = {
        "operation": "first_Kchi_homogeneous_inlet_v2",
        "group": group,
        "branch": branch,
        "field": field,
        "incoming_component": source_component,
        "chi_order": chi_order,
        "incoming_order": source_order,
        "leibniz_piece": piece,
        "source_ids": ids,
        "error_classes": classes,
        "uncertainty": bool(_meta(chi)["uncertainty"] or _meta(source)["uncertainty"] or error_only),
        "error_only": error_only,
        "symmetric_error": error_only,
        "formal_hb_power_added": added,
        "coefficient_contains_formal_hb": False,
        "baseline_K0_included": False,
        "chi_origin": _origin_descriptor(chi.get("origin")),
        "incoming_origin": _origin_descriptor(source.get("origin")),
    }
    return {
        "powers": _powers_with_hb(source, added),
        "coefficient": coefficient,
        "source_ids": ids,
        "input_ref": {"group": group, "source_ids": list(_meta(source)["source_ids"])},
        "chi_source_ids": list(_meta(chi)["source_ids"]),
        "incoming_source_ids": list(_meta(source)["source_ids"]),
        "error_classes": classes,
        "uncertainty": origin["uncertainty"],
        "origin": origin,
        "origin_texts": [text for text in (chi.get("origin_text"), source.get("origin_text")) if isinstance(text, str)],
        "group": group,
        "branch": branch,
        "field": field,
        "incoming_component": source_component,
        "chi_order": chi_order,
        "incoming_order": source_order,
        "source_order": source_order,
        "leibniz_piece": piece,
        "formal_hb_power_added": added,
        "kernel_mass": mass,
        "kernel_mass_is_coefficient_only": True,
        "coefficient_contains_formal_hb": False,
        "physical_dx_count": 1,
        "signed_source_term": not error_only,
        "error_only": error_only,
        "symmetric_error": error_only,
        "baseline_K0_included": False,
        "j1_included": True,
    }


def _first_stage(c: Any, packet: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]], chi: Mapping[int, Mapping[str, Any]], masses: Mapping[str, Any], *, group: str) -> dict[str, Any]:
    branches = _empty_branches(BRANCHES)
    branch_spec = {
        "F_from_Fc": ("F", "F", "Fc"),
        "F_from_Fxc": ("F", "Fx", "Fxc"),
        "Fx_from_Fc": ("Fx", "F", "Fc"),
        "Fx_from_Fxc": ("Fx", "Fx", "Fxc"),
    }
    for branch, (field, source_field, source_component) in branch_spec.items():
        branch_result = branches[branch]
        branch_result.update(
            {
                "field": field,
                "incoming_component": source_component,
                "operator": "-2*chi*K0*H0",
                "kernel_mass": masses[branch],
                "kernel_mass_formula": MASS_FORMULAS[branch],
                "formal_hb_power_added": FIRST_HB_POWER_ADDED[branch],
                "signed_sectors": True,
                "error_only_sectors_separate": True,
                "baseline_K0_included": False,
            }
        )
        for chi_order, source_order, output_order in ((0, 0, 0), (1, 0, 1), (0, 1, 1)):
            for source in packet[source_field][source_order]:
                output = _first_term(
                    c,
                    chi[chi_order],
                    source,
                    branch=branch,
                    field=field,
                    source_component=source_component,
                    chi_order=chi_order,
                    source_order=source_order,
                    mass=masses[branch],
                    group=group,
                )
                branch_result[f"Z{output_order}"].append(output)
                branch_result["response"][f"Z{output_order}"].append(output)
    fields = _aggregate_fields(branches, POTENTIAL_FIELDS)
    for field in ("U", "K", "P", "Ux"):
        fields[field] = {"Z0": [], "Z1": []}
    return {
        "stage": "first_potential",
        "fields": fields,
        "response": fields,
        "branches": branches,
        "structural_zero_fields": ["U", "K", "P", "Ux"],
        "first_potential_only": True,
        "baseline_K0_included": False,
        "j2_included": False,
        "ordinary_Z1_product_rule": "chi_Z*inlet_Z0 + chi*inlet_Z1",
    }


def _tail_term(
    c: Any,
    source: Mapping[str, Any],
    chi_terms: Sequence[Mapping[str, Any]],
    radius_factor: Any,
    *,
    branch: str,
    field: str,
    source_component: str,
    source_order: int,
    chi_order: int,
    group: str,
    radius_formula: str,
) -> dict[str, Any]:
    radius = _nonnegative(c, _abs_upper(c, source["coefficient"]) * radius_factor, "j>=2 tail radius")
    coefficient = _symmetric(c, radius)
    source_ids = list(_meta(source)["source_ids"])
    chi_ids = _unique([source_id for term in chi_terms for source_id in _meta(term)["source_ids"]])
    ids = _unique(source_ids + chi_ids)
    classes = _unique(list(_meta(source)["error_classes"]) + [value for term in chi_terms for value in _meta(term)["error_classes"]])
    for required in ("j2_kernel_remainder", "symmetric_error"):
        if required not in classes:
            classes.append(required)
    piece = "chi_Z*inlet_Z0" if (chi_order, source_order) == (1, 0) else "chi*inlet_Z1" if (chi_order, source_order) == (0, 1) else "chi*inlet_Z0"
    added = J2_HB_POWER_ADDED[branch]
    origin = {
        "operation": "j>=2_Kchi_homogeneous_inlet_tail_v2",
        "group": group,
        "branch": branch,
        "field": field,
        "incoming_component": source_component,
        "chi_order": chi_order,
        "incoming_order": source_order,
        "source_order": source_order,
        "leibniz_piece": piece,
        "source_ids": ids,
        "input_ref": {"group": group, "source_ids": source_ids},
        "incoming_source_ids": source_ids,
        "chi_source_ids": chi_ids,
        "error_classes": classes,
        "uncertainty": True,
        "error_only": True,
        "symmetric_error": True,
        "tail_starts_at_j": 2,
        "formal_hb_power_added": added,
        "coefficient_contains_formal_hb": False,
        "radius_formula": radius_formula,
        "chi_origins": [_origin_descriptor(term.get("origin")) for term in chi_terms],
        "incoming_origin": _origin_descriptor(source.get("origin")),
    }
    return {
        "powers": _powers_with_hb(source, added),
        "coefficient": coefficient,
        "radius": radius,
        "radius_factor": radius_factor,
        "coefficient_is_symmetric_radius": True,
        "source_ids": ids,
        "incoming_source_ids": source_ids,
        "chi_source_ids": chi_ids,
        "chi_source_classes": _unique(value for term in chi_terms for value in _meta(term)["error_classes"]),
        "error_classes": classes,
        "uncertainty": True,
        "origin": origin,
        "origin_texts": [
            text
            for term in (*chi_terms, source)
            for text in (term.get("origin_text"),)
            if isinstance(text, str)
        ],
        "group": group,
        "branch": branch,
        "field": field,
        "incoming_component": source_component,
        "source_order": source_order,
        "incoming_order": source_order,
        "chi_order": chi_order,
        "leibniz_piece": piece,
        "tail_starts_at_j": 2,
        "kernel_tail_order": "j>=2",
        "formal_hb_power_added": added,
        "coefficient_contains_formal_hb": False,
        "physical_dx_count_per_word": ">=2",
        "signed_source_term": False,
        "error_only": True,
        "symmetric_error": True,
        "baseline_K0_included": False,
        "j1_included": False,
    }


def _j2_stage(c: Any, packet: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]], chi: Mapping[int, Mapping[str, Any]], bounds: Mapping[str, Any], *, group: str, sc: Any) -> dict[str, Any]:
    branches = _empty_branches(BRANCHES)
    k0, kz, kx, kxz, T = bounds["k0"], bounds["kz"], bounds["kx"], bounds["kxz"], bounds["T"]
    base_f = k0 * k0 * T / 2
    base_fz = kz * k0 * T
    base_fx = kx * k0 * T / 2
    base_fxz = (kxz * k0 + kx * kz) * T / 2
    specs = {
        "F_from_Fc": ("F", "F", "Fc", base_f, base_fz, "B*hb^4*k0^2*T/2"),
        "F_from_Fxc": ("F", "Fx", "Fxc", sc / 4 * base_f, sc / 4 * base_fz, "B*hb^4*k0^2*T/2 with B_Fxc=(sc/4)|Fxc|"),
        "Fx_from_Fc": ("Fx", "F", "Fc", base_fx, base_fxz, "B*hb^3*kx*k0*T/2"),
        "Fx_from_Fxc": ("Fx", "Fx", "Fxc", sc / 4 * base_fx, sc / 4 * base_fxz, "B*hb^3*kx*k0*T/2 with B_Fxc=(sc/4)|Fxc|"),
    }
    for branch, (field, source_field, source_component, z0_factor, z1_chi_factor, formula) in specs.items():
        branch_result = branches[branch]
        branch_result.update(
            {
                "field": field,
                "incoming_component": source_component,
                "kernel_tail_order": "j>=2",
                "tail_starts_at_j": 2,
                "formal_hb_power_added": J2_HB_POWER_ADDED[branch],
                "signed_sectors": False,
                "all_terms_error_only": True,
                "baseline_K0_included": False,
                "j1_included": False,
                "radius_factors": {"Z0": z0_factor, "Z1_source": z0_factor, "Z1_chi": z1_chi_factor},
            }
        )
        # Z0: chi * inlet Z0.
        for source in packet[source_field][0]:
            output = _tail_term(
                c,
                source,
                (chi[0],),
                z0_factor,
                branch=branch,
                field=field,
                source_component=source_component,
                source_order=0,
                chi_order=0,
                group=group,
                radius_formula=formula,
            )
            branch_result["Z0"].append(output)
            branch_result["response"]["Z0"].append(output)
        # Z1: chi * inlet Z1.
        for source in packet[source_field][1]:
            output = _tail_term(
                c,
                source,
                (chi[0],),
                z0_factor,
                branch=branch,
                field=field,
                source_component=source_component,
                source_order=1,
                chi_order=0,
                group=group,
                radius_formula=formula.replace("B*", "BZ*") if "B*" in formula else formula,
            )
            branch_result["Z1"].append(output)
            branch_result["response"]["Z1"].append(output)
        # Z1: chi_Z * inlet Z0.  Keep both chi source IDs in provenance.
        for source in packet[source_field][0]:
            output = _tail_term(
                c,
                source,
                (chi[0], chi[1]),
                z1_chi_factor,
                branch=branch,
                field=field,
                source_component=source_component,
                source_order=0,
                chi_order=1,
                group=group,
                radius_formula=("B*" + ("hb^4*(kz*k0*T)" if field == "F" else "hb^3*(kxz*k0+kx*kz)*T/2")),
            )
            branch_result["Z1"].append(output)
            branch_result["response"]["Z1"].append(output)
    fields = _aggregate_fields(branches, POTENTIAL_FIELDS)
    for field in ("U", "K", "P", "Ux"):
        fields[field] = {"Z0": [], "Z1": []}
    return {
        "stage": "j2_tail",
        "fields": fields,
        "response": fields,
        "branches": branches,
        "structural_zero_fields": ["U", "K", "P", "Ux"],
        "j2_kernel_remainder_included": True,
        "j2_remainder_included": True,
        "symmetric_error_outputs": True,
        "baseline_K0_included": False,
        "j1_included": False,
        "kernel_bounds": bounds,
        "component_error_radii": {
            "F": {
                "B_formula": "|Fc| + hb*(sc/4)*|Fxc|",
                "BZ_formula": "|Fc_Z| + hb*(sc/4)*|Fxc_Z|",
                "R_F": "B*hb^4*k0^2*T/2",
                "R_FZ": "BZ*hb^4*k0^2*T/2 + B*hb^4*kz*k0*T",
            },
            "Fx": {
                "B_formula": "|Fc| + hb*(sc/4)*|Fxc|",
                "BZ_formula": "|Fc_Z| + hb*(sc/4)*|Fxc_Z|",
                "R_Fx": "B*hb^3*kx*k0*T/2",
                "R_FxZ": "BZ*hb^3*kx*k0*T/2 + B*hb^3*(kxz*k0+kx*kz)*T/2",
            },
        },
        "ordinary_Z1_product_rule": "chi_Z*inlet_Z0 + chi*inlet_Z1",
    }


def _group_result(
    c: Any,
    packet: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]],
    *,
    group: str,
    skip_k0: bool,
    sc: Any,
    chi: Mapping[int, Mapping[str, Any]],
    masses: Mapping[str, Any],
    tail_bounds: Mapping[str, Any],
    metadata: Mapping[str, Any],
) -> dict[str, Any]:
    k0 = _k0_stage(c, packet, group=group, sc=sc, skip=skip_k0)
    first = _first_stage(c, packet, chi, masses, group=group)
    tail = _j2_stage(c, packet, chi, tail_bounds, group=group, sc=sc)
    return {
        "group": group,
        "metadata": dict(metadata),
        "skip_K0": bool(skip_k0),
        "k0": k0,
        "first_potential": first,
        "j2_tail": tail,
        "inlet": {
            field: {
                f"Z{order}": list(packet[field][order])
                for order in ORDINARY_Z_ORDERS
            }
            for field in FIELDS
        },
        "input_error_only": any(
            term.get("error_only") is True
            for field in FIELDS
            for order in ORDINARY_Z_ORDERS
            for term in packet[field][order]
        ),
        "formal_powers_preserved": True,
        "structural_zero_fields": {
            "k0": list(FIELDS) if skip_k0 else [],
            "first_potential": ["U", "K", "P", "Ux"],
            "j2_tail": ["U", "K", "P", "Ux"],
        },
    }


def build(
    c: Any,
    incoming: Any = None,
    chi: Any = None,
    *,
    hb: Any,
    sc: Any,
    inlet_groups: Any = None,
    groups: Any = None,
    incoming_groups: Any = None,
    packets: Any = None,
    chi_Z: Any = None,
    chi_jets: Any = None,
    source_family: Any = None,
    source_atoms: Any = None,
    geometry_certificate: Any = None,
    certificate: Any = None,
    q_star: Any = 1,
    q: Any = None,
    F: Any = None,
    U: Any = None,
    K: Any = None,
    P: Any = None,
    Fx: Any = None,
    Ux: Any = None,
    Fc: Any = None,
    Uc: Any = None,
    Kc: Any = None,
    Pc: Any = None,
    Fxc: Any = None,
    Uxc: Any = None,
    **unexpected: Any,
) -> dict[str, Any]:
    """Build endpoint sectors for arbitrary caller-owned inlet groups.

Each group must provide explicit ``Z0`` and ``Z1`` rows for every field in
:data:`FIELDS`; each row is a finite list of sectors, with ``[]`` reserved for
an explicit structural zero.  Group packets may set ``skip_K0`` to retain an
    already-entered baseline while still receiving the first-potential and
    j>=2 outputs.  ``chi`` has exactly one zero-formal-power term per ordinary-Z
    order.  The endpoint is deliberately restricted to ``q=1``.
    """

    if not all(hasattr(c, key) for key in ("mpf", "make_mpf", "exp")):
        raise TypeError("c must be a directed interval context")
    if unexpected:
        raise TypeError("unexpected transport options: " + ", ".join(sorted(unexpected)))
    direct_pairs = (
        ("F", F, "Fc", Fc),
        ("U", U, "Uc", Uc),
        ("K", K, "Kc", Kc),
        ("P", P, "Pc", Pc),
        ("Fx", Fx, "Fxc", Fxc),
        ("Ux", Ux, "Uxc", Uxc),
    )
    for canonical, canonical_value, alias, alias_value in direct_pairs:
        if canonical_value is not None and alias_value is not None:
            raise ValueError(f"direct aliases {canonical} and {alias} were supplied twice")
    direct = {
        canonical: canonical_value if canonical_value is not None else alias_value
        for canonical, canonical_value, _, alias_value in direct_pairs
    }
    if any(value is not None for value in direct.values()) and any(
        value is not None for value in (incoming, inlet_groups, groups, incoming_groups, packets)
    ):
        raise ValueError("direct six-component aliases cannot be mixed with grouped inlet input")
    if chi is not None and chi_jets is not None:
        raise ValueError("chi and chi_jets were supplied twice")
    chi = chi if chi is not None else chi_jets
    if geometry_certificate is not None and certificate is not None and geometry_certificate != certificate:
        raise ValueError("geometry_certificate and certificate disagree")
    cert = _geometry(geometry_certificate if geometry_certificate is not None else certificate)
    atoms = _atoms(source_atoms)
    if q is not None:
        if q_star != 1:
            raise ValueError("q and q_star were supplied twice")
        q_star = q
    q_interval = _iv(c, q_star, "q_star")
    if _lo(c, q_interval) != 1 or _hi(c, q_interval) != 1:
        raise ValueError("q_star must be exact endpoint 1")
    hb = _strict_hb(c, hb, allow_closed=cert["strict_hb_admitted"])
    sc = _positive(c, sc, "sc")
    chi_rows = _chi(c, chi, chi_Z)
    masses = inlet_kernel_masses(c, hb=hb, sc=sc, strict_hb_admitted=cert["strict_hb_admitted"])
    tail_constants = kernel_bounds(c, sc=sc, chi=chi_rows[0], chi_Z=chi_rows[1])
    if packets is not None:
        if groups is not None or inlet_groups is not None or incoming_groups is not None or incoming is not None:
            raise ValueError("packets cannot be combined with another inlet packet argument")
        incoming_groups = packets
    parsed_groups = _groups(c, incoming, inlet_groups=inlet_groups, groups=groups, incoming_groups=incoming_groups, direct=direct)
    grouped = {
        name: _group_result(c, packet, group=name, skip_k0=skip, sc=sc, chi=chi_rows, masses=masses, tail_bounds=tail_constants, metadata=metadata)
        for name, (packet, skip, metadata) in parsed_groups.items()
    }
    names = list(grouped)
    primary = names[0]
    endpoint = c.mpf(2) * c.exp(hb * sc / 8)
    d_interval = hb * sc / 4
    geometry = {
        "q_star": q_interval,
        "q_interval": "[0,1]",
        "x_endpoint": endpoint,
        "x_formula": "2*exp(hb*sc*q/8)",
        "d": d_interval,
        "d_formula": "hb*sc/4",
        "c": sc / 8,
        "A": sc / 16,
        "hb_interval": hb,
        "sc_interval": sc,
        "admissible_hb": "0<hb<1/2",
        "strict_hb_admitted": True,
        "closed_hb_enclosure_allowed": True,
        "hb_sc_Z_independent": True,
        "rho_radius_Z_independent": True,
        "ordinary_Z_geometry_fixed": True,
        "ordinary_Z_geometry_differentiated": False,
        "formal_hb_not_materialized": True,
        "formal_source_powers_not_materialized": True,
        "midpoint_values_used": False,
        "physical_dt_count": 0,
        "physical_dx_count": 0,
        "phi_interval": _decay_bound(c, sc, c.mpf(1)),
        "phi_formula": "phi(d)=(1-exp(-d))/d in [exp(-sc/8),1]",
    }
    first_map = {name: value["first_potential"] for name, value in grouped.items()}
    tail_map = {name: value["j2_tail"] for name, value in grouped.items()}
    k0_map = {name: value["k0"] for name, value in grouped.items()}
    result = {
        "schema": SCHEMA,
        "status": "conditional_partial_six_component_qp_homogeneous_transport",
        "source_family": source_family,
        "source_atoms": atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "fields": list(FIELDS),
        "inlet_components": list(INLET_COMPONENTS),
        "ordinary_Z_orders": [0, 1],
        "groups": names,
        "inlet_groups": names,
        "group_results": grouped,
        "grouped_inlet_packets_separate": True,
        "grouped_inlet_packets_merged": False,
        "primary_inlet_group": primary,
        "k0_transport": k0_map,
        "first_potential": first_map,
        "j2_tail": tail_map,
        "chi_jets": {f"Z{order}": [chi_rows[order]] for order in ORDINARY_Z_ORDERS},
        "incoming_jets": {
            name: {
                field: {
                    f"Z{order}": list(packet[field][order])
                    for order in ORDINARY_Z_ORDERS
                }
                for field in FIELDS
            }
            for name, (packet, _, _) in parsed_groups.items()
        },
        "geometry": geometry,
        "geometry_certificate": cert,
        "certificate_passthrough": cert["raw"],
        "kernel_bounds": tail_constants,
        "constants": tail_constants,
        "kernel": {
            "K0": "six-component homogeneous transport",
            "first_potential": "-2*chi*K0*H0",
            "j2_tail": "sum_{j>=2}(-2*chi)^j K0",
            "masses": {branch: masses[branch] for branch in BRANCHES},
            "mass_formulas": dict(MASS_FORMULAS),
            "tail_bounds": tail_constants,
            "component_error_radii": {
                "F": {
                    "B_formula": "|Fc| + hb*(sc/4)*|Fxc|",
                    "BZ_formula": "|Fc_Z| + hb*(sc/4)*|Fxc_Z|",
                    "R_F": "B*hb^4*k0^2*T/2",
                    "R_FZ": "BZ*hb^4*k0^2*T/2 + B*hb^4*kz*k0*T",
                },
                "Fx": {
                    "B_formula": "|Fc| + hb*(sc/4)*|Fxc|",
                    "BZ_formula": "|Fc_Z| + hb*(sc/4)*|Fxc_Z|",
                    "R_Fx": "B*hb^3*kx*k0*T/2",
                    "R_FxZ": "BZ*hb^3*kx*k0*T/2 + B*hb^3*(kxz*k0+kx*kz)*T/2",
                },
            },
            "ordinary_Z1_product_rule": "chi_Z*inlet_Z0 + chi*inlet_Z1",
        },
        "ordinary_Z1_product_rule": "chi_Z*inlet_Z0 + chi*inlet_Z1",
        "chi_Z_required": True,
        "chi_source_sector_count_required": 1,
        "chi_formal_powers_required": (0, 0, 0, 0),
        "component_sector_count_required": "finite list per ordinary-Z row",
        "component_sector_rows_finite": True,
        "source_ids_preserved": True,
        "error_classes_preserved": True,
        "formal_powers_preserved": True,
        "formal_hb_not_materialized": True,
        "formal_source_powers_not_materialized": True,
        "signed_coefficients_retained": True,
        "uncertainty_separate": True,
        "error_only_outputs_explicit": True,
        "symmetric_error_outputs": True,
        "first_potential_only": True,
        "j2_kernel_remainder_included": True,
        "j2_remainder_included": True,
        "baseline_K0_included": any(not value["skip_K0"] for value in grouped.values()),
        "baseline_K0_reentered": False,
        "baseline_K0_skipped_groups": [name for name, value in grouped.items() if value["skip_K0"]],
        "K0_j1_j2_disjoint": True,
        "Kchi_higher_words_included": False,
        "Q_feedback_included": False,
        "P_feedback_included": False,
        "complete_QP_feedback": False,
        "allword_tail_added": False,
        "m15_computed": False,
        "complete_collar_field": False,
        "complete_collar": False,
        "complete_inlet": False,
        "finite_physical_field_claimed": False,
        "old_builder_called": False,
        "ancestor_imports": False,
        "receipt_reads": False,
        "automatic_code_execution": False,
    }
    # For a one-group caller, convenient branch aliases expose the exact same
    # unmerged objects without creating a second aggregation path.
    if len(grouped) == 1:
        one = grouped[primary]
        result.update(
            {
                "branches": one["first_potential"]["branches"],
                "response": one["first_potential"]["fields"],
            }
        )
    else:
        result["branches"] = {name: value["first_potential"]["branches"] for name, value in grouped.items()}
        result["response"] = {name: value["first_potential"]["fields"] for name, value in grouped.items()}
    return result


transport = build
build_transport = build
build_qp_transport = build
evaluate = build


