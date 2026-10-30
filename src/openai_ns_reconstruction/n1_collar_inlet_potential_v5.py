"""Independent v5 first ``Kchi`` collar-inlet potential correction.

The leaf emits only the four finite inlet branches at
``X=2*exp(hb*sc/8)``.  It retains formal sectors
``(hb,F0_anchor,Pstar_squared,Lambda)`` and adds hb powers
``F<-Fc:+2, F<-Fxc:+3, Fx<-Fc:+1, Fx<-Fxc:+2``.  Coefficient masses are
``2c**2[exp(-2A),exp(2A)]``, ``(4/3)c**3[exp(-4A),exp(2A)]``,
``2c[exp(-3A),exp(A)]``, and ``2c**2[exp(-5A),exp(A)]`` for
``c=sc/8,A=sc/16``.

Input provenance accepts both mapping and string ``origin`` metadata.  String
origins are preserved verbatim as ``origin``/``origin_text`` while source IDs,
error classes, and uncertainty flags remain in explicit metadata.  The fixed
geometry Z1 rule is ``chi_Z*incoming_Z0 + chi*incoming_Z1``.  A closed hb
enclosure is accepted only with ``strict_hb_admitted`` certificate evidence;
no midpoint, receipt, old leaf, j2 remainder, all-word tail, or m15 term is
constructed.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

SCHEMA = "n1_collar_inlet_potential_v5"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
FIELDS = ("F", "Fx")
INLET_COMPONENTS = ("Fc", "Fxc")
BRANCHES = ("F_from_Fc", "F_from_Fxc", "Fx_from_Fc", "Fx_from_Fxc")
FORMAL_HB_POWER_ADDED = {"F_from_Fc": 2, "F_from_Fxc": 3, "Fx_from_Fc": 1, "Fx_from_Fxc": 2}
MASS_FORMULAS = {"F_from_Fc": "2*c**2*[exp(-2*A),exp(2*A)]", "F_from_Fxc": "(4/3)c**3*[exp(-4A),exp(2A)]", "Fx_from_Fc": "2*c*[exp(-3A),exp(A)]", "Fx_from_Fxc": "2*c**2*[exp(-5A),exp(A)]"}
GROUP_ALIASES = {"G0g": ("G0g", "g0g"), "core_Dg": ("core_Dg", "core_dg", "Dg", "core"), "core_symmetric_tail": ("core_symmetric_tail", "core_symmetric", "symmetric_tail", "tail")}

__all__ = ["SCHEMA", "SOURCE_POWER_BASIS", "ORDINARY_Z_ORDERS", "FIELDS", "INLET_COMPONENTS", "BRANCHES", "FORMAL_HB_POWER_ADDED", "MASS_FORMULAS", "inlet_kernel_masses", "kernel_masses", "build", "build_inlet", "build_inlet_potential", "inlet_potential", "evaluate"]


def _seq(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    return [] if value is None else list(value) if _seq(value) else [value]


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    for value in values:
        if not any(value == old for old in result):
            result.append(value)
    return result


def _is_iv(value: Any) -> bool:
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
    if _is_iv(value):
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
    if _seq(value) and len(value) == 2:
        return _make(c, value[0], value[1])
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{label} is not interval-compatible") from exc


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    value = _iv(c, value)
    lo, hi = value._mpi_
    return c.make_mpf((lo, lo)), c.make_mpf((hi, hi))


def _abs_upper(c: Any, value: Any) -> Any:
    lo, hi = _bounds(c, value)
    return max(-lo, hi)


def _symmetric(c: Any, value: Any) -> Any:
    radius = _abs_upper(c, value)
    return _make(c, -radius, radius)


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lo, hi = value._mpi_
    if c.make_mpf((hi, hi)) < 0:
        raise ValueError(f"{label} is negative")
    if c.make_mpf((lo, lo)) < 0:
        lo = c.mpf(0)._mpi_[0]
    return c.make_mpf((lo, hi))


def _strict_hb(c: Any, value: Any, *, allow_closed: bool) -> Any:
    value = _iv(c, value, "hb")
    lo, hi = _bounds(c, value)
    half = c.mpf(1) / 2
    if lo < 0 or hi > half or (not allow_closed and (lo <= 0 or hi >= half)):
        raise ValueError("hb must enclose 0<hb<1/2")
    return value


def _positive(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lo, _ = _bounds(c, value)
    if lo <= 0:
        raise ValueError(f"{label} must be strictly positive")
    return value


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
    uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False))) or any(token in operation or token in kind for token in ("error", "uncertainty", "symmetric"))
    error_only = bool(value.get("error_only", value.get("symmetric_error", False)))
    for key in ("origin", "metadata", "parents"):
        child = value.get(key)
        if child is None or child is value or not isinstance(child, Mapping) and not _seq(child):
            continue
        for item in (_items(child) if _seq(child) else [child]):
            nested = _meta(item, seen=seen)
            ids.extend(nested["source_ids"])
            classes.extend(nested["error_classes"])
            uncertain = uncertain or nested["uncertainty"]
            error_only = error_only or nested["error_only"]
    return {"source_ids": _unique(ids), "error_classes": _unique(classes), "uncertainty": uncertain, "error_only": error_only}


def _origin_descriptor(value: Any) -> Any:
    if isinstance(value, str):
        return value
    if isinstance(value, Mapping):
        return dict(value)
    return None


def _powers(raw: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(raw, Mapping):
        return (0, 0, 0, 0)
    powers = raw.get("powers")
    if powers is None:
        powers = (raw.get("hb_power", 0), raw.get("F0_power", raw.get("F0_anchor_power", 0)), raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)), raw.get("Lambda_power", raw.get("lambda_power", 0)))
    if not _seq(powers) or len(powers) != 4 or any(isinstance(x, bool) or not isinstance(x, int) for x in powers):
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
    origin_meta = {"source_ids": _unique(ids), "error_classes": _unique(classes), "uncertainty": bool(info["uncertainty"] or error_only)}
    result = dict(item)
    result.update({"powers": _powers(item, label), "coefficient": _iv(c, coefficient, label + ".coefficient"), "source_ids": _unique(ids), "error_classes": _unique(classes), "uncertainty": origin_meta["uncertainty"], "origin_metadata": origin_meta})
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


def _terms(c: Any, raw: Any, label: str, default_id: str) -> list[dict[str, Any]]:
    if raw is None:
        raise TypeError(f"{label} is required")
    if isinstance(raw, Mapping) and any(key in raw for key in ("coefficient", "coefficient_interval", "value", "interval", "powers")):
        values = [raw]
    elif _seq(raw):
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
        if text in {"0", "z0", "value"}:
            return 0
        if text in {"1", "z1", "derivative", "ordinary_z1"}:
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
        return {index: _terms(c, found[index], f"{label}.Z{index}", f"{default_id}:Z{index}") for index in ORDINARY_Z_ORDERS}
    if _seq(raw) and len(raw) == 2 and not all(isinstance(item, Mapping) and any(key in item for key in ("coefficient", "powers", "value", "interval")) for item in raw):
        return {index: _terms(c, raw[index], f"{label}.Z{index}", f"{default_id}:Z{index}") for index in ORDINARY_Z_ORDERS}
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
    return {"Fc": _rows(c, fc, label + ".Fc", label + ":Fc"), "Fxc": _rows(c, fxc, label + ".Fxc", label + ":Fxc")}


def _group_name(value: Any) -> str | None:
    text = str(value).lower().replace("-", "_").replace(" ", "_")
    for canonical, aliases in GROUP_ALIASES.items():
        if text in {str(alias).lower().replace("-", "_").replace(" ", "_") for alias in aliases}:
            return canonical
    return None


def _groups(c: Any, incoming: Any, *, Fc: Any, Fxc: Any, incoming_Fc: Any, incoming_Fxc: Any, inlet_groups: Any, G0g: Any, core_Dg: Any, core_symmetric_tail: Any) -> dict[str, dict[str, dict[int, list[dict[str, Any]]]]]:
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
    if isinstance(chi, Mapping):
        for wrapper in ("ordinary_Z", "chi_jets", "jets", "rows"):
            if wrapper in chi and isinstance(chi[wrapper], (Mapping, Sequence)):
                chi = chi[wrapper]
                break
    if chi_Z is not None:
        if isinstance(chi, Mapping) and any(_order(key) == 1 for key in chi):
            raise ValueError("chi_Z supplied twice")
        z0 = chi if not isinstance(chi, Mapping) or any(key in chi for key in ("coefficient", "value", "interval", "powers")) else _lookup(chi, ("Z0", "z0", "0", "value"))
        if z0 is None:
            raise ValueError("chi Z0 is required")
        return {0: _terms(c, z0, "chi.Z0", "chi:Z0"), 1: _terms(c, chi_Z, "chi.Z1", "chi:Z1")}
    return _rows(c, chi, "chi", "chi")


def _geometry(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError("explicit geometry certificate is required")
    evidence = value.get("geometry", value)
    if not isinstance(evidence, Mapping):
        raise ValueError("geometry certificate must be a mapping")
    hb_sc = evidence.get("hb_sc_Z_independent") is True or (evidence.get("hb_Z_independent") is True and evidence.get("sc_Z_independent") is True)
    rho = evidence.get("rho_radius_Z_independent") is True or (evidence.get("rho_Z_independent") is True and evidence.get("radius_Z_independent") is True)
    strict = evidence.get("strict_hb_admitted") is True or value.get("strict_hb_admitted") is True
    proof = evidence.get("strict_hb_admission_proof", evidence.get("strict_hb_proof", value.get("strict_hb_admission_proof")))
    if proof is True:
        strict = True
    elif isinstance(proof, Mapping):
        strict = strict or proof.get("strict_hb_admitted") is True or (proof.get("lower_strict") is True and proof.get("upper_strict") is True)
    elif isinstance(proof, str):
        strict = strict or proof.replace(" ", "").lower() in {"0<hb<1/2", "0<hb<.5"}
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


def inlet_kernel_masses(c: Any, *, hb: Any, sc: Any, strict_hb_admitted: bool = False) -> dict[str, Any]:
    """Return directed coefficient masses; hb remains a formal sector atom."""

    hb = _strict_hb(c, hb, allow_closed=bool(strict_hb_admitted))
    sc = _positive(c, sc, "sc")
    phase, A = sc / 8, sc / 16

    def expbox(left: Any, right: Any) -> Any:
        return _make(c, c.exp(left).a, c.exp(right).b)

    result = {"F_from_Fc": _nonnegative(c, 2 * phase**2 * expbox(-2 * A, 2 * A), "F<-Fc mass"), "F_from_Fxc": _nonnegative(c, c.mpf(4) / 3 * phase**3 * expbox(-4 * A, 2 * A), "F<-Fxc mass"), "Fx_from_Fc": _nonnegative(c, 2 * phase * expbox(-3 * A, A), "Fx<-Fc mass"), "Fx_from_Fxc": _nonnegative(c, 2 * phase**2 * expbox(-5 * A, A), "Fx<-Fxc mass")}
    result.update({"Fc_to_F": result["F_from_Fc"], "Fxc_to_F": result["F_from_Fxc"], "Fc_to_Fx": result["Fx_from_Fc"], "Fxc_to_Fx": result["Fx_from_Fxc"]})
    return result


kernel_masses = inlet_kernel_masses


def _output(c: Any, chi_term: Mapping[str, Any], incoming_term: Mapping[str, Any], branch: str, field: str, component: str, chi_order: int, incoming_order: int, mass: Any) -> dict[str, Any]:
    base = tuple(a + b for a, b in zip(chi_term["powers"], incoming_term["powers"]))
    powers = (base[0] + FORMAL_HB_POWER_ADDED[branch], base[1], base[2], base[3])
    coefficient = -2 * chi_term["coefficient"] * incoming_term["coefficient"] * mass
    error_only = bool(_meta(chi_term)["error_only"] or _meta(incoming_term)["error_only"] or chi_term.get("symmetric_error") or incoming_term.get("symmetric_error"))
    if error_only:
        coefficient = _symmetric(c, coefficient)
    ids = _unique(_meta(chi_term)["source_ids"] + _meta(incoming_term)["source_ids"])
    classes = _unique(_meta(chi_term)["error_classes"] + _meta(incoming_term)["error_classes"])
    if error_only and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    piece = "chi_Z*incoming_Z0" if (chi_order, incoming_order) == (1, 0) else "chi*incoming_Z1" if (chi_order, incoming_order) == (0, 1) else "chi*incoming_Z0"
    uncertainty = bool(error_only or _meta(chi_term)["uncertainty"] or _meta(incoming_term)["uncertainty"])
    origin_texts = [text for text in (chi_term.get("origin_text"), incoming_term.get("origin_text")) if isinstance(text, str)]
    origin = {"operation": "first_Kchi_inlet_potential_v5", "branch": branch, "field": field, "incoming_component": component, "chi_order": chi_order, "incoming_order": incoming_order, "leibniz_piece": piece, "source_ids": ids, "error_classes": classes, "uncertainty": uncertainty, "error_only": error_only, "symmetric_error": error_only, "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch], "coefficient_contains_formal_hb": False, "baseline_K0_reentered": False, "origin_texts": origin_texts, "chi_origin": _origin_descriptor(chi_term.get("origin")), "incoming_origin": _origin_descriptor(incoming_term.get("origin"))}
    return {"powers": powers, "coefficient": coefficient, "source_ids": ids, "error_classes": classes, "uncertainty": uncertainty, "origin": origin, "origin_texts": origin_texts, "branch": branch, "field": field, "incoming_component": component, "chi_order": chi_order, "incoming_order": incoming_order, "leibniz_piece": piece, "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch], "kernel_mass": mass, "kernel_mass_is_coefficient_only": True, "coefficient_contains_formal_hb": False, "physical_dt_count": 1, "causal_positive_kernel": True, "signed_source_term": not error_only, "error_only": error_only, "symmetric_error": error_only, "baseline_K0_included": False}


def _branch_set(c: Any, incoming: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]], chi: Mapping[int, Sequence[Mapping[str, Any]]], masses: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    components = {"F_from_Fc": ("F", "Fc"), "F_from_Fxc": ("F", "Fxc"), "Fx_from_Fc": ("Fx", "Fc"), "Fx_from_Fxc": ("Fx", "Fxc")}
    for branch in BRANCHES:
        field, component = components[branch]
        rows = {"Z0": [], "Z1": []}
        for order, routes in ((0, ((0, 0),)), (1, ((1, 0), (0, 1)))):
            for chi_order, incoming_order in routes:
                for chi_term in chi[chi_order]:
                    for incoming_term in incoming[component][incoming_order]:
                        rows[f"Z{order}"].append(_output(c, chi_term, incoming_term, branch, field, component, chi_order, incoming_order, masses[branch]))
        result[branch] = {"branch": branch, "field": field, "incoming_component": component, "operator": "D=-2*chi*inlet-potential", "formula": "-2*chi*integral_2^X K0(X,t)*H0F(t)dt" if field == "F" else "-2*chi*integral_2^X (t/X)^3*H0F(t)dt", "input_basis": "H0F=Fc" if component == "Fc" else "H0F=(1-4/t**2)*Fxc", "response": rows, "Z0": rows["Z0"], "Z1": rows["Z1"], "kernel_mass": masses[branch], "kernel_mass_formula": MASS_FORMULAS[branch], "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch], "signed_sectors": True, "error_only_sectors_separate": True, "coefficient_contains_formal_hb": False}
    return result


def build(c: Any, incoming: Any = None, chi: Any = None, *, hb: Any, sc: Any, Fc: Any = None, Fxc: Any = None, incoming_Fc: Any = None, incoming_Fxc: Any = None, inlet_groups: Any = None, G0g: Any = None, core_Dg: Any = None, core_symmetric_tail: Any = None, chi_Z: Any = None, chi_jets: Any = None, source_family: Any = None, source_atoms: Any = None, geometry_certificate: Any = None, certificate: Any = None, q_star: Any = 1) -> dict[str, Any]:
    """Return finite D inlet branches for one or more distinct inlet groups."""

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
    qlo, qhi = _bounds(c, q)
    if qlo != 1 or qhi != 1:
        raise ValueError("q_star must be exact endpoint 1")
    hb = _strict_hb(c, hb, allow_closed=cert["strict_hb_admitted"])
    sc = _positive(c, sc, "sc")
    chi_rows = _chi(c, chi, chi_Z)
    masses = inlet_kernel_masses(c, hb=hb, sc=sc, strict_hb_admitted=cert["strict_hb_admitted"])
    groups = _groups(c, incoming, Fc=Fc, Fxc=Fxc, incoming_Fc=incoming_Fc, incoming_Fxc=incoming_Fxc, inlet_groups=inlet_groups, G0g=G0g, core_Dg=core_Dg, core_symmetric_tail=core_symmetric_tail)
    grouped = {name: _branch_set(c, packet, chi_rows, masses) for name, packet in groups.items()}
    primary = "core_Dg" if "core_Dg" in grouped else next(iter(grouped))
    branches = grouped[primary]
    aggregate = {"F": {"Z0": branches["F_from_Fc"]["Z0"] + branches["F_from_Fxc"]["Z0"], "Z1": branches["F_from_Fc"]["Z1"] + branches["F_from_Fxc"]["Z1"]}, "Fx": {"Z0": branches["Fx_from_Fc"]["Z0"] + branches["Fx_from_Fxc"]["Z0"], "Z1": branches["Fx_from_Fc"]["Z1"] + branches["Fx_from_Fxc"]["Z1"]}}
    response = {"F": {"Fc": branches["F_from_Fc"]["response"], "Fxc": branches["F_from_Fxc"]["response"]}, "Fx": {"Fc": branches["Fx_from_Fc"]["response"], "Fxc": branches["Fx_from_Fxc"]["response"]}}
    endpoint = c.mpf(2) * c.exp(hb * sc / 8)
    return {"schema": SCHEMA, "status": "conditional_partial_first_kchi_inlet_potential", "source_family": source_family, "source_atoms": atoms, "source_power_basis": list(SOURCE_POWER_BASIS), "fields": list(FIELDS), "inlet_components": list(INLET_COMPONENTS), "ordinary_Z_orders": [0, 1], "branches": branches, "branch_names": list(BRANCHES), "response": response, "D_inlet": response, "D": aggregate, "D_x": aggregate["Fx"], "F": aggregate["F"], "Fx": aggregate["Fx"], "inlet_groups": list(grouped), "primary_inlet_group": primary, "group_results": grouped, "grouped_inlet_packets_separate": True, "grouped_inlet_packets_merged": False, "all_branches": {f"{name}:{branch}": item for name, packet in grouped.items() for branch, item in packet.items()}, "chi_jets": {f"Z{order}": list(chi_rows[order]) for order in (0, 1)}, "incoming_jets": {component: {f"Z{order}": list(groups[primary][component][order]) for order in (0, 1)} for component in INLET_COMPONENTS}, "geometry": {"q_star": q, "q_interval": "[0,1]", "x_endpoint": endpoint, "x_formula": "2*exp(hb*sc/8)", "c": sc / 8, "A": sc / 16, "hb_interval": hb, "sc_interval": sc, "admissible_hb": "0<hb<1/2", "strict_hb_admitted": True, "closed_hb_enclosure_allowed": True, "hb_sc_Z_independent": True, "rho_radius_Z_independent": True, "ordinary_Z_geometry_fixed": True, "ordinary_Z_geometry_differentiated": False, "formal_hb_not_materialized": True, "physical_dt_attached_once": True}, "kernel": {"K0": "K0(X,t)=t*(1-t**2/X**2)/2", "F": "D=-2*chi*integral_2^X K0(X,t)*H0F(t)dt", "Fx": "D_x=-2*chi*integral_2^X (t/X)^3*H0F(t)dt", "H0F": "Fc+(1-4/t**2)*Fxc", "H0Fx": "(2/t)^3*Fxc", "masses": {branch: masses[branch] for branch in BRANCHES}, "mass_formulas": dict(MASS_FORMULAS), "causal_positive": True}, "geometry_certificate": cert, "certificate_passthrough": cert["raw"], "chi_Z_required": True, "ordinary_Z1_product_rule": "chi_Z*incoming_Z0 + chi*incoming_Z1", "ordinary_Z1_geometry_differentiated": False, "source_ids_preserved": True, "error_classes_preserved": True, "formal_powers_preserved": True, "signed_coefficients_retained": True, "uncertainty_separate": True, "error_only_outputs_explicit": True, "formal_hb_powers_added": dict(FORMAL_HB_POWER_ADDED), "finite_D_inlet_only": True, "Kchi_first_correction_only": True, "Kchi_higher_words_included": False, "j2_kernel_remainder_included": False, "j2_remainder_included": False, "baseline_K0_reentered": False, "baseline_K0_included": False, "allword_tail_added": False, "m15_computed": False, "complete_collar_field": False, "finite_physical_field_claimed": False, "old_builder_called": False, "ancestor_imports": False, "receipt_reads": False, "automatic_code_execution": False}


build_inlet = build
build_inlet_potential = build
inlet_potential = build
evaluate = build
