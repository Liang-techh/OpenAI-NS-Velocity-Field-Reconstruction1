"""Source-bound signed local second-potential word and physical m15 leaf.

The caller supplies one exact-source certificate, whole-phase ``g5`` sectors,
axis ``chi``/``chi_Z`` sectors, and leading-F sectors.  The leaf emits the
finite ``j=2`` word, its same-source directed-negative P pair, and an
independent symmetric ``j>=3`` remainder.  No old tail is replaced, summed,
or written to a ledger; no producer or legacy module is imported.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

SCHEMA = "n1_local_second_potential_m15_v2"
COMPONENT = "local_second_potential_m15"
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
SOURCE_POWER_BASIS = POWER_BASIS
ORDINARY_Z_ORDERS = (0, 1)
ORDINARY_Z_LABELS = ("Z0", "Z1")
PHYSICAL_POWER_SHIFT = (1, 0, 0, -2)
ETA_NUMERATOR = 5
ETA_DENOMINATOR = 8

__all__ = [
    "SCHEMA", "COMPONENT", "POWER_BASIS", "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS", "ORDINARY_Z_LABELS", "PHYSICAL_POWER_SHIFT",
    "build", "build_local_second_potential", "build_local_second_potential_m15", "enclose", "evaluate",
]

_MISSING = object()
_ROWS = (("Z0", "z0", 0, "0"), ("Z1", "z1", 1, "1"))


def _map(x: Any) -> bool:
    return isinstance(x, Mapping)


def _seq(x: Any) -> bool:
    return isinstance(x, Sequence) and not isinstance(x, (str, bytes, bytearray))


def _get(m: Mapping[Any, Any], *keys: Any) -> Any:
    for key in keys:
        try:
            if key in m:
                return m[key]
        except (TypeError, ValueError):
            pass
    low = {str(key).lower(): value for key, value in m.items()}
    for key in keys:
        if str(key).lower() in low:
            return low[str(key).lower()]
    return _MISSING


def _items(x: Any) -> list[Any]:
    if x is None or x is _MISSING:
        return []
    if isinstance(x, (str, bytes, bytearray)):
        return [x]
    return list(x) if _seq(x) else [x]


def _raw(x: Any, upper: bool) -> Any:
    raw = getattr(x, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    if hasattr(x, "a") and hasattr(x, "b"):
        endpoint = x.b if upper else x.a
        nested = getattr(endpoint, "_mpi_", None)
        return nested[1 if upper else 0] if nested is not None else endpoint
    return x


def _point(c: Any, raw: Any) -> Any:
    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((raw, raw))
        except Exception:
            pass
    return c.mpf(raw)


def _iv(c: Any, lo: Any, hi: Any) -> Any:
    a, b = _raw(lo, False), _raw(hi, True)
    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((a, b))
        except Exception:
            pass
    try:
        return c.mpf([_point(c, a), _point(c, b)])
    except Exception:
        return c.mpf((a, b))


def _as_iv(c: Any, x: Any, label: str = "interval") -> Any:
    if hasattr(x, "_mpi_"):
        return x
    if hasattr(x, "a") and hasattr(x, "b"):
        return _iv(c, x.a, x.b)
    if _map(x):
        nested = _get(x, "interval", "enclosure", "box", "value", "record")
        if nested is not _MISSING and nested is not x:
            return _as_iv(c, nested, label)
        lo, hi = _get(x, "lower", "lo", "left"), _get(x, "upper", "hi", "right")
        if lo is not _MISSING and hi is not _MISSING:
            return _iv(c, lo, hi)
    if _seq(x) and len(x) == 2:
        return _iv(c, x[0], x[1])
    try:
        return c.mpf(x)
    except Exception as exc:
        raise TypeError(f"{label} is not interval-compatible") from exc


def _lo(c: Any, x: Any) -> Any:
    return _point(c, _raw(x, False))


def _hi(c: Any, x: Any) -> Any:
    return _point(c, _raw(x, True))


def _zero(c: Any) -> Any:
    return c.mpf(0)


def _one(c: Any) -> Any:
    return c.mpf(1)


def _absup(c: Any, x: Any) -> Any:
    lo, hi, z = _lo(c, x), _hi(c, x), _zero(c)
    if bool(lo >= z):
        return hi
    if bool(hi <= z):
        return -lo
    return -lo if bool(-lo >= hi) else hi


def _sym(c: Any, x: Any) -> Any:
    r = _absup(c, x)
    return _iv(c, -r, r)


def _neg(c: Any, x: Any) -> Any:
    return _iv(c, -_hi(c, x), -_lo(c, x))


def _unique(xs: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    for x in xs:
        if not any(x == old for old in result):
            result.append(x)
    return result


def _walk(x: Any) -> list[Mapping[Any, Any]]:
    if not _map(x):
        return []
    result = [x]
    for child in x.values():
        if _map(child):
            result.extend(_walk(child))
        elif _seq(child):
            for item in child:
                if _map(item):
                    result.extend(_walk(item))
    return result


def _named(cert: Mapping[Any, Any], names: Sequence[str]) -> list[Any]:
    wanted = {str(n).lower() for n in names}
    return [item for m in _walk(cert) for key, item in m.items() if str(key).lower() in wanted]


def _truth(x: Any) -> bool:
    if isinstance(x, bool):
        return x
    if isinstance(x, str):
        return x.strip().lower() in {"true", "yes", "passed", "proved", "verified", "certified"}
    if _map(x):
        vals = [_get(x, key) for key in ("passed", "proved", "verified", "certified", "valid", "same")]
        vals = [value for value in vals if value is not _MISSING]
        return bool(vals) and all(value is True for value in vals)
    return False


def _flag(cert: Mapping[Any, Any], names: Sequence[str], label: str) -> bool:
    vals = _named(cert, names)
    if not vals or any(not _truth(value) for value in vals):
        raise ValueError(f"certificate is missing passed {label} proof")
    return True


def _sign(x: Any) -> int | None:
    if isinstance(x, bool):
        return 2 if x else None
    if isinstance(x, int) and not isinstance(x, bool):
        return x
    if isinstance(x, str):
        try:
            return int(x.strip().replace("+", ""))
        except ValueError:
            return None
    if _map(x):
        y = _get(x, "sign", "value", "coefficient", "P_sign")
        return None if y is _MISSING else _sign(y)
    return None


def _eta(c: Any, x: Any) -> bool:
    if x is True:
        return True
    if _map(x):
        for key in ("mass", "value", "eta", "interval", "enclosure"):
            value = _get(x, key)
            if value is not _MISSING and _eta(c, value):
                return True
        return _truth(x)
    try:
        y = _as_iv(c, x)
        e = c.mpf(5) / c.mpf(8)
        return bool(_lo(c, y) == e and _hi(c, y) == e)
    except Exception:
        return False


def _certificate(c: Any, hb: Any, sc: Any, cert: Any) -> dict[str, Any]:
    if not _map(cert):
        raise TypeError("caller certificate is required")
    declared = _get(cert, "passed", "checks_passed", "certificate_passed")
    if declared is not _MISSING and declared is not True:
        raise ValueError("caller certificate passed flag is not true")
    combined = _named(cert, ("same_source", "same_source_chi_g5", "same_chi_g5_source", "same_source_family", "same_source_chi", "same_chi_g5"))
    if combined:
        same = all(_truth(x) for x in combined)
    else:
        same = _flag(cert, ("same_g5_source", "g5_same_source", "same_g5", "same_source_g5"), "same g5 source") and _flag(
            cert, ("same_chi_source", "chi_same_source", "same_chi", "same_source_chi"), "same chi source"
        )
    if not same:
        raise ValueError("certificate must prove one exact common chi/g5 source")
    axis = _flag(cert, ("chi_is_axis_potential", "chi_axis_potential", "chi_x_independent", "chi_is_x_independent", "chi_is_x_independent_axis", "chi_is_axis", "axis_potential"), "axis chi")
    geom = _flag(cert, ("geometry_Z_independent", "geometry_z_independent", "ordinary_Z_geometry_fixed", "Z_independent_geometry"), "Z-independent geometry")
    whole = _flag(cert, ("g5_whole_phase_enclosure", "whole_phase_g5_enclosure", "g5_full_phase", "source_g5_full_q_enclosure", "full_q_source_enclosure", "whole_phase_enclosure"), "whole-phase g5")
    signs = [_sign(x) for x in _named(cert, ("P_sign", "p_sign", "pressure_sign", "P_sign_value"))]
    if not signs or any(x is not None and x != 2 for x in signs) or 2 not in signs:
        raise ValueError("certificate P_sign must be +2")
    etas = _named(cert, ("eta_mass_5_over_8", "eta_mass", "original_eta_mass", "original_cutoff_mass_5_over_8", "eta"))
    if not etas or any(not _eta(c, x) for x in etas):
        raise ValueError("certificate must prove eta mass 5/8")
    hb_i, sc_i = _as_iv(c, hb, "hb"), _as_iv(c, sc, "sc")
    if bool(_lo(c, hb_i) < _zero(c)) or bool(_hi(c, hb_i) > c.mpf(1) / c.mpf(2)):
        raise ValueError("hb closure must be contained in [0,1/2]")
    if bool(_hi(c, sc_i) <= _zero(c)):
        raise ValueError("sc must be positive")
    q = _get(cert, "q_domain", "q_bounds", "source_q_domain")
    q_i = _as_iv(c, q, "q_domain") if q is not _MISSING else _iv(c, _zero(c), _one(c))
    if bool(_lo(c, q_i) > _zero(c)) or bool(_hi(c, q_i) < _one(c)):
        raise ValueError("q domain must contain [0,1]")
    return {
        "passed": True, "same_source_chi_g5": same, "chi_is_axis_potential": axis,
        "chi_x_independent": True, "geometry_Z_independent": geom,
        "g5_whole_phase_enclosure": whole, "P_sign": 2,
        "eta_mass": c.mpf(5) / c.mpf(8), "eta_mass_formula": "5/8",
        "hb": hb_i, "hb_closure": _iv(c, _zero(c), c.mpf(1) / c.mpf(2)),
        "hb_coefficient_bound_only": True, "q_domain": q_i, "sc": sc_i,
        "raw": dict(cert),
    }


def _powers(row: Mapping[Any, Any]) -> list[int]:
    raw = _get(row, "powers", "formal_powers", "power_tuple")
    if raw is _MISSING:
        raw = [_get(row, "hb_power", "hb"), _get(row, "F0_anchor_power", "F0_power", "F0"), _get(row, "Pstar_squared_power", "Pstar2_power", "Pstar_power"), _get(row, "Lambda_power", "lambda_power", "Lambda")]
        raw = [0, 0, 0, 0] if any(x is _MISSING for x in raw) else raw
    if _map(raw):
        raw = [_get(raw, "hb", 0), _get(raw, "F0_anchor", "F0", 1), _get(raw, "Pstar_squared", "Pstar2", 2), _get(raw, "Lambda", "lambda", 3)]
    if not _seq(raw) or len(raw) != 4 or any(isinstance(x, bool) or not isinstance(x, int) for x in raw):
        raise TypeError("sector powers must contain four integers")
    return [int(x) for x in raw]


def _sector(c: Any, raw: Any, label: str, default: str) -> dict[str, Any]:
    if not _map(raw):
        raw = {"coefficient": raw}
    coef = _get(raw, "coefficient", "coefficient_interval", "interval", "enclosure", "value")
    if coef is _MISSING:
        raise KeyError(f"{label} is missing coefficient")
    origin = raw.get("origin") if _map(raw.get("origin")) else {}
    ids = _get(raw, "source_ids", "source_id", "sources", "source")
    ids = _get(origin, "source_ids", "source_id", "sources", "source") if ids is _MISSING else ids
    if ids is _MISSING:
        raise KeyError(f"{label} is missing source_ids")
    paths = _get(raw, "source_paths", "source_path", "paths", "path")
    paths = _get(origin, "source_paths", "source_path", "paths", "path") if paths is _MISSING else paths
    cls = _get(raw, "error_classes", "error_class", "uncertainty_class")
    cls = _get(origin, "error_classes", "error_class", "uncertainty_class") if cls is _MISSING else cls
    ids, paths, cls = (_items(ids), [] if paths is _MISSING else _items(paths), [] if cls is _MISSING else _items(cls))
    error = raw.get("error_only", raw.get("symmetric_error", False))
    if not isinstance(error, bool):
        raise TypeError(f"{label}.error_only must be boolean")
    result = dict(raw)
    result.update({"coefficient": _as_iv(c, coef, f"{label}.coefficient"), "coefficient_interval": _as_iv(c, coef, f"{label}.coefficient"), "powers": _powers(raw), "formal_powers": _powers(raw), "power_basis": list(POWER_BASIS), "source_ids": ids, "source_paths": paths, "error_classes": cls, "error_only": error, "symmetric_error": bool(raw.get("symmetric_error", error)), "uncertainty": bool(raw.get("uncertainty", False) or error)})
    return result


def _rows(c: Any, value: Any, label: str, default: str) -> dict[str, list[dict[str, Any]]]:
    current = value
    for _ in range(8):
        if not _map(current) or _get(current, "coefficient", "coefficient_interval", "interval", "value") is not _MISSING:
            break
        nested = _get(current, "ordinary_Z", "rows", "jets", "g5", "g5_source", "whole_phase", "F", "F1", "leading_F", "Flead", "field", "sectors", "source_sectors", "terms", "records")
        if nested is _MISSING or nested is current:
            break
        current = nested
    if _seq(current) and len(current) == 2:
        raw = {"Z0": current[0], "Z1": current[1]}
    elif _map(current):
        raw = {}
        for name, keys in zip(ORDINARY_Z_LABELS, _ROWS):
            selected = _get(current, *keys)
            if selected is _MISSING:
                raise KeyError(f"{label} is missing {name}")
            raw[name] = selected
    else:
        raise TypeError(f"{label} requires Z0 and Z1 rows")
    answer: dict[str, list[dict[str, Any]]] = {}
    for name in ORDINARY_Z_LABELS:
        values = raw[name]
        if _map(values):
            if _get(values, "coefficient", "coefficient_interval", "interval", "value") is not _MISSING:
                values = [values]
            else:
                nested_values = _get(values, "sectors", "source_sectors", "terms", "records")
                if nested_values is not _MISSING:
                    values = nested_values
                elif values and all(_map(item) for item in values.values()):
                    values = list(values.values())
                else:
                    values = [values]
        elif not _seq(values):
            values = [values]
        answer[name] = [_sector(c, item, f"{label}.{name}[{i}]", f"{default}:{name}:{i}") for i, item in enumerate(values)]
    return answer


def _chi(c: Any, chi: Any, chi_Z: Any, chi_jets: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if chi is None and chi_jets is None:
        raise TypeError("chi is required")
    if chi is not None and chi_jets is not None:
        raise ValueError("chi and chi_jets were supplied twice")
    raw = chi if chi is not None else chi_jets
    if _map(raw):
        wrapped = _get(raw, "chi", "chi_jets", "ordinary_Z", "rows", "jets")
        if wrapped is not _MISSING and wrapped is not raw:
            raw = wrapped
    if _map(raw) and _get(raw, "Z0", "z0", 0, "0") is not _MISSING:
        rows = _rows(c, raw, "chi", "chi")
        z0, z1 = rows["Z0"], rows["Z1"]
    else:
        z0, z1 = [_sector(c, raw, "chi.Z0", "chi:Z0")], []
    if chi_Z is not None:
        if z1:
            raise ValueError("chi_Z supplied twice")
        z1 = [_sector(c, chi_Z, "chi.Z1", "chi:Z1")]
    if not z1:
        raise TypeError("chi_Z is required")
    return z0, z1


def _sum_powers(*rows: Sequence[int]) -> list[int]:
    return [sum(row[i] for row in rows) for i in range(4)]


def _prov(*rows: Mapping[str, Any]) -> dict[str, Any]:
    return {"source_ids": _unique([x for row in rows for x in row.get("source_ids", [])]), "source_paths": _unique([x for row in rows for x in row.get("source_paths", [])]), "error_classes": _unique([x for row in rows for x in row.get("error_classes", [])]), "uncertainty": any(row.get("uncertainty", False) for row in rows), "error_only": any(row.get("error_only", False) for row in rows)}


def _field(c: Any, source: Mapping[str, Any], chi0: Mapping[str, Any], chiZ: Mapping[str, Any], mass: Any, route: str) -> dict[str, Any]:
    if route == "chi2_g0" or route == "chi2_gZ":
        factor, factors, piece, pos = chi0["coefficient"] ** 2, (chi0, chi0), ("chi*chi*g0" if route == "chi2_g0" else "chi*chi*g_Z"), None
    elif route in {"chiZ_outer_g0", "chiZ_inner_g0"}:
        factor, factors = chi0["coefficient"] * chiZ["coefficient"], (chi0, chiZ)
        pos, piece = ("outer" if route == "chiZ_outer_g0" else "inner"), ("chi_Z(outer)*chi*g0" if route == "chiZ_outer_g0" else "chi*chi_Z(inner)*g0")
    else:
        raise ValueError(f"unknown field route {route}")
    raw = c.mpf(4) * mass * factor * source["coefficient"]
    uses_chiZ = route in {"chiZ_outer_g0", "chiZ_inner_g0"}
    p = _prov(source, *factors)
    coef = _sym(c, raw) if p["error_only"] else raw
    powers = _sum_powers(source["powers"], *(x["powers"] for x in factors))
    powers[0] += 6
    return {"coefficient": coef, "coefficient_interval": coef, "powers": powers, "formal_powers": list(powers), "power_basis": list(POWER_BASIS), "family": "finite_j2", "route": route, "ordinary_Z_piece": piece, "ordinary_Z_order": 0 if route == "chi2_g0" else 1, "chi_Z_insertion_position": pos, "operator": "+4 K0 chi K0 chi K0 g", "operator_coefficient": 4, "kernel_mass": mass, "kernel_mass_formula": "[0,B^6/720]", "formal_hb_power_added": 6, "q0_value": _zero(c), "q0_exact_zero": True, "signed_field": not p["error_only"], "signed_uncertainty": p["uncertainty"], "finite_signed_interval": not p["error_only"], "error_only": p["error_only"], "symmetric_error": p["error_only"], "uncertainty": p["uncertainty"], "source_ids": p["source_ids"], "source_paths": p["source_paths"], "error_classes": p["error_classes"], "source_ids_by_factor": {"g5": list(source["source_ids"]), "chi": list(chi0["source_ids"]), "chi_Z": list(chiZ["source_ids"]) if uses_chiZ else []}, "source_paths_by_factor": {"g5": list(source["source_paths"]), "chi": list(chi0["source_paths"]), "chi_Z": list(chiZ["source_paths"]) if uses_chiZ else []}, "chi_Z_source_ids": list(chiZ["source_ids"]) if uses_chiZ else [], "chi_Z_source_paths": list(chiZ["source_paths"]) if uses_chiZ else [], "source_coefficient": source["coefficient"], "signed_source_coefficient_preserved": True, "point_value_selected": False, "chi_reapplied": False, "green_reapplied": False}


def _tail_field(c: Any, source: Mapping[str, Any], chi0: Mapping[str, Any], chiZ: Mapping[str, Any], C0: Any, C1: Any, route: str) -> dict[str, Any]:
    if route == "C0_g0" or route == "C0_gZ":
        factor, factors, piece, c1 = C0, (chi0,), ("chi^3*g0" if route == "C0_g0" else "chi^3*g_Z"), False
    elif route == "C1_chiZ_g0":
        factor, factors, piece, c1 = C1, (chi0, chi0, chiZ), "chi_Z*g0", True
    else:
        raise ValueError(f"unknown tail route {route}")
    radius = factor * _absup(c, source["coefficient"])
    coef, powers = _sym(c, radius), list(source["powers"])
    powers[0] += 8
    ids = _unique([x for row in (source, *factors) for x in row.get("source_ids", [])])
    paths = _unique([x for row in (source, *factors) for x in row.get("source_paths", [])])
    classes = _unique([x for row in (source, *factors) for x in row.get("error_classes", [])] + ["j>=3_symmetric_remainder", "symmetric_error"])
    return {"coefficient": coef, "coefficient_interval": coef, "radius": radius, "powers": powers, "formal_powers": list(powers), "power_basis": list(POWER_BASIS), "family": "j_ge_3", "route": route, "ordinary_Z_piece": piece, "ordinary_Z_order": 0 if route == "C0_g0" else 1, "radius_factor": factor, "radius_formula": "24*M^2*Mz*B^8*exp(M*B^2)/40320" if c1 else "(2M)^3*B^8*exp(M*B^2/2)/40320", "formal_hb_power_added": 8, "symmetric_error": True, "error_only": True, "uncertainty": True, "signed_tail_value": False, "assigned_tail_value": False, "cancellation_claimed": False, "source_ids": ids, "source_paths": paths, "error_classes": classes, "chi_source_ids": _unique([x for row in factors for x in row.get("source_ids", [])]), "chi_Z_source_ids": list(chiZ["source_ids"]) if c1 else [], "chi_Z_source_paths": list(chiZ["source_paths"]) if c1 else [], "source_ids_by_factor": {"g5": list(source["source_ids"]), "chi": list(chi0["source_ids"]), "chi_Z": list(chiZ["source_ids"]) if c1 else []}, "source_paths_by_factor": {"g5": list(source["source_paths"]), "chi": list(chi0["source_paths"]), "chi_Z": list(chiZ["source_paths"]) if c1 else []}, "source_coefficient": source["coefficient"], "signed_source_coefficient_preserved": True, "point_value_selected": False, "chi_reapplied": False, "green_reapplied": False, "majorant_drops_alternating_word_sign": True}


def _physical(c: Any, sc: Any) -> dict[str, Any]:
    upper = _hi(c, c.exp(sc / c.mpf(4)))
    env = _iv(c, _one(c), upper)
    factor = _iv(c, -c.mpf(5) * _hi(c, sc) * upper, -c.mpf(5) * _lo(c, sc))
    return {"exp_sc_over_4": env, "physical_factor": factor, "physical_factor_formula": "-5*sc*[1,exp(sc/4)]", "physical_factor_negative": bool(_hi(c, factor) < _zero(c)), "eta": c.mpf(5) / c.mpf(8), "eta_mass_formula": "5/8", "formal_power_shift": list(PHYSICAL_POWER_SHIFT)}


def _project(c: Any, physical: Mapping[str, Any], fields: Mapping[str, Sequence[Mapping[str, Any]]], leading: Mapping[str, Sequence[Mapping[str, Any]]], family: str) -> dict[str, Any]:
    out: dict[str, Any] = {"Z0": [], "Z1": []}
    for lead in leading["Z0"]:
        for field in fields["Z0"]:
            out["Z0"].append(_project_term(c, physical, lead, field, family, "Z0_leading_Z0_times_field_Z0"))
    for lead in leading["Z1"]:
        for field in fields["Z0"]:
            out["Z1"].append(_project_term(c, physical, lead, field, family, "Z1_leading_Z1_times_field_Z0"))
    for lead in leading["Z0"]:
        for field in fields["Z1"]:
            out["Z1"].append(_project_term(c, physical, lead, field, family, f"Z1_leading_Z0_times_{field['route']}"))
    out.update({"family": family, "physical_swirl": "m15", "ordinary_Z1_paths": ["leading_F_Z1*field_Z0", "leading_F_Z0*field_Z1"], "Z1_routes_separate": True, "point_values_selected": False})
    return out


def _project_term(c: Any, physical: Mapping[str, Any], lead: Mapping[str, Any], field: Mapping[str, Any], family: str, kind: str) -> dict[str, Any]:
    raw = physical["physical_factor"] * lead["coefficient"] * field["coefficient"]
    error = bool(lead.get("error_only", False) or field.get("error_only", False))
    coef = _sym(c, raw) if error else raw
    p = _prov(lead, field)
    powers = _sum_powers(lead["powers"], field["powers"], PHYSICAL_POWER_SHIFT)
    classes = _unique(p["error_classes"] + (["symmetric_error"] if error else []))
    return {"coefficient": coef, "coefficient_interval": coef, "powers": powers, "formal_powers": list(powers), "power_basis": list(POWER_BASIS), "family": family, "field_family": field.get("family"), "route": field.get("route"), "ordinary_Z_order": field.get("ordinary_Z_order"), "derivative_kind": kind, "ordinary_Z_piece": field.get("ordinary_Z_piece"), "physical_swirl": "m15", "physical_factor": physical["physical_factor"], "physical_factor_formula": physical["physical_factor_formula"], "physical_formal_power_shift": list(PHYSICAL_POWER_SHIFT), "source_ids": p["source_ids"], "source_paths": p["source_paths"], "error_classes": classes, "uncertainty": bool(p["uncertainty"] or error), "error_only": error, "symmetric_error": error, "signed_field": family == "finite_j2" and not error, "signed_uncertainty": p["uncertainty"], "signed_tail_value": family == "finite_j2" and not error, "assigned_tail_value": family == "finite_j2" and not error, "cancellation_claimed": False, "leading_source_ids": list(lead["source_ids"]), "leading_source_paths": list(lead["source_paths"]), "field_source_ids": list(field["source_ids"]), "field_source_paths": list(field["source_paths"]), "source_ids_by_factor": {"leading_F": list(lead["source_ids"]), "field": list(field["source_ids"])}, "source_paths_by_factor": {"leading_F": list(lead["source_paths"]), "field": list(field["source_paths"])}, "field_source": dict(field), "leading_source": dict(lead), "point_value_selected": False, "chi_reapplied": False, "green_reapplied": False, "chi_Z_source_ids": list(field.get("chi_Z_source_ids", [])), "chi_Z_source_paths": list(field.get("chi_Z_source_paths", []))}


def _pairs(c: Any, finite: Mapping[str, Sequence[Mapping[str, Any]]], prefix: str = "finite_j2") -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    pairs: list[dict[str, Any]] = []
    for z in ORDINARY_Z_LABELS:
        for i, term in enumerate(finite[z]):
            pair_id = f"{prefix}:{z}:{i}"
            neg = dict(term)
            neg.update({"coefficient": _neg(c, term["coefficient"]), "coefficient_interval": _neg(c, term["coefficient"]), "role": "P_compensation", "generated_from": "D2g", "correlated_with_pair_id": pair_id, "interval_operation": "exact_directed_negative", "independent_interval_sum_performed": False, "functional_identity": "P=-D2g", "source_ids": list(term["source_ids"]), "source_paths": list(term["source_paths"]), "error_classes": list(term["error_classes"]), "powers": list(term["powers"]), "formal_powers": list(term["formal_powers"])})
            rows[z].append(neg)
            parent = dict(term)
            parent.update({"role": "D2g_parent", "correlated_pair_id": pair_id, "independent_interval_sum_performed": False})
            pairs.append({"pair_id": pair_id, "ordinary_Z": z, "ordinary_Z_order": 0 if z == "Z0" else 1, "D2g": parent, "P": neg, "D2g_parent": parent, "P_compensation": neg, "formal_powers": list(term["formal_powers"]), "powers": list(term["powers"]), "source_ids": list(term["source_ids"]), "source_paths": list(term["source_paths"]), "error_classes": list(term["error_classes"]), "symbolic_sum": {"exact": "0", "exact_zero": True, "scope": "same-source correlated pair only", "independent_boxed_interval_sum_performed": False, "unconditioned_global": False}})
    return {"Z0": rows["Z0"], "Z1": rows["Z1"], "family": "P_correlated_negative", "role": "P_compensation", "correlated_pair_only": True, "same_source_required": True, "old_errors_preserved": True, "independent_box_sum": False, "point_values_selected": False}, pairs


def _copy(x: Any) -> Any:
    if _map(x):
        return {key: _copy(value) for key, value in x.items()}
    if _seq(x):
        return [_copy(value) for value in x]
    return x


def build(c: Any, hb: Any = None, sc: Any = None, chi: Any = None, chi_Z: Any = None, g5: Any = None, leading_F: Any = None, certificate: Any = None, *, chi_jets: Any = None, g5_source: Any = None, source: Any = None, leading: Any = None, leading_source: Any = None, old_j2_tail: Any = None, prior_j2_tail: Any = None, retained_j2_tail: Any = None, source_family: Any = None, source_atoms: Any = None, geometry_certificate: Any = None, source_certificate: Any = None, **aliases: Any) -> dict[str, Any]:
    """Build finite ``j=2``, correlated P, and independent ``j>=3`` packets."""
    if not all(hasattr(c, key) for key in ("mpf", "exp")):
        raise TypeError("c must provide directed interval arithmetic")
    hb = hb if hb is not None else aliases.pop("hb_interval", aliases.pop("h_b", None))
    sc = sc if sc is not None else aliases.pop("sc_interval", aliases.pop("s_c", None))
    if hb is None or sc is None:
        raise TypeError("hb and sc are required")
    candidates = [value for value in (g5, g5_source, source) if value is not None]
    if len(candidates) > 1:
        raise ValueError("g5 aliases supplied more than once")
    g5 = candidates[0] if candidates else aliases.pop("whole_phase_g5", aliases.pop("g5_whole_phase", aliases.pop("whole_phase_source", None)))
    if g5 is None:
        raise TypeError("whole-phase g5 sectors are required")
    lead_candidates = [value for value in (leading_F, leading, leading_source) if value is not None]
    if len(lead_candidates) > 1:
        raise ValueError("leading F aliases supplied more than once")
    leading_F = lead_candidates[0] if lead_candidates else aliases.pop("F_leading", aliases.pop("leading_F_source", None))
    if leading_F is None:
        raise TypeError("leading F sectors are required")
    certs = [value for value in (certificate, geometry_certificate, source_certificate) if value is not None]
    if len(certs) > 1 and any(value != certs[0] for value in certs[1:]):
        raise ValueError("certificate aliases disagree")
    certificate = certs[0] if certs else aliases.pop("proof", aliases.pop("caller_certificate", None))
    if certificate is None:
        raise TypeError("caller certificate is required")
    if chi is None:
        chi = aliases.pop("chi_source", None)
    if chi_Z is None:
        chi_Z = aliases.pop("chi_Z_source", aliases.pop("chi_z_source", None))
    if retained_j2_tail is None and prior_j2_tail is None and old_j2_tail is None:
        old_j2_tail = aliases.pop("old_tail", aliases.pop("old_j2_error_bound", aliases.pop("j2_tail", None)))
    if aliases:
        raise TypeError("unexpected local second-potential arguments: " + ", ".join(sorted(aliases)))
    cert = _certificate(c, hb, sc, certificate)
    chi0_rows, chiZ_rows = _chi(c, chi, chi_Z, chi_jets)
    g5_rows, lead_rows = _rows(c, g5, "g5", "g5"), _rows(c, leading_F, "leading_F", "leading_F")
    chi0, chiZ = chi0_rows[0], chiZ_rows[0]
    M, Mz = _iv(c, _zero(c), max((_absup(c, row["coefficient"]) for row in chi0_rows), default=_zero(c))), _iv(c, _zero(c), max((_absup(c, row["coefficient"]) for row in chiZ_rows), default=_zero(c)))
    B = (cert["sc"] / c.mpf(4)) * c.exp(cert["sc"] / c.mpf(16))
    mass = _iv(c, _zero(c), B**6 / c.mpf(720))
    C0 = ((c.mpf(2) * M) ** 3) * B**8 * c.exp(M * B**2 / c.mpf(2)) / c.mpf(40320)
    C1 = c.mpf(24) * M**2 * Mz * B**8 * c.exp(M * B**2) / c.mpf(40320)
    physical = _physical(c, cert["sc"])
    finite_field: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    finite_field["Z0"] = [_field(c, row, chi0, chiZ, mass, "chi2_g0") for row in g5_rows["Z0"]]
    finite_field["Z1"] = [_field(c, row, chi0, chiZ, mass, "chi2_gZ") for row in g5_rows["Z1"]]
    finite_field["Z1"] += [_field(c, row, chi0, chiZ, mass, route) for route in ("chiZ_outer_g0", "chiZ_inner_g0") for row in g5_rows["Z0"]]
    tail_field: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    tail_field["Z0"] = [_tail_field(c, row, chi0, chiZ, C0, C1, "C0_g0") for row in g5_rows["Z0"]]
    tail_field["Z1"] = [_tail_field(c, row, chi0, chiZ, C0, C1, "C0_gZ") for row in g5_rows["Z1"]]
    tail_field["Z1"] += [_tail_field(c, row, chi0, chiZ, C0, C1, "C1_chiZ_g0") for row in g5_rows["Z0"]]
    finite = _project(c, physical, finite_field, lead_rows, "finite_j2")
    tail = _project(c, physical, tail_field, lead_rows, "j_ge_3")
    p_group, pairs = _pairs(c, finite)
    p_field_group, field_pairs = _pairs(c, finite_field, prefix="field_finite_j2")
    inherited = retained_j2_tail if retained_j2_tail is not None else prior_j2_tail if prior_j2_tail is not None else old_j2_tail
    atoms = list(POWER_BASIS) if source_atoms is None else list(source_atoms)
    if set(atoms) != set(POWER_BASIS):
        raise ValueError("source_atoms must contain the four formal atoms")
    geometry = {**physical, "hb": cert["hb"], "hb_closure": cert["hb_closure"], "hb_coefficient_bound_only": True, "hb_formal_not_materialized": True, "sc": cert["sc"], "B": B, "B_formula": "(sc/4)*exp(sc/16)", "K0_bound": "K0(x,t)<=x-t", "L_bound": "L=x-2<=hb*B", "j2_positive_mass_envelope": mass, "j2_positive_mass_formula": "[0,B^6/720]*hb^6", "j2_formal_hb_power": 6, "M": M, "Mz": Mz, "M_formula": "sup|chi|", "Mz_formula": "sup|chi_Z|", "j3_C0": C0, "j3_C1": C1, "j3_C0_formula": "(2M)^3*B^8*exp(M*B^2/2)/40320", "j3_C1_formula": "24*M^2*Mz*B^8*exp(M*B^2)/40320", "j3_formal_hb_power": 8, "j3_physical_formal_hb_power": 9, "j3_no_M_division": True, "series_proof": {"j_equals_k_plus_3": True, "j_bound": "j=k+3<=3*2^k", "factorial_bound": "(2k+8)!>=8!*k!", "M_division_used": False}}
    finite_group = {**finite, "field": finite_field, "field_rows": finite_field, "physical_m15": finite, "operator": "+4 K0 chi K0 chi K0 g", "finite_signed_word": True}
    p_group.update({"field": {"Z0": list(p_field_group["Z0"]), "Z1": list(p_field_group["Z1"])}, "field_correlated_pairs": field_pairs, "physical_m15": {"Z0": list(p_group["Z0"]), "Z1": list(p_group["Z1"])}, "DPW0": "-D2g", "old_errors_preserved": True})
    tail_group = {**tail, "field": tail_field, "field_rows": tail_field, "physical_m15": tail, "tail_is_symmetric_error_only": True, "old_j2_tail_remains_retained": True, "cancellation_claimed": False}
    return {"schema": SCHEMA, "component": COMPONENT, "source_only": True, "caller_owned_context": True, "source_family": source_family, "source_atoms": atoms, "source_power_basis": list(POWER_BASIS), "power_basis": list(POWER_BASIS), "ordinary_Z_orders": list(ORDINARY_Z_ORDERS), "ordinary_Z_labels": list(ORDINARY_Z_LABELS), "certificate": cert, "caller_certificate_verified": True, "geometry": geometry, "chi": {"Z0": chi0_rows, "Z1": chiZ_rows}, "chi_jets": {"Z0": chi0_rows, "Z1": chiZ_rows}, "M": M, "Mz": Mz, "g5_jets": g5_rows, "leading_F_jets": lead_rows, "finite_j2": finite_group, "finite_word": finite_group, "j2": finite_group, "D2g": finite_group, "D²g": finite_group, "P_correlated_negative": p_group, "P_correlated": p_group, "P": p_group, "DPW0": p_group, "j_ge_3": tail_group, "j>=3": tail_group, "j3_remainder": tail_group, "tail": tail_group, "groups": {"finite_j2": finite_group, "P_correlated_negative": p_group, "j_ge_3": tail_group}, "groups_separate": True, "finite_field": finite_field, "finite_physical_m15": finite, "P_correlated_pairs": pairs, "P_field_correlated_pairs": field_pairs, "correlated_pairs": pairs, "old_j2_tail": _copy(inherited) if inherited is not None else [], "retained_old_j2_tail": _copy(inherited) if inherited is not None else [], "old_j2_tail_input_present": inherited is not None, "old_j2_tail_retained": True, "old_j2_tail_remains_retained": True, "old_j2_tail_overwritten": False, "old_j2_tail_cancelled": False, "old_j2_tail_independent_box_added": False, "old_errors_preserved": True, "all_error_classes_preserved": True, "all_formal_powers_preserved": True, "chi_Z_atom_ids_preserved": True, "ordinary_Z1_product_rule": "leading_F_Z1*field_Z0 + leading_F_Z0*field_Z1", "ordinary_Z1_field_product_rule": "chi*chi*g_Z + chi_Z(outer)*chi*g0 + chi*chi_Z(inner)*g0", "chi_Z_insertion_positions": ["outer", "inner"], "ordinary_Z1_leading_routes": ["leading_F_Z1*field_Z0", "leading_F_Z0*field_Z1"], "ordinary_Z1_routes_separate": True, "P_is_only_correlated_pair": True, "P_sign": 2, "eta_mass": c.mpf(5) / c.mpf(8), "independent_boxed_interval_sum_performed": False, "unconditioned_global_zero": False, "tail_cancellation_claimed": False, "j_ge_3_assigned_value": False, "j_ge_3_symmetric_error_only": True, "ledger_written": False, "writes_ledger": False, "no_ledger": True, "ledger_join_claimed": False, "exact_global_equality_claim": False, "old_j2_error_bound_remains": True, "old_j2_error_bound_cancelled": False, "ancestor_imports": False, "point_values_selected": False, "midpoint_values_used": False, "physical_swirl": "m15", "physical_factor": physical["physical_factor"], "physical_factor_formula": physical["physical_factor_formula"], "formal_physical_shift": list(PHYSICAL_POWER_SHIFT), "finite_j2_formal_hb_power": 6, "j_ge_3_formal_hb_power": 8, "j_ge_3_physical_hb_power": 9, "q_domain": cert["q_domain"], "q0_value": _zero(c), "q0_exact_zero": True}


def build_local_second_potential_m15(c: Any, *args: Any, **kwargs: Any) -> dict[str, Any]:
    return build(c, *args, **kwargs)


build_local_second_potential = build
enclose = build
evaluate = build
