"""Independent local-g5 potential projection leaf for the physical m15 term.

This module is intentionally a data-only leaf.  It accepts caller-owned
directed interval arithmetic and already saved source sectors; it does not
import or call an ancestor, replay a source provider, run a ledger, or merge
the pressure/local-G0g results.  The public entry point is::

    build(c, sc, local_g5_j1, j2_tail, leading_F, certificate)

Each source argument is a mapping with ``Z0`` and ``Z1`` sector lists (an
optional outer ``F1`` key is accepted).  A sector carries ``coefficient``,
``powers`` in the order ``[hb, F0, Pstar^2, Lambda]``, ``source_ids``,
``source_paths``, ``error_classes``, and ``error_only``.  ``certificate``
must carry the caller's ``hb`` interval and explicit proofs that the source
chi/g enclosure is valid on all ``q in [0, 1]`` and that one exponent bound is
uniform there.

The projected m15 integrand is

    -2 Lambda^-2 integral eta rho Flead F1 d rho,
    rho = 4 exp(hb sc q/4),  d rho = hb (sc/4) rho d q,
    integral_0^1 eta(q) dq = 5/8.

The volume bound is ``rho^2 in [16, 16 exp(sc/4)]``.  After the factor 16 is
factored out, the physical interval is ``-5 sc [1, exp(sc/4)]`` and the
formal power shift is ``[+1, 0, 0, -2]``.  J1 endpoint sectors are finite
signed intervals and are extended over the phase by ``hull(0, endpoint)``;
J2 sectors are pre-existing whole-collar symmetric errors and are copied as
such.  No point representative is selected.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


__all__ = [
    "FIELDS",
    "POWER_BASIS",
    "ETA_NUMERATOR",
    "ETA_DENOMINATOR",
    "build",
    "build_local_g5_m15",
]


FIELDS = ("F1",)
POWER_BASIS = ("hb", "F0", "Pstar^2", "Lambda")
ETA_NUMERATOR = 5
ETA_DENOMINATOR = 8
POWER_SHIFT = (1, 0, 0, -2)

_MISSING = object()
_Z_KEYS = (("Z0", "z0", 0, "0"), ("Z1", "z1", 1, "1"))


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _lookup(mapping: Mapping[Any, Any], *keys: Any) -> Any:
    for key in keys:
        try:
            if key in mapping:
                return mapping[key]
        except (TypeError, ValueError):
            continue
    return _MISSING


def _require(mapping: Mapping[Any, Any], label: str, *keys: Any) -> Any:
    value = _lookup(mapping, *keys)
    if value is _MISSING or value is None:
        raise KeyError(f"{label} is missing one of {keys!r}")
    return value


def _point_from_raw(c: Any, raw: Any) -> Any:
    """Make a point interval from an mpmath raw endpoint when possible."""

    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((raw, raw))
        except Exception:
            pass
    for pair in ((raw, raw), [raw, raw]):
        try:
            return c.mpf(pair)
        except Exception:
            pass
    return c.mpf(raw)


def _endpoint(c: Any, value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return _point_from_raw(c, raw[1 if upper else 0])
    if hasattr(value, "a") and hasattr(value, "b"):
        candidate = value.b if upper else value.a
        nested_raw = getattr(candidate, "_mpi_", None)
        if nested_raw is not None:
            return _point_from_raw(c, nested_raw[1 if upper else 0])
        return c.mpf(candidate)
    return c.mpf(value)


def _interval(c: Any, lower: Any, upper: Any) -> Any:
    """Construct an interval while retaining both directed endpoints."""

    if _is_interval(lower):
        lower = _endpoint(c, lower, False)
    if _is_interval(upper):
        upper = _endpoint(c, upper, True)
    for pair in ([lower, upper], (lower, upper)):
        try:
            return c.mpf(pair)
        except Exception:
            pass
    raise TypeError("directed context cannot construct an interval")


def _as_interval(c: Any, value: Any, label: str) -> Any:
    if _is_interval(value):
        return value
    if _is_mapping(value):
        nested = _lookup(value, "interval", "enclosure", "record")
        if nested is not _MISSING:
            return _as_interval(c, nested, label)
        lower = _lookup(value, "lower", "lo")
        upper = _lookup(value, "upper", "hi")
        if lower is not _MISSING and upper is not _MISSING:
            return _interval(c, lower, upper)
        raise TypeError(f"{label} is not an interval record")
    if _is_sequence(value) and len(value) == 2:
        return _interval(c, value[0], value[1])
    if not hasattr(c, "mpf"):
        raise TypeError("c must provide directed interval mpf arithmetic")
    try:
        return c.mpf(value)
    except Exception as exc:
        raise TypeError(f"{label} cannot be converted by the caller context") from exc


def _zero(c: Any) -> Any:
    return c.mpf(0)


def _one(c: Any) -> Any:
    return c.mpf(1)


def _half(c: Any) -> Any:
    return _one(c) / c.mpf(2)


def _strictly_positive(c: Any, value: Any) -> bool:
    return bool(_endpoint(c, value, False) > _zero(c))


def _nonnegative(c: Any, value: Any) -> bool:
    return bool(_endpoint(c, value, False) >= _zero(c))


def _strictly_less(c: Any, value: Any, bound: Any) -> bool:
    return bool(_endpoint(c, value, True) < bound)


def _exp(c: Any, value: Any) -> Any:
    return c.exp(value)


def _symmetrize(c: Any, value: Any) -> Any:
    """Return the symmetric error hull of a directed interval."""

    lo = _endpoint(c, value, False)
    hi = _endpoint(c, value, True)
    z = _zero(c)
    alo = -lo if bool(lo < z) else lo
    ahi = -hi if bool(hi < z) else hi
    radius = alo if bool(alo >= ahi) else ahi
    return _interval(c, -radius, radius)


def _is_symmetric(c: Any, value: Any) -> bool:
    lo = _endpoint(c, value, False)
    hi = _endpoint(c, value, True)
    return bool(lo == -hi)


def _hull_zero(c: Any, value: Any) -> Any:
    """Directed ``hull({0}, value)``; no midpoint or point value is used."""

    lo = _endpoint(c, value, False)
    hi = _endpoint(c, value, True)
    z = _zero(c)
    lower = lo if bool(lo <= z) else z
    upper = hi if bool(hi >= z) else z
    return _interval(c, lower, upper)


def _source_list(value: Any, label: str) -> list[Mapping[str, Any]]:
    if _is_mapping(value):
        nested = _lookup(value, "sectors", "terms", "source_sectors", "records")
        if nested is not _MISSING:
            return _source_list(nested, label)
        if _lookup(value, "coefficient", "coefficient_interval") is not _MISSING:
            return [value]
        values = list(value.values())
        if values and all(_is_mapping(item) for item in values):
            return values
        raise TypeError(f"{label} is not a sector or sector sequence")
    if _is_sequence(value):
        result = list(value)
        if not all(_is_mapping(item) for item in result):
            raise TypeError(f"{label} must contain mappings")
        return result
    raise TypeError(f"{label} must be a mapping or a sequence of mappings")


def _z_sources(value: Any, label: str) -> dict[str, list[Mapping[str, Any]]]:
    current = value
    if _is_mapping(current):
        nested = _lookup(current, "F1", "f1", "field", "field_F1", "F1_sources")
        if nested is not _MISSING:
            current = nested
    if _is_sequence(current) and len(current) == 2:
        is_one_sector_list = all(
            _is_mapping(item)
            and _lookup(item, "coefficient", "coefficient_interval") is not _MISSING
            for item in current
        )
        if not is_one_sector_list:
            return {"Z0": _source_list(current[0], f"{label}.Z0"),
                    "Z1": _source_list(current[1], f"{label}.Z1")}
    if not _is_mapping(current):
        raise TypeError(f"{label} must contain Z0 and Z1 sectors")
    result: dict[str, list[Mapping[str, Any]]] = {}
    for canonical, keys in zip(("Z0", "Z1"), _Z_KEYS):
        selected = _lookup(current, *keys)
        if selected is _MISSING:
            raise KeyError(f"{label} is missing {canonical}")
        result[canonical] = _source_list(selected, f"{label}.{canonical}")
    return result


def _powers(record: Mapping[str, Any], label: str) -> tuple[int, int, int, int]:
    raw = _lookup(record, "powers", "formal_powers", "power_tuple")
    if raw is _MISSING:
        raise KeyError(f"{label} is missing powers")
    if _is_mapping(raw):
        raw = [
            _require(raw, label, "hb", 0),
            _require(raw, label, "F0", 1),
            _require(raw, label, "Pstar^2", "Pstar2", 2),
            _require(raw, label, "Lambda", 3),
        ]
    if not _is_sequence(raw) or len(raw) != 4:
        raise TypeError(f"{label}.powers must have four entries")
    result = []
    for index, item in enumerate(raw):
        if not isinstance(item, int) or isinstance(item, bool):
            raise TypeError(f"{label}.powers[{index}] must be an integer")
        result.append(item)
    return tuple(result)  # type: ignore[return-value]


def _labels(record: Mapping[str, Any], label: str) -> dict[str, Any]:
    origin = record.get("origin") if _is_mapping(record.get("origin")) else {}
    ids = _lookup(record, "source_ids", "source_id")
    if ids is _MISSING:
        ids = _lookup(origin, "source_ids", "source_id")
    if ids is _MISSING:
        raise KeyError(f"{label} is missing source_ids")
    if isinstance(ids, (str, bytes, bytearray)) or not _is_sequence(ids):
        ids = [ids]
    else:
        ids = list(ids)
    paths = _lookup(record, "source_paths", "source_path")
    if paths is _MISSING:
        paths = _lookup(origin, "source_paths", "source_path")
    if paths is _MISSING:
        paths = []
    elif isinstance(paths, (str, bytes, bytearray)) or not _is_sequence(paths):
        paths = [paths]
    else:
        paths = list(paths)
    classes = _lookup(record, "error_classes", "error_class")
    if classes is _MISSING:
        classes = _lookup(origin, "error_classes", "error_class")
    if classes is _MISSING or classes is None:
        classes = []
    elif isinstance(classes, (str, bytes, bytearray)) or not _is_sequence(classes):
        classes = [classes]
    else:
        classes = list(classes)
    error_only = record.get("error_only", False)
    if not isinstance(error_only, bool):
        raise TypeError(f"{label}.error_only must be boolean")
    return {
        "source_ids": list(ids),
        "source_paths": list(paths),
        "error_classes": list(classes),
        "error_only": error_only,
        "uncertainty": bool(record.get("uncertainty", origin.get("uncertainty", False)) or error_only),
    }


def _sector(c: Any, record: Mapping[str, Any], label: str) -> dict[str, Any]:
    coefficient = _lookup(record, "coefficient", "coefficient_interval")
    if coefficient is _MISSING:
        raise KeyError(f"{label} is missing coefficient")
    result = dict(record)
    result["coefficient"] = _as_interval(c, coefficient, f"{label}.coefficient")
    result["coefficient_interval"] = result["coefficient"]
    result["powers"] = list(_powers(record, label))
    result.update(_labels(record, label))
    return result


def _certificate_flag(certificate: Mapping[str, Any], names: Sequence[str]) -> bool:
    nested = _lookup(certificate, "source_certificate", "proof", "proofs", "uniformity")
    candidates: list[Mapping[str, Any]] = [certificate]
    if _is_mapping(nested):
        candidates.append(nested)
    for current in candidates:
        for name in names:
            value = _lookup(current, name)
            if value is _MISSING:
                continue
            if value is True:
                return True
            if _is_mapping(value) and value.get("proved") is True:
                return True
    return False


def _validate_certificate(c: Any, sc: Any, certificate: Any) -> tuple[Any, dict[str, Any]]:
    if not _is_mapping(certificate):
        raise TypeError("certificate must be a mapping")
    hb_raw = _lookup(certificate, "hb", "hb_interval", "h_b")
    if hb_raw is _MISSING:
        nested = _lookup(certificate, "phase", "source_certificate")
        if _is_mapping(nested):
            hb_raw = _lookup(nested, "hb", "hb_interval", "h_b")
    if hb_raw is _MISSING:
        raise KeyError("certificate is missing hb")
    hb = _as_interval(c, hb_raw, "certificate.hb")
    strict = _strictly_positive(c, hb) and _strictly_less(c, hb, _half(c))
    closed = (bool(_endpoint(c, hb, False) >= _zero(c)) and bool(_endpoint(c, hb, True) <= _half(c)) and certificate.get("strict_hb_admitted") is True)
    if not (strict or closed):
        raise ValueError("certificate must prove 0 < hb < 1/2 (closed hull requires strict_hb_admitted)")
    if certificate.get("geometry_Z_independent") is not True or certificate.get("original_cutoff_mass_5_over_8") is not True:
        raise ValueError("fixed ordinary-Z geometry and original cutoff mass are required")
    if not _certificate_flag(
        certificate,
        ("source_chi_g_full_q_enclosure", "chi_g_full_q_enclosure",
         "full_q_source_enclosure", "source_full_q_enclosure"),
    ):
        raise ValueError("certificate must explicitly prove source chi/g enclosure for all q in [0,1]")
    if not _certificate_flag(
        certificate,
        ("unified_exponent_bound", "uniform_exponent_bound",
         "source_exponent_bound", "exponent_bound_uniform_on_q"),
    ):
        raise ValueError("certificate must explicitly prove one uniform exponent bound on q in [0,1]")
    sc_i = _as_interval(c, sc, "sc")
    if not _strictly_positive(c, sc_i):
        raise ValueError("sc must be strictly positive")
    q_domain = _lookup(certificate, "q_domain", "q_bounds", "source_q_domain")
    if q_domain is not _MISSING:
        q_interval = _as_interval(c, q_domain, "certificate.q_domain")
        if bool(_endpoint(c, q_interval, False) > _zero(c)) or bool(
            _endpoint(c, q_interval, True) < _one(c)
        ):
            raise ValueError("certificate q domain must contain [0,1]")
    eta = _lookup(certificate, "eta", "eta_mass", "original_eta_mass")
    if eta is not _MISSING:
        eta_interval = _as_interval(c, eta, "certificate.eta")
        expected = c.mpf(ETA_NUMERATOR) / c.mpf(ETA_DENOMINATOR)
        if not (bool(_endpoint(c, eta_interval, False) == expected)
                and bool(_endpoint(c, eta_interval, True) == expected)):
            raise ValueError("certificate eta must be the original mass 5/8")
    return hb, {
        "hb": hb,
        "sc": sc_i,
        "eta": c.mpf(ETA_NUMERATOR) / c.mpf(ETA_DENOMINATOR),
        "q_domain": _interval(c, _zero(c), _one(c)),
        "source_chi_g_full_q_enclosure": True,
        "unified_exponent_bound": True,
        "x_definition": "rho=4*exp(hb*sc*q/4)",
        "rho_derivative_definition": "d_rho/dq=hb*(sc/4)*rho",
    }


def _normalise_j1(c: Any, value: Any, label: str) -> dict[str, list[dict[str, Any]]]:
    raw = _z_sources(value, label)
    result: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    for z_order in ("Z0", "Z1"):
        for index, record in enumerate(raw[z_order]):
            item = _sector(c, record, f"{label}.{z_order}[{index}]")
            q = _lookup(record, "q", "q_endpoint", "phase_q")
            if q is not _MISSING:
                q_i = _as_interval(c, q, f"{label}.{z_order}[{index}].q")
                if not (bool(_endpoint(c, q_i, False) == _one(c))
                        and bool(_endpoint(c, q_i, True) == _one(c))):
                    raise ValueError(f"{label}.{z_order}[{index}] is not a q=1 endpoint sector")
                item["endpoint_q"] = q_i
            item["source_kind"] = "J1_endpoint_finite_signed"
            result[z_order].append(item)
    return result


def _normalise_j2(c: Any, value: Any, label: str) -> dict[str, list[dict[str, Any]]]:
    raw = _z_sources(value, label)
    result: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    for z_order in ("Z0", "Z1"):
        for index, record in enumerate(raw[z_order]):
            item = _sector(c, record, f"{label}.{z_order}[{index}]")
            if not item["error_only"]:
                raise ValueError(f"{label}.{z_order}[{index}] must be error_only for J2")
            if not _is_symmetric(c, item["coefficient"]):
                raise ValueError(f"{label}.{z_order}[{index}] must be a symmetric whole-collar error")
            item["source_kind"] = "J2_whole_collar_symmetric_error"
            result[z_order].append(item)
    return result


def _normalise_leading(c: Any, value: Any, label: str) -> dict[str, list[dict[str, Any]]]:
    raw = _z_sources(value, label)
    result: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    for z_order in ("Z0", "Z1"):
        for index, record in enumerate(raw[z_order]):
            item = _sector(c, record, f"{label}.{z_order}[{index}]")
            item["source_kind"] = "same_source_leading_F"
            result[z_order].append(item)
    return result


def _formal_power_product(left: Sequence[int], right: Sequence[int]) -> list[int]:
    return [left[i] + right[i] + POWER_SHIFT[i] for i in range(4)]


def _source_metadata(left: Mapping[str, Any], right: Mapping[str, Any]) -> dict[str, Any]:
    left_ids = list(left["source_ids"])
    right_ids = list(right["source_ids"])
    left_paths = list(left["source_paths"])
    right_paths = list(right["source_paths"])
    left_classes = list(left["error_classes"])
    right_classes = list(right["error_classes"])
    return {
        "source_ids": left_ids + right_ids,
        "source_paths": left_paths + right_paths,
        "error_classes": left_classes + right_classes,
        "source_ids_by_factor": {"leading_F": left_ids, "F1": right_ids},
        "source_paths_by_factor": {"leading_F": left_paths, "F1": right_paths},
        "error_classes_by_factor": {"leading_F": left_classes, "F1": right_classes},
    }


def _product_sector(
    c: Any,
    physical_factor: Any,
    leading: Mapping[str, Any],
    field: Mapping[str, Any],
    *,
    family: str,
    derivative_kind: str,
    q_hull_applied: bool,
) -> dict[str, Any]:
    lead_coefficient = leading["coefficient"]
    field_coefficient = field["coefficient"]
    if q_hull_applied:
        field_coefficient = _hull_zero(c, field_coefficient)
    raw = physical_factor * lead_coefficient * field_coefficient
    error_only = family == "J2" or leading["error_only"] or field["error_only"]
    coefficient = _symmetrize(c, raw) if error_only else raw
    powers = _formal_power_product(leading["powers"], field["powers"])
    item = {
        "coefficient": coefficient,
        "coefficient_interval": coefficient,
        "powers": powers,
        "formal_powers": powers,
        "power_basis": list(POWER_BASIS),
        "projection_power_shift": list(POWER_SHIFT),
        "leading_F_coefficient": lead_coefficient,
        "field_F1_coefficient": field_coefficient,
        "field_F1_endpoint_coefficient": field.get("coefficient") if q_hull_applied else None,
        "physical_factor": physical_factor,
        "family": family,
        "derivative_kind": derivative_kind,
        "ordinary_Z_order": 0 if derivative_kind == "Z0" else 1,
        "error_only": error_only,
        "symmetric_error": error_only,
        "q_hull_applied": q_hull_applied,
        "finite_signed_interval": family == "J1" and not error_only,
        "uncertainty": bool(leading["uncertainty"] or field["uncertainty"] or error_only),
        "point_value_selected": False,
        "chi_reapplied": False,
        "green_reapplied": False,
        "source_equalities_claimed": False,
        "source_components": {
            "leading_F": dict(leading),
            "F1": dict(field),
        },
    }
    item.update(_source_metadata(leading, field))
    return item


def _projection_geometry(c: Any, sc: Any, hb: Any) -> dict[str, Any]:
    one = _one(c)
    four = c.mpf(4)
    sixteen = c.mpf(16)
    quarter_sc = sc / four
    exp_sc_quarter = _exp(c, quarter_sc)
    exp_upper = _endpoint(c, exp_sc_quarter, True)
    exp_sc_quarter_envelope = _interval(c, one, exp_upper)
    rho2 = _interval(c, sixteen, sixteen * exp_upper)
    physical_factor = -c.mpf(5) * sc * exp_sc_quarter_envelope
    return {
        "rho_definition": "4*exp(hb*sc*q/4)",
        "rho_derivative": "hb*(sc/4)*rho",
        "rho2_interval": rho2,
        "rho2_formula": "[16,16*exp(sc/4)]",
        "exp_sc_over_4": exp_sc_quarter_envelope,
        "physical_factor": physical_factor,
        "physical_factor_formula": "-5*sc*[1,exp(sc/4)]",
        "physical_factor_negative": True,
        "eta": c.mpf(ETA_NUMERATOR) / c.mpf(ETA_DENOMINATOR),
        "eta_mass_formula": "5/8",
        "formal_power_shift": list(POWER_SHIFT),
        "formal_power_shift_formula": "hb^(+1)*Lambda^(-2)",
        "hb_bound_used_for_rho2": "0<hb<1/2",
        "hb": hb,
    }


def _family_output(
    c: Any,
    family: str,
    field_sources: Mapping[str, list[Mapping[str, Any]]],
    leading_sources: Mapping[str, list[Mapping[str, Any]]],
    physical_factor: Any,
) -> dict[str, list[dict[str, Any]]]:
    output: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    # Z0 = Flead_Z0 * field_Z0.
    for leading in leading_sources["Z0"]:
        for field in field_sources["Z0"]:
            output["Z0"].append(_product_sector(
                c, physical_factor, leading, field, family=family,
                derivative_kind="Z0", q_hull_applied=(family == "J1")))
    # Z1 = Flead_Z1 * field_Z0 + Flead_Z0 * field_Z1.  The saved field Z1
    # already contains chi_Z*g0 and chi*g1; this leaf performs no second
    # chi or Green action.
    for leading_z1 in leading_sources["Z1"]:
        for field_z0 in field_sources["Z0"]:
            output["Z1"].append(_product_sector(
                c, physical_factor, leading_z1, field_z0, family=family,
                derivative_kind="Z1_leading_Z1_times_field_Z0",
                q_hull_applied=(family == "J1")))
    for leading_z0 in leading_sources["Z0"]:
        for field_z1 in field_sources["Z1"]:
            output["Z1"].append(_product_sector(
                c, physical_factor, leading_z0, field_z1, family=family,
                derivative_kind="Z1_leading_Z0_times_field_Z1",
                q_hull_applied=(family == "J1")))
    return output


def build(
    c: Any,
    sc: Any,
    local_g5_j1: Any,
    j2_tail: Any,
    leading_F: Any,
    certificate: Any,
) -> dict[str, Any]:
    """Project saved local-g5 J1/J2 F1 sectors into the physical m15 term.

    ``local_g5_j1`` is the finite signed q=1 endpoint source.  ``j2_tail``
    is an already whole-collar symmetric error source.  ``leading_F`` is the
    same-source leading F sector.  All three have independent Z0/Z1 lists;
    no list is collapsed and no source id is deduplicated.
    """

    hb, certificate_record = _validate_certificate(c, sc, certificate)
    sc_i = certificate_record["sc"]
    j1 = _normalise_j1(c, local_g5_j1, "local_g5_j1")
    j2 = _normalise_j2(c, j2_tail, "j2_tail")
    leading = _normalise_leading(c, leading_F, "leading_F")
    geometry = _projection_geometry(c, sc_i, hb)
    physical_factor = geometry["physical_factor"]
    j1_output = _family_output(c, "J1", j1, leading, physical_factor)
    j2_output = _family_output(c, "J2", j2, leading, physical_factor)

    result = {
        "schema": "n1_local_g5_potential_physical_m15_v1",
        "component": "local_g5_potential_physical_m15",
        "source_only": True,
        "caller_owned_context": True,
        "certificate": certificate_record,
        "geometry": geometry,
        "phase": {
            "q_domain": _interval(c, _zero(c), _one(c)),
            "q_zero": _zero(c),
            "q_one": _interval(c, _one(c), _one(c)),
            "J1_q_zero_source_value": "zero endpoint included by hull(0, endpoint_interval)",
            "J1_q_one_source_value": "caller-provided finite endpoint interval",
        },
        "J1": j1_output,
        "J2": j2_output,
        "j1": j1_output,
        "j2": j2_output,
        "groups": {"J1": j1_output, "J2": j2_output},
        "groups_separate": True,
        "z_derivative_identity": "Z1=Flead_Z*field_Z0 + Flead*field_Z1",
        "z1_terms_are_not_merged": True,
        "field_Z1_chi_already_bound": True,
        "chi_reapplied": False,
        "green_reapplied": False,
        "pressure_reference_only": True,
        "local_G0g_reference_only": True,
        "pressure_and_local_G0g_merged": False,
        "ledger_212_written": False,
        "old_prefix_nonoverlap_claim": False,
        "same_source_equality_claim": False,
        "ancestor_replayed": False,
        "point_values_selected": False,
        "J1_contract": {
            "q_endpoint": _interval(c, _one(c), _one(c)),
            "phase_q_domain": _interval(c, _zero(c), _one(c)),
            "phase_q_zero": _zero(c),
            "whole_phase_coefficient": "hull(0,endpoint_interval)",
            "q_zero_is_in_hull": True,
            "positive_kernel_mass": "q^4 times a positive monotone radius factor",
            "positive_kernel_mass_monotone_in_q": True,
            "unified_exponent_envelope": "exp([-4*A,+4*A])",
            "unified_exponent_envelope_required_from_certificate": True,
            "positive_q4_mass_monotone": True,
            "finite_signed_interval": True,
            "source_chi_g_full_q_enclosure_consumed": True,
            "unified_exponent_bound_consumed": True,
        },
        "J2_contract": {
            "whole_collar": True,
            "symmetric_error_only": True,
            "input_interval_copied_before_projection": True,
            "source_coefficient_unchanged": True,
        },
        "limits": [
            "Only the physical m15 projection is emitted.",
            "Pressure and local G0g results are references only; they are not merged.",
            "No chi or Green factor is applied to the already saved field Z1 rows.",
            "J1 hulls and all coefficients remain directed intervals, never point values.",
        ],
    }
    return result


def build_local_g5_m15(
    c: Any,
    sc: Any,
    local_g5_j1: Any,
    j2_tail: Any,
    leading_F: Any,
    certificate: Any,
) -> dict[str, Any]:
    """Named alias for :func:`build`, kept explicit for callers."""

    return build(c, sc, local_g5_j1, j2_tail, leading_F, certificate)

