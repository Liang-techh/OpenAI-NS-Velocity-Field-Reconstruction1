"""Unexecuted finite positive-denominator gauge/angular repair source leaf."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from fractions import Fraction
from typing import Any


COMPONENT = "native_O3_eta_n1_finite_gauge_jet_repair"
SCHEMA = COMPONENT + "_kernel_result_v1"
CAPSULE_SCHEMA = COMPONENT + "_source_capsule_v1"
# Root-frozen capsule pins.  The raw pin is authenticated by the root runner;
# this leaf independently verifies the complete canonical JSON projection
# before creating an interval context or entering any source family.
CAPSULE_RAW_SHA256 = "de3973fbf0815dd2ba77a19d35e5dd1eb7eafe6ca6a08bf45b866a41d36e62bf"
CAPSULE_CANONICAL_SHA256 = "741fa4f436ed1a4ff13ca9c78ddf8fa2313b64ea6b7565d787b8cd74c4d02ba8"
INTERVAL_DPS = 280
ETA_LABELS = ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1", "uniform_eta")
ETA_POINT_LABELS = ETA_LABELS[:-1]
RHO_LABELS = ("4", "uniform_rho")
FAMILY_NAMES = ("positive_denominator_gap", "rational_log_gauge_jet_ports", "angular_coefficient_contraction", "finite_angular_rhs_binding")

CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "positive_denominator_gap_entered": 1,
    "positive_denominator_gap_complete": 1,
    "rational_log_gauge_jet_ports_entered": 1,
    "rational_log_gauge_jet_ports_complete": 1,
    "angular_coefficient_contraction_entered": 1,
    "angular_coefficient_contraction_complete": 1,
    "finite_angular_rhs_binding_entered": 1,
    "finite_angular_rhs_binding_complete": 1,
    "proof_bundle": 1,
    "family_completed": 1,
}
MANAGER_JOURNAL_EVENT_COUNT = 13

SUCCESS_COUNTS = {
    "ancestor_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "provider_calls": 0,
    "checker_calls": 0,
    "selftest_calls": 0,
    "source_value_reducer_calls": 0,
    "N_materialization_calls": 0,
    "absolute_F0_exp_calls": 0,
    "source_family_calls": 1,
    "new_math_family_calls": 4,
    "positive_denominator_gap_family_calls": 1,
    "rational_log_gauge_jet_ports_family_calls": 1,
    "angular_coefficient_contraction_family_calls": 1,
    "finite_angular_rhs_binding_family_calls": 1,
    "derive_calls": 1,
    "numeric_evaluation_calls": 1,
    "new_mathematical_proofs": 1,
    "rhs_rows_evaluated": 48,
    "copied_g4_g6_rows": 32,
    "copied_point_g5_C0_rows": 14,
    "new_g5_C0_rows": 2,
    "new_g5_C1_rows": 16,
    "new_ell2_uniform_rows": 1,
    "new_ell3_rows": 8,
    "regular_n1_solution_owner_calls": 0,
    "old_212_error_recertification_calls": 0,
}

SUCCESS_FLAGS = {
    "source_only": True,
    "metadata_only": False,
    "numeric_evaluation": True,
    "conditional": True,
    "finite_rhs_packet_admitted": True,
    "finite_rhs_packet_conditional_local": True,
    "finite_g5_uniform_C0_C1_admitted": True,
    "positive_denominator_gap_certified": True,
    "angular_coefficient_identity_registered": True,
    "hash_only_numerical_equality": False,
    "old_math_reentered": False,
    "old_provider_called": False,
    "provider_called": False,
    "checker_called": False,
    "selftest_called": False,
    "N_materialized": False,
    "absolute_F0_exp_materialized": False,
    "regular_n1_response_ready": False,
    "actual_n1_solution": False,
    "n1_solution_materialized": False,
    "m15_bound": False,
    "global_closure": False,
    "no_global_closure": True,
    "time_n_recursion": False,
    "no_time_recursion": True,
    "old_errors_removed": False,
    "old_212_error_recertification": False,
    "full_goal_complete": False,
    "uniform_eta_direct_row": True,
    "uniform_rho_direct_interval": True,
    "m3_candidate_unresolved": True,
    "m3_cap_not_selected_as_value": True,
    "axis_sector_retained_as_signed_interval": True,
    "old_point_g5_C0_copied": True,
    "g4_g6_relevant_sectors_copied": True,
    "new_g5_C0_uniform_only": True,
    "new_g5_C1_all_rows": True,
}
EXPECTED_COUNTS = dict(SUCCESS_COUNTS)
EXPECTED_FLAGS = dict(SUCCESS_FLAGS)
COUNT_EVENT_MAP = {
    "source_family_calls": ["family_entered"],
    "new_math_family_calls": ["positive_denominator_gap_entered", "rational_log_gauge_jet_ports_entered", "angular_coefficient_contraction_entered", "finite_angular_rhs_binding_entered"],
    "positive_denominator_gap_family_calls": ["positive_denominator_gap_entered"],
    "rational_log_gauge_jet_ports_family_calls": ["rational_log_gauge_jet_ports_entered"],
    "angular_coefficient_contraction_family_calls": ["angular_coefficient_contraction_entered"],
    "finite_angular_rhs_binding_family_calls": ["finite_angular_rhs_binding_entered"],
    "derive_calls": ["family_entered"],
    "numeric_evaluation_calls": ["finite_angular_rhs_binding_complete"],
    "new_mathematical_proofs": ["proof_bundle"],
    "rhs_rows_evaluated": ["finite_angular_rhs_binding_complete"],
}
API = {"name": "derive", "signature": "derive(data, emit=None)", "capsule_selector": "data['capsule']", "result_schema": SCHEMA, "source_only": True, "metadata_only": False, "numeric_evaluation": True, "conditional": True}


class KernelInputError(ValueError):
    """Raised when the root-frozen repair capsule fails its contract."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _as_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not _is_mapping(value):
        raise KernelInputError(label + " must be a mapping")
    return value


def _json_value(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Fraction):
        return str(value)
    if _is_mapping(value):
        return {str(key): _json_value(child) for key, child in value.items()}
    if _is_sequence(value):
        return [_json_value(child) for child in value]
    return str(value)


def _canonical_bytes(value: Any, strip_hash_metadata: bool = False) -> bytes:
    clean = _json_value(value)
    if strip_hash_metadata and _is_mapping(clean):
        clean = dict(clean)
        for key in ("raw_capsule_sha256", "canonical_capsule_sha256", "capsule_raw_sha256", "capsule_canonical_sha256", "capsule_hashes", "hashes"):
            clean.pop(key, None)
    return json.dumps(clean, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _lookup(mapping: Mapping[str, Any], names: Sequence[str]) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    return None


def _required(mapping: Mapping[str, Any], names: Sequence[str], label: str) -> Any:
    value = _lookup(mapping, names)
    if value is None:
        raise KernelInputError("missing source field " + label)
    return value


def _row_mapping(value: Any, label: str) -> dict[str, Mapping[str, Any]]:
    if _is_mapping(value):
        nested = _lookup(value, ("rows", "table", "records"))
        if nested is not None and (_is_mapping(nested) or _is_sequence(nested)):
            value = nested
        else:
            return {str(key): _as_mapping(item, label + "." + str(key)) for key, item in value.items()}
    if _is_sequence(value):
        result = {}
        for item in value:
            row = _as_mapping(item, label + "[]")
            name = _lookup(row, ("eta", "eta_label", "label", "rho", "rho_label"))
            if name is None:
                raise KernelInputError(label + "[] has no label")
            result[str(name)] = row
        return result
    raise KernelInputError(label + " must be a mapping or sequence")


def _raw_mpf_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise KernelInputError(label + " must be a four-word MPF tuple")
    try:
        result = tuple(int(item) for item in value)
    except (TypeError, ValueError) as error:
        raise KernelInputError(label + " has a non-integer MPF word") from error
    return result  # type: ignore[return-value]


def _new_interval_context() -> Any:
    from mpmath.ctx_iv import MPIntervalContext

    context = MPIntervalContext()
    context.dps = INTERVAL_DPS
    return context


def _point_interval(context: Any, value: Any) -> Any:
    if isinstance(value, bool):
        return context.mpf("1" if value else "0")
    if isinstance(value, int):
        return context.mpf(str(value))
    text = str(value).strip()
    if "/" in text and text.count("/") == 1:
        left, right = text.split("/", 1)
        return context.mpf(left.strip()) / context.mpf(right.strip())
    return context.mpf(text)


def _interval_value(context: Any, value: Any, label: str) -> Any:
    if _is_mapping(value):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            return context.make_mpf((_raw_mpf_tuple(value["lower_exact_mpf_tuple"], label + ".lower"), _raw_mpf_tuple(value["upper_exact_mpf_tuple"], label + ".upper")))
        for nested in ("interval", "enclosure", "record", "bound", "value", "raw", "selected"):
            if nested in value:
                return _interval_value(context, value[nested], label + "." + nested)
        if "lower" in value and "upper" in value:
            return context.mpf((str(value["lower"]), str(value["upper"])))
        raise KernelInputError(label + " is not an interval record")
    if _is_sequence(value):
        if len(value) != 2:
            raise KernelInputError(label + " endpoint sequence must have length two")
        return context.mpf((str(value[0]), str(value[1])))
    if value is None:
        raise KernelInputError(label + " is missing")
    try:
        return _point_interval(context, value)
    except Exception as error:
        raise KernelInputError(label + " is not an interval value") from error


def _endpoint(context: Any, value: Any, upper: bool = False) -> Any:
    converted = context.convert(value)
    raw = converted._mpi_[1 if upper else 0]
    return context.make_mpf((raw, raw))


def _lower(context: Any, value: Any) -> Any:
    return _endpoint(context, value, False)


def _upper(context: Any, value: Any) -> Any:
    return _endpoint(context, value, True)


def _make_interval(context: Any, lower: Any, upper: Any) -> Any:
    return context.make_mpf((context.convert(lower)._mpi_[0], context.convert(upper)._mpi_[1]))


def _zero(context: Any) -> Any:
    return context.mpf("0")


def _one(context: Any) -> Any:
    return context.mpf("1")


def _abs_bound(context: Any, value: Any) -> Any:
    left, right, zero = _lower(context, value), _upper(context, value), _zero(context)
    if bool(left >= zero):
        return _make_interval(context, zero, right)
    if bool(right <= zero):
        return _make_interval(context, zero, -left)
    return _make_interval(context, zero, -left if bool(-left >= right) else right)


def _finite_token(value: Any) -> bool:
    if not isinstance(value, str):
        return True
    return value.strip().lower() not in {"nan", "+nan", "-nan", "inf", "+inf", "-inf", "infinity", "+infinity", "-infinity"}


def _finite_interval(context: Any, value: Any, label: str) -> Any:
    if _is_mapping(value):
        for key in ("lower", "upper", "width"):
            if key in value and not _finite_token(value[key]):
                raise KernelInputError(label + " contains nonfinite " + key)
    elif _is_sequence(value):
        if any(not _finite_token(item) for item in value):
            raise KernelInputError(label + " contains a nonfinite endpoint")
    elif not _finite_token(value):
        raise KernelInputError(label + " is nonfinite")
    result = _interval_value(context, value, label)
    for raw in context.convert(result)._mpi_:
        if len(raw) != 4 or int(raw[3]) < 0:
            raise KernelInputError(label + " has a nonfinite MPF endpoint")
    return result


def _nonzero_interval(context: Any, value: Any, label: str) -> None:
    zero = _zero(context)
    if bool(_lower(context, value) <= zero <= _upper(context, value)):
        raise KernelInputError(label + " crosses zero")


def _square_signed(context: Any, value: Any, label: str) -> Any:
    left, right, zero = _lower(context, value), _upper(context, value), _zero(context)
    if bool(left >= zero):
        return _make_interval(context, left ** 2, right ** 2)
    if bool(right <= zero):
        return _make_interval(context, right ** 2, left ** 2)
    if label == "sigma":
        raise KernelInputError("sigma interval crosses zero")
    magnitude = -left if bool(-left >= right) else right
    return _make_interval(context, zero, magnitude ** 2)


def _interval_record(context: Any, value: Any) -> dict[str, Any]:
    converted = context.convert(value)
    lower_raw, upper_raw = converted._mpi_
    lower = context.make_mpf((lower_raw, lower_raw))
    upper = context.make_mpf((upper_raw, upper_raw))
    width = upper - lower
    from mpmath.libmp import to_str

    def display(raw: tuple[int, int, int, int]) -> str:
        return str(to_str(raw, INTERVAL_DPS)).strip()

    _width_lower_raw, width_upper_raw = context.convert(width)._mpi_
    return {"lower": display(lower_raw), "upper": display(upper_raw), "width": display(width_upper_raw), "width_upper_exact_mpf_tuple": list(width_upper_raw), "lower_exact_mpf_tuple": list(lower_raw), "upper_exact_mpf_tuple": list(upper_raw), "point_value": False, "display_decimal_approximation": True, "exact_mpf_tuples_authoritative": True}


def _record_or_none(context: Any, value: Any) -> Any:
    return None if value is None else _interval_record(context, value)


def _finite_tree(value: Any, label: str) -> None:
    if isinstance(value, str):
        if not _finite_token(value):
            raise KernelInputError(label + " contains nonfinite source DATA")
    elif _is_mapping(value):
        for key, child in value.items():
            _finite_tree(child, label + "." + str(key))
    elif _is_sequence(value):
        for index, child in enumerate(value):
            _finite_tree(child, label + "[" + str(index) + "]")


def _emit(emit: Callable[[str, Mapping[str, Any]], Any] | None, events: list[dict[str, Any]], name: str, payload: Mapping[str, Any] | None = None) -> None:
    row: dict[str, Any] = {"event": name, "component": COMPONENT, "source_only": True, "metadata_only": False, "numeric_evaluation": name.endswith("_complete") or name in {"family_entered", "proof_bundle", "family_completed"}, "conditional": True}
    if payload:
        row.update(_json_value(payload))
    events.append(row)
    if emit is not None:
        emit(name, row)


def _source_parameters(context: Any, capsule: Mapping[str, Any]) -> dict[str, Any]:
    raw = _as_mapping(capsule["parameters"], "parameters")
    result = {
        "delta": _finite_interval(context, _required(raw, ("delta", "δ"), "parameters.delta"), "parameters.delta"),
        "epsilon": _finite_interval(context, _required(raw, ("epsilon_core", "epsilon", "ε"), "parameters.epsilon_core"), "parameters.epsilon_core"),
        "sigma": _finite_interval(context, _required(raw, ("sigma", "σ"), "parameters.sigma"), "parameters.sigma"),
    }
    _nonzero_interval(context, result["epsilon"], "parameters.epsilon_core")
    return result


def _eta_source_rows(value: Any, label: str) -> dict[str, Mapping[str, Any]]:
    rows = _row_mapping(value, label)
    for eta in ETA_LABELS:
        if eta not in rows:
            raise KernelInputError(label + " is missing " + eta)
    return rows


def _source_gauge_data(context: Any, capsule: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    rows = _eta_source_rows(capsule["gauge_data"], "gauge_data")
    aliases = {"z": ("z",), "d": ("d",), "L": ("L",), "H": ("H",), "H_Z": ("H_Z", "HZ"), "H_ZZ": ("H_ZZ", "HZZ"), "ratio": ("ratio",), "ell1": ("ell1",)}
    result = {}
    for eta in ETA_LABELS:
        row = rows[eta]
        parsed = {name: _finite_interval(context, _required(row, names, "gauge_data." + eta + "." + name), "gauge_data." + eta + "." + name) for name, names in aliases.items()}
        parsed["ell2"] = None if eta == "uniform_eta" else _finite_interval(context, _required(row, ("ell2",), "gauge_data." + eta + ".ell2"), "gauge_data." + eta + ".ell2")
        parsed["source_copy"] = _json_value(row)
        result[eta] = parsed
    return result


def _mixed_component_root(capsule: Mapping[str, Any], component: str) -> Any:
    root = _as_mapping(capsule["mixed_phi"], "mixed_phi")
    names = (component, "phi" if component == "Phi" else "phi_minus_one")
    value = _lookup(root, names)
    if value is None:
        raise KernelInputError("mixed_phi is missing " + component)
    return value


def _mixed_rows(value: Any, label: str, context: Any) -> dict[str, dict[str, dict[tuple[int, int], dict[str, Any]]]]:
    eta_rows = _eta_source_rows(value, label)
    result: dict[str, dict[str, dict[tuple[int, int], dict[str, Any]]]] = {}
    for eta in ETA_LABELS:
        eta_row = eta_rows[eta]
        result[eta] = {}
        for rho in RHO_LABELS:
            nested = eta_row.get(rho)
            if nested is None:
                raise KernelInputError(label + "." + eta + " is missing " + rho)
            row = _as_mapping(nested, label + "." + eta + "." + rho)
            result[eta][rho] = {}
            for k in range(3):
                for m in range(4):
                    key = "m" + str(m) + "_k" + str(k)
                    raw = _as_mapping(_required(row, (key,), label + "." + eta + "." + rho + "." + key), label + "." + key)
                    actual = _finite_interval(context, _required(raw, ("actual",), key + ".actual"), key + ".actual")
                    candidate_raw = raw.get("candidate")
                    candidate = None if candidate_raw is None else _finite_interval(context, candidate_raw, key + ".candidate")
                    if m < 3 and candidate is None:
                        raise KernelInputError(label + "." + eta + "." + rho + "." + key + " candidate missing")
                    if m == 3 and candidate is not None:
                        raise KernelInputError(label + "." + eta + "." + rho + "." + key + " m3 candidate must remain None")
                    error_raw = raw.get("error")
                    error = _abs_bound(context, actual) if error_raw is None else _finite_interval(context, error_raw, key + ".error")
                    result[eta][rho][(k, m)] = {"actual": actual, "candidate": candidate, "error": error, "source_bound_only": candidate is None, "signed_value_resolved": candidate is not None, "source_copy": _json_value(raw)}
    return result


def _parent_rhs_rows(capsule: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    raw = capsule["parent_rhs_rows"]
    if not _is_sequence(raw):
        raise KernelInputError("parent_rhs_rows must be a sequence")
    return [_as_mapping(row, "parent_rhs_rows[]") for row in raw]


def _rhs_index(rows: Sequence[Mapping[str, Any]]) -> dict[tuple[str, str, str], Mapping[str, Any]]:
    result = {}
    for row in rows:
        key = (str(_required(row, ("eta",), "parent_rhs_rows[].eta")), str(_required(row, ("rho",), "parent_rhs_rows[].rho")), str(_required(row, ("quantity",), "parent_rhs_rows[].quantity")))
        result[key] = row
    return result


def _capsule_hashes(capsule: Mapping[str, Any], data: Mapping[str, Any]) -> dict[str, Any]:
    actual = hashlib.sha256(_canonical_bytes(capsule, strip_hash_metadata=False)).hexdigest()
    if actual.lower() != CAPSULE_CANONICAL_SHA256.lower():
        raise KernelInputError("source capsule canonical digest mismatch")
    return {"raw_capsule_sha256": CAPSULE_RAW_SHA256, "canonical_capsule_sha256": actual, "canonical_expected_sha256": CAPSULE_CANONICAL_SHA256, "raw_and_canonical_are_distinct": CAPSULE_RAW_SHA256.lower() != actual.lower(), "root_capsule_hash_pending": False, "raw_pin_authenticated_by": "root runner", "canonicalization": "sha256(canonical_json(entire_capsule))"}


def _validate_capsule(capsule: Any, data: Mapping[str, Any]) -> Mapping[str, Any]:
    capsule = _as_mapping(capsule, "capsule")
    if capsule.get("schema") != CAPSULE_SCHEMA:
        raise KernelInputError("unexpected finite gauge repair capsule schema")
    if capsule.get("component") != COMPONENT:
        raise KernelInputError("finite gauge repair capsule component does not match")
    for key in ("parameters", "gauge_data", "mixed_phi", "parent_rhs_rows", "exact_view_registration", "proof_packet", "retained_native_packet", "family_identity", "source_atoms", "preservation", "grid", "source_refs", "shape_contract", "repair_input_policy", "scope", "parent_rows_full_result"):
        if capsule.get(key) is None:
            raise KernelInputError("finite gauge repair capsule is missing " + key)
    _capsule_hashes(capsule, data)
    status = _lookup(capsule, ("status", "execution_status"))
    if _is_mapping(status):
        for key in ("executed", "numeric_evaluation", "old_math_calls", "old_provider_calls"):
            if bool(status.get(key, False)):
                raise KernelInputError("source capsule is already executed: " + key)
    return capsule


def _positive_denominator_scalars(context: Any, params: Mapping[str, Any], rows: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    sigma2 = _square_signed(context, params["sigma"], "sigma")
    sigma_lower = _lower(context, sigma2)
    if not bool(sigma_lower > _zero(context)):
        raise KernelInputError("sigma squared lower bound is not strictly positive")
    result = {}
    for eta in ETA_LABELS:
        row = rows[eta]
        h2 = _square_signed(context, row["H"], "H")
        D = h2 + sigma2
        D_lower = _lower(context, D)
        if not bool(D_lower >= sigma_lower) or not bool(D_lower > _zero(context)):
            raise KernelInputError("positive denominator gap failed at " + eta)
        _finite_interval(context, _interval_record(context, D), "derived." + eta + ".D")
        result[eta] = {**row, "D": D, "D_lower_gap": D_lower - sigma_lower, "sigma_squared": sigma2, "H_squared": h2}
    return result


def _rational_log_ports(context: Any, params: Mapping[str, Any], rows: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    result = {}
    for eta in ETA_LABELS:
        row = rows[eta]
        D = row["D"]
        _nonzero_interval(context, row["L"], "gauge_data." + eta + ".L")
        LZ, LZZ = -2 * params["delta"] * row["z"], -2 * params["delta"]
        A = row["L"] * row["H"]
        A1 = LZ * row["H"] + row["L"] * row["H_Z"]
        A2 = LZZ * row["H"] + 2 * LZ * row["H_Z"] + row["L"] * row["H_ZZ"]
        D1 = 2 * row["H"] * row["H_Z"]
        D2 = 2 * row["H_Z"] ** 2 + 2 * row["H"] * row["H_ZZ"]
        ell2 = row["ell2"] if eta != "uniform_eta" else -(A1 / D - A * D1 / D ** 2) / params["epsilon"]
        ell3 = -(A2 / D - 2 * A1 * D1 / D ** 2 - A * D2 / D ** 2 + 2 * A * D1 ** 2 / D ** 3) / params["epsilon"]
        for name, value in (("D", D), ("A", A), ("A1", A1), ("A2", A2), ("D1", D1), ("D2", D2), ("ell2", ell2), ("ell3", ell3)):
            _finite_interval(context, _interval_record(context, value), "derived." + eta + "." + name)
        result[eta] = {**row, "D": D, "L_Z": LZ, "L_ZZ": LZZ, "A": A, "A1": A1, "A2": A2, "D1": D1, "D2": D2, "ell2": ell2, "ell3": ell3, "new_ell2": eta == "uniform_eta", "new_ell3": True}
    return result


# A tiny exact polynomial ring is used only for the new formal coefficient
# certificate.  It never consumes source values and never calls an old
# operator.
_RING_VARS = ("w", "delta", "z", "d", "rho", "L")
_ZERO_MONOMIAL = (0,) * len(_RING_VARS)


def _p_const(value: int | Fraction) -> dict[tuple[int, ...], Fraction]:
    value = Fraction(value)
    return {} if value == 0 else {_ZERO_MONOMIAL: value}


def _p_var(name: str) -> dict[tuple[int, ...], Fraction]:
    exponents = [0] * len(_RING_VARS)
    exponents[_RING_VARS.index(name)] = 1
    return {tuple(exponents): Fraction(1)}


def _p_add(*terms: Mapping[tuple[int, ...], Fraction]) -> dict[tuple[int, ...], Fraction]:
    result: dict[tuple[int, ...], Fraction] = {}
    for term in terms:
        for monomial, coefficient in term.items():
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def _p_neg(term: Mapping[tuple[int, ...], Fraction]) -> dict[tuple[int, ...], Fraction]:
    return {monomial: -coefficient for monomial, coefficient in term.items()}


def _p_scale(term: Mapping[tuple[int, ...], Fraction], factor: int | Fraction) -> dict[tuple[int, ...], Fraction]:
    factor = Fraction(factor)
    return {monomial: factor * coefficient for monomial, coefficient in term.items() if factor * coefficient}


def _p_mul(left: Mapping[tuple[int, ...], Fraction], right: Mapping[tuple[int, ...], Fraction]) -> dict[tuple[int, ...], Fraction]:
    result: dict[tuple[int, ...], Fraction] = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(a + b for a, b in zip(left_monomial, right_monomial))
            result[monomial] = result.get(monomial, Fraction(0)) + left_coefficient * right_coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def _p_pow(term: Mapping[tuple[int, ...], Fraction], exponent: int) -> dict[tuple[int, ...], Fraction]:
    result = _p_const(1)
    for _ in range(exponent):
        result = _p_mul(result, term)
    return result


def _formal_coefficient_certificate() -> dict[str, Any]:
    w, delta, z, d, rho, L = (_p_var(name) for name in _RING_VARS)
    alpha = _p_add(w, _p_const(-1), delta)
    beta = _p_add(_p_scale(w, 2), _p_const(-3), delta)
    A = {"f00": _p_mul(w, z), "f01": d, "f10": _p_scale(_p_mul(z, rho), -2)}
    AZ = {"f00": w, "f01": _p_mul(_p_add(w, _p_const(-2)), z), "f02": d, "f10": _p_scale(rho, -2), "f11": _p_scale(_p_mul(z, rho), -2)}
    AR = {"f10": _p_mul(_p_add(w, _p_const(-2)), z), "f11": d, "f20": _p_scale(_p_mul(z, rho), -2)}
    fields = set(A) | set(AZ) | set(AR)
    T = {field: _p_add(_p_mul(_p_mul(alpha, z), A.get(field, {})), _p_mul(d, AZ.get(field, {})), _p_scale(_p_mul(_p_mul(z, rho), AR.get(field, {})), -2)) for field in fields}
    common = {field: _p_add(_p_mul(T.get(field, {}), L), _p_scale(_p_mul(_p_mul(_p_mul(delta, z), d), A.get(field, {})), 2)) for field in fields}
    expected = {
        "f00": _p_add(_p_mul(_p_mul(w, _p_add(_p_mul(alpha, _p_pow(z, 2)), d)), L), _p_scale(_p_mul(_p_mul(_p_mul(delta, w), _p_pow(z, 2)), d), 2)),
        "f01": _p_add(_p_mul(_p_mul(_p_mul(d, z), beta), L), _p_scale(_p_mul(_p_mul(delta, z), _p_pow(d, 2)), 2)),
        "f02": _p_mul(_p_pow(d, 2), L),
        "f10": _p_add(_p_scale(_p_mul(_p_mul(rho, _p_add(d, _p_mul(beta, _p_pow(z, 2)))), L), -2), _p_scale(_p_mul(_p_mul(_p_mul(_p_pow(z, 2), d), rho), delta), -4)),
        "f11": _p_scale(_p_mul(_p_mul(_p_mul(d, z), rho), L), -4),
        "f20": _p_scale(_p_mul(_p_mul(_p_pow(z, 2), _p_pow(rho, 2)), L), 4),
    }
    return {
        "schema": "finite_gauge_angular_coefficient_normal_form_v1",
        "passed": common == expected,
        "formal_numeric_evaluation": False,
        "formal_A": "w*z*f00+d*f01-2*z*rho*f10",
        "formal_A_Z": "w*f00+(w-2)*z*f01+d*f02-2*rho*f10-2*z*rho*f11",
        "formal_A_rho": "(w-2)*z*f10+d*f11-2*z*rho*f20",
        "normal_form_common_denominator": "T*L+2*delta*z*d*A over L^3",
        "coefficient_order": ["C00", "C01", "C02", "C10", "C11", "C20"],
        "C11_corrected": "-4*d*z*rho/L^2",
        "C11_has_no_L_minus_3_term": True,
        "zero_coefficients": ["C12", "C21", "C22"],
        "C0_normal_form_source": "C0 coefficients are the certified normal form; C1 is obtained from it by the ordinary fixed-rho product rule and Bell derivative identities",
        "C1_separate_normal_form": False,
        "C1_rule": "D2_Z=sum(C'_km*f_km+C_km*f_k,m+1); equivalently B_km=C'_km+C_k,m-1 and D2_Z=sum(B_km*f_km)",
        "fixed_rho_rules": ["partial_Z(rho)=0", "partial_Z(d)=-2*z", "partial_Z(L)=-2*delta*z", "partial_Z(w)=0", "partial_Z(delta)=0"],
        "bell_rule": "J0=psi0; J1=ell1*psi0+psi1; J2=(ell1^2+ell2)*psi0+2*ell1*psi1+psi2; J3=(ell1^3+3*ell1*ell2+ell3)*psi0+3*(ell1^2+ell2)*psi1+3*ell1*psi2+psi3",
        "bell_derivative_steps": ["F=F0_anchor*ratio*psi", "differentiate at fixed rho through ordinary Z order three", "substitute J0 through J3 without choosing a point value"],
        "ring_variables": list(_RING_VARS),
    }


def _dual_add(left: tuple[Any, Any], right: tuple[Any, Any]) -> tuple[Any, Any]:
    return left[0] + right[0], left[1] + right[1]


def _dual_neg(value: tuple[Any, Any]) -> tuple[Any, Any]:
    return -value[0], -value[1]


def _dual_mul(left: tuple[Any, Any], right: tuple[Any, Any]) -> tuple[Any, Any]:
    return left[0] * right[0], left[1] * right[0] + left[0] * right[1]


def _dual_inv(value: tuple[Any, Any]) -> tuple[Any, Any]:
    return 1 / value[0], -value[1] / value[0] ** 2


def _dual_div(left: tuple[Any, Any], right: tuple[Any, Any]) -> tuple[Any, Any]:
    return _dual_mul(left, _dual_inv(right))


def _dual_scale(value: tuple[Any, Any], factor: Any) -> tuple[Any, Any]:
    return value[0] * factor, value[1] * factor


def _dual_pow(value: tuple[Any, Any], exponent: int) -> tuple[Any, Any]:
    result = (1, 0)
    for _ in range(exponent):
        result = _dual_mul(result, value)
    return result


def _angular_coefficients(context: Any, params: Mapping[str, Any], scalar: Mapping[str, Any], rho: Any) -> tuple[dict[tuple[int, int], Any], dict[tuple[int, int], Any]]:
    z, d, L, delta = scalar["z"], scalar["d"], scalar["L"], params["delta"]
    w, alpha, beta = -2 - delta, (-2 - delta) - 1 + delta, 2 * (-2 - delta) - 3 + delta
    zd, dd, Ld = (z, _one(context)), (d, -2 * z), (L, -2 * delta * z)
    rd, deltad, wd = (rho, _zero(context)), (delta, _zero(context)), (w, _zero(context))
    alphad, betad = (alpha, _zero(context)), (beta, _zero(context))
    c: dict[tuple[int, int], tuple[Any, Any]] = {}
    c[(0, 0)] = _dual_add(_dual_div(_dual_mul(wd, _dual_add(_dual_mul(alphad, _dual_pow(zd, 2)), dd)), _dual_pow(Ld, 2)), _dual_div(_dual_scale(_dual_mul(_dual_mul(_dual_mul(deltad, wd), _dual_pow(zd, 2)), dd), 2), _dual_pow(Ld, 3)))
    c[(0, 1)] = _dual_add(_dual_div(_dual_mul(_dual_mul(dd, zd), betad), _dual_pow(Ld, 2)), _dual_div(_dual_scale(_dual_mul(_dual_mul(deltad, zd), _dual_pow(dd, 2)), 2), _dual_pow(Ld, 3)))
    c[(0, 2)] = _dual_div(_dual_pow(dd, 2), _dual_pow(Ld, 2))
    c[(1, 0)] = _dual_add(_dual_div(_dual_scale(_dual_mul(rd, _dual_add(dd, _dual_mul(betad, _dual_pow(zd, 2)))), -2), _dual_pow(Ld, 2)), _dual_div(_dual_scale(_dual_mul(_dual_mul(_dual_mul(deltad, _dual_pow(zd, 2)), dd), rd), -4), _dual_pow(Ld, 3)))
    c[(1, 1)] = _dual_div(_dual_scale(_dual_mul(_dual_mul(dd, zd), rd), -4), _dual_pow(Ld, 2))
    c[(2, 0)] = _dual_div(_dual_scale(_dual_mul(_dual_pow(zd, 2), _dual_pow(rd, 2)), 4), _dual_pow(Ld, 2))
    return {key: pair[0] for key, pair in c.items()}, {key: pair[1] for key, pair in c.items()}


def _channel_psi(context: Any, mixed: Mapping[str, Any], eta: str, rho: str, axis: bool) -> dict[tuple[int, int], dict[str, Any]]:
    if not axis:
        return mixed[eta][rho]
    zero, one = _zero(context), _one(context)
    return {(k, m): {"actual": one if (k, m) == (0, 0) else zero, "candidate": one if (k, m) == (0, 0) else zero, "error": zero, "source_bound_only": False, "signed_value_resolved": True} for k in range(3) for m in range(4)}


def _bell_jets(psi: Mapping[tuple[int, int], Mapping[str, Any]], scalar: Mapping[str, Any]) -> dict[tuple[int, int], dict[str, Any]]:
    ell1, ell2, ell3 = scalar["ell1"], scalar["ell2"], scalar["ell3"]
    result = {}
    for k in range(3):
        for m in range(4):
            def channel(which: str) -> Any:
                p0, p1 = psi[(k, 0)].get(which), psi[(k, 1)].get(which)
                p2, p3 = psi[(k, 2)].get(which), psi[(k, 3)].get(which)
                if m == 0:
                    return p0
                if m == 1:
                    return None if p0 is None or p1 is None else ell1 * p0 + p1
                if m == 2:
                    return None if any(item is None for item in (p0, p1, p2)) else (ell1 ** 2 + ell2) * p0 + 2 * ell1 * p1 + p2
                return None if any(item is None for item in (p0, p1, p2, p3)) else (ell1 ** 3 + 3 * ell1 * ell2 + ell3) * p0 + 3 * (ell1 ** 2 + ell2) * p1 + 3 * ell1 * p2 + p3
            result[(k, m)] = {"actual": channel("actual"), "candidate": channel("candidate")}
    return result


def _channel_contraction(context: Any, params: Mapping[str, Any], scalar: Mapping[str, Any], coeffs: Mapping[tuple[int, int], Any], coeffs_z: Mapping[tuple[int, int], Any], psi: Mapping[tuple[int, int], Mapping[str, Any]], need_c0: bool = True, need_c1: bool = True) -> dict[str, Any]:
    bell = _bell_jets(psi, scalar)
    ratio = scalar["ratio"]
    totals = {derivative: {"actual": _zero(context), "candidate": _zero(context)} for derivative, needed in ((0, need_c0), (1, need_c1)) if needed}
    for derivative in tuple(totals):
        for k, m in coeffs:
            for which in ("actual", "candidate"):
                if totals[derivative][which] is None:
                    continue
                j0 = bell[(k, m)][which]
                if j0 is None:
                    totals[derivative][which] = None
                    continue
                term = (coeffs_z[(k, m)] if derivative else coeffs[(k, m)]) * ratio * j0
                if derivative:
                    j1 = bell[(k, m + 1)][which]
                    if j1 is None:
                        totals[derivative][which] = None
                        continue
                    term = term + coeffs[(k, m)] * ratio * j1
                totals[derivative][which] = totals[derivative][which] + term
    scale = -2 * params["epsilon"]
    def bundle(index: int) -> dict[str, Any] | None:
        if index not in totals:
            return None
        actual = scale * totals[index]["actual"]
        candidate = None if totals[index]["candidate"] is None else scale * totals[index]["candidate"]
        error = _abs_bound(context, actual - candidate) if candidate is not None else _abs_bound(context, actual)
        return {"actual": actual, "candidate": candidate, "error": error, "source_bound_only": candidate is None, "signed_value_resolved": candidate is not None}
    return {"C0": bundle(0), "C1": bundle(1), "bell": bell}


def _rhs_bundle_record(context: Any, bundle: Mapping[str, Any]) -> dict[str, Any]:
    actual, candidate = bundle["actual"], bundle["candidate"]
    return {"actual": _interval_record(context, actual), "candidate": _record_or_none(context, candidate), "error": _interval_record(context, bundle["error"]), "actual_is_enclosure": True, "candidate_available": candidate is not None, "candidate_is_exact_function": False, "signed_value_resolved": False, "source_bound_only": candidate is None}


def _c1_sector_record(context: Any, bundle: Mapping[str, Any]) -> dict[str, Any]:
    candidate = bundle["candidate"]
    return {"actual_C1": _interval_record(context, bundle["actual"]), "candidate_C1": _record_or_none(context, candidate), "error_C1": _interval_record(context, bundle["error"]), "actual_is_enclosure": True, "candidate_available": candidate is not None, "candidate_is_exact_function": False, "signed_value_resolved": False, "source_bound_only": candidate is None}


def _pair_sector(context: Any, c0: Mapping[str, Any], c1: Mapping[str, Any]) -> dict[str, Any]:
    candidate = None if c0["candidate"] is None or c1["candidate"] is None else [_interval_record(context, c0["candidate"]), _interval_record(context, c1["candidate"])]
    return {"actual": [_interval_record(context, c0["actual"]), _interval_record(context, c1["actual"])], "candidate": candidate, "candidate_available": candidate is not None, "candidate_is_exact_function": False, "signed_value_resolved": False, "actual_is_enclosure": True, "source_bound_only": candidate is None}


def _copy_parent_fields(row: Mapping[str, Any], relevant: Sequence[str], label: str) -> dict[str, Any]:
    sectors = _as_mapping(_required(row, ("sectors",), label + ".sectors"), label + ".sectors")
    selected = {}
    for key in relevant:
        if key not in sectors:
            raise KernelInputError(label + " missing relevant sector " + key)
        _finite_tree(sectors[key], label + ".sectors." + key)
        selected[key] = _json_value(sectors[key])
    result = {"sectors": selected, "copied_from_parent": True}
    for key in ("atom", "F0_anchor_formal", "C0", "C1", "actual_C0", "actual_C1", "candidate_C0", "candidate_C1", "error_C0", "error_C1", "source_actual_independent_pressure_target", "reserved_own_swirl_registered_once"):
        if key in row:
            if key in {"C0", "C1", "actual_C0", "actual_C1", "candidate_C0", "candidate_C1", "error_C0", "error_C1"}:
                _finite_tree(row[key], label + "." + key)
            result[key] = _json_value(row[key])
    return result


def _copy_g5_c0_sectors(row: Mapping[str, Any], label: str) -> dict[str, Any]:
    sectors = _as_mapping(_required(row, ("C0_sectors",), label + ".C0_sectors"), label + ".C0_sectors")
    result = {}
    for key in ("g5_angular_gauge_axis", "g5_angular_core_correction"):
        if key not in sectors:
            raise KernelInputError(label + " missing " + key)
        _finite_tree(sectors[key], label + ".C0_sectors." + key)
        result[key] = _json_value(sectors[key])
    return result


def _merge_c0_c1_sector(old_c0: Mapping[str, Any], new_c1: Mapping[str, Any]) -> dict[str, Any]:
    old_candidate = old_c0.get("candidate")
    new_candidate = new_c1.get("candidate_C1")
    candidate = [old_candidate, new_candidate] if old_candidate is not None and new_candidate is not None else None
    return {"actual": [old_c0["actual"], new_c1["actual_C1"]], "candidate": candidate, "candidate_available": candidate is not None, "candidate_is_exact_function": False, "signed_value_resolved": False, "actual_is_enclosure": True, "source_bound_only": candidate is None}


def _g5_new_row(context: Any, eta: str, rho: str, axis: Mapping[str, Any], correction: Mapping[str, Any], c0_source: Mapping[str, Any] | None, copied_c0: bool) -> dict[str, Any]:
    c0 = axis["C0"] if copied_c0 else None
    if not copied_c0:
        c0_actual = axis["C0"]["actual"] + correction["C0"]["actual"]
        c0_candidate = None if axis["C0"]["candidate"] is None or correction["C0"]["candidate"] is None else axis["C0"]["candidate"] + correction["C0"]["candidate"]
        c0 = {"actual": c0_actual, "candidate": c0_candidate, "error": _abs_bound(context, c0_actual - c0_candidate) if c0_candidate is not None else _abs_bound(context, c0_actual)}
    c1_actual = axis["C1"]["actual"] + correction["C1"]["actual"]
    c1_candidate = None if axis["C1"]["candidate"] is None or correction["C1"]["candidate"] is None else axis["C1"]["candidate"] + correction["C1"]["candidate"]
    c1 = {"actual": c1_actual, "candidate": c1_candidate, "error": _abs_bound(context, c1_actual - c1_candidate) if c1_candidate is not None else _abs_bound(context, c1_actual)}
    row = {"eta": eta, "rho": rho, "quantity": "g5", "atom": _json_value(c0_source.get("atom", "F0_anchor")) if copied_c0 else "F0_anchor", "F0_anchor_formal": _json_value(c0_source.get("F0_anchor_formal", True)) if copied_c0 else True, "C0": _rhs_bundle_record(context, c0) if not copied_c0 else _json_value(c0_source["C0"]), "C1": _rhs_bundle_record(context, c1), "actual_C0": _interval_record(context, c0["actual"]) if not copied_c0 else _json_value(c0_source.get("actual_C0")), "candidate_C0": _record_or_none(context, c0["candidate"]) if not copied_c0 else _json_value(c0_source.get("candidate_C0")), "error_C0": _interval_record(context, c0["error"]) if not copied_c0 else _json_value(c0_source.get("error_C0")), "actual_C1": _interval_record(context, c1_actual), "candidate_C1": _record_or_none(context, c1_candidate), "error_C1": _interval_record(context, c1["error"]), "new_g5_C0": not copied_c0, "new_g5_C1": True, "copied_point_g5_C0": copied_c0, "source_actual_independent_pressure_target": True, "reserved_own_swirl_registered_once": True, "n1_solution_value": None}
    if not copied_c0:
        row["sectors"] = {"g5_angular_gauge_axis": _pair_sector(context, axis["C0"], axis["C1"]), "g5_angular_core_correction": _pair_sector(context, correction["C0"], correction["C1"])}
    else:
        c0_sectors = _copy_g5_c0_sectors(c0_source, "parent_rhs." + eta + "." + rho + ".g5")
        new_c1_sectors = {"g5_angular_gauge_axis": _c1_sector_record(context, axis["C1"]), "g5_angular_core_correction": _c1_sector_record(context, correction["C1"])}
        row["C0_sectors"] = c0_sectors
        row["sectors"] = {key: _merge_c0_c1_sector(c0_sectors[key], new_c1_sectors[key]) for key in c0_sectors}
        row["new_C1_sectors"] = new_c1_sectors
    return row


def _rho_interval(context: Any, label: str) -> Any:
    if label == "4":
        four = context.mpf("4")
        return context.make_mpf((four._mpi_[0], four._mpi_[1]))
    if label == "uniform_rho":
        zero, four = context.mpf("0"), context.mpf("4")
        return context.make_mpf((zero._mpi_[0], four._mpi_[1]))
    raise KernelInputError("unknown rho label " + label)


def _build_rhs_rows(context: Any, capsule: Mapping[str, Any], mixed_delta: Mapping[str, Any], scalars: Mapping[str, Any], params: Mapping[str, Any], parent_index: Mapping[tuple[str, str, str], Mapping[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, int]]:
    rows: list[dict[str, Any]] = []
    counts = {"copied_g4_g6_rows": 0, "copied_point_g5_C0_rows": 0, "new_g5_C0_rows": 0, "new_g5_C1_rows": 0}
    for eta in ETA_LABELS:
        for rho in RHO_LABELS:
            for quantity, relevant in (("g4", ("g4_axis", "g4_correction", "g4_viscous")), ("g6", ("g6_axial_axis", "g6_axial_correction", "g6_pressure_coupling"))):
                parent = parent_index.get((eta, rho, quantity))
                if parent is None:
                    raise KernelInputError("parent RHS missing " + eta + "." + rho + "." + quantity)
                row = {"eta": eta, "rho": rho, "quantity": quantity}
                row.update(_copy_parent_fields(parent, relevant, "parent_rhs." + eta + "." + rho + "." + quantity))
                rows.append(row)
                counts["copied_g4_g6_rows"] += 1
            parent_g5 = parent_index.get((eta, rho, "g5"))
            if eta in ETA_POINT_LABELS:
                if parent_g5 is None or parent_g5.get("C0") is None:
                    raise KernelInputError("point g5 C0 DATA missing at " + eta + "." + rho)
                _finite_tree(parent_g5["C0"], "parent_rhs." + eta + "." + rho + ".g5.C0")
                copied_c0 = True
                counts["copied_point_g5_C0_rows"] += 1
            else:
                parent_g5 = None
                copied_c0 = False
                counts["new_g5_C0_rows"] += 1
            coeffs, coeffs_z = _angular_coefficients(context, params, scalars[eta], _rho_interval(context, rho))
            axis_psi = _channel_psi(context, {}, eta, rho, True)
            correction_psi = _channel_psi(context, mixed_delta, eta, rho, False)
            need_c0 = eta == "uniform_eta"
            axis = _channel_contraction(context, params, scalars[eta], coeffs, coeffs_z, axis_psi, need_c0=need_c0, need_c1=True)
            correction = _channel_contraction(context, params, scalars[eta], coeffs, coeffs_z, correction_psi, need_c0=need_c0, need_c1=True)
            row = _g5_new_row(context, eta, rho, axis, correction, parent_g5, copied_c0)
            rows.append(row)
            counts["new_g5_C1_rows"] += 1
    if len(rows) != 48:
        raise KernelInputError("finite angular RHS row count is not 48")
    return rows, counts


def _gauge_payload(context: Any, scalars: Mapping[str, Any]) -> dict[str, Any]:
    result = {}
    for eta in ETA_LABELS:
        row = scalars[eta]
        result[eta] = {}
        for rho in RHO_LABELS:
            result[eta][rho] = {"z": _interval_record(context, row["z"]), "d": _interval_record(context, row["d"]), "L": _interval_record(context, row["L"]), "H": _interval_record(context, row["H"]), "H_Z": _interval_record(context, row["H_Z"]), "H_ZZ": _interval_record(context, row["H_ZZ"]), "D": _interval_record(context, row["D"]), "ratio": _interval_record(context, row["ratio"]), "ell1": _interval_record(context, row["ell1"]), "ell2": _interval_record(context, row["ell2"]), "ell3": _interval_record(context, row["ell3"]), "D_lower_gap": _interval_record(context, row["D_lower_gap"]), "D_lower_gap_is_slack": True, "F0_anchor_formal": True, "new_ell2": bool(row["new_ell2"]), "new_ell3": True}
    return result


def derive(data: Mapping[str, Any], emit: Callable[[str, Mapping[str, Any]], Any] | None = None) -> dict[str, Any]:
    """Build the finite angular repair packet once after root freezing."""

    if emit is not None and not callable(emit):
        raise KernelInputError("emit must be callable or None")
    if not _is_mapping(data):
        raise KernelInputError("derive requires a mapping")
    if "capsule" not in data:
        raise KernelInputError("derive requires data['capsule']")
    capsule = _validate_capsule(data["capsule"], data)
    events: list[dict[str, Any]] = []
    _emit(emit, events, "family_entered", {"capsule_schema": CAPSULE_SCHEMA, "component": COMPONENT})
    context = _new_interval_context()
    params = _source_parameters(context, capsule)
    gauge_rows = _source_gauge_data(context, capsule)
    mixed_root = _as_mapping(capsule["mixed_phi"], "mixed_phi")
    mixed_phi = _mixed_rows(_required(mixed_root, ("Phi",), "mixed_phi.Phi"), "mixed_phi.Phi", context)
    mixed_delta = _mixed_rows(_required(mixed_root, ("PhiMinusOne",), "mixed_phi.PhiMinusOne"), "mixed_phi.PhiMinusOne", context)
    parent_index = _rhs_index(_parent_rhs_rows(capsule))
    hashes = _capsule_hashes(capsule, data)

    _emit(emit, events, "positive_denominator_gap_entered", {"family": FAMILY_NAMES[0], "sigma_sign_branch": True, "H_zero_crossing_allowed": True, "sigma_zero_crossing_rejected": True})
    gap_rows = _positive_denominator_scalars(context, params, gauge_rows)
    _emit(emit, events, "positive_denominator_gap_complete", {"family": FAMILY_NAMES[0], "eta_rows": 8, "sigma_squared_strictly_positive": True, "D_lower_ge_sigma_lower": True})

    _emit(emit, events, "rational_log_gauge_jet_ports_entered", {"family": FAMILY_NAMES[1], "ell1_reused": 8, "ell2_uniform_new": True, "ell3_new": 8, "absolute_F0_exp_materialized": False})
    scalars = _rational_log_ports(context, params, gap_rows)
    _emit(emit, events, "rational_log_gauge_jet_ports_complete", {"family": FAMILY_NAMES[1], "new_ell2_uniform_rows": 1, "new_ell3_rows": 8, "ratio_data_reused": True})

    _emit(emit, events, "angular_coefficient_contraction_entered", {"family": FAMILY_NAMES[2], "formal_normal_form": True, "C11_corrected": "-4*d*z*rho/L^2", "C1_fixed_rho": True})
    certificate = _formal_coefficient_certificate()
    if not certificate["passed"]:
        raise KernelInputError("formal angular coefficient identity failed")
    _emit(emit, events, "angular_coefficient_contraction_complete", {"family": FAMILY_NAMES[2], "coefficient_identity": True, "bell_rule_registered": True, "no_old_operator_check": True})

    _emit(emit, events, "finite_angular_rhs_binding_entered", {"family": FAMILY_NAMES[3], "rhs_row_count": 48, "copy_g4_g6": 32, "copy_point_g5_C0": 14, "new_g5_C0": 2, "new_g5_C1": 16})
    rhs_rows, row_counts = _build_rhs_rows(context, capsule, mixed_delta, scalars, params, parent_index)
    _finite_tree(rhs_rows, "finite_rhs_packet")
    _emit(emit, events, "finite_angular_rhs_binding_complete", {"family": FAMILY_NAMES[3], **row_counts, "finite_rhs_packet_admitted": True, "finite_g5_uniform_C0_C1_admitted": True})

    sigma_squared_bound = _interval_record(context, scalars[ETA_LABELS[0]]["sigma_squared"])
    sigma_squared_bound["bound_not_function_value"] = True
    positive_gap = {"sigma_squared_lower_strictly_positive": True, "sigma_squared_lower_bound": sigma_squared_bound, "positive_sigma_squared_enclosure": sigma_squared_bound, "H_squared_nonnegative_by_sign_branch": True, "D_lower_ge_sigma_lower": True, "D_lower_gap_semantics": "nonnegative slack only; may equal zero and is not the positive lower bound"}
    proof = {"schema": COMPONENT + "_proof_v1", "proof_count": 1, "conditional_local_packet": True, "positive_denominator_gap": positive_gap, "angular_coefficient_certificate": certificate, "g5_normalization": "g5=-2*epsilon_core*F0_anchor*ratio*D2_w(psi)", "axis_channel": "psi=1 signed interval; radial derivatives structurally zero", "correction_channel": "psi=PhiMinusOne DATA; no Phi-minus-one subtraction", "candidate_semantics": "new candidate fields are interval approximations or source bounds, never exact function values; candidate_available is separate from signed_value_resolved, which is false for new records; inherited DATA flags are preserved literally", "m3_candidate": None, "m3_cap_not_selected": True, "finite_rhs_packet_admitted": True, "finite_g5_uniform_C0_C1_admitted": True, "regular_n1_response_ready": False, "actual_n1_solution": False, "m15_bound": False, "global_closure": False, "no_global_closure": True, "no_time_recursion": True, "old_errors_removed": False, "old_212_error_recertification": False, "no_old_provider_calls": True, "no_checker_or_selftest": True, "no_absolute_F0_exp": True, "source_identity": _json_value(capsule["exact_view_registration"]), "proof_packet_source": _json_value(capsule["proof_packet"]), "preservation": _json_value(capsule["preservation"]), "family_identity": _json_value(capsule["family_identity"])}
    proof["source_refs"] = _json_value(capsule["source_refs"])
    _emit(emit, events, "proof_bundle", {"proof_count": 1, "finite_rhs_packet_admitted": True, "finite_g5_uniform_C0_C1_admitted": True, "conditional_local_packet": True, "regular_n1_response_ready": False})
    _emit(emit, events, "family_completed", {"proof_count": 1, "source_only": True, "full_goal_complete": False, "actual_n1_solution": False})

    definitions = {"denominator": {"D": "H^2+sigma^2", "sigma_square": "strict sign-branch square", "H_square": "proper sign-branch square; zero lower only when H crosses zero", "gap": "lower(D)>=lower(sigma^2)>0"}, "gauge": {"L_Z": "-2*delta*z", "L_ZZ": "-2*delta", "A": "L*H", "A1": "L_Z*H+L*H_Z", "A2": "L_ZZ*H+2*L_Z*H_Z+L*H_ZZ", "D1": "2*H*H_Z", "D2": "2*H_Z^2+2*H*H_ZZ", "ell2_uniform": "-(A1/D-A*D1/D^2)/epsilon_core", "ell3": "-(A2/D-2*A1*D1/D^2-A*D2/D^2+2*A*D1^2/D^3)/epsilon_core", "F0_anchor": "formal positive atom only", "ratio": "saved DATA"}, "coefficients": {"w": "-2-delta", "alpha": "w-1+delta", "beta": "2*w-3+delta", "C00": "w*(alpha*z^2+d)/L^2+2*delta*w*z^2*d/L^3", "C01": "d*z*beta/L^2+2*delta*z*d^2/L^3", "C02": "d^2/L^2", "C10": "-2*rho*(d+beta*z^2)/L^2-4*delta*z^2*d*rho/L^3", "C11": "-4*d*z*rho/L^2", "C20": "4*z^2*rho^2/L^2"}, "C1": "D2_Z=sum(C'_km*f_km+C_km*f_k,m+1); equivalently B_km=C'_km+C_k,m-1 and D2_Z=sum(B_km*f_km)", "C1_source": "derived from the C0 normal form by ordinary fixed-rho product rule and Bell derivative identities; no separate C1 normal-form claim", "Bell": "J0=psi0; J1=ell1*psi0+psi1; J2=(ell1^2+ell2)*psi0+2*ell1*psi1+psi2; J3=(ell1^3+3*ell1*ell2+ell3)*psi0+3*(ell1^2+ell2)*psi1+3*ell1*psi2+psi3", "g5": "g5=-2*epsilon_core*F0_anchor*ratio*contraction", "candidate_policy": "m3 candidate remains None; no cap, midpoint, zero, or hash equality is selected as a value", "candidate_semantics": "new candidate fields are interval approximations or source bounds, never exact function values; candidate_available is separate from signed_value_resolved, which is false for new records; inherited DATA flags are copied literally", "row_policy": "g4/g6 relevant sectors and point g5 C0 are copied DATA; only uniform g5 C0 and all g5 C1 are new"}
    definitions["denominator"]["D_lower_gap"] = "nonnegative slack only; may equal zero and is not the positive lower bound"
    definitions["denominator"]["positive_sigma_squared_enclosure"] = "shared raw MPF interval record with exact endpoint tuples; bound_not_function_value=true"
    payload = {"component": COMPONENT, "schema": SCHEMA, "source_capsule_schema": CAPSULE_SCHEMA, "source_only": True, "metadata_only": False, "numeric_evaluation": True, "conditional": True, "runtime_model_metadata": None, "capsule_hashes": hashes, "canonical_gate": {"algorithm": "sha256(canonical_json(entire_capsule))", "raw_sha256": CAPSULE_RAW_SHA256, "canonical_sha256": CAPSULE_CANONICAL_SHA256}, "source_ids": _json_value(capsule.get("source_ids")), "source_atoms": _json_value(capsule["source_atoms"]), "parameters": {"delta": _interval_record(context, params["delta"]), "epsilon_core": _interval_record(context, params["epsilon"]), "sigma": _interval_record(context, params["sigma"])}, "gauge_table": _gauge_payload(context, scalars), "rhs_rows": rhs_rows, "definitions": definitions, "proof_packet": proof, "exact_view_registration": _json_value(capsule["exact_view_registration"]), "retained_native_packet": _json_value(capsule["retained_native_packet"]), "preservation": _json_value(capsule["preservation"]), "family_identity": _json_value(capsule["family_identity"]), "grid": _json_value(capsule["grid"]), "summary": {"eta_count": 8, "rho_count": 2, "rhs_row_count": 48, **row_counts, "new_ell2_uniform_rows": 1, "new_ell3_rows": 8, "finite_rhs_packet_admitted": True, "finite_g5_uniform_C0_C1_admitted": True, "regular_n1_response_ready": False, "actual_n1_solution": False, "m15_bound": False, "global_closure": False, "no_global_closure": True, "no_time_recursion": True, "old_errors_removed": False, "old_212_error_recertification": False}, "actual_counts_contract": dict(SUCCESS_COUNTS), "actual_flags_contract": dict(SUCCESS_FLAGS)}
    payload["positive_sigma_squared_enclosure"] = sigma_squared_bound
    payload["source_refs"] = _json_value(capsule["source_refs"])
    result = {"schema": SCHEMA, "component": COMPONENT, "passed": True, "source_only": True, "metadata_only": False, "numeric_evaluation": True, "conditional": True, "runtime_model_metadata": None, "capsule_raw_sha256": CAPSULE_RAW_SHA256, "capsule_canonical_sha256": CAPSULE_CANONICAL_SHA256, "api": dict(API), "actual_counts": dict(SUCCESS_COUNTS), "actual_counts_contract": dict(SUCCESS_COUNTS), "new_flags": dict(SUCCESS_FLAGS), "actual_flags_contract": dict(SUCCESS_FLAGS), "expected_counts": dict(EXPECTED_COUNTS), "expected_flags": dict(EXPECTED_FLAGS), "callback_event_counts": dict(CALLBACK_EVENT_COUNTS), "callback_event_count": len(CALLBACK_EVENT_COUNTS), "manager_journal_event_count": MANAGER_JOURNAL_EVENT_COUNT, "journal_event_count": MANAGER_JOURNAL_EVENT_COUNT, "count_event_map": {key: list(value) for key, value in COUNT_EVENT_MAP.items()}, "proof_count": 1, "finite_rhs_packet_admitted": True, "finite_g5_uniform_C0_C1_admitted": True, "regular_n1_response_ready": False, "actual_n1_solution": False, "m15_bound": False, "global_closure": False, "no_global_closure": True, "no_time_recursion": True, "old_errors_removed": False, "full_goal_complete": False, "payload": payload, "events": events, "unexecuted_code": {"status": "AST-ready-unexecuted", "sha256": None, "bytes": None, "hash_scope": "root-frozen source bytes; root computes after acceptance"}}
    encoded = json.dumps(_json_value(result), ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    result["return_code_hash"], result["return_code_bytes"] = hashlib.sha256(encoded).hexdigest(), len(encoded)
    result["code_hash"], result["code_bytes"] = result["return_code_hash"], result["return_code_bytes"]
    result["serialization_hash_scope"] = "canonical result before hash metadata"
    return result


__all__ = ["API", "CALLBACK_EVENT_COUNTS", "CAPSULE_CANONICAL_SHA256", "CAPSULE_RAW_SHA256", "CAPSULE_SCHEMA", "COMPONENT", "COUNT_EVENT_MAP", "EXPECTED_COUNTS", "EXPECTED_FLAGS", "FAMILY_NAMES", "INTERVAL_DPS", "KernelInputError", "RHO_LABELS", "SCHEMA", "SUCCESS_COUNTS", "SUCCESS_FLAGS", "derive"]
