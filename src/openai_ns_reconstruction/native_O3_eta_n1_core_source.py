"""Unexecuted native O3 eta-chart n=1 source leaf.

This worker is intentionally self contained.  It consumes a root-frozen
``native_O3_eta_n1_core_source_source_capsule_v1`` packet and emits a typed local source packet;
it does not import a project module, an old provider, a checker, or an
ancestor leaf.  The root runner is the only caller which may execute this
file after freezing the source capsule.

All interval calculations use a private ``MPIntervalContext``.  A scalar
endpoint is rendered from its raw MPF tuple, rather than through
``ctx.nstr(singleton_interval)``; consequently a lower or upper display can
never accidentally contain an interval literal.  The endpoint tuples, not
their decimal displays, are authoritative.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "native_O3_eta_n1_core_source_kernel_result_v1"
COMPONENT = "native_O3_eta_n1_core_source"
CAPSULE_SCHEMA = "native_O3_eta_n1_core_source_source_capsule_v1"
# Root-frozen bytes supplied for this leaf; the canonical JSON digest is
# computed separately and is deliberately not substituted for this raw pin.
CAPSULE_RAW_SHA256 = "f95b25404c4a5f4ffdc46c740353ecc5ec3130baa17da3b20fc28e7a36fcbeed"
ROOT_FREEZE_SHA = "b0df4a01383d69933d7ca2070c9358d9df13cd6e8b8dc79e65f1d0898c369ff6"
INTERVAL_DPS = 280

ETA_LABELS = ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1", "uniform_eta")
ETA_POINT_LABELS = ETA_LABELS[:-1]
RHO_LABELS = ("4", "uniform_rho")
QUANTITIES = ("g4", "g5", "g6")
ORDINARY_Z_ORDERS = (0, 1, 2, 3)
RADIAL_ORDERS = (0, 1, 2)
FAMILY_NAMES = (
    "exact_view_registration",
    "mixed_core_jets",
    "gauge_third_jet",
    "n1_differential_operators",
    "canonical_rhs_assembly",
)

CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "exact_view_registration_entered": 1,
    "exact_view_registration_complete": 1,
    "mixed_core_jets_entered": 1,
    "mixed_core_jets_complete": 1,
    "gauge_third_jet_entered": 1,
    "gauge_third_jet_complete": 1,
    "n1_differential_operators_entered": 1,
    "n1_differential_operators_complete": 1,
    "canonical_rhs_assembly_entered": 1,
    "canonical_rhs_assembly_complete": 1,
    "proof_bundle": 1,
    "family_completed": 1,
}
MANAGER_JOURNAL_EVENT_COUNT = 15

# Literal audit contracts.  The zero entries are intentional: this leaf does
# not materialize N, old/provider functions, source reducers, checkers, or
# self-tests.
SUCCESS_COUNTS = {
    "N_materialization_calls": 0,
    "ancestor_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "provider_calls": 0,
    "checker_calls": 0,
    "selftest_calls": 0,
    "source_value_reducer_calls": 0,
    "F0_absolute_exp_calls": 0,
    "source_family_calls": 1,
    "new_math_family_calls": 5,
    "exact_view_registration_family_calls": 1,
    "mixed_core_jets_family_calls": 1,
    "gauge_third_jet_family_calls": 1,
    "n1_differential_operators_family_calls": 1,
    "canonical_rhs_assembly_family_calls": 1,
    "derive_calls": 1,
    "numeric_evaluation_calls": 1,
    "new_mathematical_proofs": 1,
    "rhs_rows_evaluated": 48,
    "regular_n1_solution_owner_calls": 0,
    "old_212_error_recertification_calls": 0,
}

SUCCESS_FLAGS = {
    "source_only": True,
    "metadata_only": False,
    "numeric_evaluation": True,
    "conditional": True,
    "exact_view_registration": True,
    "conditional_leading_identity_registered": True,
    "leading_identity_uses_source_definition": True,
    "leading_identity_uses_shared_theorem": True,
    "hash_only_numerical_equality": False,
    "regular_n1_view": False,
    "regular_n1_solution_owner": False,
    "n1_solution_materialized": False,
    "n1_coefficients_solved": False,
    "global_closure": False,
    "time_n_recursion": False,
    "full_physical": False,
    "NS_claim": False,
    "full_goal": False,
    "full_goal_complete": False,
    "N_materialized": False,
    "old_math_reentered": False,
    "old_provider_called": False,
    "provider_called": False,
    "checker_called": False,
    "selftest_called": False,
    "source_value_point_selected": False,
    "source_parameter_point_selected": False,
    "F0_absolute_exp_materialized": False,
    "m3_actual_jet_candidate": False,
    "m3_source_bound_only": True,
    "uniform_eta_direct_row": True,
    "uniform_eta_sample_hull": False,
    "uniform_rho_is_real_interval": True,
    "rho_4_is_exact_real_label": True,
    "paired_mpf_interval": True,
    "raw_canonical_sha_distinct": True,
    "epsilon_applied_to_DU_once": True,
    "DU_data_already_physical": True,
    "epsilon_applied_only_to_DU_norm_bounds": True,
    "mean_primitives_consumed_without_reintegration": True,
    "no_extra_Pstar": True,
    "no_extra_half": True,
    "no_extra_sqrt2": True,
    "old_212_error_recertification": False,
    "leading_moments_appended_to_212": False,
    "independent_pressure_target_retained": True,
    "reserved_own_swirl_retained_once": True,
}

COUNT_EVENT_MAP = {
    "source_family_calls": ["family_entered"],
    "new_math_family_calls": [
        "exact_view_registration_entered",
        "mixed_core_jets_entered",
        "gauge_third_jet_entered",
        "n1_differential_operators_entered",
        "canonical_rhs_assembly_entered",
    ],
    "exact_view_registration_family_calls": ["exact_view_registration_entered"],
    "mixed_core_jets_family_calls": ["mixed_core_jets_entered"],
    "gauge_third_jet_family_calls": ["gauge_third_jet_entered"],
    "n1_differential_operators_family_calls": ["n1_differential_operators_entered"],
    "canonical_rhs_assembly_family_calls": ["canonical_rhs_assembly_entered"],
    "derive_calls": ["family_entered"],
    "numeric_evaluation_calls": ["canonical_rhs_assembly_complete"],
    "new_mathematical_proofs": ["proof_bundle"],
    "rhs_rows_evaluated": ["canonical_rhs_assembly_complete"],
}

EXPECTED_COUNTS = dict(SUCCESS_COUNTS)
EXPECTED_FLAGS = dict(SUCCESS_FLAGS)
API = {
    "name": "derive",
    "signature": "derive(data, emit=None)",
    "capsule_selector": "data['capsule']",
    "result_schema": SCHEMA,
    "source_only": True,
    "metadata_only": False,
    "numeric_evaluation": True,
    "conditional": True,
}


class KernelInputError(ValueError):
    """Raised when the root-frozen source capsule fails this leaf contract."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _as_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not _is_mapping(value):
        raise KernelInputError(label + " must be a mapping")
    return value


def _json_value(value: Any) -> Any:
    """Copy JSON-shaped source data without evaluating opaque objects."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if _is_mapping(value):
        return {str(key): _json_value(child) for key, child in value.items()}
    if _is_sequence(value):
        return [_json_value(child) for child in value]
    return str(value)


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        _json_value(value),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


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


def _emit(
    emit: Callable[[str, Mapping[str, Any]], Any] | None,
    events: list[dict[str, Any]],
    name: str,
    payload: Mapping[str, Any] | None = None,
) -> None:
    row: dict[str, Any] = {
        "event": name,
        "component": COMPONENT,
        "source_only": True,
        "metadata_only": False,
        "numeric_evaluation": name.endswith("_complete")
        or name in {"family_entered", "proof_bundle", "family_completed"},
        "conditional": True,
    }
    if payload:
        row.update(_json_value(payload))
    events.append(row)
    if emit is not None:
        # Callback ABI is exactly two positional arguments.
        emit(name, row)


def _raw_mpf_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise KernelInputError(label + " must be a four-word MPF tuple")
    try:
        words = tuple(int(item) for item in value)
    except (TypeError, ValueError) as error:
        raise KernelInputError(label + " contains a non-integer MPF word") from error
    return words  # type: ignore[return-value]


def _new_interval_context() -> Any:
    # No ambient mpmath context is used, and no project module is imported.
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
        numerator, denominator = text.split("/", 1)
        return context.mpf(numerator.strip()) / context.mpf(denominator.strip())
    return context.mpf(text)


def _interval_value(context: Any, value: Any, label: str) -> Any:
    """Read an interval, preferring paired exact MPF tuples."""

    if _is_mapping(value):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            lower = _raw_mpf_tuple(value["lower_exact_mpf_tuple"], label + ".lower")
            upper = _raw_mpf_tuple(value["upper_exact_mpf_tuple"], label + ".upper")
            return context.make_mpf((lower, upper))
        for nested in (
            "interval",
            "enclosure",
            "record",
            "bound",
            "data",
            "value",
            "raw",
            "selected",
        ):
            if nested in value:
                try:
                    return _interval_value(context, value[nested], label + "." + nested)
                except KernelInputError:
                    pass
        if "lower" in value and "upper" in value:
            return context.mpf((str(value["lower"]), str(value["upper"])))
        raise KernelInputError(label + " is not an interval record")
    if _is_sequence(value):
        if len(value) != 2:
            raise KernelInputError(label + " endpoint sequence must have length two")
        if _is_mapping(value[0]) or _is_mapping(value[1]):
            left = _interval_value(context, value[0], label + ".lower")
            right = _interval_value(context, value[1], label + ".upper")
            return _make_interval(context, _lower(context, left), _upper(context, right))
        return context.mpf((str(value[0]), str(value[1])))
    if value is None:
        raise KernelInputError(label + " is missing")
    try:
        return _point_interval(context, value)
    except Exception as error:  # pragma: no cover - defensive boundary
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
    left = context.convert(lower)._mpi_[0]
    right = context.convert(upper)._mpi_[1]
    return context.make_mpf((left, right))


def _zero(context: Any) -> Any:
    return context.mpf("0")


def _one(context: Any) -> Any:
    return context.mpf("1")


def _abs_bound(context: Any, value: Any) -> Any:
    left = _lower(context, value)
    right = _upper(context, value)
    zero = _zero(context)
    if bool(left >= zero):
        return _make_interval(context, zero, right)
    if bool(right <= zero):
        return _make_interval(context, zero, -left)
    return _make_interval(context, zero, -left if bool(-left >= right) else right)


def _square_nonnegative(context: Any, value: Any) -> Any:
    bound = _abs_bound(context, value)
    return _make_interval(context, _zero(context), _upper(context, bound) ** 2)


def _widen(context: Any, value: Any, error: Any) -> Any:
    amount = _upper(context, error)
    return _make_interval(context, _lower(context, value) - amount, _upper(context, value) + amount)


def _symmetric(context: Any, bound: Any) -> Any:
    amount = _upper(context, bound)
    return _make_interval(context, -amount, amount)


def _display_scalar(raw: tuple[int, int, int, int]) -> str:
    # ``to_str`` consumes a raw MPF tuple and cannot produce ``[a,b]``.
    from mpmath.libmp import to_str

    return str(to_str(raw, INTERVAL_DPS)).strip()


def _display_endpoint(context: Any, value: Any, upper: bool = False) -> str:
    raw = context.convert(value)._mpi_[1 if upper else 0]
    return _display_scalar(raw)


def _interval_record(context: Any, value: Any) -> dict[str, Any]:
    interval = context.convert(value)
    lower_raw, upper_raw = interval._mpi_
    lower = context.make_mpf((lower_raw, lower_raw))
    upper = context.make_mpf((upper_raw, upper_raw))
    width = upper - lower
    return {
        "lower": _display_endpoint(context, lower),
        "upper": _display_endpoint(context, upper),
        "width": _display_endpoint(context, width),
        "lower_exact_mpf_tuple": list(lower_raw),
        "upper_exact_mpf_tuple": list(upper_raw),
        "point_value": False,
    }


def _record_or_none(context: Any, value: Any) -> Any:
    return None if value is None else _interval_record(context, value)


def _row_mapping(value: Any, label: str) -> dict[str, Mapping[str, Any]]:
    if _is_mapping(value):
        nested = _lookup(value, ("rows", "table", "records"))
        if nested is not None and (_is_mapping(nested) or _is_sequence(nested)):
            value = nested
        else:
            out: dict[str, Mapping[str, Any]] = {}
            for key, item in value.items():
                alias = "uniform_eta" if str(key) in {"uniform", "eta_uniform"} else str(key)
                out[alias] = _as_mapping(item, label + "." + alias)
            return out
    if _is_sequence(value):
        out = {}
        for item in value:
            row = _as_mapping(item, label + "[]")
            name = _lookup(row, ("eta", "eta_label", "label", "rho", "rho_label"))
            if name is None:
                raise KernelInputError(label + "[] has no label")
            alias = str(name)
            if alias in {"uniform", "eta_uniform"}:
                alias = "uniform_eta"
            out[alias] = row
        return out
    raise KernelInputError(label + " must be a mapping or sequence")


def _parse_polynomial(context: Any, raw: Any, label: str) -> dict[int, Any]:
    if _is_mapping(raw):
        nested = _lookup(raw, ("coefficients", "values", "records", "items"))
        if nested is not None:
            raw = nested
        else:
            out: dict[int, Any] = {}
            for key, value in raw.items():
                try:
                    degree = int(key)
                except (TypeError, ValueError):
                    continue
                out[degree] = _interval_value(context, value, label + "." + str(degree))
            if out:
                return {degree: out[degree] for degree in sorted(out)}
    if not _is_sequence(raw):
        raise KernelInputError(label + " must be a coefficient list")
    out = {}
    for index, item in enumerate(raw):
        if _is_mapping(item):
            degree = _lookup(item, ("degree", "power", "n"))
            if degree is None:
                degree = index
            coefficient = _lookup(item, ("interval", "value", "coefficient", "record", "enclosure"))
            if coefficient is None:
                coefficient = item
        elif _is_sequence(item) and len(item) == 2:
            degree, coefficient = item
        else:
            degree, coefficient = index, item
        out[int(degree)] = _interval_value(context, coefficient, label + "[" + str(degree) + "]")
    return {degree: out[degree] for degree in sorted(out)}


def _poly_add(context: Any, *polynomials: Mapping[int, Any]) -> dict[int, Any]:
    result: dict[int, Any] = {}
    for polynomial in polynomials:
        for degree, value in polynomial.items():
            result[degree] = result.get(degree, _zero(context)) + value
    return {degree: result[degree] for degree in sorted(result)}


def _poly_scale(polynomial: Mapping[int, Any], factor: Any) -> dict[int, Any]:
    return {degree: factor * value for degree, value in polynomial.items()}


def _poly_diff(polynomial: Mapping[int, Any], order: int) -> dict[int, Any]:
    result = dict(polynomial)
    for _ in range(order):
        result = {degree - 1: degree * value for degree, value in result.items() if degree > 0}
    return {degree: result[degree] for degree in sorted(result)}


def _poly_eval(context: Any, polynomial: Mapping[int, Any], rho: Any) -> Any:
    value = _zero(context)
    for degree, coefficient in polynomial.items():
        value = value + coefficient * rho ** degree
    return value


def _poly_average(polynomial: Mapping[int, Any]) -> dict[int, Any]:
    return {degree: value / (degree + 1) for degree, value in polynomial.items()}


def _poly_drop_constant(polynomial: Mapping[int, Any]) -> dict[int, Any]:
    """Retain Phi-1 without subtracting interval-valued function records."""

    return {degree: value for degree, value in polynomial.items() if degree >= 1}


def _source_parameters(context: Any, capsule: Mapping[str, Any]) -> dict[str, Any]:
    raw = _as_mapping(capsule["parameters"], "parameters")
    aliases = {
        "delta": ("delta", "δ"),
        "epsilon": ("epsilon", "epsilon_core", "ε"),
        "h": ("h",),
        "sigma": ("sigma", "σ"),
        "j": ("j",),
        "s": ("s",),
        "a": ("a",),
        "U0a": ("U0a", "U0", "U0_a"),
        "d_a": ("d_a", "da", "dA"),
        "L_a": ("L_a", "La", "L0a"),
        "H1a": ("H1a", "H_1a", "H1_a"),
        "H2a": ("H2a", "H_2a", "H2_a"),
    }
    result = {}
    for name, names in aliases.items():
        result[name] = _interval_value(context, _required(raw, names, "parameters." + name), "parameters." + name)
    return result


def _source_polynomials(context: Any, capsule: Mapping[str, Any]) -> dict[str, dict[int, Any]]:
    raw = _as_mapping(capsule["anchor_polynomials"], "anchor_polynomials")
    aliases = {
        "Phi1": ("Phi1", "Phi", "phi", "phi1"),
        "DU1": ("DU1", "DU", "DeltaU1", "deltaU1", "du"),
        "PhiZ": ("PhiZ", "Phi1_Z", "phiZ", "phi1_Z"),
        "DUZ": ("DUZ", "DU1_Z", "DeltaU1_Z", "duZ"),
        "PhiZZ": ("PhiZZ", "Phi1_ZZ", "phiZZ", "phi1_ZZ"),
        "DUZZ": ("DUZZ", "DU1_ZZ", "DeltaU1_ZZ", "duZZ"),
    }
    result = {}
    for name, names in aliases.items():
        result[name] = _parse_polynomial(context, _required(raw, names, "anchor_polynomials." + name), "anchor_polynomials." + name)
    return result


def _source_norm_data(context: Any, capsule: Mapping[str, Any]) -> dict[str, Any]:
    raw = _as_mapping(capsule["norm_data"], "norm_data")
    aliases = {
        "E_norm": ("E_norm", "Enorm", "error_norm"),
        "Mphi": ("Mphi", "M_Phi", "Mphi_norm"),
        "Mpsi": ("Mpsi", "M_Psi", "Mpsi_norm"),
        "PhiZZ_error_existing": ("PhiZZ_error_existing", "Phi1ZZ_error", "E_PhiZZ"),
        "DUZZ_error_existing": ("DUZZ_error_existing", "DeltaU1ZZ_error", "E_DUZZ"),
    }
    result = {}
    for name, names in aliases.items():
        value = _lookup(raw, names)
        result[name] = None if value is None else _interval_value(context, value, "norm_data." + name)
    for name in ("E_norm", "Mphi", "Mpsi"):
        if result[name] is None:
            raise KernelInputError("norm_data." + name + " is missing")
    return result


def _primitive_to_average(polynomial: Mapping[int, Any]) -> dict[int, Any]:
    """Turn a saved JD primitive into Mavg by reindexing, never dividing by rho."""

    return {degree - 1: value for degree, value in polynomial.items() if degree >= 1}


def _source_mean_primitives(context: Any, capsule: Mapping[str, Any]) -> dict[str, dict[int, dict[str, Any]]]:
    raw = capsule.get("mean_primitives")
    if raw is None:
        return {}
    rows = _row_mapping(raw, "mean_primitives")
    result: dict[str, dict[int, dict[str, Any]]] = {}
    for eta, row in rows.items():
        result[eta] = {}
        for m, names in ((0, ("JD_candidate", "JD")), (1, ("JDZ_candidate", "JD_Z", "JDz"))):
            candidate_raw = _lookup(row, names)
            enclosed_raw = _lookup(row, ("JD_enclosed", "JD_actual", "JD_bound") if m == 0 else ("JDZ_enclosed", "JDZ_actual", "JDZ_bound"))
            if candidate_raw is None:
                continue
            candidate_poly = _primitive_to_average(_parse_polynomial(context, candidate_raw, "mean_primitives." + eta + ".candidate"))
            actual_poly = candidate_poly if enclosed_raw is None else _primitive_to_average(_parse_polynomial(context, enclosed_raw, "mean_primitives." + eta + ".enclosed"))
            # JD/JDZ are already physical DeltaU primitives in the frozen
            # moment DATA.  Reindex only; do not apply epsilon_core again.
            result[eta][m] = {
                "candidate": candidate_poly,
                "actual": actual_poly,
            }
    return result


def _source_eta_rows(capsule: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    rows = _row_mapping(capsule["eta_rows"], "eta_rows")
    for eta in ETA_LABELS:
        if eta not in rows:
            raise KernelInputError("eta_rows is missing " + eta)
    return rows


def _source_row_atoms(context: Any, row: Mapping[str, Any], label: str) -> dict[str, Any]:
    aliases = {
        "shift": ("shift", "eta_shift"),
        "Ubase": ("Ubase", "U_base", "u_base"),
        "ratio": ("ratio", "F0_ratio", "relative_profile"),
        "ell1": ("ell1", "ell_1", "logF0_Z"),
        "ell2_known": ("ell2_known", "ell2", "ell_2"),
    }
    result = {}
    for name, names in aliases.items():
        raw = _lookup(row, names)
        if raw is None:
            if name == "ell2_known":
                result[name] = None
                continue
            raise KernelInputError(label + " is missing " + name)
        result[name] = _interval_value(context, raw, label + "." + name)
    return result


def _validate_capsule(capsule: Any) -> Mapping[str, Any]:
    capsule = _as_mapping(capsule, "capsule")
    if capsule.get("schema") != CAPSULE_SCHEMA:
        raise KernelInputError("unexpected source capsule schema")
    if capsule.get("component") != COMPONENT:
        raise KernelInputError("source capsule component does not match")
    if len(ROOT_FREEZE_SHA) != 64 or any(char not in "0123456789abcdef" for char in ROOT_FREEZE_SHA.lower()):
        raise KernelInputError("ROOT_FREEZE_SHA is not frozen")
    actual_canonical_sha = hashlib.sha256(_canonical_bytes(capsule)).hexdigest()
    if actual_canonical_sha.lower() != ROOT_FREEZE_SHA.lower():
        raise KernelInputError("source capsule canonical digest mismatch")
    required = ("parameters", "anchor_polynomials", "norm_data", "eta_rows", "value_rows", "exact_view_contract", "source_ids", "grid")
    for key in required:
        if capsule.get(key) is None:
            raise KernelInputError("source capsule is missing " + key)
    status = _lookup(capsule, ("status", "execution_status"))
    if _is_mapping(status):
        for key in ("executed", "numeric_evaluation", "new_math_calls", "old_math_calls"):
            if bool(status.get(key, False)):
                raise KernelInputError("source capsule is already executed: " + key)
    return capsule


def _rho_interval(context: Any, label: str) -> Any:
    if label == "4":
        four = context.mpf("4")
        return context.make_mpf((four._mpi_[0], four._mpi_[1]))
    if label == "uniform_rho":
        zero = context.mpf("0")
        four = context.mpf("4")
        # Direct real [0,4], never a seven-point sample hull.
        return context.make_mpf((zero._mpi_[0], four._mpi_[1]))
    raise KernelInputError("unknown rho label " + label)


def _bound_B(context: Any, norm_value: Any, h: Any, m: int, k: int) -> Any:
    rising = 1
    for index in range(k):
        rising *= m + 1 + index
    factorial = 1
    for index in range(2, m + 1):
        factorial *= index
    t = context.mpf("1") - context.mpf("4") / context.mpf("20")
    numerator = _upper(context, norm_value) * factorial * rising
    denominator = (m + 1) ** 2 * h ** m * context.mpf("20") ** k * t ** (m + k + 1)
    return _make_interval(context, _zero(context), numerator / denominator)


def _error_transport(context: Any, norm: Mapping[str, Any], h: Any, shift: Any, component: str, m: int, k: int, epsilon: Any, existing: Any = None) -> Any:
    """Transport source errors, with physical DeltaU scaling exactly once."""

    is_du = component == "DU"
    norm_key = "Mpsi" if is_du else "Mphi"
    epsilon_abs = _abs_bound(context, epsilon)
    t3_raw = _bound_B(context, norm[norm_key], h, 3, k)
    abs_shift = _abs_bound(context, shift)
    if m == 2 and k == 0 and existing is not None:
        # Existing PhiZZ/DeltaUZZ anchor errors are already in their source
        # units (DeltaUZZ is physical).  Add the full shifted T3 tail rather
        # than treating the anchor enclosure as the whole eta error.
        t3 = epsilon_abs * t3_raw if is_du else t3_raw
        return _abs_bound(context, existing) + t3 * abs_shift
    e0 = _bound_B(context, norm["E_norm"], h, 0, k)
    e1 = _bound_B(context, norm["E_norm"], h, 1, k)
    e2 = _bound_B(context, norm["E_norm"], h, 2, k)
    shift2 = _square_nonnegative(context, shift)
    if m == 0:
        base = e0 + abs_shift * e1 + shift2 * e2 / 2 + t3_raw * abs_shift ** 3 / 6
    elif m == 1:
        base = e1 + abs_shift * e2 + t3_raw * shift2 / 2
    else:
        base = e2 + t3_raw * abs_shift
    return epsilon_abs * base if is_du else base


def _looks_triplet(value: Any) -> bool:
    if not _is_mapping(value):
        return False
    keys = {str(key).lower() for key in value}
    return bool(keys & {"actual", "candidate", "error", "enclosure", "actual_value"})


def _extract_triplet(row: Mapping[str, Any], component: str, m: int) -> Mapping[str, Any] | None:
    comp_aliases = {
        "Phi": ("phi", "Phi", "Phi1", "phi1"),
        "PhiMinusOne": ("phi_minus_one", "PhiMinusOne", "phiMinusOne", "deltaPhi"),
        "DU": ("du", "DU", "DU1", "DeltaU", "DeltaU1", "deltaU1"),
    }[component]
    channel_aliases = ("C0", "c0", "m0") if m == 0 else ("C1", "c1", "m1")
    for comp_key in comp_aliases:
        item = row.get(comp_key)
        if _is_mapping(item):
            # Current-moment DATA stores one direct component record together
            # with its ordinary-Z record as candidate_Z/error_Z/actual_Z.
            # Preserve those records directly; do not reevaluate old fields.
            if m == 1 and any(key in item for key in ("candidate_Z", "error_Z", "actual_Z")):
                return {
                    "candidate": item.get("candidate_Z"),
                    "error": item.get("error_Z"),
                    "actual": item.get("actual_Z"),
                    "source_projection": "candidate_Z/error_Z/actual_Z",
                }
            for channel in channel_aliases:
                child = item.get(channel)
                if _looks_triplet(child):
                    return child
            if _looks_triplet(item):
                return item
    for channel in channel_aliases:
        item = row.get(channel)
        if _is_mapping(item):
            for comp_key in comp_aliases:
                child = item.get(comp_key)
                if _looks_triplet(child):
                    return child
    for comp_key in comp_aliases:
        for channel in channel_aliases:
            for name in (comp_key + "_" + channel, comp_key + channel, channel + "_" + comp_key):
                item = row.get(name)
                if _looks_triplet(item):
                    return item
    item = row.get(comp_aliases[0])
    return item if _looks_triplet(item) else None


def _triplet_interval(context: Any, raw: Any, label: str) -> Any:
    if raw is None:
        return None
    try:
        return _interval_value(context, raw, label)
    except KernelInputError:
        if _is_mapping(raw) and raw.get("source_record") is not None:
            value = _interval_value(context, raw["source_record"], label + ".source_record")
            # Source DATA records, including descriptor-wrapped records, are
            # frozen physical DeltaU/Phi records.  No epsilon rescaling is
            # permitted on this input path.
            return value
        raise


def _triplet_jet(context: Any, raw: Mapping[str, Any], label: str) -> dict[str, Any]:
    actual = _triplet_interval(context, _lookup(raw, ("actual", "actual_value", "enclosure", "value")), label + ".actual")
    candidate = _triplet_interval(context, _lookup(raw, ("candidate", "candidate_value", "approximation")), label + ".candidate")
    error = _triplet_interval(context, _lookup(raw, ("error", "remainder", "bound")), label + ".error")
    if actual is None and candidate is not None:
        actual = _widen(context, candidate, error if error is not None else _zero(context))
    if actual is None:
        raise KernelInputError(label + " has no actual/enclosure record")
    if error is None:
        error = _abs_bound(context, actual - candidate) if candidate is not None else _abs_bound(context, actual)
    return {"candidate": candidate, "actual": actual, "error": _abs_bound(context, error), "source_copy": _json_value(raw), "signed_value_resolved": candidate is not None}


def _jet(context: Any, candidate: Any, actual: Any, error: Any, source_bound_only: bool = False) -> dict[str, Any]:
    if actual is None and candidate is not None:
        actual = _widen(context, candidate, error if error is not None else _zero(context))
    if actual is None:
        actual = _symmetric(context, error if error is not None else _zero(context))
    if error is None:
        error = _abs_bound(context, actual - candidate) if candidate is not None else _abs_bound(context, actual)
    return {"candidate": candidate, "actual": actual, "error": _abs_bound(context, error), "source_bound_only": source_bound_only, "signed_value_resolved": candidate is not None}


def _jet_add(context: Any, left: Mapping[str, Any], right: Mapping[str, Any]) -> dict[str, Any]:
    candidate = left["candidate"] + right["candidate"] if left["candidate"] is not None and right["candidate"] is not None else None
    return _jet(context, candidate, left["actual"] + right["actual"], left["error"] + right["error"], left.get("source_bound_only", False) or right.get("source_bound_only", False))


def _source_mixed_jets(context: Any, capsule: Mapping[str, Any], params: Mapping[str, Any], polys: Mapping[str, Mapping[int, Any]], norm: Mapping[str, Any], rows: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    value_rows = _row_mapping(capsule["value_rows"], "value_rows")
    result: dict[str, Any] = {"Phi": {}, "PhiMinusOne": {}, "DU": {}}
    half = context.mpf("1") / 2
    for eta in ETA_LABELS:
        row_atoms = _source_row_atoms(context, rows[eta], "eta_rows." + eta)
        result["Phi"][eta] = {}
        result["PhiMinusOne"][eta] = {}
        result["DU"][eta] = {}
        for rho_label in RHO_LABELS:
            rho = _rho_interval(context, rho_label)
            value_row = value_rows.get(eta, {}).get(rho_label, value_rows.get(eta, {}))
            value_row = _as_mapping(value_row, "value_rows." + eta + "." + rho_label)
            for component, base_name in (("Phi", "Phi1"), ("PhiMinusOne", "Phi1"), ("DU", "DU1")):
                result[component][eta][rho_label] = {}
                if component == "PhiMinusOne":
                    z0 = _poly_drop_constant(polys[base_name])
                    z1 = _poly_drop_constant(polys["PhiZ"])
                    z2 = _poly_drop_constant(polys["PhiZZ"])
                else:
                    z0 = polys[base_name]
                    z1 = polys["PhiZ" if component == "Phi" else "DUZ"]
                    z2 = polys["PhiZZ" if component == "Phi" else "DUZZ"]
                for k in RADIAL_ORDERS:
                    p0, p1, p2 = _poly_diff(z0, k), _poly_diff(z1, k), _poly_diff(z2, k)
                    candidate_polys = {
                        0: _poly_add(context, p0, _poly_scale(p1, row_atoms["shift"]), _poly_scale(p2, half * row_atoms["shift"] ** 2)),
                        1: _poly_add(context, p1, _poly_scale(p2, row_atoms["shift"])),
                        2: p2,
                    }
                    for m in ORDINARY_Z_ORDERS:
                        source = _extract_triplet(value_row, component, m) if k == 0 and m in (0, 1) else None
                        if source is not None:
                            source_jet = _triplet_jet(context, source, "value_rows." + eta + "." + rho_label + "." + component + ".m" + str(m))
                            # Frozen value_rows.DU are already physical
                            # DeltaU=epsilon_core*Psi records; preserve them
                            # directly, just as Phi and Phi-1 DATA are reused.
                            result[component][eta][rho_label][(k, m)] = source_jet
                            continue
                        if m == 3:
                            bound = _bound_B(context, norm["Mpsi" if component == "DU" else "Mphi"], params["h"], 3, k)
                            if component == "DU":
                                bound *= _abs_bound(context, params["epsilon"])
                            result[component][eta][rho_label][(k, m)] = _jet(context, None, _symmetric(context, bound), bound, True)
                            continue
                        existing = norm["DUZZ_error_existing" if component == "DU" else "PhiZZ_error_existing"] if m == 2 and k == 0 else None
                        error = _error_transport(context, norm, params["h"], row_atoms["shift"], component, m, k, params["epsilon"], existing)
                        candidate = _poly_eval(context, candidate_polys[m], rho)
                        result[component][eta][rho_label][(k, m)] = _jet(context, candidate, None, error)
    return result


def _build_fields(context: Any, capsule: Mapping[str, Any], params: Mapping[str, Any], polys: Mapping[str, Mapping[int, Any]], mixed: Mapping[str, Any], rows: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    fields: dict[str, Any] = {"U": {}, "DU": {}, "M": {}, "DM": {}}
    mean_primitives = _source_mean_primitives(context, capsule)
    half = context.mpf("1") / 2
    for eta in ETA_LABELS:
        atoms = _source_row_atoms(context, rows[eta], "eta_rows." + eta)
        shift = atoms["shift"]
        z_du = _poly_add(context, polys["DU1"], _poly_scale(polys["DUZ"], shift), _poly_scale(polys["DUZZ"], half * shift ** 2))
        z_duz = _poly_add(context, polys["DUZ"], _poly_scale(polys["DUZZ"], shift))
        du_polys: dict[int, dict[int, Any] | None] = {0: z_du, 1: z_duz, 2: polys["DUZZ"], 3: None}
        ubase_polys: dict[int, dict[int, Any]] = {0: {0: atoms["Ubase"]}, 1: {0: context.mpf("4")}, 2: {}, 3: {}}
        for key in fields:
            fields[key][eta] = {}
        for rho_label in RHO_LABELS:
            rho = _rho_interval(context, rho_label)
            for key in fields:
                fields[key][eta][rho_label] = {}
            for k in RADIAL_ORDERS:
                for m in ORDINARY_Z_ORDERS:
                    du = mixed["DU"][eta][rho_label][(k, m)]
                    base_candidate = _poly_eval(context, _poly_diff(ubase_polys[m], k), rho) if m <= 2 else _zero(context)
                    base = _jet(context, base_candidate, base_candidate, _zero(context))
                    fields["DU"][eta][rho_label][(k, m)] = du
                    fields["U"][eta][rho_label][(k, m)] = _jet_add(context, base, du)
                    primitive = mean_primitives.get(eta, {}).get(m)
                    if primitive is not None and m in (0, 1) and du_polys[m] is not None:
                        # For k=0, reuse the saved pointwise JD/JDZ box after
                        # symbolic division by rho (reindexing).  For k>0,
                        # differentiate only the candidate polynomial and use
                        # a fresh mixed-jet bound; enclosed/error polynomials
                        # are never treated as differentiable functions.
                        mean_candidate_poly = primitive["candidate"]
                        dm_candidate = _poly_eval(context, _poly_diff(mean_candidate_poly, k), rho)
                        if k == 0:
                            dm_actual = _poly_eval(context, primitive["actual"], rho)
                            dm_error = _abs_bound(context, dm_actual - dm_candidate)
                        else:
                            dm_error = du["error"] / (k + 1)
                            dm_actual = _widen(context, dm_candidate, dm_error)
                    else:
                        # Anchor DUZZ is already physical DeltaUZZ; average
                        # and radial-differentiate without another epsilon.
                        dm_candidate = None if du_polys[m] is None else _poly_eval(context, _poly_diff(_poly_average(du_polys[m]), k), rho)
                        dm_error = du["error"] / (k + 1)
                        dm_actual = _symmetric(context, dm_error) if dm_candidate is None else _widen(context, dm_candidate, dm_error)
                    dm = _jet(context, dm_candidate, dm_actual, dm_error, m == 3)
                    fields["DM"][eta][rho_label][(k, m)] = dm
                    base_m_candidate = _poly_eval(context, _poly_diff(_poly_average(ubase_polys[m]), k), rho) if m <= 2 else _zero(context)
                    fields["M"][eta][rho_label][(k, m)] = _jet_add(context, _jet(context, base_m_candidate, base_m_candidate, _zero(context)), dm)
    return fields


def _q0_jets(context: Any, params: Mapping[str, Any], atoms: Mapping[str, Any]) -> dict[int, dict[str, Any]]:
    z = params["a"] + atoms["shift"]
    ub = atoms["Ubase"]
    delta = params["delta"]
    d = _one(context) - z ** 2
    L = _one(context) - delta * z ** 2
    Lz, Lzz = -2 * delta * z, -2 * delta
    n = (_one(context) + delta) * z * ub - 4 * d
    nz = (_one(context) + delta) * (ub + 4 * z) + 8 * z
    nzz = 8 * (_one(context) + delta) + 8
    q0 = n / L
    q0z = nz / L - n * Lz / L ** 2
    q0zz = nzz / L - 2 * nz * Lz / L ** 2 - n * Lzz / L ** 2 + 2 * n * Lz ** 2 / L ** 3
    return {0: _jet(context, q0, q0, _zero(context)), 1: _jet(context, q0z, q0z, _zero(context)), 2: _jet(context, q0zz, q0zz, _zero(context))}


def _q_increment(
    context: Any,
    params: Mapping[str, Any],
    z: Any,
    d: Any,
    L: Any,
    du: Mapping[int, Mapping[str, Any]],
    dm: Mapping[int, Mapping[str, Any]],
    dmz: Mapping[int, Mapping[str, Any]],
    which: str,
) -> dict[int, Any]:
    """Evaluate the Q numerator and its ordinary Z derivatives through C2.

    ``du`` contains DeltaU jets, ``dm`` contains mean-DeltaU jets, and ``dmz``
    is the shifted mean table ``dm[m + 1]``.  Each ``which`` channel is
    propagated independently; in particular, a missing m3 candidate makes
    only the C2 candidate unresolved rather than manufacturing a value from
    the actual enclosure.  This is ordinary differentiation of
    ``Q = N/L`` with ``Z=(Z,1,0)``, ``d=(d,-2Z,-2)``, and
    ``L=(L,-2*delta*Z,-2*delta)``.
    """

    def value(table: Mapping[int, Mapping[str, Any]], order: int) -> Any:
        return table[order][which]

    delta = params["delta"]
    L1, L2 = -2 * delta * z, -2 * delta

    n: dict[int, Any] = {0: None, 1: None, 2: None}
    if not any(value(table, order) is None for table, order in ((du, 0), (dm, 0), (dmz, 1))):
        du0, dm0, dmz1 = value(du, 0), value(dm, 0), value(dmz, 1)
        n[0] = 2 * z * du0 - (1 - delta) * z * dm0 - d * dmz1
    if not any(value(table, order) is None for table, order in ((du, 0), (du, 1), (dm, 0), (dm, 1), (dmz, 1), (dmz, 2))):
        du0, du1 = value(du, 0), value(du, 1)
        dm0, dm1 = value(dm, 0), value(dm, 1)
        dmz1, dmz2 = value(dmz, 1), value(dmz, 2)
        n[1] = (
            2 * du0
            + 2 * z * du1
            - (1 - delta) * (dm0 + z * dmz1)
            + 2 * z * dmz1
            - d * dmz2
        )
    if not any(value(table, order) is None for table, order in ((du, 1), (du, 2), (dm, 1), (dm, 2), (dmz, 1), (dmz, 2), (dmz, 3))):
        du1, du2 = value(du, 1), value(du, 2)
        dm1, dm2 = value(dm, 1), value(dm, 2)
        dmz1, dmz2, dmz3 = value(dmz, 1), value(dmz, 2), value(dmz, 3)
        n[2] = (
            4 * du1
            + 2 * z * du2
            - (1 - delta) * (2 * dm1 + z * dmz2)
            + 2 * dmz1
            + 4 * z * dmz2
            - d * dmz3
        )

    q: dict[int, Any] = {0: None, 1: None, 2: None}
    if n[0] is not None:
        q[0] = n[0] / L
    if n[1] is not None and n[0] is not None:
        q[1] = n[1] / L - n[0] * L1 / L ** 2
    if n[2] is not None and n[1] is not None and n[0] is not None:
        q[2] = n[2] / L - 2 * n[1] * L1 / L ** 2 - n[0] * L2 / L ** 2 + 2 * n[0] * L1 ** 2 / L ** 3
    return q


def _q_bundle(context: Any, params: Mapping[str, Any], rows: Mapping[str, Mapping[str, Any]], fields: Mapping[str, Any]) -> dict[str, Any]:
    result = {}
    for eta in ETA_LABELS:
        atoms = _source_row_atoms(context, rows[eta], "eta_rows." + eta)
        result[eta] = {}
        for rho_label in RHO_LABELS:
            z = params["a"] + atoms["shift"]
            d = _one(context) - z ** 2
            L = _one(context) - params["delta"] * z ** 2
            q0 = _q0_jets(context, params, atoms)
            result[eta][rho_label] = {}
            result[eta][rho_label]["delta_q"] = {}
            for k in RADIAL_ORDERS:
                du_jets = {m: fields["DU"][eta][rho_label][(k, m)] for m in (0, 1, 2)}
                dm_jets = {m: fields["DM"][eta][rho_label][(k, m)] for m in (0, 1, 2)}
                dmz_jets = {m: fields["DM"][eta][rho_label][(k, m)] for m in (1, 2, 3)}
                actual_q = _q_increment(context, params, z, d, L, du_jets, dm_jets, dmz_jets, "actual")
                candidate_q = _q_increment(context, params, z, d, L, du_jets, dm_jets, dmz_jets, "candidate")
                for m in (0, 1, 2):
                    actual = actual_q[m]
                    candidate = candidate_q[m]
                    if actual is None:
                        raise KernelInputError("Q actual C" + str(m) + " jet is unresolved")
                    error = _abs_bound(context, actual - candidate) if candidate is not None else _abs_bound(context, actual)
                    increment = _jet(context, candidate, actual, error, candidate is None)
                    base = q0[m] if k == 0 else _jet(context, _zero(context), _zero(context), _zero(context))
                    result[eta][rho_label]["delta_q"][(k, m)] = increment
                    result[eta][rho_label][(k, m)] = _jet_add(context, base, increment)
            result[eta][rho_label]["q0"] = q0
    return result


def _sflow(context: Any, q: Mapping[tuple[int, int], Mapping[str, Any]], u: Mapping[tuple[int, int], Mapping[str, Any]], params: Mapping[str, Any], z: Any, d: Any, L: Any, rho: Any, which: str) -> tuple[Any, Any]:
    keys = ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 1))
    if any(q[key][which] is None for key in keys) or u[(0, 0)][which] is None or u[(0, 1)][which] is None:
        return None, None
    qv, qz, qzz = q[(0, 0)][which], q[(0, 1)][which], q[(0, 2)][which]
    qr, qrz, qrr, qrrz = q[(1, 0)][which], q[(1, 1)][which], q[(2, 0)][which], q[(2, 1)][which]
    uv, uz = u[(0, 0)][which], u[(0, 1)][which]
    h = _one(context) - params["delta"]
    B = qv + h * z * qz / 2 + rho * qr
    Bz = qz + h * (qz + z * qzz) / 2 + rho * qrz
    C = -2 * z * qv + d * qz - 2 * z * rho * qr
    Cz = -2 * qv - 4 * z * qz + d * qzz - 2 * rho * qr - 2 * z * rho * qrz
    S = B / L + qv * (qv / 2 + rho * qr) + uv * C / L
    Lz = -2 * params["delta"] * z
    Sz = Bz / L - B * Lz / L ** 2 + qv * qz + rho * (qz * qr + qv * qrz) + (uz * C + uv * Cz) / L - uv * C * Lz / L ** 2
    return S, Sz


def _d_add(a: tuple[Any, Any], b: tuple[Any, Any]) -> tuple[Any, Any]:
    return a[0] + b[0], a[1] + b[1]


def _d_scale(a: tuple[Any, Any], factor: Any) -> tuple[Any, Any]:
    return a[0] * factor, a[1] * factor


def _d_mul(a: tuple[Any, Any], b: tuple[Any, Any]) -> tuple[Any, Any]:
    return a[0] * b[0], a[1] * b[0] + a[0] * b[1]


def _d_inv(a: tuple[Any, Any]) -> tuple[Any, Any]:
    return 1 / a[0], -a[1] / a[0] ** 2


def _d_div(a: tuple[Any, Any], b: tuple[Any, Any]) -> tuple[Any, Any]:
    return _d_mul(a, _d_inv(b))


def _pair_scale(pair: tuple[Any, Any], factor: Any) -> tuple[Any, Any]:
    return pair[0] * factor, pair[1] * factor


def _pair_add(left: tuple[Any, Any], right: tuple[Any, Any]) -> tuple[Any, Any]:
    return left[0] + right[0], left[1] + right[1]


def _pair_bundle(context: Any, actual: tuple[Any, Any], candidate: tuple[Any, Any] | None) -> dict[str, Any]:
    return {
        "actual": [_interval_record(context, actual[0]), _interval_record(context, actual[1])],
        "candidate": None if candidate is None else [_interval_record(context, candidate[0]), _interval_record(context, candidate[1])],
        "candidate_available": candidate is not None,
        "source_bound_only": candidate is None,
    }


def _coupling_pair(context: Any, params: Mapping[str, Any], z: Any, L: Any, rho: Any, pair: tuple[Any, Any]) -> tuple[Any, Any]:
    Lz = -2 * params["delta"] * z
    return (
        -4 * params["epsilon"] * z * rho * pair[0] / L,
        -4 * params["epsilon"] * rho * ((pair[0] + z * pair[1]) / L - z * pair[0] * Lz / L ** 2),
    )


def _sector_s_pairs(context: Any, params: Mapping[str, Any], atoms: Mapping[str, Any], rho: Any, q: Mapping[Any, Any], fields: Mapping[str, Any], which: str) -> tuple[tuple[Any, Any], tuple[Any, Any], tuple[Any, Any]]:
    z = params["a"] + atoms["shift"]
    d = _one(context) - z ** 2
    L = _one(context) - params["delta"] * z ** 2
    h = _one(context) - params["delta"]
    q0 = q["q0"][0][which]
    q0z = q["q0"][1][which]
    q0zz = q["q0"][2][which]
    dq = q["delta_q"][(0, 0)][which]
    dqz = q["delta_q"][(0, 1)][which]
    dqzz = q["delta_q"][(0, 2)][which]
    qr, qrz, qrr, qrrz = q["delta_q"][(1, 0)][which], q["delta_q"][(1, 1)][which], q["delta_q"][(2, 0)][which], q["delta_q"][(2, 1)][which]
    ub = atoms["Ubase"]
    duv = fields["DU_value"][which]
    duz = fields["DU_Z_value"][which]
    ub_pair, z_pair = (ub, context.mpf("4")), (z, context.mpf("1"))
    d_pair, L_pair = (d, -2 * z), (L, -2 * params["delta"] * z)
    q0_pair, q0z_pair = (q0, q0z), (q0z, q0zz)
    dq_pair, dqz_pair = (dq, dqz), (dqz, dqzz)
    qy_pair = _d_scale((qr, qrz), rho)
    H0 = _d_add(_d_scale(_d_mul(z_pair, q0_pair), -2), _d_mul(d_pair, q0z_pair))
    c = h / 2
    axis = _d_add(_d_div(_d_add(q0_pair, _d_scale(_d_mul(z_pair, q0z_pair), c)), L_pair), _d_scale(_d_mul(q0_pair, q0_pair), context.mpf("1") / 2))
    axis = _d_add(axis, _d_div(_d_mul(ub_pair, H0), L_pair))
    HD = _d_add(_d_add(_d_scale(_d_mul(z_pair, dq_pair), -2), _d_mul(d_pair, dqz_pair)), _d_scale(_d_mul(z_pair, qy_pair), -2))
    # The first Q/L contribution is divided by L; the quadratic Q terms are
    # outside that quotient, exactly as in S_corr.
    corr = _d_div(
        _d_add(
            _d_add(dq_pair, _d_scale(_d_mul(z_pair, dqz_pair), c)),
            qy_pair,
        ),
        L_pair,
    )
    corr = _d_add(
        corr,
        _d_add(
            _d_add(_d_mul(q0_pair, dq_pair), _d_scale(_d_mul(dq_pair, dq_pair), context.mpf("1") / 2)),
            _d_add(_d_mul(q0_pair, qy_pair), _d_mul(dq_pair, qy_pair)),
        ),
    )
    corr = _d_add(corr, _d_div(_d_add(_d_add(_d_mul(ub_pair, HD), _d_mul((duv, duz), H0)), _d_mul((duv, duz), HD)), L_pair))
    viscous = (2 * qr + rho * qrr, 2 * qrz + rho * qrrz)
    return axis, corr, viscous


def _d2(context: Any, field: Mapping[tuple[int, int], Mapping[str, Any]], w: Any, params: Mapping[str, Any], z: Any, d: Any, L: Any, rho: Any, which: str, derivative: bool = False) -> Any:
    needed = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 1)]
    if derivative:
        needed += [(0, 3), (1, 2)]
    if any(field[(k, m)][which] is None for k, m in needed):
        return None
    f = lambda k, m: field[(k, m)][which]
    A = w * z * f(0, 0) + d * f(0, 1) - 2 * z * rho * f(1, 0)
    AZ = w * f(0, 0) + (w - 2) * z * f(0, 1) + d * f(0, 2) - 2 * rho * f(1, 0) - 2 * z * rho * f(1, 1)
    AR = (w - 2) * z * f(1, 0) + d * f(1, 1) - 2 * z * rho * f(2, 0)
    B = (w - 1 + params["delta"]) * z * A + d * AZ - 2 * z * rho * AR
    if not derivative:
        return B / L ** 2 + 2 * params["delta"] * z * d * A / L ** 3
    AZZ = (2 * w - 2) * f(0, 1) + (w - 4) * z * f(0, 2) + d * f(0, 3) - 4 * rho * f(1, 1) - 2 * z * rho * f(1, 2)
    ARZ = (w - 2) * f(1, 0) + (w - 4) * z * f(1, 1) + d * f(1, 2) - 2 * rho * f(2, 0) - 2 * z * rho * f(2, 1)
    BZ = (w - 1 + params["delta"]) * (A + z * AZ) - 2 * z * AZ + d * AZZ - 2 * rho * AR - 2 * z * rho * ARZ
    C = 2 * params["delta"] * z * d * A
    CZ = 2 * params["delta"] * ((d - 2 * z ** 2) * A + z * d * AZ)
    LZ = -2 * params["delta"] * z
    return BZ / L ** 2 - 2 * B * LZ / L ** 3 + CZ / L ** 3 - 3 * C * LZ / L ** 4


def _gauge_scalars(context: Any, params: Mapping[str, Any], atoms: Mapping[str, Any], eta: str) -> dict[str, Any]:
    """Build eta-only gauge scalars once; rho rows reuse these exact records."""

    shift, z = atoms["shift"], params["a"] + atoms["shift"]
    shift2 = shift ** 2
    d = params["d_a"] - 2 * params["a"] * shift - shift2
    L = params["L_a"] - 2 * params["delta"] * params["a"] * shift - params["delta"] * shift2
    H = params["H1a"] * shift + params["H2a"] * shift2 / 2 - 4 * shift ** 3
    HZ, HZZ = params["H1a"] + params["H2a"] * shift - 12 * shift2, params["H2a"] - 24 * shift
    D = _square_nonnegative(context, H) + _square_nonnegative(context, params["sigma"])
    LZ, LZZ = -2 * params["delta"] * z, -2 * params["delta"]
    A, A1, A2 = L * H, LZ * H + L * HZ, LZZ * H + 2 * LZ * HZ + L * HZZ
    D1, D2 = 2 * H * HZ, 2 * _square_nonnegative(context, HZ) + 2 * H * HZZ
    # ell1 and the seven point-eta ell2 values are direct DATA records.  Do
    # not recompute those old values; only the true uniform-eta row receives
    # the new ell2 formula.  ell3 is the one new third-jet value per eta.
    ell1_formula = None
    ell2_formula = None if eta != "uniform_eta" else -(A1 / D - A * D1 / D ** 2) / params["epsilon"]
    ell3 = -(A2 / D - 2 * A1 * D1 / D ** 2 - A * D2 / D ** 2 + 2 * A * D1 ** 2 / D ** 3) / params["epsilon"]
    ell1, ell2 = atoms["ell1"], ell2_formula if eta == "uniform_eta" else atoms["ell2_known"]
    return {"atoms": atoms, "shift": shift, "z": z, "d": d, "L": L, "H": H, "H_Z": HZ, "H_ZZ": HZZ, "D": D, "A": A, "A1": A1, "A2": A2, "D1": D1, "D2": D2, "ell1": ell1, "ell2": ell2, "ell3": ell3, "ell1_formula": ell1_formula, "ell2_formula": ell2_formula}


def _gauge(context: Any, params: Mapping[str, Any], rows: Mapping[str, Mapping[str, Any]], mixed: Mapping[str, Any], rho_label: str, eta: str, scalar_data: Mapping[str, Any] | None = None) -> tuple[dict[str, Any], dict[tuple[int, int], dict[str, Any]]]:
    if scalar_data is None:
        scalar_data = _gauge_scalars(context, params, _source_row_atoms(context, rows[eta], "eta_rows." + eta), eta)
    atoms = scalar_data["atoms"]
    shift = scalar_data["shift"]
    z, d, L = scalar_data["z"], scalar_data["d"], scalar_data["L"]
    H, HZ, HZZ = scalar_data["H"], scalar_data["H_Z"], scalar_data["H_ZZ"]
    D = scalar_data["D"]
    A, A1, A2 = scalar_data["A"], scalar_data["A1"], scalar_data["A2"]
    D1, D2 = scalar_data["D1"], scalar_data["D2"]
    ell1, ell2, ell3 = scalar_data["ell1"], scalar_data["ell2"], scalar_data["ell3"]
    ell1_formula, ell2_formula = scalar_data["ell1_formula"], scalar_data["ell2_formula"]
    gauge = {"z": z, "d": d, "L": L, "H": H, "H_Z": HZ, "H_ZZ": HZZ, "D": D, "A": A, "A1": A1, "A2": A2, "D1": D1, "D2": D2, "ell1": ell1, "ell2": ell2, "ell3": ell3, "ell1_formula": ell1_formula, "ell2_formula": ell2_formula, "ratio": atoms["ratio"], "F0_anchor_formal": True, "H_at_anchor_exact_zero": True}
    fjets, faxis, fcorr = {}, {}, {}
    ratio = atoms["ratio"]
    for k in RADIAL_ORDERS:
        for m in ORDINARY_Z_ORDERS:
            p0, p1, p2 = mixed["Phi"][eta][rho_label][(k, 0)], mixed["Phi"][eta][rho_label][(k, 1)], mixed["Phi"][eta][rho_label][(k, 2)]
            p3 = mixed["Phi"][eta][rho_label][(k, 3)]
            def formula(which: str) -> Any:
                if m == 0:
                    return ratio * p0[which]
                if m == 1:
                    return None if p1[which] is None or p0[which] is None else ratio * (ell1 * p0[which] + p1[which])
                if m == 2:
                    return None if any(p[which] is None for p in (p0, p1, p2)) else ratio * ((ell1 ** 2 + ell2) * p0[which] + 2 * ell1 * p1[which] + p2[which])
                return None if any(p[which] is None for p in (p0, p1, p2, p3)) else ratio * ((ell1 ** 3 + 3 * ell1 * ell2 + ell3) * p0[which] + 3 * (ell1 ** 2 + ell2) * p1[which] + 3 * ell1 * p2[which] + p3[which])
            candidate, actual = formula("candidate"), formula("actual")
            fjets[(k, m)] = _jet(context, candidate, actual, None, m == 3)
            # The gauge-axis coefficient is the ratio times the Bell
            # polynomial in ell; all radial derivatives of the axis sector
            # are structurally zero.  The correction uses the separately
            # retained Phi-1 mixed jets and never subtracts Phi from one.
            if k > 0:
                axis_candidate = axis_actual = _zero(context)
            elif m == 0:
                axis_candidate = axis_actual = ratio
            elif m == 1:
                axis_candidate = axis_actual = ratio * ell1
            elif m == 2:
                axis_candidate = axis_actual = ratio * (ell1 ** 2 + ell2)
            else:
                axis_candidate = axis_actual = ratio * (ell1 ** 3 + 3 * ell1 * ell2 + ell3)
            faxis[(k, m)] = _jet(context, axis_candidate, axis_actual, _zero(context))
            d0, d1, d2 = mixed["PhiMinusOne"][eta][rho_label][(k, 0)], mixed["PhiMinusOne"][eta][rho_label][(k, 1)], mixed["PhiMinusOne"][eta][rho_label][(k, 2)]
            d3 = mixed["PhiMinusOne"][eta][rho_label][(k, 3)]
            def correction(which: str) -> Any:
                if m == 0:
                    return ratio * d0[which]
                if m == 1:
                    return None if d0[which] is None or d1[which] is None else ratio * (ell1 * d0[which] + d1[which])
                if m == 2:
                    return None if any(part[which] is None for part in (d0, d1, d2)) else ratio * ((ell1 ** 2 + ell2) * d0[which] + 2 * ell1 * d1[which] + d2[which])
                return None if any(part[which] is None for part in (d0, d1, d2, d3)) else ratio * ((ell1 ** 3 + 3 * ell1 * ell2 + ell3) * d0[which] + 3 * (ell1 ** 2 + ell2) * d1[which] + 3 * ell1 * d2[which] + d3[which])
            corr_candidate, corr_actual = correction("candidate"), correction("actual")
            fcorr[(k, m)] = _jet(context, corr_candidate, corr_actual, None, m == 3)
    gauge["F_jets"], gauge["F_axis_jets"], gauge["F_correction_jets"] = fjets, faxis, fcorr
    return gauge, fjets


def _record_bundle(context: Any, item: Mapping[str, Any]) -> dict[str, Any]:
    return {"actual": _interval_record(context, item["actual"]), "candidate": _record_or_none(context, item["candidate"]), "error": _interval_record(context, item["error"]), "source_bound_only": bool(item.get("source_bound_only", False)), "signed_value_resolved": bool(item.get("signed_value_resolved", False))}


def _candidate_expr(context: Any, actual: Any, candidate: Any) -> dict[str, Any]:
    return {"actual": actual, "candidate": candidate, "error": _abs_bound(context, actual - candidate) if candidate is not None else _abs_bound(context, actual), "source_bound_only": candidate is None, "signed_value_resolved": candidate is not None}


def _proof_view(capsule: Mapping[str, Any]) -> dict[str, Any]:
    contract = _as_mapping(capsule["exact_view_contract"], "exact_view_contract")
    owner = _lookup(contract, ("old_owner", "owner_record", "old_source_owner"))
    owner = _as_mapping(owner, "exact_view_contract.old_owner") if _is_mapping(owner) else {}
    definitions = _lookup(contract, ("formal_definitions", "exact_functions", "definitions", "exact_field_definitions"))
    theorem = _lookup(contract, ("accepted_theorem_facts", "accepted_common_theorem", "accepted_common_source_theorem", "shared_theorem", "common_unique_theorem", "theorem_facts"))
    old_constructor = _lookup(contract, ("old_component_constructor_ast", "old_constructor_ast", "old_component_constructor"))
    if old_constructor is None:
        old_constructor = _lookup(owner, ("constructor_source", "constructor_ast", "old_component_constructor_ast"))
    old_component = _lookup(contract, ("old_component_source_ast", "component_source", "old_component_source"))
    if old_component is None:
        old_component = _lookup(owner, ("component_source", "component_ast", "old_component_source_ast"))
    pressure_identity = _lookup(contract, ("pressure_scale_identity", "pressure_identity", "P0_identity"))
    old_p0 = _lookup(contract, ("old_P0_assignment_ast", "old_p0_assignment_ast", "old_P0_assignment"))
    if old_p0 is None and _is_mapping(pressure_identity):
        old_p0 = _lookup(pressure_identity, ("old_P0_assignment_ast", "old_p0_assignment_ast", "axis_old_P0_provider", "old_P0_provider"))
    source_pins = _lookup(contract, ("sourcepins", "source_pins", "old_source_pins", "projection_pins"))
    if source_pins is None:
        source_pins = capsule.get("projection_pins")
    source_definition_present = (_is_mapping(definitions) or bool(_lookup(contract, ("source_function_definition", "exact_source_definition", "source_function_interpretation")))
        or bool(old_component)) and bool(old_constructor)
    p0_definition_present = bool(old_p0) or bool(pressure_identity)
    theorem_present = theorem is not None or bool(_lookup(contract, ("conditional_common_unique_fixed_point", "common_theorem_consumed")))
    if not source_definition_present or not p0_definition_present or not theorem_present:
        raise KernelInputError("exact_view_contract lacks source definitions, pressure identity, or shared theorem")
    if theorem is not None and _is_mapping(theorem):
        for key in ("passed", "accepted", "unique_actual_analytic_fixed_point", "true_H1_datum", "fourteen_stage_source"):
            if key in theorem and theorem[key] is False:
                raise KernelInputError("accepted common theorem flag is false: " + key)
    if definitions is not None and _is_mapping(definitions):
        for key in ("Paxis_eta_exact", "Ptotal_eta_exact", "oldPphysical", "Paxis_scaled_exact", "Ptotal_scaled_exact", "Pphysical_old"):
            if key in definitions and "epsilon" not in str(definitions[key]):
                raise KernelInputError("exact pressure definition lacks epsilon scale: " + key)
    return {"schema": "exact_view_registration_v1", "owner_uri": _json_value(_lookup(contract, ("old_owner_uri", "owner_uri", "owner"))) or _json_value(_lookup(owner, ("path", "uri", "owner_uri"))), "provider_hashes": _json_value(_lookup(contract, ("provider_hashes", "source_provider_hashes"))), "old_component_constructor_ast": _json_value(old_constructor), "old_component_source_ast": _json_value(old_component), "old_P0_assignment_ast": _json_value(old_p0), "old_raw_system_source": _json_value(_lookup(contract, ("old_raw_system_source", "raw_system_source", "axial_source"))), "source_pins": _json_value(source_pins), "formal_definitions": _json_value(definitions), "accepted_theorem_facts": _json_value(theorem), "old_source_uri_aliases": {"Phi": "old current_original_unique_analytic_core_component.Phi", "Psi": "old current_original_unique_analytic_core_component.Psi", "meanPsi": "integral_t_Psi", "U": "old U0+epsilon_core*Psi", "M": "radial average of U", "F0": "old F0(Z)", "Paxis": "epsilon_core*old P0(Z)", "Ptotal": "epsilon_core*oldPphysical(Z)", "Q": "radial_Q(delta,U,M)"}, "conditional_leading_identity_registered": True, "registration_basis": ["source_function_definition", "shared_common_unique_theorem"], "hash_only_numerical_equality": False, "candidate_polynomial_is_exact": False, "first_Picard_is_exact": False, "interval_overlap_is_not_identity": True, "regular_n1_view": False}


def _compact_native_packet(value: Any) -> Any:
    if _is_sequence(value):
        if len(value) > 32:
            return {"sequence_count": len(value), "retained_without_duplication": True}
        return [_compact_native_packet(child) if _is_mapping(child) or _is_sequence(child) else _json_value(child) for child in value]
    if not _is_mapping(value):
        return _json_value(value)
    omit = {"rows", "values", "samples", "payload", "result", "raw_terms", "full_table"}
    out = {}
    for key, child in value.items():
        if str(key) in omit:
            out[str(key) + "_retained_without_duplication"] = True
        else:
            out[str(key)] = _compact_native_packet(child) if _is_mapping(child) or _is_sequence(child) else _json_value(child)
    return out


def _rhs_record(context: Any, bundle: Mapping[str, Any]) -> dict[str, Any]:
    actual, candidate = bundle["actual"], bundle["candidate"]
    error = _abs_bound(context, actual - candidate) if candidate is not None else _abs_bound(context, actual)
    return {"actual": _interval_record(context, actual), "candidate": _record_or_none(context, candidate), "error": _interval_record(context, error), "actual_is_enclosure": True, "signed_value_resolved": candidate is not None, "source_bound_only": candidate is None}


def _source_hashes(capsule: Mapping[str, Any]) -> dict[str, Any]:
    source_ids = _as_mapping(capsule["source_ids"], "source_ids")
    raw = _lookup(capsule, ("raw_capsule_sha256", "raw_sha256", "capsule_raw_sha256")) or _lookup(source_ids, ("raw_capsule_sha256", "raw_sha256", "capsule_raw_sha256")) or CAPSULE_RAW_SHA256
    canonical = hashlib.sha256(_canonical_bytes(capsule)).hexdigest()
    return {"raw_capsule_sha256": raw, "canonical_capsule_sha256": canonical, "root_freeze_sha256": ROOT_FREEZE_SHA, "raw_and_canonical_are_distinct": raw is None or str(raw).lower() != canonical.lower(), "root_capsule_hash_pending": False}


def _serialize_result(value: Mapping[str, Any]) -> tuple[str, str, int]:
    text = json.dumps(_json_value(value), ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)
    encoded = text.encode("utf-8")
    return text, hashlib.sha256(encoded).hexdigest(), len(encoded)


def derive(data: Mapping[str, Any], emit: Callable[[str, Mapping[str, Any]], Any] | None = None) -> dict[str, Any]:
    """Build the local eta n=1 source packet once, after root freezing."""

    if emit is not None and not callable(emit):
        raise KernelInputError("emit must be callable or None")
    if not _is_mapping(data):
        raise KernelInputError("derive requires a mapping")
    capsule = data.get("capsule", data)
    capsule = _validate_capsule(capsule)
    events: list[dict[str, Any]] = []
    _emit(emit, events, "family_entered", {"capsule_schema": CAPSULE_SCHEMA, "precision_dps": INTERVAL_DPS, "N_materialized": False})
    context = _new_interval_context()
    params, polys = _source_parameters(context, capsule), _source_polynomials(context, capsule)
    norm, rows = _source_norm_data(context, capsule), _source_eta_rows(capsule)
    source_ids, grid, hashes = _json_value(capsule["source_ids"]), _json_value(capsule["grid"]), _source_hashes(capsule)

    _emit(emit, events, "exact_view_registration_entered", {"family": FAMILY_NAMES[0], "before_math": True})
    view = _proof_view(capsule)
    _emit(emit, events, "exact_view_registration_complete", {"family": FAMILY_NAMES[0], "conditional": True, "hash_only_identity": False})

    _emit(emit, events, "mixed_core_jets_entered", {"family": FAMILY_NAMES[1], "before_math": True, "eta_count": 8, "rho_labels": list(RHO_LABELS)})
    mixed = _source_mixed_jets(context, capsule, params, polys, norm, rows)
    fields = _build_fields(context, capsule, params, polys, mixed, rows)
    _emit(emit, events, "mixed_core_jets_complete", {"family": FAMILY_NAMES[1], "ordinary_Z_orders": list(ORDINARY_Z_ORDERS), "radial_orders": list(RADIAL_ORDERS), "m3_candidate": None, "uniform_eta_direct_row": True})

    _emit(emit, events, "gauge_third_jet_entered", {"family": FAMILY_NAMES[2], "before_math": True, "F0_anchor_formal": True})
    gauges, gauge_fields = {}, {}
    gauge_scalars = {}
    for eta in ETA_LABELS:
        gauge_scalars[eta] = _gauge_scalars(context, params, _source_row_atoms(context, rows[eta], "eta_rows." + eta), eta)
        gauges[eta], gauge_fields[eta] = {}, {}
        for rho_label in RHO_LABELS:
            gauges[eta][rho_label], gauge_fields[eta][rho_label] = _gauge(context, params, rows, mixed, rho_label, eta, gauge_scalars[eta])
    _emit(emit, events, "gauge_third_jet_complete", {"family": FAMILY_NAMES[2], "ell1_reused": 8, "ell2_recomputed_uniform_only": True, "ell3_rows": 8, "gauge_scalar_rows": 8, "gauge_scalars_reused_for_rho_rows": 16, "H_anchor_zero_exact": True})

    _emit(emit, events, "n1_differential_operators_entered", {"family": FAMILY_NAMES[3], "before_math": True, "ordinary_Z": True})
    qbundle = _q_bundle(context, params, rows, fields)
    _emit(emit, events, "n1_differential_operators_complete", {"family": FAMILY_NAMES[3], "D2_orders": [0, 1], "Q_radial_orders": [0, 1, 2], "C1_ordinary_derivative": True, "Q_C2_ordinary_derivative": True, "candidate_channels_independent": True})

    _emit(emit, events, "canonical_rhs_assembly_entered", {"family": FAMILY_NAMES[4], "before_math": True, "row_count": 48})
    rhs_rows: list[dict[str, Any]] = []
    for eta in ETA_LABELS:
        atoms = _source_row_atoms(context, rows[eta], "eta_rows." + eta)
        for rho_label in RHO_LABELS:
            rho = _rho_interval(context, rho_label)
            z, d, L = params["a"] + atoms["shift"], _one(context) - (params["a"] + atoms["shift"]) ** 2, _one(context) - params["delta"] * (params["a"] + atoms["shift"]) ** 2
            qslice, uslice = qbundle[eta][rho_label], fields["U"][eta][rho_label]
            q = {(k, m): qslice[(k, m)] for k in RADIAL_ORDERS for m in (0, 1, 2)}
            u = uslice
            flow_actual, flowz_actual = _sflow(context, q, u, params, z, d, L, rho, "actual")
            flow_candidate, flowz_candidate = _sflow(context, q, u, params, z, d, L, rho, "candidate")
            p_actual, pz_actual = -params["epsilon"] * flow_actual / 2 + 2 * q[(1, 0)]["actual"] + rho * q[(2, 0)]["actual"], -params["epsilon"] * flowz_actual / 2 + 2 * q[(1, 1)]["actual"] + rho * q[(2, 1)]["actual"]
            p_candidate = None if flow_candidate is None else -params["epsilon"] * flow_candidate / 2 + 2 * q[(1, 0)]["candidate"] + rho * q[(2, 0)]["candidate"]
            pz_candidate = None if flowz_candidate is None else -params["epsilon"] * flowz_candidate / 2 + 2 * q[(1, 1)]["candidate"] + rho * q[(2, 1)]["candidate"]
            F, U = gauge_fields[eta][rho_label], u
            zsrc = {"actual": -params["epsilon"] * _d2(context, U, -1 - params["delta"], params, z, d, L, rho, "actual"), "candidate": None}
            zsrc["candidate"] = None if _d2(context, U, -1 - params["delta"], params, z, d, L, rho, "candidate") is None else -params["epsilon"] * _d2(context, U, -1 - params["delta"], params, z, d, L, rho, "candidate")
            ts = {"actual": -params["epsilon"] * _d2(context, F, -2 - params["delta"], params, z, d, L, rho, "actual"), "candidate": None}
            ts["candidate"] = None if _d2(context, F, -2 - params["delta"], params, z, d, L, rho, "candidate") is None else -params["epsilon"] * _d2(context, F, -2 - params["delta"], params, z, d, L, rho, "candidate")
            zsz = {"actual": -params["epsilon"] * _d2(context, U, -1 - params["delta"], params, z, d, L, rho, "actual", True), "candidate": None}
            zsz["candidate"] = None if _d2(context, U, -1 - params["delta"], params, z, d, L, rho, "candidate", True) is None else -params["epsilon"] * _d2(context, U, -1 - params["delta"], params, z, d, L, rho, "candidate", True)
            tsz = {"actual": -params["epsilon"] * _d2(context, F, -2 - params["delta"], params, z, d, L, rho, "actual", True), "candidate": None}
            tsz["candidate"] = None if _d2(context, F, -2 - params["delta"], params, z, d, L, rho, "candidate", True) is None else -params["epsilon"] * _d2(context, F, -2 - params["delta"], params, z, d, L, rho, "candidate", True)
            g4a, g4za, g4c, g4zc = 2 * context.sqrt(rho) * p_actual, 2 * context.sqrt(rho) * pz_actual, (None if p_candidate is None else 2 * context.sqrt(rho) * p_candidate), (None if pz_candidate is None else 2 * context.sqrt(rho) * pz_candidate)
            g5a, g5za, g5c, g5zc = 2 * ts["actual"], 2 * tsz["actual"], (None if ts["candidate"] is None else 2 * ts["candidate"]), (None if tsz["candidate"] is None else 2 * tsz["candidate"])
            Lz = -2 * params["delta"] * z
            g6a = 2 * zsrc["actual"] - 4 * params["epsilon"] * z * rho * p_actual / L
            g6za = 2 * zsz["actual"] - 4 * params["epsilon"] * rho * ((p_actual + z * pz_actual) / L - z * p_actual * Lz / L ** 2)
            g6c = None if zsrc["candidate"] is None or p_candidate is None else 2 * zsrc["candidate"] - 4 * params["epsilon"] * z * rho * p_candidate / L
            g6zc = None if zsz["candidate"] is None or pz_candidate is None else 2 * zsz["candidate"] - 4 * params["epsilon"] * rho * ((p_candidate + z * pz_candidate) / L - z * p_candidate * Lz / L ** 2)
            # Keep axis/correction/viscous sectors separate.  No full-minus-axis subtraction is used.
            base_U, corr_U = {}, {}
            for k in RADIAL_ORDERS:
                for m in ORDINARY_Z_ORDERS:
                    if k == 0 and m == 0:
                        base_value = atoms["Ubase"]
                    elif k == 0 and m == 1:
                        base_value = context.mpf("4")
                    else:
                        base_value = _zero(context)
                    base_U[(k, m)] = _jet(context, base_value, base_value, _zero(context))
                    corr_U[(k, m)] = fields["DU"][eta][rho_label][(k, m)]
            du_sector_fields = {"DU_value": {"actual": corr_U[(0, 0)]["actual"], "candidate": corr_U[(0, 0)]["candidate"]}, "DU_Z_value": {"actual": corr_U[(0, 1)]["actual"], "candidate": corr_U[(0, 1)]["candidate"]}}
            actual_sector_pairs = _sector_s_pairs(context, params, atoms, rho, qslice, du_sector_fields, "actual")
            candidate_sector_pairs = None
            axis_candidate_pairs = None
            sector_required = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 1)]
            zero_jet = _jet(context, _zero(context), _zero(context), _zero(context))
            qaxis_candidate = {"q0": qslice["q0"], "delta_q": {key: zero_jet for key in sector_required}}
            for k in RADIAL_ORDERS:
                for m in (0, 1, 2):
                    qaxis_candidate[(k, m)] = qslice["q0"][m] if k == 0 else zero_jet
            zero_du_fields = {"DU_value": {"actual": _zero(context), "candidate": _zero(context)}, "DU_Z_value": {"actual": _zero(context), "candidate": _zero(context)}}
            axis_candidate_pairs = _sector_s_pairs(context, params, atoms, rho, qaxis_candidate, zero_du_fields, "candidate")
            if all(qslice["delta_q"][key]["candidate"] is not None for key in sector_required) and all(du_sector_fields[name]["candidate"] is not None for name in ("DU_value", "DU_Z_value")):
                candidate_sector_pairs = _sector_s_pairs(context, params, atoms, rho, qslice, du_sector_fields, "candidate")
            p_axis_actual = _pair_scale(actual_sector_pairs[0], -params["epsilon"] / 2)
            p_corr_actual = _pair_scale(actual_sector_pairs[1], -params["epsilon"] / 2)
            p_visc_actual = actual_sector_pairs[2]
            p_axis_candidate = _pair_scale(axis_candidate_pairs[0], -params["epsilon"] / 2)
            p_corr_candidate = None if candidate_sector_pairs is None else _pair_scale(candidate_sector_pairs[1], -params["epsilon"] / 2)
            p_visc_candidate = None if any(qslice["delta_q"][key]["candidate"] is None for key in ((1, 0), (1, 1), (2, 0), (2, 1))) else (2 * qslice["delta_q"][(1, 0)]["candidate"] + rho * qslice["delta_q"][(2, 0)]["candidate"], 2 * qslice["delta_q"][(1, 1)]["candidate"] + rho * qslice["delta_q"][(2, 1)]["candidate"])
            F_axis, F_corr = gauges[eta][rho_label]["F_axis_jets"], gauges[eta][rho_label]["F_correction_jets"]
            axis_t_actual = -params["epsilon"] * _d2(context, F_axis, -2 - params["delta"], params, z, d, L, rho, "actual")
            corr_t_actual = -params["epsilon"] * _d2(context, F_corr, -2 - params["delta"], params, z, d, L, rho, "actual")
            axis_t_candidate_raw = _d2(context, F_axis, -2 - params["delta"], params, z, d, L, rho, "candidate")
            corr_t_candidate_raw = _d2(context, F_corr, -2 - params["delta"], params, z, d, L, rho, "candidate")
            axis_t_candidate = None if axis_t_candidate_raw is None else -params["epsilon"] * axis_t_candidate_raw
            corr_t_candidate = None if corr_t_candidate_raw is None else -params["epsilon"] * corr_t_candidate_raw
            axis_tz_actual = -params["epsilon"] * _d2(context, F_axis, -2 - params["delta"], params, z, d, L, rho, "actual", True)
            corr_tz_actual = -params["epsilon"] * _d2(context, F_corr, -2 - params["delta"], params, z, d, L, rho, "actual", True)
            axis_tz_candidate_raw = _d2(context, F_axis, -2 - params["delta"], params, z, d, L, rho, "candidate", True)
            corr_tz_candidate_raw = _d2(context, F_corr, -2 - params["delta"], params, z, d, L, rho, "candidate", True)
            axis_tz_candidate = None if axis_tz_candidate_raw is None else -params["epsilon"] * axis_tz_candidate_raw
            corr_tz_candidate = None if corr_tz_candidate_raw is None else -params["epsilon"] * corr_tz_candidate_raw
            axis_z_actual = -params["epsilon"] * _d2(context, base_U, -1 - params["delta"], params, z, d, L, rho, "actual")
            corr_z_actual = -params["epsilon"] * _d2(context, corr_U, -1 - params["delta"], params, z, d, L, rho, "actual")
            axis_z_candidate_raw = _d2(context, base_U, -1 - params["delta"], params, z, d, L, rho, "candidate")
            corr_z_candidate_raw = _d2(context, corr_U, -1 - params["delta"], params, z, d, L, rho, "candidate")
            axis_z_candidate = None if axis_z_candidate_raw is None else -params["epsilon"] * axis_z_candidate_raw
            corr_z_candidate = None if corr_z_candidate_raw is None else -params["epsilon"] * corr_z_candidate_raw
            p_axis_g4_actual, p_corr_g4_actual, p_visc_g4_actual = _pair_scale(p_axis_actual, 2 * context.sqrt(rho)), _pair_scale(p_corr_actual, 2 * context.sqrt(rho)), _pair_scale(p_visc_actual, 2 * context.sqrt(rho))
            p_axis_g4_candidate = None if p_axis_candidate is None else _pair_scale(p_axis_candidate, 2 * context.sqrt(rho))
            p_corr_g4_candidate = None if p_corr_candidate is None else _pair_scale(p_corr_candidate, 2 * context.sqrt(rho))
            p_visc_g4_candidate = None if p_visc_candidate is None else _pair_scale(p_visc_candidate, 2 * context.sqrt(rho))
            g6_axis_actual = _pair_add((2 * axis_z_actual, 2 * (-params["epsilon"] * _d2(context, base_U, -1 - params["delta"], params, z, d, L, rho, "actual", True))), _coupling_pair(context, params, z, L, rho, p_axis_actual))
            g6_corr_actual = _pair_add((2 * corr_z_actual, 2 * (-params["epsilon"] * _d2(context, corr_U, -1 - params["delta"], params, z, d, L, rho, "actual", True))), _coupling_pair(context, params, z, L, rho, p_corr_actual))
            g6_visc_actual = _coupling_pair(context, params, z, L, rho, p_visc_actual)
            g6_axis_candidate = None if axis_z_candidate is None or p_axis_candidate is None or _d2(context, base_U, -1 - params["delta"], params, z, d, L, rho, "candidate", True) is None else _pair_add((2 * axis_z_candidate, 2 * (-params["epsilon"] * _d2(context, base_U, -1 - params["delta"], params, z, d, L, rho, "candidate", True))), _coupling_pair(context, params, z, L, rho, p_axis_candidate))
            g6_corr_candidate = None if corr_z_candidate is None or p_corr_candidate is None or _d2(context, corr_U, -1 - params["delta"], params, z, d, L, rho, "candidate", True) is None else _pair_add((2 * corr_z_candidate, 2 * (-params["epsilon"] * _d2(context, corr_U, -1 - params["delta"], params, z, d, L, rho, "candidate", True))), _coupling_pair(context, params, z, L, rho, p_corr_candidate))
            g6_visc_candidate = None if p_visc_candidate is None else _coupling_pair(context, params, z, L, rho, p_visc_candidate)
            sectors = {"pressure_flow_axis": _pair_bundle(context, p_axis_actual, p_axis_candidate), "pressure_flow_correction": _pair_bundle(context, p_corr_actual, p_corr_candidate), "pressure_viscous": _pair_bundle(context, p_visc_actual, p_visc_candidate), "g4_axis": _pair_bundle(context, p_axis_g4_actual, p_axis_g4_candidate), "g4_correction": _pair_bundle(context, p_corr_g4_actual, p_corr_g4_candidate), "g4_viscous": _pair_bundle(context, p_visc_g4_actual, p_visc_g4_candidate), "g5_angular_gauge_axis": _pair_bundle(context, (2 * axis_t_actual, 2 * axis_tz_actual), None if axis_t_candidate is None or axis_tz_candidate is None else (2 * axis_t_candidate, 2 * axis_tz_candidate)), "g5_angular_core_correction": _pair_bundle(context, (2 * corr_t_actual, 2 * corr_tz_actual), None if corr_t_candidate is None or corr_tz_candidate is None else (2 * corr_t_candidate, 2 * corr_tz_candidate)), "g6_axial_axis": _pair_bundle(context, g6_axis_actual, g6_axis_candidate), "g6_axial_correction": _pair_bundle(context, g6_corr_actual, g6_corr_candidate), "g6_pressure_coupling": _pair_bundle(context, g6_visc_actual, g6_visc_candidate)}
            bundles = {"g4": ((g4a, g4c), (g4za, g4zc), "one"), "g5": ((g5a, g5c), (g5za, g5zc), "F0_anchor"), "g6": ((g6a, g6c), (g6za, g6zc), "one")}
            for quantity in QUANTITIES:
                c0, c1, atom = bundles[quantity]
                rhs0 = _rhs_record(context, {"actual": c0[0], "candidate": c0[1]})
                rhs1 = _rhs_record(context, {"actual": c1[0], "candidate": c1[1]})
                rhs_rows.append({"eta": eta, "rho": rho_label, "quantity": quantity, "atom": atom, "F0_anchor_formal": atom == "F0_anchor", "C0": rhs0, "C1": rhs1, "actual_C0": _interval_record(context, c0[0]), "actual_C1": _interval_record(context, c1[0]), "candidate_C0": _record_or_none(context, c0[1]), "candidate_C1": _record_or_none(context, c1[1]), "error_C0": rhs0["error"], "error_C1": rhs1["error"], "sectors": sectors, "source_actual_independent_pressure_target": True, "reserved_own_swirl_registered_once": True, "n1_solution_value": None})
    _emit(emit, events, "canonical_rhs_assembly_complete", {"family": FAMILY_NAMES[4], "rhs_rows": len(rhs_rows), "C0_C1": True, "n1_solution_materialized": False})

    mixed_payload = {component: {eta: {rho: {"m" + str(m) + "_k" + str(k): _record_bundle(context, mixed[component][eta][rho][(k, m)]) for k in RADIAL_ORDERS for m in ORDINARY_Z_ORDERS} for rho in RHO_LABELS} for eta in ETA_LABELS} for component in ("Phi", "PhiMinusOne", "DU")}
    gauge_payload = {eta: {rho: {"z": _interval_record(context, gauges[eta][rho]["z"]), "d": _interval_record(context, gauges[eta][rho]["d"]), "L": _interval_record(context, gauges[eta][rho]["L"]), "H": _interval_record(context, gauges[eta][rho]["H"]), "H_Z": _interval_record(context, gauges[eta][rho]["H_Z"]), "H_ZZ": _interval_record(context, gauges[eta][rho]["H_ZZ"]), "ratio": _interval_record(context, gauges[eta][rho]["ratio"]), "ell1": _interval_record(context, gauges[eta][rho]["ell1"]), "ell2": _interval_record(context, gauges[eta][rho]["ell2"]), "ell3": _interval_record(context, gauges[eta][rho]["ell3"]), "F0_anchor_formal": True, "F_jets": {"m" + str(m) + "_k" + str(k): _record_bundle(context, gauges[eta][rho]["F_jets"][(k, m)]) for k in RADIAL_ORDERS for m in ORDINARY_Z_ORDERS}, "F_axis_jets": {"m" + str(m) + "_k" + str(k): _record_bundle(context, gauges[eta][rho]["F_axis_jets"][(k, m)]) for k in RADIAL_ORDERS for m in ORDINARY_Z_ORDERS}, "F_correction_jets": {"m" + str(m) + "_k" + str(k): _record_bundle(context, gauges[eta][rho]["F_correction_jets"][(k, m)]) for k in RADIAL_ORDERS for m in ORDINARY_Z_ORDERS}} for rho in RHO_LABELS} for eta in ETA_LABELS}

    _emit(emit, events, "proof_bundle", {"proof_count": 1, "conditional_exact_view": True, "m3_candidate": None, "regular_n1_view": False, "full_goal_complete": False})
    _emit(emit, events, "family_completed", {"proof_count": 1, "source_only": True, "n1_solution_materialized": False, "full_physical": False, "global_closure": False})
    definitions = {"exact_view": {"coordinate": "Z_eta=a+s*eta", "Phi_eta_exact": "Phi_old(rho,Z_eta)", "Psi_eta_exact": "Psi_old(rho,Z_eta)", "meanPsi_eta_exact": "integral_0^1 Psi_old(t*rho,Z_eta) dt", "U_eta_exact": "4*Z_eta+j+epsilon_core*Psi_eta_exact", "F0_eta_exact": "F0_old(Z_eta)=exp(-logC-Lambda*G(Z_eta))", "Paxis_eta_exact": "epsilon_core*P0_old(Z_eta)", "Ptotal_eta_exact": "epsilon_core*P0_old(Z_eta)+epsilon_core^2*F0_eta_exact^2*V(Phi_eta_exact^2)", "oldPphysical": "P0_old(Z_eta)+epsilon_core*F0_eta_exact^2*V(Phi_eta_exact^2)", "scaled_pressure_identity": "Ptotal_eta_exact=epsilon_core*oldPphysical", "candidate_polynomials_are_envelopes": True}, "mixed_jets": {"shift": "s*eta", "ordinary_Z": "fixed rho", "B_mk": "M*m!*(m+1)_k/[(m+1)^2*h^m*20^k*(1-4/20)^(m+k+1)]", "m0": "E0k+|shift|E1k+shift^2*E2k/2+T3k*|shift|^3/6", "m1": "E1k+|shift|E2k+T3k*shift^2/2", "m2": "E2k+T3k*|shift|", "m3": "candidate=None; source enclosure=[-T3k,+T3k]", "Mavg": "average polynomial coefficients by (n+1), then radial differentiate", "epsilon_rule": "DeltaU=epsilon_core*DeltaU1 exactly once"}, "gauge": {"H": "H1a*shift+H2a*shift^2/2-4*shift^3", "L": "L_a-2*delta*a*shift-delta*shift^2", "D": "H^2+sigma^2>0", "ell1": "-(A/D)/epsilon", "ell2": "-(A1/D-A*D1/D^2)/epsilon", "ell3": "-(A2/D-2*A1*D1/D^2-A*D2/D^2+2*A*D1^2/D^3)/epsilon", "F_jets": "F0_anchor*ratio times ordinary product-rule coefficients", "F0_anchor_materialized": False}, "operators": {"A": "w*Z*f+d*fZ-2*Z*rho*fr", "A_Z": "w*f+(w-2)*Z*fZ+d*fZZ-2*rho*fr-2*Z*rho*frZ", "A_rho": "(w-2)*Z*fr+d*frZ-2*Z*rho*frr", "D2": "[(w-1+delta)*Z*A+d*A_Z-2*Z*rho*A_rho]/L^2+2*delta*Z*d*A/L^3", "Nt": "-D2_{-2-delta}[F]", "Nz": "-D2_{-1-delta}[U]", "Sflow": "(Q+(1-delta)*Z*QZ/2+rho*Qr)/L+Q*(Q/2+rho*Qr)+U*(-2*Z*Q+d*QZ-2*Z*rho*Qr)/L", "Psrc": "-epsilon*Sflow/2+2*Qr+rho*Qrr", "g4": "2*sqrt(rho)*Psrc", "g5": "2*Tsrc", "g6": "2*Zsrc-4*epsilon*Z*rho*Psrc/L", "C1": "ordinary derivative at fixed rho via product/quotient rule"}}
    definitions["mixed_jets"]["epsilon_rule"] = "anchor/value_rows/mean-primitives DU already physical DeltaU=epsilon_core*Psi; epsilon_core is applied only to norm-derived DU bounds"
    definitions["units"] = {"DU_anchor_polynomials": "physical DeltaU", "DU_value_rows": "physical DeltaU", "JD_JDZ": "physical DeltaU primitives", "norm_Mpsi": "Psi norm before physical epsilon scaling", "DU_norm_bound_scaling": "multiply epsilon_core exactly once"}
    definitions["operators"]["Q_C2_ordinary_rules"] = {
        "Z_jet": "(Z,1,0)",
        "d_jet": "(d,-2*Z,-2)",
        "L_jet": "(L,-2*delta*Z,-2*delta)",
        "N0": "2*Z*DU0-(1-delta)*Z*DM0-d*DMZ1",
        "N1": "2*DU0+2*Z*DU1-(1-delta)*(DM0+Z*DMZ1)+2*Z*DMZ1-d*DMZ2",
        "N2": "4*DU1+2*Z*DU2-(1-delta)*(2*DM1+Z*DMZ2)+2*DMZ1+4*Z*DMZ2-d*DMZ3",
        "Q0": "N0/L",
        "Q1": "N1/L-N0*L1/L^2",
        "Q2": "N2/L-2*N1*L1/L^2-N0*L2/L^2+2*N0*L1^2/L^3",
        "candidate_rule": "candidate and actual channels are evaluated independently; unresolved candidate orders remain None",
    }
    definitions["gauge"]["scalar_evaluation"] = "ell1 and saved ell2 are reused; one new ell3 and one uniform-eta ell2 are computed per eta then shared by both rho rows"
    definitions["pressure_ports"] = {"P0_old": "old leading axis pressure provider", "Np": "Psrc coefficient inhomogeneity", "P1": "not solved or identified with Np", "m15": "not solved or identified with P0/Np/P1"}
    proof = {"kind": "component_native_O3_eta_n1_core_source_proof_v1", "proof_count": 1, "conditional_exact_view": view, "domain": {"rho": "[0,4]", "eta": "[-1,1]", "uniform_rho_direct": True, "extension_4_1_used": False}, "ordinary_Z_C1": True, "paired_exact_mpf_endpoints_authoritative": True, "m3_actual_jet_unresolved": True, "m3_candidate": None, "regular_n1_solution_owner": False, "old_212_error_recertification": False, "leading_moments_appended_to_212": False, "actual_independent_pressure_target_retained": True, "reserved_own_swirl_retained_once": True, "no_N_materialization": True, "no_old_provider_calls": True, "no_checker_or_selftest": True, "no_global_closure": True, "no_time_recursion": True, "no_full_physical_field": True, "no_NS_claim": True}
    parameter_payload = {("epsilon_core" if name == "epsilon" else name): _interval_record(context, value) for name, value in params.items()}
    payload = {"component": COMPONENT, "schema": SCHEMA, "source_capsule_schema": CAPSULE_SCHEMA, "source_only": True, "metadata_only": False, "numeric_evaluation": True, "conditional": True, "runtime_model_metadata": None, "source_ids": source_ids, "source_hashes": hashes, "canonical_gate": {"algorithm": "sha256(canonical_json(entire_capsule))", "actual_sha256": hashes["canonical_capsule_sha256"], "root_freeze_sha256": ROOT_FREEZE_SHA, "raw_sha256": CAPSULE_RAW_SHA256}, "grid": grid, "parameters": parameter_payload, "exact_view_registration": view, "mixed_jet_table": mixed_payload, "gauge_table": gauge_payload, "rhs_rows": rhs_rows, "definitions": definitions, "proof_packet": proof, "retained_native_packet": _compact_native_packet(capsule.get("retained_native_packet")), "source_data_reuse": {"value_rows": "direct k0 m0/m1 records consumed once; full source table omitted", "mean_primitives": "JD/JDZ candidate and pointwise enclosed polynomials consumed without reintegration; full source table omitted"}, "summary": {"eta_count": 8, "rho_count": 2, "rhs_row_count": len(rhs_rows), "rhs_expected_row_count": 48, "quantities": list(QUANTITIES), "ordinary_Z_orders": list(ORDINARY_Z_ORDERS), "radial_orders": list(RADIAL_ORDERS), "m3_candidate": None, "m3_source_bound_only": True, "uniform_eta_direct_row": True, "uniform_rho": "real [0,4] direct interval", "mean_primitives_consumed_without_reintegration": bool(capsule.get("mean_primitives") is not None), "F0_anchor": "formal source atom only", "regular_n1_view": False, "n1_solution_materialized": False, "global_closure": False, "full_physical": False, "NS_claim": False}, "actual_counts_contract": dict(SUCCESS_COUNTS), "actual_flags_contract": dict(SUCCESS_FLAGS)}
    result = {"schema": SCHEMA, "component": COMPONENT, "passed": True, "metadata_only": False, "source_only": True, "numeric_evaluation": True, "conditional": True, "runtime_model_metadata": None, "root_freeze_sha256": ROOT_FREEZE_SHA, "capsule_raw_sha256": CAPSULE_RAW_SHA256, "api": dict(API), "actual_counts": dict(SUCCESS_COUNTS), "actual_counts_contract": dict(SUCCESS_COUNTS), "new_flags": dict(SUCCESS_FLAGS), "actual_flags_contract": dict(SUCCESS_FLAGS), "expected_counts": dict(EXPECTED_COUNTS), "expected_flags": dict(EXPECTED_FLAGS), "callback_event_counts": dict(CALLBACK_EVENT_COUNTS), "callback_event_count": len(CALLBACK_EVENT_COUNTS), "manager_journal_event_count": MANAGER_JOURNAL_EVENT_COUNT, "journal_event_count": MANAGER_JOURNAL_EVENT_COUNT, "count_event_map": {key: list(value) for key, value in COUNT_EVENT_MAP.items()}, "proof_count": 1, "regular_n1_view": False, "n1_solution_materialized": False, "global_closure": False, "full_physical": False, "full_goal_complete": False, "payload": payload, "events": events, "unexecuted_code": {"status": "AST-ready-unexecuted", "sha256": None, "bytes": None, "hash_scope": "root-frozen source bytes; root computes after this leaf is accepted"}}
    _, result_hash, result_bytes = _serialize_result(result)
    result["return_code_hash"], result["return_code_bytes"] = result_hash, result_bytes
    result["code_hash"], result["code_bytes"] = result_hash, result_bytes
    result["serialization_hash_scope"] = "canonical result before hash metadata"
    return result


__all__ = ["API", "CALLBACK_EVENT_COUNTS", "CAPSULE_RAW_SHA256", "CAPSULE_SCHEMA", "COMPONENT", "COUNT_EVENT_MAP", "EXPECTED_COUNTS", "EXPECTED_FLAGS", "FAMILY_NAMES", "INTERVAL_DPS", "KernelInputError", "RHO_LABELS", "ROOT_FREEZE_SHA", "SCHEMA", "SUCCESS_COUNTS", "SUCCESS_FLAGS", "derive"]
