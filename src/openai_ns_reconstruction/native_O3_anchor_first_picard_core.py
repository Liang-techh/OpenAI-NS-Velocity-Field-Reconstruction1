"""Fresh exact-anchor first-Picard kernel for the native O3 core.

This file is intentionally a leaf.  It consumes one root-frozen capsule and
evaluates the first Picard image of the original coupled operator at the
exact implicit anchor.  It does not import an entered kernel, reconstruct a
Taylor recurrence, select an interval midpoint, or materialize the absolute
``F0`` atom.  The root owns capsule freezing, the single invocation, and
result persistence.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "native_O3_anchor_first_picard_core_kernel_result_v1"
CAPSULE_SCHEMA = "native_O3_anchor_first_picard_core_source_capsule_v1"
COMPONENT = "native_O3_anchor_first_picard_core"
INTERVAL_DPS = 280
SOURCE_PRECISION_BITS = 867

# The root replaces this exact literal after it freezes the complete capsule.
# A placeholder is deliberately rejected by _validate_capsule.
ROOT_FREEZE_SHA = "acde8e01042a2d5ee2049b0502e81cac8f382123d63c93343324c458af0456ff"

SOURCE_IDS = {
    "anchor_axis_result_sha256": (
        "c5a6f8273105a8f8710a44455de072ca93b5cd02b60f382249b8a0e5a3b016e8"
    ),
    "common_fixed_point_capsule_sha256": (
        "012511ca881a0207e04b7f494ce8ca275d7ea8a0a04d2914a824b5f846d64e0e"
    ),
}
FAMILY_IDENTITY = {
    "actual_five_defect_family_sha256": (
        "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894"
    ),
    "datum_enclosure_sha256": (
        "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010"
    ),
    "implicit_source_sha256": (
        "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae"
    ),
}

FAMILY_NAMES = (
    "anchored_anchor_seed",
    "anchored_picard_polynomial",
    "anchored_residual_certificate",
    "anchored_radial_value_ports",
    "anchored_moment_pressure_ports",
)
CALLBACK_EVENTS = (
    "family_entered",
    "anchored_anchor_seed_entered",
    "anchored_anchor_seed_complete",
    "anchored_picard_polynomial_entered",
    "anchored_picard_polynomial_complete",
    "anchored_residual_certificate_entered",
    "anchored_residual_certificate_complete",
    "anchored_radial_value_ports_entered",
    "anchored_radial_value_ports_complete",
    "anchored_moment_pressure_ports_entered",
    "anchored_moment_pressure_ports_complete",
    "proof_bundle",
    "family_completed",
)
CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "anchored_anchor_seed_entered": 1,
    "anchored_anchor_seed_complete": 1,
    "anchored_picard_polynomial_entered": 1,
    "anchored_picard_polynomial_complete": 1,
    "anchored_residual_certificate_entered": 1,
    "anchored_residual_certificate_complete": 1,
    "anchored_radial_value_ports_entered": 1,
    "anchored_radial_value_ports_complete": 1,
    "anchored_moment_pressure_ports_entered": 1,
    "anchored_moment_pressure_ports_complete": 1,
    "proof_bundle": 1,
    "family_completed": 1,
}
MANAGER_JOURNAL_EVENT_COUNT = 15

# This is one source entry and one numerical invocation containing five fresh
# mathematical families.  It is not five provider or ancestor invocations.
SUCCESS_COUNTS = {
    "N_materialization_calls": 0,
    "ancestor_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "provider_calls": 0,
    "checker_calls": 0,
    "selftest_calls": 0,
    "source_family_calls": 1,
    "new_math_family_calls": 5,
    "anchored_anchor_seed_family_calls": 1,
    "anchored_picard_polynomial_family_calls": 1,
    "anchored_residual_certificate_family_calls": 1,
    "anchored_radial_value_ports_family_calls": 1,
    "anchored_moment_pressure_ports_family_calls": 1,
    "derive_calls": 1,
    "numeric_evaluation_calls": 1,
    "new_mathematical_proofs": 1,
    "source_value_reducer_calls": 0,
    "F0_absolute_exp_calls": 0,
    "exp_calls": 0,
}

SUCCESS_FLAGS = {
    "source_only": True,
    "metadata_only": False,
    "numeric_evaluation": True,
    "conditional": True,
    "local_actual_fixed_point_anchor": True,
    "same_source_complete_pair": True,
    "exact_anchor_atom": True,
    "anchor_midpoint_selected": False,
    "absolute_F0_exp_materialized": False,
    "N_materialized": False,
    "old_math_reentered": False,
    "old_provider_called": False,
    "provider_called": False,
    "checker_called": False,
    "selftest_called": False,
    "source_value_point_selected": False,
    "source_parameter_point_selected": False,
    "full_physical": False,
    "physical_layer_completed": False,
    "global_five_moment": False,
    "time_n_recursion": False,
    "full_goal": False,
    "full_goal_complete": False,
}

COUNT_EVENT_MAP = {
    "source_family_calls": ["family_entered"],
    "new_math_family_calls": [
        "anchored_anchor_seed_entered",
        "anchored_picard_polynomial_entered",
        "anchored_residual_certificate_entered",
        "anchored_radial_value_ports_entered",
        "anchored_moment_pressure_ports_entered",
    ],
    "anchored_anchor_seed_family_calls": ["anchored_anchor_seed_entered"],
    "anchored_picard_polynomial_family_calls": [
        "anchored_picard_polynomial_entered"
    ],
    "anchored_residual_certificate_family_calls": [
        "anchored_residual_certificate_entered"
    ],
    "anchored_radial_value_ports_family_calls": [
        "anchored_radial_value_ports_entered"
    ],
    "anchored_moment_pressure_ports_family_calls": [
        "anchored_moment_pressure_ports_entered"
    ],
    "derive_calls": ["family_entered"],
    "numeric_evaluation_calls": ["family_entered"],
    "new_mathematical_proofs": ["proof_bundle"],
}

# Static aliases are part of the runner contract.  Keep them literal and
# separate from the mutable event journal returned by derive.
EXPECTED_COUNTS = {
    "N_materialization_calls": 0,
    "ancestor_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "provider_calls": 0,
    "checker_calls": 0,
    "selftest_calls": 0,
    "source_family_calls": 1,
    "new_math_family_calls": 5,
    "anchored_anchor_seed_family_calls": 1,
    "anchored_picard_polynomial_family_calls": 1,
    "anchored_residual_certificate_family_calls": 1,
    "anchored_radial_value_ports_family_calls": 1,
    "anchored_moment_pressure_ports_family_calls": 1,
    "derive_calls": 1,
    "numeric_evaluation_calls": 1,
    "new_mathematical_proofs": 1,
    "source_value_reducer_calls": 0,
    "F0_absolute_exp_calls": 0,
    "exp_calls": 0,
}
EXPECTED_FLAGS = {
    "source_only": True,
    "metadata_only": False,
    "numeric_evaluation": True,
    "conditional": True,
    "local_actual_fixed_point_anchor": True,
    "same_source_complete_pair": True,
    "exact_anchor_atom": True,
    "anchor_midpoint_selected": False,
    "absolute_F0_exp_materialized": False,
    "N_materialized": False,
    "old_math_reentered": False,
    "old_provider_called": False,
    "provider_called": False,
    "checker_called": False,
    "selftest_called": False,
    "source_value_point_selected": False,
    "source_parameter_point_selected": False,
    "full_physical": False,
    "physical_layer_completed": False,
    "global_five_moment": False,
    "time_n_recursion": False,
    "full_goal": False,
    "full_goal_complete": False,
}

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
API_METADATA = {
    "function": "derive",
    "signature": "derive(data, emit=None)",
    "input": "data['capsule']",
    "metadata_only": False,
    "source_only": True,
    "numeric_evaluation": True,
    "conditional": True,
    "runtime_model_metadata": None,
}


class KernelInputError(ValueError):
    """Raised when the root capsule is missing or unauthenticated."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _as_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not _is_mapping(value):
        raise KernelInputError(label + " must be a mapping")
    return value


def _json_value(value: Any) -> Any:
    """Copy JSON-shaped metadata without interpreting source objects."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if _is_mapping(value):
        return {str(key): _json_value(child) for key, child in value.items()}
    if _is_sequence(value):
        return [_json_value(child) for child in value]
    return str(value)


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
        emit(name, row)


def _canonical_bytes(capsule: Mapping[str, Any]) -> bytes:
    body = _json_value(capsule)
    return json.dumps(
        body,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _valid_sha(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(char in "0123456789abcdefABCDEF" for char in value)
    )


def _lookup(mapping: Mapping[str, Any], names: Sequence[str]) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    return None


def _raw_mpf_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise KernelInputError(label + " must be a four-word MPF tuple")
    try:
        words = tuple(int(item) for item in value)
    except (TypeError, ValueError) as exc:
        raise KernelInputError(label + " contains a non-integer MPF word") from exc
    return words  # type: ignore[return-value]


def _interval_value(context: Any, value: Any, label: str) -> Any:
    """Read a directed exact record and construct it with a tuple pair."""

    if not _is_mapping(value):
        raise KernelInputError(label + " must be an exact interval record")
    for nested in ("interval", "enclosure", "record", "value"):
        if nested in value and not (
            "lower_exact_mpf_tuple" in value
            or "upper_exact_mpf_tuple" in value
        ):
            return _interval_value(context, value[nested], label + "." + nested)
    if "lower_exact_mpf_tuple" not in value or "upper_exact_mpf_tuple" not in value:
        raise KernelInputError(
            label + " requires lower_exact_mpf_tuple and upper_exact_mpf_tuple"
        )
    lower = _raw_mpf_tuple(
        value["lower_exact_mpf_tuple"], label + ".lower_exact_mpf_tuple"
    )
    upper = _raw_mpf_tuple(
        value["upper_exact_mpf_tuple"], label + ".upper_exact_mpf_tuple"
    )
    # MPIntervalContext.make_mpf requires ((s,m,e,b),(s,m,e,b)); a raw tuple
    # by itself is intentionally never passed to it.
    return context.make_mpf((lower, upper))


def _raw_endpoint(context: Any, value: Any, upper: bool = False) -> Any:
    interval = context.convert(value)
    raw = interval._mpi_[1 if upper else 0]
    return context.make_mpf((raw, raw))


def _lower(context: Any, value: Any) -> Any:
    return _raw_endpoint(context, value, False)


def _upper(context: Any, value: Any) -> Any:
    return _raw_endpoint(context, value, True)


def _make_interval(context: Any, lower: Any, upper: Any) -> Any:
    """Construct an interval from endpoint values using the required pair."""

    lo = context.convert(lower)
    hi = context.convert(upper)
    return context.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _display_endpoint(context: Any, value: Any) -> str:
    text = str(context.nstr(context.convert(value), n=INTERVAL_DPS)).strip()
    if text.startswith("[") and text.endswith("]") and "," in text:
        text = text[1:-1].split(",", 1)[0].strip()
    return text


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


def _strict_positive(context: Any, value: Any) -> bool:
    return bool(_lower(context, value) > context.mpf("0"))


def _exact_zero(context: Any, value: Any) -> bool:
    zero = context.mpf("0")
    return bool(_lower(context, value) == zero and _upper(context, value) == zero)


def _abs_upper(context: Any, value: Any) -> Any:
    lo = _lower(context, value)
    hi = _upper(context, value)
    zero = context.mpf("0")
    if bool(lo >= zero):
        return hi
    if bool(hi <= zero):
        return -lo
    return -lo if bool(-lo >= hi) else hi


def _widen(context: Any, value: Any, error: Any) -> Any:
    err = _upper(context, error)
    return _make_interval(
        context,
        _lower(context, value) - err,
        _upper(context, value) + err,
    )


def _parameter_record(
    context: Any,
    parameters: Mapping[str, Any],
    name: str,
) -> tuple[Any, Any]:
    raw = parameters.get(name)
    if raw is None:
        raise KernelInputError("missing source parameter " + name)
    return raw, _interval_value(context, raw, "parameters." + name)


def _capsule_record(
    context: Any,
    mapping: Mapping[str, Any],
    name: str,
    label: str,
) -> tuple[Any, Any]:
    raw = mapping.get(name)
    if raw is None:
        raise KernelInputError("missing source record " + label)
    return raw, _interval_value(context, raw, label)


def _validate_capsule(
    capsule: Any,
) -> tuple[
    Mapping[str, Any],
    Mapping[str, Any],
    Mapping[str, Any],
    Mapping[str, Any],
    Mapping[str, Any],
]:
    capsule = _as_mapping(capsule, "data.capsule")
    if capsule.get("schema") != CAPSULE_SCHEMA:
        raise KernelInputError("unexpected first-Picard source capsule schema")
    if capsule.get("component") != COMPONENT:
        raise KernelInputError("source capsule component does not match")
    if not _valid_sha(ROOT_FREEZE_SHA) or ROOT_FREEZE_SHA == "0" * 64:
        raise KernelInputError("ROOT_FREEZE_SHA is not frozen")
    actual = hashlib.sha256(_canonical_bytes(capsule)).hexdigest()
    if actual.lower() != ROOT_FREEZE_SHA.lower():
        raise KernelInputError("source capsule canonical digest mismatch")

    source_ids = _as_mapping(capsule.get("source_ids"), "capsule.source_ids")
    for name, expected in SOURCE_IDS.items():
        if source_ids.get(name) != expected:
            raise KernelInputError("source id mismatch: " + name)

    parameters = _as_mapping(capsule.get("parameters"), "capsule.parameters")
    anchor = _as_mapping(capsule.get("anchor"), "capsule.anchor")
    pressure = _as_mapping(
        capsule.get("pressure_scaled"), "capsule.pressure_scaled"
    )
    banach_container = _as_mapping(capsule.get("banach"), "capsule.banach")
    if all(
        name in banach_container
        for name in ("correction_bound", "pair_lipschitz_bound", "self_map_bound")
    ):
        banach = banach_container
    else:
        banach = _lookup(
            banach_container,
            ("complete_source_proof_DATA", "complete_source_proof_data"),
        )
        if banach is None:
            banach = banach_container
    banach = _as_mapping(banach, "capsule.banach.complete_source_proof_DATA")
    grid = _as_mapping(capsule.get("grid"), "capsule.grid")
    if list(grid.get("rho", ())) != ["0", "1/4", "1", "2", "4.1"]:
        raise KernelInputError("capsule rho grid does not match the frozen grid")
    if list(grid.get("uniform_rho", ())) != ["0", "4.1"]:
        raise KernelInputError("capsule uniform rho grid does not match the frozen grid")
    for name in ("P0", "P0_Z", "P0_ZZ"):
        if name not in pressure:
            raise KernelInputError("pressure_scaled is missing " + name)
    for name in ("correction_bound", "pair_lipschitz_bound", "self_map_bound"):
        if name not in banach:
            raise KernelInputError("banach is missing " + name)
    if banach.get("same_source_complete_pair") is not True:
        raise KernelInputError("same-source complete-pair guard is required")
    if list(banach.get("uniform_real_centers", ())) not in (
        [-1, 1],
        ["-1", "1"],
    ):
        raise KernelInputError("uniform real centers must be [-1,1]")
    if str(banach.get("rho_reference")) != "20":
        raise KernelInputError("rho_reference must be the literal 20")
    if not str(banach.get("norm_definition", "")).strip():
        raise KernelInputError("banach norm definition is required")

    family_identity = _lookup(
        capsule,
        ("family_identity", "family_identity_hashes", "family_identity3hashes"),
    )
    if family_identity is not None:
        family_identity = _as_mapping(family_identity, "capsule.family_identity")
        for name, expected in FAMILY_IDENTITY.items():
            value = family_identity.get(name)
            if not _valid_sha(value) or value.lower() != expected.lower():
                raise KernelInputError("family identity hash is invalid: " + name)
    else:
        raise KernelInputError("family identity hashes are required")
    guards = _lookup(capsule, ("preservation_guards", "guards", "preservation"))
    if guards is None:
        raise KernelInputError("preservation guards are required")
    guards = _as_mapping(guards, "capsule.preservation_guards")
    if guards is not None:
        for name in (
            "global_five_moment",
            "global_five_moment_closure",
            "time_n_recursion",
            "physical_layer_completed",
            "full_goal",
            "full_goal_complete",
        ):
            if name in guards and guards[name] is not False:
                raise KernelInputError("preservation guard is not false: " + name)
        for name in (
            "N_materialization_calls",
            "absolute_F0_exponential_materialization_calls",
            "old_math_calls",
            "old_provider_calls",
        ):
            if name in guards and guards[name] != 0:
                raise KernelInputError("preservation count is not zero: " + name)
    definition = str(
        _lookup(anchor, ("exact_definition", "definition", "root_definition"))
        or ""
    )
    normalized_definition = definition.replace(" ", "")
    if "H(a)=0" not in normalized_definition or "[-j,0]" not in normalized_definition:
        raise KernelInputError("exact unique H(a)=0 anchor definition is missing")
    if anchor.get("midpoint_selected", False) is True or anchor.get(
        "root_interval_midpoint_selected", False
    ) is True:
        raise KernelInputError("anchor midpoint selection is forbidden")
    return capsule, parameters, anchor, pressure, banach


def _new_interval_context() -> Any:
    # This is the only runtime numerical import.  The root invokes derive once.
    from mpmath.ctx_iv import MPIntervalContext

    context = MPIntervalContext()
    context.dps = INTERVAL_DPS
    return context


def _poly_eval(context: Any, coefficients: Sequence[Any], rho: Any) -> Any:
    value = context.mpf("0")
    for coefficient in reversed(tuple(coefficients)):
        value = coefficient + rho * value
    return value


def _exact_text(context: Any, text: str) -> Any:
    if "/" in text:
        numerator, denominator = text.split("/", 1)
        return context.mpf(int(numerator)) / context.mpf(int(denominator))
    return context.mpf(text)


def _poly_mul(context: Any, left: Sequence[Any], right: Sequence[Any]) -> list[Any]:
    result = [context.mpf("0") for _ in range(len(left) + len(right) - 1)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] = result[i + j] + left_value * right_value
    return result


def _poly_integral(context: Any, coefficients: Sequence[Any]) -> list[Any]:
    result = [context.mpf("0")]
    for index, coefficient in enumerate(coefficients):
        result.append(coefficient / context.mpf(index + 1))
    return result


def _raw_parameter_packet(
    context: Any,
    parameters: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    raw_packet: dict[str, Any] = {}
    packet: dict[str, Any] = {}
    for name in ("delta", "epsilon_core", "j", "logC"):
        raw, value = _parameter_record(context, parameters, name)
        raw_packet[name] = _json_value(raw)
        packet[name] = value
    if not _strict_positive(context, packet["epsilon_core"]):
        raise KernelInputError("epsilon_core must be strictly positive")
    if not _strict_positive(context, packet["j"]):
        raise KernelInputError("j must be strictly positive")
    if not _strict_positive(context, packet["delta"]):
        raise KernelInputError("delta must be strictly positive")
    if not bool(_upper(context, packet["delta"]) < context.mpf("1")):
        raise KernelInputError("delta must be below one")
    if not _strict_positive(context, packet["logC"]):
        raise KernelInputError("logC must be strictly positive")
    return raw_packet, packet


def _source_scale_bound(
    context: Any,
    epsilon: Any,
    log_c: Any,
    epsilon_raw: Any,
) -> tuple[Any, dict[str, Any]]:
    """Bound S without evaluating either exp or log."""

    lower_record = epsilon_raw
    while _is_mapping(lower_record) and "lower_exact_mpf_tuple" not in lower_record:
        nested = _lookup(lower_record, ("interval", "enclosure", "record", "value"))
        if nested is None:
            break
        lower_record = nested
    if not _is_mapping(lower_record):
        raise KernelInputError("epsilon_core exact lower tuple is missing")
    lower_tuple = _raw_mpf_tuple(
        lower_record.get("lower_exact_mpf_tuple"),
        "parameters.epsilon_core.lower_exact_mpf_tuple",
    )
    exponent, bitcount = lower_tuple[2], lower_tuple[3]
    exponent_floor = exponent + bitcount - 1
    if exponent_floor >= 0:
        raise KernelInputError("epsilon_core is not in the small-source branch")
    E = -exponent_floor
    if not bool(_upper(context, epsilon) < context.mpf("1")):
        raise KernelInputError("epsilon_core must be below one for the source bound")
    if not bool(2 * _lower(context, log_c) >= 8 * context.mpf(E)):
        raise KernelInputError("logC does not certify the source-scale bound")
    epsilon_upper = _upper(context, epsilon)
    source_upper = epsilon_upper**8
    source_interval = _make_interval(context, context.mpf("0"), source_upper)
    proof = {
        "epsilon_lower_power_two_exponent": exponent_floor,
        "E": E,
        "ln2_lt_one_used": True,
        "inequality": "2*logC.lower >= 8*E",
        "S_upper_formula": "epsilon_core.upper**8",
        "S_bound_not_amplitude_choice": True,
        "strictly_positive_formal_atom": True,
        "exp_materialized": False,
        "log_materialized": False,
    }
    return source_interval, proof


def _anchor_packet(
    context: Any,
    anchor: Mapping[str, Any],
) -> tuple[Any, dict[str, Any]]:
    raw = anchor.get("record")
    if raw is None:
        raw = anchor.get("interval")
    if raw is None:
        raise KernelInputError("anchor.record is required")
    a = _interval_value(context, raw, "anchor.record")
    return a, {
        "record": _interval_record(context, a),
        "exact_definition": str(anchor.get("exact_definition", "")),
        "root_tag": str(anchor.get("root_tag", "exact implicit root a of H(a)=0")),
        "midpoint_selected": False,
        "Z_atom": "a; exact interval atom, never its midpoint",
    }


def _first_picard_data(
    context: Any,
    packet: Mapping[str, Any],
    source_interval: Any,
    anchor_value: Any,
    pressure: Mapping[str, Any],
    rho_grid: Sequence[str],
) -> dict[str, Any]:
    delta = packet["delta"]
    epsilon = packet["epsilon_core"]
    j = packet["j"]
    a = anchor_value
    d = 1 - a**2
    L = 1 - delta * a**2
    U0 = -((1 - delta) * a) / (2 * d)
    U0_Z = context.mpf("4")
    H1 = (context.mpf("9") - delta) / 2 - 12 * a**2 - 2 * j * a
    A = (1 + delta) / 2 * (1 - 2 * a * U0)
    A_Z = -(1 + delta) * (U0 + 4 * a)
    _, p = _capsule_record(context, pressure, "P0", "pressure_scaled.P0")
    _, p_Z = _capsule_record(context, pressure, "P0_Z", "pressure_scaled.P0_Z")
    _, p_ZZ = _capsule_record(context, pressure, "P0_ZZ", "pressure_scaled.P0_ZZ")
    H0_term = context.mpf("0")
    N = epsilon * A * U0 + 2 * (1 + delta) * a * p - d * p_Z
    N_Z = (
        epsilon * (A_Z * U0 + 4 * A - 4 * H1)
        + 2 * (1 + delta) * p
        + 2 * (2 + delta) * a * p_Z
        - d * p_ZZ
    )
    C = N / (2 * L)
    C_Z = N_Z / (2 * L) + (2 * delta * a / L) * C
    c = C / epsilon
    c_Z = C_Z / epsilon
    beta = (3 - delta / 2 + j * a) / L
    K = (3 + delta) / 2 - (3 + delta) * a * U0

    S = source_interval
    B = epsilon**2 * S
    phi_coefficients = (
        context.mpf("1"),
        -epsilon * beta / 4,
        -epsilon * ((1 + delta) * a * C + d * C_Z) / (24 * L),
    )
    delta_u_coefficients = (
        C,
        epsilon * K * C / (8 * L) - B * (2 + delta) * a / (4 * L),
        epsilon * (d * C * C_Z - (3 + delta) * a * C**2) / (36 * L),
    )
    psi_coefficients = tuple(value / epsilon for value in delta_u_coefficients)
    phi_square_coefficients = _poly_mul(context, phi_coefficients, phi_coefficients)
    pressure_integral_coefficients = _poly_integral(
        context, phi_square_coefficients
    )
    return {
        "anchor": a,
        "d": d,
        "L": L,
        "U0": U0,
        "U0_Z": U0_Z,
        "H1": H1,
        "A": A,
        "A_Z": A_Z,
        "p": p,
        "p_Z": p_Z,
        "p_ZZ": p_ZZ,
        "H0_term": H0_term,
        "N": N,
        "N_Z": N_Z,
        "C": C,
        "C_Z": C_Z,
        "c": c,
        "c_Z": c_Z,
        "beta": beta,
        "K": K,
        "S": S,
        "B": B,
        "phi_coefficients": phi_coefficients,
        "delta_u_coefficients": delta_u_coefficients,
        "psi_coefficients": psi_coefficients,
        "phi_square_coefficients": phi_square_coefficients,
        "pressure_integral_coefficients": pressure_integral_coefficients,
        "rho_grid": tuple(rho_grid),
        "definitions": {
            "Z": "a exact interval atom, not midpoint",
            "d": "1-a^2",
            "L": "1-delta*a^2",
            "U0": "-((1-delta)*a)/(2*d), from H(a)=0",
            "U0_affine_identity": "U0=4*a+j, exact H(a)=0 correlation",
            "U0_Z": "4",
            "H1": "(9-delta)/2-12*a^2-2*j*a",
            "A": "(1+delta)/2*(1-2*a*U0)",
            "A_Z": "-(1+delta)*(U0+4*a)",
            "p": "pressure_scaled.P0 = epsilon*physical_P0 at anchor",
            "N": "epsilon*A*U0+2*(1+delta)*a*p-d*p_Z; H0 term exact zero",
            "N_Z": (
                "epsilon*(A_Z*U0+4*A-4*H1)+2*(1+delta)*p+"
                "2*(2+delta)*a*p_Z-d*p_ZZ"
            ),
            "C": "N/(2*L)=epsilon*c",
            "C_Z": "N_Z/(2*L)+(2*delta*a/L)*C",
            "beta": "(3-delta/2+j*a)/L",
            "K": "(3+delta)/2-(3+delta)*a*U0",
            "operator": (
                "x0=(R_chi*1,c*rho); T(X)=x0+epsilon/2*"
                "(R_chi*J2*Etheta,J1*Ez); R_chi=(I+J2*chi/2)^-1"
            ),
        },
    }


def _radial_error(context: Any, error_norm: Any, order: int) -> Any:
    factorial = (1, 1, 2)[order]
    radius = context.mpf("4.1")
    ratio = 1 - radius / 20
    return error_norm * factorial / (20**order * ratio ** (order + 1))


def _actual_row(
    context: Any,
    first: Mapping[str, Any],
    rho: Any,
    errors: Mapping[int, Any],
    range_row: bool = False,
) -> dict[str, Any]:
    phi_coefficients = first["phi_coefficients"]
    du_coefficients = first["delta_u_coefficients"]
    psi_coefficients = first["psi_coefficients"]
    phi = _poly_eval(context, phi_coefficients, rho)
    # DeltaU1 and Psi1 coefficient packets begin at rho^1.  The leading
    # coefficient is C (respectively C/epsilon), not a degree-zero value.
    delta_u = rho * _poly_eval(context, du_coefficients, rho)
    psi = rho * _poly_eval(context, psi_coefficients, rho)
    delta_phi = rho * (
        phi_coefficients[1] + phi_coefficients[2] * rho
    )
    if _exact_zero(context, rho):
        phi_actual = context.mpf("1")
        delta_phi_actual = context.mpf("0")
        delta_actual = context.mpf("0")
        psi_actual = context.mpf("0")
        axis_exact = True
    else:
        delta_phi_actual = _widen(context, delta_phi, errors[0])
        # Keep the small nonlinear correction visible instead of subtracting
        # two rounded full values near one.  This is the same E0 enclosure.
        phi_actual = _make_interval(
            context,
            1 + _lower(context, delta_phi_actual),
            1 + _upper(context, delta_phi_actual),
        )
        delta_actual = _widen(
            context, delta_u, first["epsilon_value"] * errors[0]
        )
        psi_actual = _widen(context, psi, errors[0])
        axis_exact = False

    phi_rho = phi_coefficients[1] + 2 * phi_coefficients[2] * rho
    phi_rhorho = 2 * phi_coefficients[2]
    return {
        "rho": _interval_record(context, rho),
        "range_row": range_row,
        "Phi1_candidate": _interval_record(context, phi),
        "DeltaPhi1_candidate": _interval_record(context, delta_phi),
        "DeltaPhi_actual": _interval_record(context, delta_phi_actual),
        "Phi_actual": _interval_record(context, phi_actual),
        "Phi_actual_definition": "1+DeltaPhi_actual",
        "DeltaU1_candidate": _interval_record(context, delta_u),
        "DeltaU_actual": _interval_record(context, delta_actual),
        "U1_actual": _interval_record(
            context, first["U0"] + delta_actual
        ),
        "Psi1_candidate": _interval_record(context, psi),
        "Psi_actual": _interval_record(context, psi_actual),
        "Phi_radial_derivatives": {
            "k1_candidate": _interval_record(context, phi_rho),
            "k1_actual": _interval_record(
                context, _widen(context, phi_rho, errors[1])
            ),
            "k2_candidate": _interval_record(context, phi_rhorho),
            "k2_actual": _interval_record(
                context, _widen(context, phi_rhorho, errors[2])
            ),
        },
        "ordinary_Z_derivative_claimed": False,
        "structural_axis_error_zero": axis_exact,
        "point_value": False,
    }


def _moment_pressure_row(
    context: Any,
    first: Mapping[str, Any],
    rho: Any,
    errors: Mapping[int, Any],
    uniform_abs_phi: Any,
) -> dict[str, Any]:
    """Compute the fifth-family moment and pressure ports for one rho row."""

    du_coefficients = first["delta_u_coefficients"]
    moment_over_epsilon = (
        rho * first["U0"]
        + du_coefficients[0] * rho**2 / 2
        + du_coefficients[1] * rho**3 / 3
        + du_coefficients[2] * rho**4 / 4
    )
    moment_physical = first["epsilon_value"] * moment_over_epsilon
    pressure_integral = _poly_eval(
        context, first["pressure_integral_coefficients"], rho
    )
    rho_upper = _upper(context, rho)
    I_error = rho_upper * _upper(context, errors[0]) * (
        2 * uniform_abs_phi + _upper(context, errors[0])
    )
    I_actual = _widen(context, pressure_integral, I_error)
    if _exact_zero(context, rho):
        moment_over_epsilon_actual = context.mpf("0")
        moment_physical_actual = context.mpf("0")
        I_actual = context.mpf("0")
        scaled_pressure = first["p"]
        physical_pressure = first["p"] / first["epsilon_value"]
    else:
        moment_over_epsilon_error = rho_upper * _upper(
            context, first["epsilon_value"] * errors[0]
        )
        moment_over_epsilon_actual = _widen(
            context, moment_over_epsilon, moment_over_epsilon_error
        )
        moment_physical_actual = _widen(
            context,
            moment_physical,
            first["epsilon_value"] * moment_over_epsilon_error,
        )
        scaled_pressure = first["p"] + first["B"] * I_actual
        physical_pressure = (
            first["p"] / first["epsilon_value"]
            + first["epsilon_value"] * first["S"] * I_actual
        )
    return {
        "M_z_over_epsilon_candidate": _interval_record(
            context, moment_over_epsilon
        ),
        "M_z_over_epsilon_actual": _interval_record(
            context, moment_over_epsilon_actual
        ),
        "M_z_physical_candidate": _interval_record(context, moment_physical),
        "M_z_physical_actual": _interval_record(context, moment_physical_actual),
        "I_candidate": _interval_record(context, pressure_integral),
        "I_actual": _interval_record(context, I_actual),
        "pressure_scaled_actual": _interval_record(context, scaled_pressure),
        "pressure_physical_formal_actual": _interval_record(
            context, physical_pressure
        ),
    }


def _serialize_result(value: Mapping[str, Any]) -> tuple[str, str, int]:
    text = json.dumps(
        _json_value(value),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    encoded = text.encode("utf-8")
    return text, hashlib.sha256(encoded).hexdigest(), len(encoded)


def derive(
    data: Mapping[str, Any],
    emit: Callable[[str, Mapping[str, Any]], Any] | None = None,
) -> dict[str, Any]:
    """Evaluate one bounded fresh first-Picard source leaf."""

    if not _is_mapping(data) or "capsule" not in data:
        raise KernelInputError("derive requires data['capsule']")
    if emit is not None and not callable(emit):
        raise KernelInputError("emit must be callable or None")

    capsule, parameters, anchor_raw, pressure_raw, banach_raw = _validate_capsule(
        data["capsule"]
    )
    events: list[dict[str, Any]] = []
    _emit(
        emit,
        events,
        "family_entered",
        {
            "capsule_schema": CAPSULE_SCHEMA,
            "precision_dps": INTERVAL_DPS,
            "precision_bits": SOURCE_PRECISION_BITS,
            "source_ids": _json_value(capsule["source_ids"]),
            "exact_anchor_atom": True,
            "old_or_provider_calls": 0,
        },
    )

    context = _new_interval_context()
    rho_texts = ("0", "1/4", "1", "2", "4.1")

    _emit(
        emit,
        events,
        "anchored_anchor_seed_entered",
        {"family": "anchored_anchor_seed", "before_math": True},
    )
    raw_parameters, parameter_values = _raw_parameter_packet(context, parameters)
    raw_epsilon = parameters["epsilon_core"]
    source_scale, source_scale_proof = _source_scale_bound(
        context,
        parameter_values["epsilon_core"],
        parameter_values["logC"],
        raw_epsilon,
    )
    anchor_value, anchor_packet = _anchor_packet(context, anchor_raw)
    _emit(
        emit,
        events,
        "anchored_anchor_seed_complete",
        {
            "family": "anchored_anchor_seed",
            "exact_H_anchor_identity": True,
            "midpoint_selected": False,
        },
    )

    _emit(
        emit,
        events,
        "anchored_picard_polynomial_entered",
        {"family": "anchored_picard_polynomial", "before_math": True},
    )
    first = _first_picard_data(
        context,
        parameter_values,
        source_scale,
        anchor_value,
        pressure_raw,
        rho_texts,
    )
    first["epsilon_value"] = parameter_values["epsilon_core"]
    _emit(
        emit,
        events,
        "anchored_picard_polynomial_complete",
        {
            "family": "anchored_picard_polynomial",
            "degree_Phi": 2,
            "degree_DeltaU": 3,
            "degree_pressure_integral": 5,
            "absolute_F0_exp_materialized": False,
        },
    )

    _emit(
        emit,
        events,
        "anchored_residual_certificate_entered",
        {"family": "anchored_residual_certificate", "before_math": True},
    )
    raw_D, D = _capsule_record(context, banach_raw, "correction_bound", "banach.correction_bound")
    raw_q, q = _capsule_record(
        context, banach_raw, "pair_lipschitz_bound", "banach.pair_lipschitz_bound"
    )
    raw_self, self_map = _capsule_record(
        context, banach_raw, "self_map_bound", "banach.self_map_bound"
    )
    D_U = _upper(context, D)
    q_U = _upper(context, q)
    self_map_U = _upper(context, self_map)
    half = context.mpf("1") / 2
    if not bool(q_U < half):
        raise KernelInputError("pair Lipschitz upper endpoint must be below 1/2")
    if not bool(D_U < half) or not bool(self_map_U < half):
        raise KernelInputError("correction/self-map upper endpoint must be below 1/2")
    residual_bound = q_U * D_U
    fixed_error = residual_bound / (1 - q_U)
    errors = {order: _radial_error(context, fixed_error, order) for order in range(3)}
    _emit(
        emit,
        events,
        "anchored_residual_certificate_complete",
        {
            "family": "anchored_residual_certificate",
            "D_upper_endpoint_only": _display_endpoint(context, D_U),
            "q_upper_endpoint_only": _display_endpoint(context, q_U),
            "self_map_upper_endpoint_only": _display_endpoint(context, self_map_U),
            "residual_is_q_times_D": True,
            "fixed_point_error_is_qD_over_1_minus_q": True,
            "same_source_complete_pair": True,
        },
    )

    _emit(
        emit,
        events,
        "anchored_radial_value_ports_entered",
        {"family": "anchored_radial_value_ports", "before_math": True},
    )
    phi_abs_bound = context.mpf("0")
    uniform_rho = context.mpf("4.1")
    for index, coefficient in enumerate(first["phi_coefficients"]):
        phi_abs_bound = phi_abs_bound + _abs_upper(context, coefficient) * uniform_rho**index
    radial_rows: dict[str, Any] = {}
    for rho_text in rho_texts:
        rho = _exact_text(context, rho_text)
        radial_rows[rho_text] = _actual_row(context, first, rho, errors)
    rho_range = _make_interval(context, context.mpf("0"), uniform_rho)
    radial_rows["uniform_rho"] = _actual_row(
        context, first, rho_range, errors, range_row=True
    )
    _emit(
        emit,
        events,
        "anchored_radial_value_ports_complete",
        {
            "family": "anchored_radial_value_ports",
            "uniform_rho": ["0", "4.1"],
            "radial_embedding_orders": [0, 1, 2],
            "anchor_axis_structural_error_zero": True,
        },
    )

    _emit(
        emit,
        events,
        "anchored_moment_pressure_ports_entered",
        {"family": "anchored_moment_pressure_ports", "before_math": True},
    )
    for rho_text in rho_texts:
        rho = _exact_text(context, rho_text)
        radial_rows[rho_text].update(
            _moment_pressure_row(context, first, rho, errors, phi_abs_bound)
        )
    radial_rows["uniform_rho"].update(
        _moment_pressure_row(context, first, rho_range, errors, phi_abs_bound)
    )
    _emit(
        emit,
        events,
        "anchored_moment_pressure_ports_complete",
        {
            "family": "anchored_moment_pressure_ports",
            "M_z_over_epsilon_and_physical_distinct": True,
            "pressure_integral_direct_finite_convolution": True,
            "absolute_F0_exp_materialized": False,
        },
    )

    proof = {
        "kind": "native_O3_exact_anchor_first_picard_proof_v1",
        "proof_count": 1,
        "local_actual_fixed_point_anchor": True,
        "global_first_picard_image": True,
        "same_source_complete_pair": True,
        "source_operator": first["definitions"]["operator"],
        "x0_at_anchor": {
            "Phi0": "1",
            "Phi0_Z": "0",
            "Psi0": "c*rho",
            "Psi0_Z": "c_Z*rho",
        },
        "x1_is_T_of_x0": True,
        "x1_in_ball": bool(D_U < half),
        "D_and_q_are_upper_endpoint_norm_constants": True,
        "residual_norm_bound": _interval_record(context, residual_bound),
        "actual_fixed_point_error_bound": _interval_record(context, fixed_error),
        "radial_embedding_error_formula": (
            "E_k=errnorm*k!/(20^k*(1-4.1/20)^(k+1)), k=0,1,2"
        ),
        "anchor_Z_derivative_claimed": False,
        "interval_overlap_argument": False,
        "point_value_selection": False,
    }
    _emit(emit, events, "proof_bundle", proof)
    _emit(
        emit,
        events,
        "family_completed",
        {
            "proof_count": 1,
            "source_only": True,
            "metadata_only": False,
            "numeric_evaluation": True,
            "conditional": True,
            "physical_layer_completed": False,
            "global_five_moment": False,
            "time_n_recursion": False,
            "full_goal": False,
        },
    )

    parameter_data = {
        name: {
            "raw": raw_parameters[name],
            "interval": _interval_record(context, parameter_values[name]),
        }
        for name in ("delta", "epsilon_core", "j", "logC")
    }
    source_packet = {
        "formal_definition": "S=exp(-2*logC), F0(a)>0 formal atom",
        "strictly_positive_formal_atom": True,
        "materialized": False,
        "numerical_enclosure": _interval_record(context, source_scale),
        "bound_proof": source_scale_proof,
        "B_definition": "B=epsilon_core^2*S",
        "B_numerical_enclosure": _interval_record(context, first["B"]),
        "absolute_F0_value_available": False,
        "absolute_F0_exp_materialized": False,
    }
    polynomial_packet = {
        "anchor": anchor_packet,
        "d": _interval_record(context, first["d"]),
        "L": _interval_record(context, first["L"]),
        "U0": _interval_record(context, first["U0"]),
        "U0_Z": _interval_record(context, first["U0_Z"]),
        "H1": _interval_record(context, first["H1"]),
        "A": _interval_record(context, first["A"]),
        "A_Z": _interval_record(context, first["A_Z"]),
        "p": _interval_record(context, first["p"]),
        "p_Z": _interval_record(context, first["p_Z"]),
        "p_ZZ": _interval_record(context, first["p_ZZ"]),
        "H0_term": _interval_record(context, first["H0_term"]),
        "N": _interval_record(context, first["N"]),
        "N_Z": _interval_record(context, first["N_Z"]),
        "C": _interval_record(context, first["C"]),
        "C_Z": _interval_record(context, first["C_Z"]),
        "c": _interval_record(context, first["c"]),
        "c_Z": _interval_record(context, first["c_Z"]),
        "beta": _interval_record(context, first["beta"]),
        "K": _interval_record(context, first["K"]),
        "Phi1_coefficients": [
            _interval_record(context, value) for value in first["phi_coefficients"]
        ],
        "DeltaU1_coefficients": [
            _interval_record(context, value)
            for value in first["delta_u_coefficients"]
        ],
        "Psi1_coefficients": [
            _interval_record(context, value) for value in first["psi_coefficients"]
        ],
        "pressure_integral_coefficients": [
            _interval_record(context, value)
            for value in first["pressure_integral_coefficients"]
        ],
        "formal_Phi1": (
            "1-epsilon*beta*rho/4-"
            "epsilon*((1+delta)*a*C+d*C_Z)*rho^2/(24*L)"
        ),
        "formal_DeltaU1": (
            "C*rho+epsilon*K*C*rho^2/(8*L)+"
            "epsilon*(d*C*C_Z-(3+delta)*a*C^2)*rho^3/(36*L)-"
            "B*(2+delta)*a*rho^2/(4*L)"
        ),
        "formal_Psi1": "DeltaU1/epsilon",
        "formal_I": "V(Phi1^2), direct degree-5 finite convolution",
        "X0_anchor": {
            "Phi0": _interval_record(context, context.mpf("1")),
            "Phi0_Z": _interval_record(context, context.mpf("0")),
            "Psi0": "c*rho",
            "Psi0_Z": "c_Z*rho",
        },
        "X1_anchor": {
            "Phi1": "Phi1 polynomial above",
            "DeltaU1": "DeltaU1 polynomial above",
            "Psi1": "Psi1=DeltaU1/epsilon",
            "radial_projection_tail": "zero for exact first-Picard restriction",
        },
    }
    banach_packet = {
        "correction_bound": {
            "raw": _json_value(raw_D),
            "interval": _interval_record(context, D),
            "upper_endpoint_used": _display_endpoint(context, D_U),
        },
        "pair_lipschitz_bound": {
            "raw": _json_value(raw_q),
            "interval": _interval_record(context, q),
            "upper_endpoint_used": _display_endpoint(context, q_U),
        },
        "self_map_bound": {
            "raw": _json_value(raw_self),
            "interval": _interval_record(context, self_map),
            "upper_endpoint_used": _display_endpoint(context, self_map_U),
        },
        "same_source_complete_pair": True,
        "q_upper_below_half": True,
        "D_upper_below_half": True,
        "self_map_upper_below_half": True,
        "uniform_real_centers": [-1, 1],
        "rho_reference": "20",
        "norm_definition": _json_value(banach_raw.get("norm_definition")),
        "residual_norm_bound": _interval_record(context, residual_bound),
        "fixed_point_error_bound": _interval_record(context, fixed_error),
    }
    family_identity_value = _lookup(
        capsule,
        ("family_identity", "family_identity_hashes", "family_identity3hashes"),
    )
    payload = {
        "component": COMPONENT,
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "conditional": True,
        "full_physical": False,
        "global_five_moment": False,
        "time_n_recursion": False,
        "full_goal": False,
        "source_capsule_schema": CAPSULE_SCHEMA,
        "runtime_model_metadata": None,
        "api": dict(API),
        "api_metadata": dict(API_METADATA),
        "canonical_gate": {
            "algorithm": "sha256(canonical_json(entire_capsule))",
            "actual_sha256": hashlib.sha256(_canonical_bytes(capsule)).hexdigest(),
            "root_constant_used": ROOT_FREEZE_SHA,
            "canonical_field_in_capsule": False,
        },
        "source_ids": _json_value(capsule["source_ids"]),
        "parameters": parameter_data,
        "anchor": anchor_packet,
        "pressure_scaled": {
            "definition": "epsilon*physical_P0 at exact anchor; not normalized pressure",
            "P0": _interval_record(context, first["p"]),
            "P0_Z": _interval_record(context, first["p_Z"]),
            "P0_ZZ": _interval_record(context, first["p_ZZ"]),
        },
        "banach": banach_packet,
        "family_identity": _json_value(family_identity_value),
        "preservation_guards": _json_value(
            capsule.get("preservation_guards", capsule.get("guards", capsule.get("preservation")))
        ),
        "grid": {"rho": list(rho_texts), "uniform_rho": ["0", "4.1"]},
        "source_scale": source_packet,
        "anchored_polynomial": polynomial_packet,
        "actual_same_fixed_point_value_ports": {
            "rows": radial_rows,
            "error_norm": _interval_record(context, fixed_error),
            "radial_errors": {
                str(order): _interval_record(context, errors[order])
                for order in range(3)
            },
            "uniform_abs_Phi1_bound": _interval_record(context, phi_abs_bound),
            "I_actual_error_formula": (
                "rho*E0*(2*uniform_abs_Phi1+E0)"
            ),
            "ordinary_Z_derivatives_at_anchor_claimed": False,
        },
        "moment_pressure_ports": {
            "M_z_over_epsilon_is_not_ordinary_Z_derivative": True,
            "M_z_physical_definition": (
                "epsilon*(rho*U0+integral_0^rho DeltaU_actual(s) ds)"
            ),
            "M_z_error_formula": "epsilon*rho*epsilon*E0",
            "physical_M_z_rows": {
                key: value["M_z_physical_actual"] for key, value in radial_rows.items()
            },
            "M_z_over_epsilon_rows": {
                key: value["M_z_over_epsilon_actual"]
                for key, value in radial_rows.items()
            },
            "I_actual_rows": {
                key: value["I_actual"] for key, value in radial_rows.items()
            },
            "scaled_pressure_rows": {
                key: value["pressure_scaled_actual"]
                for key, value in radial_rows.items()
            },
            "formal_physical_pressure_rows": {
                key: value["pressure_physical_formal_actual"]
                for key, value in radial_rows.items()
            },
        },
        "proof_packet": proof,
        "events": events,
        "actual_counts_contract": dict(SUCCESS_COUNTS),
        "actual_flags_contract": dict(SUCCESS_FLAGS),
        "expected_counts": dict(EXPECTED_COUNTS),
        "expected_flags": dict(EXPECTED_FLAGS),
        "callback_event_counts": dict(CALLBACK_EVENT_COUNTS),
        "callback_event_count": len(CALLBACK_EVENTS),
        "manager_journal_event_count": MANAGER_JOURNAL_EVENT_COUNT,
        "journal_event_count": MANAGER_JOURNAL_EVENT_COUNT,
        "count_event_map": {
            key: list(value) for key, value in COUNT_EVENT_MAP.items()
        },
        "preservation": {
            "physical_layer_completed": False,
            "global_five_moment": False,
            "time_n_recursion": False,
            "full_goal": False,
            "absolute_F0_exp_materialized": False,
            "N_materialized": False,
        },
        "static_limitations": [
            "anchor-local first Picard image only; no full physical field",
            "global five-moment closure and time-n recursion remain false",
            "absolute F0 and N remain formal and unmaterialized",
            "pressure output is a formal physical port plus scaled numerical port",
            "no ordinary Z derivative at the exact anchor is claimed",
        ],
    }
    result_out = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "passed": True,
        "proof_count": 1,
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "conditional": True,
        "full_physical": False,
        "runtime_model_metadata": None,
        "api": dict(API),
        "api_metadata": dict(API_METADATA),
        "actual_counts": dict(SUCCESS_COUNTS),
        "new_flags": dict(SUCCESS_FLAGS),
        "expected_counts": dict(EXPECTED_COUNTS),
        "expected_flags": dict(EXPECTED_FLAGS),
        "callback_event_counts": dict(CALLBACK_EVENT_COUNTS),
        "callback_event_count": len(CALLBACK_EVENTS),
        "manager_journal_event_count": MANAGER_JOURNAL_EVENT_COUNT,
        "journal_event_count": MANAGER_JOURNAL_EVENT_COUNT,
        "count_event_map": {
            key: list(value) for key, value in COUNT_EVENT_MAP.items()
        },
        "physical_layer_completed": False,
        "global_five_moment": False,
        "time_n_recursion": False,
        "full_goal": False,
        "full_goal_complete": False,
        "payload": payload,
    }
    serialized, return_hash, return_bytes = _serialize_result(result_out)
    result_out["return_code_hash"] = return_hash
    result_out["return_code_bytes"] = return_bytes
    result_out["code_hash"] = return_hash
    result_out["code_bytes"] = return_bytes
    result_out["serialization_hash_scope"] = "canonical result before hash metadata"
    # Serialize the complete returned JSON inside this one derive call.  The
    # hash/byte pair above identifies the pre-metadata canonical return code.
    _serialize_result(result_out)
    return result_out


__all__ = [
    "API",
    "API_METADATA",
    "CALLBACK_EVENT_COUNTS",
    "CAPSULE_SCHEMA",
    "COMPONENT",
    "COUNT_EVENT_MAP",
    "EXPECTED_COUNTS",
    "EXPECTED_FLAGS",
    "FAMILY_IDENTITY",
    "INTERVAL_DPS",
    "KernelInputError",
    "MANAGER_JOURNAL_EVENT_COUNT",
    "ROOT_FREEZE_SHA",
    "SCHEMA",
    "SOURCE_IDS",
    "SOURCE_PRECISION_BITS",
    "SUCCESS_COUNTS",
    "SUCCESS_FLAGS",
    "derive",
]
