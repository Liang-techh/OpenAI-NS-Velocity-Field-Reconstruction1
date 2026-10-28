"""Exact-anchor axial Z-jet, Q-source, and physical coefficient leaf.

This module is deliberately a new leaf.  It accepts one root-frozen capsule,
reads the already accepted first-Picard source records, differentiates the
original source formula at the exact anchor, and exposes only local typed
ports.  It never imports or calls an ancestor/provider kernel, selects an
interval midpoint, materializes the absolute ``F0`` atom, or claims a global
physical field.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "native_O3_anchor_axial_Q_velocity_kernel_result_v1"
CAPSULE_SCHEMA = "native_O3_anchor_axial_Q_velocity_source_capsule_v1"
COMPONENT = "native_O3_anchor_axial_Q_velocity"
INTERVAL_DPS = 280
SOURCE_PRECISION_BITS = 867

# The root replaces this literal after freezing the complete capsule.
ROOT_FREEZE_SHA = "f16beb22b6aafb3166e0f4778df8d0fac8f8789d85d97a3b67bce52a7c041fea"

SOURCE_IDS = {
    "first_picard_result_sha256": (
        "ce516bcdd0d268509992688b3a64655f9a4a4309c337ad0445dc3e50e91108b3"
    ),
    "anchor_axis_result_sha256": (
        "c5a6f8273105a8f8710a44455de072ca93b5cd02b60f382249b8a0e5a3b016e8"
    ),
    "shared_tube_sha256": (
        "1e0de95d87056a3841c753877b51a07655d2d04420e8dc0e7a9d42ba59a608d6"
    ),
}

FAMILY_NAMES = (
    "axial_pressure_seed",
    "differentiated_picard",
    "axial_rows",
    "average_Q_pressure_gradient",
    "physical_velocity_coefficients",
)
CALLBACK_EVENTS = (
    "family_entered",
    "axial_pressure_seed_entered",
    "axial_pressure_seed_complete",
    "differentiated_picard_entered",
    "differentiated_picard_complete",
    "axial_rows_entered",
    "axial_rows_complete",
    "average_Q_pressure_gradient_entered",
    "average_Q_pressure_gradient_complete",
    "physical_velocity_coefficients_entered",
    "physical_velocity_coefficients_complete",
    "proof_bundle",
    "family_completed",
)
CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "axial_pressure_seed_entered": 1,
    "axial_pressure_seed_complete": 1,
    "differentiated_picard_entered": 1,
    "differentiated_picard_complete": 1,
    "axial_rows_entered": 1,
    "axial_rows_complete": 1,
    "average_Q_pressure_gradient_entered": 1,
    "average_Q_pressure_gradient_complete": 1,
    "physical_velocity_coefficients_entered": 1,
    "physical_velocity_coefficients_complete": 1,
    "proof_bundle": 1,
    "family_completed": 1,
}
MANAGER_JOURNAL_EVENT_COUNT = 15

# One source entry and one derive/numeric invocation contain five new
# mathematical families.  The interval-power implementation may use an
# internal exponential; no generic exp-call count is claimed here.
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
    "axial_pressure_seed_family_calls": 1,
    "differentiated_picard_family_calls": 1,
    "axial_rows_family_calls": 1,
    "average_Q_pressure_gradient_family_calls": 1,
    "physical_velocity_coefficients_family_calls": 1,
    "derive_calls": 1,
    "numeric_evaluation_calls": 1,
    "new_mathematical_proofs": 1,
    "source_value_reducer_calls": 0,
    "F0_absolute_exp_calls": 0,
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
    "local_ordinary_Z_jet": True,
    "local_Q_source": True,
    "local_velocity_coefficients": True,
    "local_ordinary_Z_jet_enclosed": True,
    "local_Q_source_enclosed": True,
    "local_velocity_coefficients_enclosed": True,
    "global_actual_source_jet_admission": False,
    "global_source_jet_admission": False,
    "Q_source_global_claimed": False,
    "physical_materialized": False,
    "physical_field_materialized": False,
    "physical_layer_completed": False,
    "global_five_moment": False,
    "time_n_recursion": False,
    "full_goal": False,
    "full_goal_complete": False,
    "generic_exp_calls_counted": False,
    "interval_power_internal_exp_allowed": True,
}

COUNT_EVENT_MAP = {
    "source_family_calls": ["family_entered"],
    "new_math_family_calls": [
        "axial_pressure_seed_entered",
        "differentiated_picard_entered",
        "axial_rows_entered",
        "average_Q_pressure_gradient_entered",
        "physical_velocity_coefficients_entered",
    ],
    "axial_pressure_seed_family_calls": ["axial_pressure_seed_entered"],
    "differentiated_picard_family_calls": ["differentiated_picard_entered"],
    "axial_rows_family_calls": ["axial_rows_entered"],
    "average_Q_pressure_gradient_family_calls": [
        "average_Q_pressure_gradient_entered"
    ],
    "physical_velocity_coefficients_family_calls": [
        "physical_velocity_coefficients_entered"
    ],
    "derive_calls": ["family_entered"],
    "numeric_evaluation_calls": ["family_entered"],
    "new_mathematical_proofs": ["proof_bundle"],
}

# Keep expected contracts literal and separate from mutable event state.
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
    "axial_pressure_seed_family_calls": 1,
    "differentiated_picard_family_calls": 1,
    "axial_rows_family_calls": 1,
    "average_Q_pressure_gradient_family_calls": 1,
    "physical_velocity_coefficients_family_calls": 1,
    "derive_calls": 1,
    "numeric_evaluation_calls": 1,
    "new_mathematical_proofs": 1,
    "source_value_reducer_calls": 0,
    "F0_absolute_exp_calls": 0,
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
    "local_ordinary_Z_jet": True,
    "local_Q_source": True,
    "local_velocity_coefficients": True,
    "local_ordinary_Z_jet_enclosed": True,
    "local_Q_source_enclosed": True,
    "local_velocity_coefficients_enclosed": True,
    "global_actual_source_jet_admission": False,
    "global_source_jet_admission": False,
    "Q_source_global_claimed": False,
    "physical_materialized": False,
    "physical_field_materialized": False,
    "physical_layer_completed": False,
    "global_five_moment": False,
    "time_n_recursion": False,
    "full_goal": False,
    "full_goal_complete": False,
    "generic_exp_calls_counted": False,
    "interval_power_internal_exp_allowed": True,
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
    """Raised when the root capsule is absent, malformed, or unauthenticated."""


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


def _contains_text(value: Any, needles: Sequence[str]) -> bool:
    if isinstance(value, str):
        compact = value.replace(" ", "")
        return any(needle.replace(" ", "") in compact for needle in needles)
    if _is_mapping(value):
        return any(_contains_text(child, needles) for child in value.values())
    if _is_sequence(value):
        return any(_contains_text(child, needles) for child in value)
    return False


def _raw_mpf_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise KernelInputError(label + " must be a four-word MPF tuple")
    try:
        words = tuple(int(item) for item in value)
    except (TypeError, ValueError) as exc:
        raise KernelInputError(label + " contains a non-integer MPF word") from exc
    return words  # type: ignore[return-value]


def _interval_value(context: Any, value: Any, label: str) -> Any:
    """Read a directed exact record and construct it from one tuple pair."""

    if not _is_mapping(value):
        raise KernelInputError(label + " must be an exact interval record")
    for nested in ("interval", "enclosure", "record", "value", "bound"):
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


def _new_interval_context() -> Any:
    # Private runtime numerical import.  The root invokes derive once.
    from mpmath.ctx_iv import MPIntervalContext

    context = MPIntervalContext()
    context.dps = INTERVAL_DPS
    return context


def _exact_text(context: Any, text: str) -> Any:
    if "/" in text:
        numerator, denominator = text.split("/", 1)
        return context.mpf(int(numerator)) / context.mpf(int(denominator))
    return context.mpf(text)


def _poly_eval(context: Any, coefficients: Sequence[Any], rho: Any) -> Any:
    value = context.mpf("0")
    for coefficient in reversed(tuple(coefficients)):
        value = coefficient + rho * value
    return value


def _finite_power_abs_bound(
    context: Any,
    coefficients: Sequence[Any],
    radius: Any,
    first_degree: int = 0,
) -> Any:
    bound = context.mpf("0")
    for offset, coefficient in enumerate(coefficients):
        degree = first_degree + offset
        bound = bound + _abs_upper(context, coefficient) * radius**degree
    return _make_interval(context, bound, bound)


def _coefficient_values(
    context: Any,
    mapping: Mapping[str, Any],
    name: str,
    label: str,
    count: int,
) -> list[Any]:
    raw = mapping.get(name)
    if raw is None:
        raise KernelInputError("missing coefficient packet " + label)
    if _is_mapping(raw):
        nested = _lookup(raw, ("coefficients", "values", "records", "items"))
        if nested is not None:
            raw = nested
        else:
            keyed: list[Any] = []
            for degree in range(count):
                candidates = (str(degree), degree)
                selected = next((raw[key] for key in candidates if key in raw), None)
                if selected is None:
                    raise KernelInputError(
                        label + " is missing coefficient degree " + str(degree)
                    )
                keyed.append(selected)
            raw = keyed
    if not _is_sequence(raw) or len(raw) != count:
        raise KernelInputError(label + " must contain exactly " + str(count) + " records")
    return [
        _interval_value(context, value, label + "[" + str(index) + "]")
        for index, value in enumerate(raw)
    ]


def _row_value(
    context: Any,
    row: Mapping[str, Any],
    names: Sequence[str],
    label: str,
) -> Any:
    selected = _lookup(row, names)
    if selected is None:
        nested = _lookup(row, ("actual_same_fixed_point_value_ports", "row"))
        if _is_mapping(nested):
            selected = _lookup(nested, names)
    if selected is None:
        raise KernelInputError(label + " is missing from parent row")
    return _interval_value(context, selected, label)


def _record_or_interval(
    context: Any,
    mapping: Mapping[str, Any],
    names: Sequence[str],
    label: str,
) -> Any:
    selected = _lookup(mapping, names)
    if selected is None:
        raise KernelInputError("missing source record " + label)
    return _interval_value(context, selected, label)


def _poly_convolution(
    context: Any,
    left: Sequence[Any],
    left_first_degree: int,
    right: Sequence[Any],
    right_first_degree: int,
) -> dict[int, Any]:
    result: dict[int, Any] = {}
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            degree = left_first_degree + i + right_first_degree + j
            if degree not in result:
                result[degree] = context.mpf("0")
            result[degree] = result[degree] + left_value * right_value
    return result


def _poly_at_degrees(context: Any, coefficients: Mapping[int, Any], rho: Any) -> Any:
    value = context.mpf("0")
    for degree in sorted(coefficients):
        value = value + coefficients[degree] * rho**degree
    return value


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


def _required_mapping(capsule: Mapping[str, Any], names: Sequence[str], label: str) -> Mapping[str, Any]:
    value = _lookup(capsule, names)
    if value is None:
        raise KernelInputError("missing source mapping " + label)
    return _as_mapping(value, label)


def _validate_capsule(capsule: Any) -> Mapping[str, Any]:
    capsule = _as_mapping(capsule, "data.capsule")
    if capsule.get("schema") != CAPSULE_SCHEMA:
        raise KernelInputError("unexpected axial Q/velocity source capsule schema")
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

    parameters = _required_mapping(capsule, ("parameters",), "capsule.parameters")
    for name in ("delta", "epsilon_core", "j", "sigma", "h"):
        if parameters.get(name) is None:
            raise KernelInputError("parameters is missing " + name)

    seed = _required_mapping(capsule, ("seed_inputs", "seed"), "capsule.seed_inputs")
    for name in (
        "a",
        "d",
        "L",
        "U0",
        "U0_Z",
        "H1",
        "A_Z",
        "C",
        "C_Z",
        "beta",
        "K",
        "p_Z",
        "p_ZZ",
        "Q0",
    ):
        if seed.get(name) is None:
            raise KernelInputError("seed_inputs is missing " + name)

    parent_polynomials = _required_mapping(
        capsule, ("parent_polynomials",), "capsule.parent_polynomials"
    )
    for name in ("Phi1_coefficients", "DeltaU1_coefficients"):
        if parent_polynomials.get(name) is None:
            raise KernelInputError("parent_polynomials is missing " + name)
    parent_rows = _required_mapping(capsule, ("parent_rows",), "capsule.parent_rows")
    for name in ("0", "1/4", "1", "2", "4.1", "uniform_rho"):
        if name not in parent_rows:
            raise KernelInputError("parent_rows is missing " + name)

    pressure_measure = _required_mapping(
        capsule, ("pressure_measure",), "capsule.pressure_measure"
    )
    for name in ("B1", "core_seed_scale"):
        if pressure_measure.get(name) is None:
            raise KernelInputError("pressure_measure is missing " + name)
    source_function = str(
        _lookup(pressure_measure, ("source_function", "source_definition", "definition"))
        or ""
    ).strip()
    if not source_function:
        raise KernelInputError("pressure_measure source function is required")
    beta_guard = _lookup(
        pressure_measure,
        ("beta_guard", "beta_range", "beta_bounds", "beta", "beta_support"),
    )
    if beta_guard is None:
        raise KernelInputError("pressure_measure beta [0,2] guard is required")
    if _is_sequence(beta_guard) and len(beta_guard) == 2:
        beta_lo, beta_hi = beta_guard[0], beta_guard[1]
    elif _is_mapping(beta_guard):
        beta_lo = _lookup(beta_guard, ("lower", "min", "lower_bound"))
        beta_hi = _lookup(beta_guard, ("upper", "max", "upper_bound"))
    else:
        beta_lo = beta_hi = None
    if str(beta_lo) not in ("0", "0.0") or str(beta_hi) not in ("2", "2.0"):
        raise KernelInputError("pressure_measure beta guard must be literal [0,2]")

    remainder = _required_mapping(capsule, ("remainder",), "capsule.remainder")
    for name in ("fixed_error_norm", "E0"):
        if remainder.get(name) is None:
            raise KernelInputError("remainder is missing " + name)
    source_scale = _required_mapping(capsule, ("source_scale",), "capsule.source_scale")
    for name in ("S_bound", "B_bound"):
        if source_scale.get(name) is None:
            raise KernelInputError("source_scale is missing " + name)
    formal_s = str(
        _lookup(
            source_scale,
            ("formalS", "formal_definition", "formal_definition_S", "S_definition"),
        )
        or ""
    ).replace(" ", "")
    if "exp(-2logC)" not in formal_s and "exp(-2*logC)" not in formal_s:
        raise KernelInputError("source_scale formal S definition is required")
    if source_scale.get("materialized") is True:
        raise KernelInputError("formal S must remain unmaterialized")

    family_identity = _lookup(
        capsule, ("family_identity", "family_identity_hashes", "family_identity3hashes")
    )
    if not _is_mapping(family_identity) or not family_identity:
        raise KernelInputError("family identity hashes are required")
    for name, value in family_identity.items():
        if not _valid_sha(value):
            raise KernelInputError("family identity hash is invalid: " + str(name))

    grid = _required_mapping(capsule, ("grid",), "capsule.grid")
    if list(grid.get("rho", ())) != ["0", "1/4", "1", "2", "4.1"]:
        raise KernelInputError("capsule rho grid does not match the frozen grid")
    if list(grid.get("uniform_rho", ())) != ["0", "4.1"]:
        raise KernelInputError("capsule uniform rho grid does not match the frozen grid")
    if list(grid.get("physical_tau", ())) != ["1", "1/4", "1/16", "1/64"]:
        raise KernelInputError("capsule physical tau grid does not match the frozen grid")
    if list(grid.get("physical_azimuth_labels", ())) != [
        "0",
        "pi/2",
        "pi",
        "3pi/2",
    ]:
        raise KernelInputError("capsule physical azimuth grid does not match the frozen grid")

    guards = _lookup(capsule, ("preservation_guards", "guards", "preservation"))
    if guards is None:
        raise KernelInputError("preservation guards are required")
    guards = _as_mapping(guards, "capsule.preservation_guards")
    false_names = (
        "global_actual_source_jet_admission",
        "global_source_jet_admission",
        "physical_materialized",
        "physical_field_materialized",
        "global_five_moment",
        "global_five_moment_closure",
        "time_n_recursion",
        "full_goal",
        "full_goal_complete",
    )
    for name in false_names:
        if name in guards and guards[name] is not False:
            raise KernelInputError("preservation guard is not false: " + name)
    for name in (
        "N_materialization_calls",
        "absolute_F0_exponential_materialization_calls",
        "F0_absolute_exp_calls",
        "ancestor_calls",
        "provider_calls",
        "checker_calls",
        "selftest_calls",
    ):
        if name in guards and guards[name] != 0:
            raise KernelInputError("preservation count is not zero: " + name)

    definitions = _lookup(
        capsule, ("definition_guards", "definition", "definitions")
    )
    if definitions is not None and _is_mapping(definitions):
        for name in ("ordinary_Z_frame", "source_theta", "physical_azimuth"):
            if name in definitions and not str(definitions[name]).strip():
                raise KernelInputError("empty definition guard: " + name)
    recipe = _lookup(
        capsule,
        (
            "source_recipe",
            "source_recipe_proof_metadata",
            "source_recipe_literal_identity",
        ),
    )
    if recipe is None or not _is_mapping(recipe):
        raise KernelInputError("source recipe proof metadata is required")

    anchor = _lookup(capsule, ("anchor", "anchor_definition"))
    anchor_definition_seen = False
    if anchor is not None and _is_mapping(anchor):
        definition = str(
            _lookup(anchor, ("exact_definition", "definition", "root_definition")) or ""
        ).replace(" ", "")
        if definition and "H(a)=0" not in definition and "H0(a)=0" not in definition:
            raise KernelInputError("exact H(a)=0 anchor definition is missing")
        anchor_definition_seen = "H(a)=0" in definition or "H0(a)=0" in definition
        if anchor.get("midpoint_selected", False) is True or anchor.get(
            "root_interval_midpoint_selected", False
        ) is True:
            raise KernelInputError("anchor midpoint selection is forbidden")
    if seed.get("midpoint_selected", False) is True:
        raise KernelInputError("seed midpoint selection is forbidden")
    if not anchor_definition_seen and not _contains_text(
        (definitions, recipe), ("H(a)=0", "H0(a)=0")
    ):
        raise KernelInputError("exact H(a)=0 anchor definition is required")
    if _contains_text((definitions, recipe), ("midpoint_selected=true", "midpoint selected")):
        raise KernelInputError("anchor midpoint selection is forbidden")
    return capsule


def _parameter_packet(
    context: Any, parameters: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    raw_packet: dict[str, Any] = {}
    packet: dict[str, Any] = {}
    for name in ("delta", "epsilon_core", "j", "sigma", "h"):
        raw, value = _parameter_record(context, parameters, name)
        raw_packet[name] = _json_value(raw)
        packet[name] = value
    for name in ("delta", "epsilon_core", "j", "sigma", "h"):
        if not _strict_positive(context, packet[name]):
            raise KernelInputError(name + " must be strictly positive")
    if not bool(_upper(context, packet["delta"]) < context.mpf("1")):
        raise KernelInputError("delta must be below one")
    return raw_packet, packet


def _seed_packet(
    context: Any, seed: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    names = (
        "a",
        "d",
        "L",
        "U0",
        "U0_Z",
        "H1",
        "A_Z",
        "C",
        "C_Z",
        "beta",
        "K",
        "p_Z",
        "p_ZZ",
        "Q0",
    )
    raw_packet: dict[str, Any] = {}
    packet: dict[str, Any] = {}
    for name in names:
        raw, value = _capsule_record(context, seed, name, "seed_inputs." + name)
        raw_packet[name] = _json_value(raw)
        packet[name] = value
    if not bool(_upper(context, packet["a"]) < context.mpf("0")):
        raise KernelInputError("exact anchor must be strictly negative")
    for name in ("d", "L"):
        if not _strict_positive(context, packet[name]):
            raise KernelInputError("seed " + name + " must be strictly positive")
    return raw_packet, packet


def _pressure_seed(
    context: Any,
    parameters: Mapping[str, Any],
    seed: Mapping[str, Any],
    pressure_measure: Mapping[str, Any],
) -> dict[str, Any]:
    """Build C_ZZ from the accepted pressure third-derivative enclosure."""

    delta = parameters["delta"]
    epsilon = parameters["epsilon_core"]
    j = parameters["j"]
    a = seed["a"]
    d = seed["d"]
    L = seed["L"]
    U0 = seed["U0"]
    C = seed["C"]
    C_Z = seed["C_Z"]
    p_Z = seed["p_Z"]
    p_ZZ = seed["p_ZZ"]
    _, B1 = _capsule_record(
        context, pressure_measure, "B1", "pressure_measure.B1"
    )
    _, core_seed_scale = _capsule_record(
        context,
        pressure_measure,
        "core_seed_scale",
        "pressure_measure.core_seed_scale",
    )
    qP = 1 + a**2
    qP_lower = _lower(context, qP)
    qP_upper = _upper(context, qP)
    I2 = _make_interval(
        context,
        _lower(context, B1) / qP_upper**4,
        3 * _upper(context, B1) / qP_lower**2,
    )
    I3 = _make_interval(
        context,
        2 * _lower(context, B1) / qP_upper**5,
        12 * _upper(context, B1) / qP_lower**3,
    )
    Pnorm3 = -12 * a * I2 + 8 * a**3 * I3
    p3 = core_seed_scale * Pnorm3
    H2 = -24 * a - 2 * j
    A_ZZ = -8 * (1 + delta)
    N_ZZ = (
        epsilon * (A_ZZ * U0 + 8 * seed["A_Z"] - 4 * H2)
        + (6 + 4 * delta) * p_Z
        + (6 + 2 * delta) * a * p_ZZ
        - d * p3
    )
    C_ZZ = N_ZZ / (2 * L) + (4 * delta * a / L) * C_Z + (2 * delta / L) * C
    return {
        "B1": B1,
        "core_seed_scale": core_seed_scale,
        "qP": qP,
        "I2": I2,
        "I3": I3,
        "Pnorm3": Pnorm3,
        "p3": p3,
        "H2": H2,
        "A_ZZ": A_ZZ,
        "N_ZZ": N_ZZ,
        "C_ZZ": C_ZZ,
    }


def _differentiate_picard(
    context: Any,
    parameters: Mapping[str, Any],
    seed: Mapping[str, Any],
    pressure_seed: Mapping[str, Any],
    source_scale: Mapping[str, Any],
) -> dict[str, Any]:
    """Differentiate the original coupled source before restricting H(a)=0."""

    delta = parameters["delta"]
    epsilon = parameters["epsilon_core"]
    sigma = parameters["sigma"]
    a = seed["a"]
    d = seed["d"]
    L = seed["L"]
    U0 = seed["U0"]
    H1 = seed["H1"]
    C = seed["C"]
    C_Z = seed["C_Z"]
    C_ZZ = pressure_seed["C_ZZ"]
    beta = seed["beta"]
    K = seed["K"]
    B = _record_or_interval(
        context,
        source_scale,
        ("B_bound", "B_numerical_enclosure", "B"),
        "source_scale.B_bound",
    )
    L_Z = -2 * delta * a
    beta_Z = (parameters["j"] + 2 * delta * a * beta) / L
    K_Z = -(3 + delta) * (U0 + 4 * a)
    chi_ZZ = 2 * H1**2 / sigma**2
    Qtheta = (1 + delta) * a * C + d * C_Z
    Qhat = (1 + delta) * C - (1 - delta) * a * C_Z + d * C_ZZ
    D = d * C * C_Z - (3 + delta) * a * C**2
    Dstar = (
        d * (C * C_ZZ + C_Z**2)
        - 2 * (4 + delta) * a * C * C_Z
        - (3 + delta) * C**2
    )
    fz1 = -epsilon * beta_Z / 4
    fz2 = (
        epsilon * L_Z * Qtheta / (24 * L**2)
        - epsilon * Qhat / (24 * L)
        - d * H1 * C / (12 * sigma**2)
    )
    fz3 = -epsilon * d * C * chi_ZZ / (96 * L)
    uz1 = C_Z
    uz2 = (
        epsilon
        * (K_Z * C + (K + H1) * C_Z - (L_Z / L) * K * C)
        / (8 * L)
        + B * (2 + delta) * (a * L_Z - L) / (4 * L**2)
        - B * d * H1 / (4 * epsilon * sigma**2)
    )
    uz3 = (
        epsilon * (Dstar - (L_Z / L) * D) / (36 * L)
        - B * d * chi_ZZ / (72 * L)
    )
    S_ZZ_over_S = -2 * L * H1 / (epsilon * sigma**2)
    return {
        "L_Z": L_Z,
        "beta_Z": beta_Z,
        "K_Z": K_Z,
        "chi_ZZ": chi_ZZ,
        "Qtheta": Qtheta,
        "Qhat": Qhat,
        "D": D,
        "Dstar": Dstar,
        "B": B,
        "S_ZZ_over_S": S_ZZ_over_S,
        "Phi1_Z_coefficients": (fz1, fz2, fz3),
        "DeltaU1_Z_coefficients": (uz1, uz2, uz3),
        "definitions": {
            "source_frame": "fixed rho ordinary Z derivative before H(a)=0 restriction",
            "Qtheta": "(1+delta)*a*C+d*C_Z; not Q_source",
            "Qhat": "(1+delta)*C-(1-delta)*a*C_Z+d*C_ZZ",
            "D": "d*C*C_Z-(3+delta)*a*C^2",
            "Dstar": "d*(C*C_ZZ+C_Z^2)-2*(4+delta)*a*C*C_Z-(3+delta)*C^2",
            "source_Phi0_ZZ": "-chiZZ*rho/4",
            "source_S_ZZ_over_S": "-2*L*H1/(epsilon*sigma^2)",
            "Phi1_Z": (
                "fz1*rho+fz2*rho^2+fz3*rho^3; "
                "coefficients are not divided by F0"
            ),
            "DeltaU1_Z": "uz1*rho+uz2*rho^2+uz3*rho^3",
        },
    }


def _axial_row(
    context: Any,
    parent_row: Mapping[str, Any],
    rho: Any,
    differentiated: Mapping[str, Any],
    epsilon: Any,
    ez_bound: Any,
    range_row: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    fz = differentiated["Phi1_Z_coefficients"]
    uz = differentiated["DeltaU1_Z_coefficients"]
    fz_candidate = rho * _poly_eval(context, fz, rho)
    uz_candidate = rho * _poly_eval(context, uz, rho)
    if _exact_zero(context, rho):
        phi_z_candidate = context.mpf("0")
        phi_z_actual = context.mpf("0")
        delta_u_z_candidate = context.mpf("0")
        delta_u_z_actual = context.mpf("0")
        u_z_actual = context.mpf("4")
        axis_exact = True
    else:
        phi_z_candidate = fz_candidate
        phi_z_actual = _widen(context, fz_candidate, ez_bound)
        delta_u_z_candidate = uz_candidate
        delta_u_z_actual = _widen(context, uz_candidate, epsilon * ez_bound)
        u_z_actual = 4 + delta_u_z_actual
        axis_exact = False
    output = dict(_json_value(parent_row))
    output.update(
        {
            "Phi1_Z_candidate": _interval_record(context, phi_z_candidate),
            "Phi1_Z_actual": _interval_record(context, phi_z_actual),
            "Phi_Z_actual": _interval_record(context, phi_z_actual),
            "DeltaU1_Z_candidate": _interval_record(context, delta_u_z_candidate),
            "DeltaU1_Z_actual": _interval_record(context, delta_u_z_actual),
            "DeltaU_Z_actual": _interval_record(context, delta_u_z_actual),
            "U_Z_actual": _interval_record(context, u_z_actual),
            "Phi_Z_error_bound": _interval_record(context, ez_bound),
            "DeltaU_Z_error_bound": _interval_record(context, epsilon * ez_bound),
            "F_Z_formal": {
                "formal_expression": "F0*Phi_Z at anchor, with F0_Z=0",
                "Phi_Z_coefficient": _interval_record(context, phi_z_actual),
                "absolute_F0_materialized": False,
                "point_value": False,
            },
            "rho": _interval_record(context, rho),
            "range_row": range_row,
            "ordinary_Z_frame": "fixed rho",
            "structural_axis_error_zero": axis_exact,
            "point_value": False,
        }
    )
    return output, {
        "rho": rho,
        "parent_row": parent_row,
        "Phi_actual": _row_value(
            context, parent_row, ("Phi_actual",), "parent_row.Phi_actual"
        ),
        "DeltaPhi_actual": _row_value(
            context,
            parent_row,
            ("DeltaPhi_actual", "DeltaPhi1_actual"),
            "parent_row.DeltaPhi_actual",
        ),
        "DeltaU_actual": _row_value(
            context, parent_row, ("DeltaU_actual",), "parent_row.DeltaU_actual"
        ),
        "U_actual": _row_value(
            context,
            parent_row,
            ("U1_actual", "U_actual"),
            "parent_row.U1_actual",
        ),
        "Phi_Z_actual": phi_z_actual,
        "DeltaU_Z_actual": delta_u_z_actual,
        "U_Z_actual": u_z_actual,
        "range_row": range_row,
    }


def _average_candidate(
    context: Any, coefficients: Sequence[Any], rho: Any
) -> Any:
    """Integrate a degree 1..3 coefficient packet divided by k+1."""

    return sum(
        (coefficient * rho ** (index + 1) / (index + 2)
         for index, coefficient in enumerate(coefficients)),
        context.mpf("0"),
    )


def _affine_F0_reference(
    expression: str,
    base_field: str,
    base_sign: int,
    f0_field: str,
    f0_sign: int,
) -> dict[str, Any]:
    """Typed affine port referring to one canonical physical row."""

    return {
        "formal_expression": expression,
        "constant_term_ref": {"field": base_field, "sign": base_sign},
        "F0_coefficient_ref": {"field": f0_field, "sign": f0_sign},
        "canonical_numeric_coefficients": True,
        "absolute_F0_materialized": False,
        "point_value": False,
    }


def _average_Q_pressure(
    context: Any,
    parameters: Mapping[str, Any],
    seed: Mapping[str, Any],
    source_scale: Mapping[str, Any],
    parent_coefficients: Sequence[Any],
    parent_polynomials: Mapping[str, Any],
    axial_rows: Mapping[str, Mapping[str, Any]],
    differentiated: Mapping[str, Any],
    e0: Any,
    ez_bound: Any,
) -> tuple[dict[str, Any], dict[int, Any], Any]:
    """Expose local radial-average Q and the differentiated pressure port."""

    delta = parameters["delta"]
    epsilon = parameters["epsilon_core"]
    a = seed["a"]
    d = seed["d"]
    L = seed["L"]
    U0 = seed["U0"]
    Q0 = seed["Q0"]
    B = _record_or_interval(
        context,
        source_scale,
        ("B_bound", "B_numerical_enclosure", "B"),
        "source_scale.B_bound",
    )
    # The source Phi polynomial is accepted data.  Only its finite
    # convolution with the new Z coefficients is formed here.
    old_phi = _coefficient_values(
        context,
        parent_polynomials,
        "Phi1_coefficients",
        "parent_polynomials.Phi1_coefficients",
        3,
    )
    new_phi_z = differentiated["Phi1_Z_coefficients"]
    convolution = _poly_convolution(context, old_phi, 0, new_phi_z, 1)
    I_Z_coefficients = {
        degree + 1: 2 * coefficient / (degree + 1)
        for degree, coefficient in convolution.items()
    }
    radius = _exact_text(context, "4.1")
    uniform_abs_phi_z = _finite_power_abs_bound(context, new_phi_z, radius, 1)
    parent_abs = _lookup(
        parent_polynomials,
        ("uniform_abs_Phi1", "uniform_abs_Phi1_bound"),
    )
    if parent_abs is None:
        raise KernelInputError(
            "accepted parent uniform_abs_Phi1_bound is required; "
            "old Phi1 bound recomputation is forbidden"
        )
    uniform_abs_phi = _interval_value(
        context, parent_abs, "parent_polynomials.uniform_abs_Phi1_bound"
    )
    rows_out: dict[str, Any] = {}
    for rho_text, row_state in axial_rows.items():
        if rho_text.startswith("__"):
            continue
        rho = row_state["rho"]
        parent_row = row_state["parent_row"]
        delta_u_actual = row_state["DeltaU_actual"]
        delta_m_candidate = _average_candidate(context, parent_coefficients, rho)
        delta_m_error = epsilon * e0
        if _exact_zero(context, rho):
            delta_m_actual = context.mpf("0")
            m_average = U0
        else:
            delta_m_actual = _widen(context, delta_m_candidate, delta_m_error)
            m_average = U0 + delta_m_actual
        new_uz_coeff = differentiated["DeltaU1_Z_coefficients"]
        delta_m_z_candidate = _average_candidate(context, new_uz_coeff, rho)
        delta_m_z_error = epsilon * ez_bound
        if _exact_zero(context, rho):
            delta_m_z_actual = context.mpf("0")
            m_average_z = context.mpf("4")
        else:
            delta_m_z_actual = _widen(
                context, delta_m_z_candidate, delta_m_z_error
            )
            m_average_z = 4 + delta_m_z_actual
        if _exact_zero(context, rho):
            delta_q = context.mpf("0")
            q_source = Q0
        else:
            delta_q = (
                2 * a * delta_u_actual
                - (1 - delta) * a * delta_m_actual
                - d * delta_m_z_actual
            ) / L
            q_source = Q0 + delta_q

        iz_candidate = _poly_at_degrees(context, I_Z_coefficients, rho)
        rho_upper = _upper(context, rho)
        iz_error_upper = 2 * rho_upper * (
            _upper(context, e0) * _upper(context, uniform_abs_phi_z)
            + _upper(context, ez_bound) * _upper(context, uniform_abs_phi)
            + _upper(context, e0) * _upper(context, ez_bound)
        )
        iz_error = _make_interval(context, context.mpf("0"), iz_error_upper)
        if _exact_zero(context, rho):
            iz_actual = context.mpf("0")
        else:
            iz_actual = _widen(context, iz_candidate, iz_error)
        epsilon_p_z = seed["p_Z"] + B * iz_actual
        rows_out[rho_text] = {
            "rho": _interval_record(context, rho),
            "parent_base_reused": True,
            "DeltaM_average_candidate": _interval_record(context, delta_m_candidate),
            "DeltaM_average_actual": _interval_record(context, delta_m_actual),
            "M_average": _interval_record(context, m_average),
            "DeltaM_average_Z_candidate": _interval_record(
                context, delta_m_z_candidate
            ),
            "DeltaM_average_Z_actual": _interval_record(
                context, delta_m_z_actual
            ),
            "M_average_Z": _interval_record(context, m_average_z),
            "DeltaQ": _interval_record(context, delta_q),
            "Q_source": _interval_record(context, q_source),
            "I1_Z_candidate": _interval_record(context, iz_candidate),
            "I1_Z_actual": _interval_record(context, iz_actual),
            "I1_Z_error_bound": _interval_record(context, iz_error),
            "epsilonP_Z": _interval_record(context, epsilon_p_z),
            "formal_incrementP_Z": {
                "formal_expression": "epsilon_core*S*I1_Z_actual",
                "S_formal": "exp(-2*logC)",
                "source_S_Z_at_anchor": "0",
                "materialized": False,
                "point_value": False,
            },
            "DeltaM_average_error_bound": _interval_record(
                context, delta_m_error
            ),
            "DeltaM_average_Z_error_bound": _interval_record(
                context, delta_m_z_error
            ),
            "ordinary_Q_definition": "Q_source=[2*Z*U-(1-delta)*Z*M-d*M_Z]/L",
            "M_definition": "U0+(1/rho)*integral_0^rho DeltaU(s,Z) ds",
            "M_Z_definition": "4+(1/rho)*integral_0^rho DeltaU_Z(s,Z) ds",
            "axis_Q_source_not_forced_zero": True,
            "axis_Q_source_exact_Q0": _exact_zero(context, rho),
            "point_value": False,
        }
    return (
        rows_out,
        I_Z_coefficients,
        uniform_abs_phi_z,
    )


def _physical_velocity_ports(
    context: Any,
    parameters: Mapping[str, Any],
    seed: Mapping[str, Any],
    axial_states: Mapping[str, Mapping[str, Any]],
    average_rows: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Create typed physical coefficient ports without materializing F0."""

    delta = parameters["delta"]
    epsilon = parameters["epsilon_core"]
    d = seed["d"]
    a = seed["a"]
    tau_labels = ("1", "1/4", "1/16", "1/64")
    rho_labels = ("0", "1/4", "1", "2", "4.1")
    angle_labels = ("0", "pi/2", "pi", "3pi/2")
    rotations = {
        "0": ("1", "0"),
        "pi/2": ("0", "1"),
        "pi": ("-1", "0"),
        "3pi/2": ("0", "-1"),
    }
    time_slices: dict[str, Any] = {}
    for tau_label in tau_labels:
        tau = _exact_text(context, tau_label)
        if not bool(_lower(context, tau) > context.mpf("0")):
            raise KernelInputError("physical tau must be positive")
        t = 1 - tau
        lambda_phys = context.sqrt(tau / d)
        kappa = lambda_phys ** (1 - delta)
        rho_slices: dict[str, Any] = {}
        for rho_label in rho_labels + ("uniform_rho",):
            state = axial_states[rho_label]
            avg = average_rows[rho_label]
            rho = state["rho"]
            R = epsilon * rho
            r_phys = lambda_phys * context.sqrt(2 * R)
            z_phys = kappa * a
            source_Q = _interval_value(
                context, avg["Q_source"], "average_Q_pressure_gradient.Q_source"
            )
            M_average = _interval_value(
                context, avg["M_average"], "average_Q_pressure_gradient.M_average"
            )
            parent_U = state["U_actual"]
            parent_delta_U = state["DeltaU_actual"]
            parent_phi = state["Phi_actual"]
            parent_delta_phi = state["DeltaPhi_actual"]
            delta_Q = _interval_value(
                context, avg["DeltaQ"], "average_Q_pressure_gradient.DeltaQ"
            )
            if _exact_zero(context, rho):
                ur = context.mpf("0")
                delta_ur = context.mpf("0")
                swirl = context.mpf("0")
                delta_swirl = context.mpf("0")
                source_U = seed["U0"]
                parent_U = seed["U0"]
                parent_delta_U = context.mpf("0")
                parent_phi = context.mpf("1")
                parent_delta_phi = context.mpf("0")
            else:
                ur = lambda_phys**(-1) * context.sqrt(R / 2) * source_Q
                delta_ur = lambda_phys**(-1) * context.sqrt(R / 2) * delta_Q
                source_U = parent_U
                swirl = (
                    lambda_phys ** (-1 - delta)
                    * context.sqrt(2 * R)
                    * parent_phi
                )
                delta_swirl = (
                    lambda_phys ** (-1 - delta)
                    * context.sqrt(2 * R)
                    * parent_delta_phi
                )
            uz = lambda_phys ** (-1 - delta) * source_U
            delta_uz = lambda_phys ** (-1 - delta) * parent_delta_U
            psi = kappa * R * M_average
            angle_ports: dict[str, Any] = {}
            for angle_label in angle_labels:
                cos_text, sin_text = rotations[angle_label]
                cos_theta = _exact_text(context, cos_text)
                sin_theta = _exact_text(context, sin_text)
                base_u_sign = int(cos_text)
                base_v_sign = int(sin_text)
                f0_u_sign = -int(sin_text)
                f0_v_sign = int(cos_text)
                angle_ports[angle_label] = {
                    "theta_frame": "physical_azimuth",
                    "source_phase_is_not_physical_azimuth": True,
                    "rotation_coefficients": {
                        "cos_theta": _interval_record(context, cos_theta),
                        "sin_theta": _interval_record(context, sin_theta),
                    },
                    "u": _affine_F0_reference(
                        "ur*cos(theta_phys)-F0*swirlcoef*sin(theta_phys)",
                        "ur",
                        base_u_sign,
                        "swirlcoef",
                        f0_u_sign,
                    ),
                    "v": _affine_F0_reference(
                        "ur*sin(theta_phys)+F0*swirlcoef*cos(theta_phys)",
                        "ur",
                        base_v_sign,
                        "swirlcoef",
                        f0_v_sign,
                    ),
                    "u_const": {"field": "ur", "sign": base_u_sign},
                    "u_F0coef": {"field": "swirlcoef", "sign": f0_u_sign},
                    "v_const": {"field": "ur", "sign": base_v_sign},
                    "v_F0coef": {"field": "swirlcoef", "sign": f0_v_sign},
                    "w": {
                        "canonical_field_ref": "uz",
                        "point_value": False,
                    },
                    "point_value": False,
                }
            rho_slices[rho_label] = {
                "rho": _interval_record(context, rho),
                "R": _interval_record(context, R),
                "r_phys": _interval_record(context, r_phys),
                "z_phys": _interval_record(context, z_phys),
                "source_Q": _interval_record(context, source_Q),
                "source_U": _interval_record(context, source_U),
                "ur": _interval_record(context, ur),
                "Deltaur": _interval_record(context, delta_ur),
                "uz": _interval_record(context, uz),
                "DeltaUz": _interval_record(context, delta_uz),
                "swirlcoef": _interval_record(context, swirl),
                "delta_swirlcoef": _interval_record(context, delta_swirl),
                "actualutheta": {
                    "formal_expression": "F0*swirlcoef",
                    "F0_coefficient": _interval_record(context, swirl),
                    "absolute_F0_materialized": False,
                    "positive_formal_factor": True,
                    "point_value": False,
                },
                "streamfunction": {
                    "psi": _interval_record(context, psi),
                    "formal_expression": "psi=kappa*R*M_average",
                    "ur_pullback": "ur=-psi_z/r_phys",
                    "uz_pullback": "uz=psi_r/r_phys",
                    "mixed_partials_divergence_zero": True,
                    "numeric_cartesian_NS_check": False,
                    "point_value": False,
                },
                "angles": angle_ports,
                "axis_regularization": {
                    "r_zero": _exact_zero(context, rho),
                    "ur_zero": _exact_zero(context, rho),
                    "swirlcoef_zero": _exact_zero(context, rho),
                    "axis_source_U_is_U0": _exact_zero(context, rho),
                    "axis_Q_source_not_forced_zero": True,
                },
                "point_value": False,
            }
        time_slices[tau_label] = {
            "tau": _interval_record(context, tau),
            "t": _interval_record(context, t),
            "lambda_phys": _interval_record(context, lambda_phys),
            "kappa": _interval_record(context, kappa),
            "rho_slices": rho_slices,
            "positive_tau": True,
            "positive_power_operations_allowed": True,
            "point_value": False,
        }
    return {
        "domain": {
            "Z": "a exact anchor atom only",
            "rho": ["0", "1/4", "1", "2", "4.1"],
            "uniform_rho": ["0", "4.1"],
            "physical_tau": list(tau_labels),
            "physical_azimuth_labels": list(angle_labels),
            "physical_theta_is_free_azimuth": True,
            "source_phase_is_not_physical_azimuth": True,
            "source_profile_scale": "Lambda_source=1/epsilon_core",
            "physical_inverse_scale": "lambda_phys=sqrt(tau/d); distinct from Lambda_source and loop q",
            "pressure_q_distinct_from_Q_source": True,
            "loop_q_distinct_from_Q_source": True,
        },
        "time_slices": time_slices,
        "conditional_streamfunction_divergence": {
            "frame": "physical cylindrical",
            "streamfunction": "psi=kappa*R*M_average",
            "regularity": "C2",
            "theta_mode": "physical_azimuth_independent",
            "phase_mode": "ordinary_Z",
            "certificate": "mixed partial cancellation; structural conditional only",
            "numeric_cartesian_NS_claim": False,
            "full_physical_field_materialized": False,
        },
        "viscosity": {"paper_unit": 1, "residual_defined_force": False},
        "point_value": False,
    }


def derive(
    data: Mapping[str, Any],
    emit: Callable[[str, Mapping[str, Any]], Any] | None = None,
) -> dict[str, Any]:
    """Evaluate one bounded exact-anchor axial/Q/velocity source leaf."""

    if not _is_mapping(data) or "capsule" not in data:
        raise KernelInputError("derive requires data['capsule']")
    if emit is not None and not callable(emit):
        raise KernelInputError("emit must be callable or None")

    capsule = _validate_capsule(data["capsule"])
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
            "generic_exp_calls_counted": False,
        },
    )

    context = _new_interval_context()
    rho_texts = ("0", "1/4", "1", "2", "4.1")
    grid = _as_mapping(capsule["grid"], "capsule.grid")
    parameters = _as_mapping(capsule["parameters"], "capsule.parameters")
    seed_source = _as_mapping(
        _lookup(capsule, ("seed_inputs", "seed")), "capsule.seed_inputs"
    )
    parent_polynomials = _as_mapping(
        capsule["parent_polynomials"], "capsule.parent_polynomials"
    )
    parent_rows_source = _as_mapping(capsule["parent_rows"], "capsule.parent_rows")
    pressure_measure = _as_mapping(
        capsule["pressure_measure"], "capsule.pressure_measure"
    )
    remainder = _as_mapping(capsule["remainder"], "capsule.remainder")
    source_scale_source = _as_mapping(capsule["source_scale"], "capsule.source_scale")

    _emit(
        emit,
        events,
        "axial_pressure_seed_entered",
        {"family": "axial_pressure_seed", "before_math": True},
    )
    raw_parameters, parameter_values = _parameter_packet(context, parameters)
    raw_seed, seed_values = _seed_packet(context, seed_source)
    raw_B1, B1 = _capsule_record(
        context, pressure_measure, "B1", "pressure_measure.B1"
    )
    raw_core_seed, core_seed_scale = _capsule_record(
        context,
        pressure_measure,
        "core_seed_scale",
        "pressure_measure.core_seed_scale",
    )
    pressure_seed = _pressure_seed(
        context, parameter_values, seed_values, pressure_measure
    )
    _emit(
        emit,
        events,
        "axial_pressure_seed_complete",
        {
            "family": "axial_pressure_seed",
            "exact_anchor_negative": True,
            "midpoint_selected": False,
            "pressure_derivative_order": 3,
            "p3_is_core_seed_scaled": True,
        },
    )

    _emit(
        emit,
        events,
        "differentiated_picard_entered",
        {"family": "differentiated_picard", "before_math": True},
    )
    source_scale_records = {
        "S_bound": source_scale_source["S_bound"],
        "B_bound": source_scale_source["B_bound"],
    }
    differentiated = _differentiate_picard(
        context,
        parameter_values,
        seed_values,
        pressure_seed,
        source_scale_records,
    )
    _emit(
        emit,
        events,
        "differentiated_picard_complete",
        {
            "family": "differentiated_picard",
            "source_derivative_before_anchor_restriction": True,
            "ordinary_Z_frame": "fixed rho",
            "Phi_Z_degree_range": [1, 3],
            "DeltaU_Z_degree_range": [1, 3],
            "absolute_F0_exp_materialized": False,
        },
    )

    _emit(
        emit,
        events,
        "axial_rows_entered",
        {"family": "axial_rows", "before_math": True},
    )
    _, fixed_error_norm = _capsule_record(
        context, remainder, "fixed_error_norm", "remainder.fixed_error_norm"
    )
    _, E0 = _capsule_record(context, remainder, "E0", "remainder.E0")
    _, h = _parameter_record(context, parameters, "h")
    ratio = 1 - _exact_text(context, "4.1") / 20
    ez_upper = _upper(context, fixed_error_norm) / (4 * _lower(context, h) * ratio)
    EzBound = _make_interval(context, context.mpf("0"), _upper(context, ez_upper))
    parent_rows: dict[str, Mapping[str, Any]] = {}
    for key in rho_texts + ("uniform_rho",):
        parent_rows[key] = _as_mapping(parent_rows_source[key], "parent_rows." + key)
    axial_rows_payload: dict[str, Any] = {}
    axial_states: dict[str, Mapping[str, Any]] = {}
    for rho_text in rho_texts:
        rho = _exact_text(context, rho_text)
        payload_row, state = _axial_row(
            context,
            parent_rows[rho_text],
            rho,
            differentiated,
            parameter_values["epsilon_core"],
            EzBound,
            False,
        )
        axial_rows_payload[rho_text] = payload_row
        axial_states[rho_text] = state
    rho_range = _make_interval(
        context, context.mpf("0"), _exact_text(context, "4.1")
    )
    payload_row, state = _axial_row(
        context,
        parent_rows["uniform_rho"],
        rho_range,
        differentiated,
        parameter_values["epsilon_core"],
        EzBound,
        True,
    )
    axial_rows_payload["uniform_rho"] = payload_row
    axial_states["uniform_rho"] = state
    _emit(
        emit,
        events,
        "axial_rows_complete",
        {
            "family": "axial_rows",
            "rows": list(rho_texts) + ["uniform_rho"],
            "axis_exact": True,
            "new_error_bound_from_weighted_norm": True,
            "Phi_Z_coefficient_not_divided_by_F0": True,
        },
    )

    _emit(
        emit,
        events,
        "average_Q_pressure_gradient_entered",
        {"family": "average_Q_pressure_gradient", "before_math": True},
    )
    parent_delta_u_coefficients = _coefficient_values(
        context,
        parent_polynomials,
        "DeltaU1_coefficients",
        "parent_polynomials.DeltaU1_coefficients",
        3,
    )
    parent_polynomials_for_calc: Mapping[str, Any] = parent_polynomials
    if (
        "uniform_abs_Phi1" not in parent_polynomials
        and "uniform_abs_Phi1_bound" not in parent_polynomials
        and (
            capsule.get("uniform_abs_Phi1_bound") is not None
            or remainder.get("uniform_abs_Phi1_bound") is not None
        )
    ):
        parent_polynomials_for_calc = dict(parent_polynomials)
        parent_polynomials_for_calc["uniform_abs_Phi1_bound"] = _lookup(
            capsule, ("uniform_abs_Phi1_bound",)
        ) or remainder["uniform_abs_Phi1_bound"]
    average_rows, I_Z_coefficients, uniform_abs_phi_z = _average_Q_pressure(
        context,
        parameter_values,
        seed_values,
        source_scale_records,
        parent_delta_u_coefficients,
        parent_polynomials_for_calc,
        axial_states,
        differentiated,
        E0,
        EzBound,
    )
    _emit(
        emit,
        events,
        "average_Q_pressure_gradient_complete",
        {
            "family": "average_Q_pressure_gradient",
            "ordinary_Q_source_anchor_only": True,
            "axis_DeltaQ_exact_zero": True,
            "axis_Q0_retained": True,
            "pressure_gradient_is_new_finite_convolution": True,
            "global_Q_source_claimed": False,
        },
    )

    _emit(
        emit,
        events,
        "physical_velocity_coefficients_entered",
        {"family": "physical_velocity_coefficients", "before_math": True},
    )
    physical_ports = _physical_velocity_ports(
        context,
        parameter_values,
        seed_values,
        axial_states,
        average_rows,
    )
    _emit(
        emit,
        events,
        "physical_velocity_coefficients_complete",
        {
            "family": "physical_velocity_coefficients",
            "tau_slices": ["1", "1/4", "1/16", "1/64"],
            "rho_slices": list(rho_texts) + ["uniform_rho"],
            "physical_azimuth_labels": ["0", "pi/2", "pi", "3pi/2"],
            "formal_F0_products_retained": True,
            "physical_field_materialized": False,
        },
    )

    proof = {
        "kind": "native_O3_exact_anchor_axial_Q_velocity_proof_v1",
        "proof_count": 1,
        "local_actual_fixed_point_anchor": True,
        "same_source_complete_pair": True,
        "ordinary_Z_frame": "fixed rho before anchor restriction",
        "anchor_is_exact_negative_H_root": True,
        "anchor_midpoint_selected": False,
        "CZZ_seed_uses_pressure_p3": True,
        "p3_source": "core_seed_scale*Pnorm3; pressure P0_ZZZ datum, not P1 or m15",
        "weighted_norm_derivative_bound": (
            "EzBound=upper(fixed_error_norm)/(4*h.lower*(1-4.1/20))"
        ),
        "axis_derivative_data": "Phi_Z(0)=0 and U_Z(0)=4 exactly",
        "Q_source_definition": "[2*a*U-(1-delta)*a*M-d*M_Z]/L at anchor only",
        "Q_axis_limit": "Q_source(0)=Q0; it is not forced to zero",
        "pressure_gradient": "I1_Z=2*V(Phi1*Phi1_Z), with new finite convolution",
        "physical_map": "tau>0, lambda=sqrt(tau/d), kappa=lambda^(1-delta)",
        "streamfunction_divergence": "conditional mixed-partial certificate only",
        "absolute_F0_exp_materialized": False,
        "N_materialized": False,
        "generic_exp_calls_counted": False,
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
            "local_ordinary_Z_jet_enclosed": True,
            "local_Q_source_enclosed": True,
            "local_velocity_coefficients_enclosed": True,
            "global_actual_source_jet_admission": False,
            "physical_materialized": False,
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
        for name in ("delta", "epsilon_core", "j", "sigma", "h")
    }
    accepted_seed_records = {
        name: _interval_record(context, seed_values[name]) for name in raw_seed
    }
    pressure_seed_records = {
        name: _interval_record(context, pressure_seed[name])
        for name in ("B1", "core_seed_scale", "qP", "I2", "I3", "Pnorm3", "p3", "H2", "A_ZZ", "N_ZZ", "C_ZZ")
    }
    axial_polynomial = {
        "accepted_parent_polynomials": _json_value(parent_polynomials),
        "Phi1_Z_coefficients": [
            {"degree": index + 1, "interval": _interval_record(context, value)}
            for index, value in enumerate(differentiated["Phi1_Z_coefficients"])
        ],
        "DeltaU1_Z_coefficients": [
            {"degree": index + 1, "interval": _interval_record(context, value)}
            for index, value in enumerate(differentiated["DeltaU1_Z_coefficients"])
        ],
        "source_scale_B": _interval_record(context, differentiated["B"]),
        "derived": {
            name: _interval_record(context, differentiated[name])
            for name in (
                "L_Z",
                "beta_Z",
                "K_Z",
                "chi_ZZ",
                "Qtheta",
                "Qhat",
                "D",
                "Dstar",
                "S_ZZ_over_S",
            )
        },
        "formal_definitions": _json_value(differentiated["definitions"]),
        "F_Z_definition": "F_Z=F0*Phi_Z at anchor and F0_Z=0; Phi_Z coefficient retained",
        "parent_base_reused": True,
    }
    axial_rows_packet = {
        "rows": axial_rows_payload,
        "error_norm": _interval_record(context, fixed_error_norm),
        "E0": _interval_record(context, E0),
        "h": _interval_record(context, h),
        "EzBound": _interval_record(context, EzBound),
        "embedding_formula": "upper(E_norm)/(4*h.lower*(1-4.1/20))",
        "axis_exact": {"Phi_Z": True, "U_Z": 4, "DeltaU_Z": 0},
        "parent_rows_copied_without_recomputation": True,
    }
    average_packet = {
        "rows": average_rows,
        "DeltaU1_parent_coefficients": [
            {"degree": index + 1, "interval": _interval_record(context, value)}
            for index, value in enumerate(parent_delta_u_coefficients)
        ],
        "I1_Z_coefficients": [
            {"degree": degree, "interval": _interval_record(context, value)}
            for degree, value in sorted(I_Z_coefficients.items())
        ],
        "uniform_abs_Phi1_Z": _interval_record(context, uniform_abs_phi_z),
        "uniform_abs_Phi1": _json_value(
            _lookup(parent_polynomials_for_calc, ("uniform_abs_Phi1", "uniform_abs_Phi1_bound"))
        ),
        "ordinary_Z_definition": "fixed rho source derivative; anchor values only",
        "M_definition": "U0+(1/rho)*integral_0^rho DeltaU(s,Z) ds",
        "M_Z_definition": "4+(1/rho)*integral_0^rho DeltaU_Z(s,Z) ds",
        "Q_source_definition": "[2*a*U-(1-delta)*a*M-d*M_Z]/L",
        "Q_source_is_not_loop_q": True,
        "Q_source_is_not_pressure_q": True,
        "pressure_gradient_definition": "I1_Z=2*V(Phi1*Phi1_Z)",
        "source_S_Z_at_anchor": "0",
        "parent_Phi1_reused": True,
    }
    source_scale_packet = {
        "formalS": str(
            _lookup(
                source_scale_source,
                (
                    "formalS",
                    "formal_definition",
                    "formal_definition_S",
                    "S_definition",
                ),
            )
            or "S=exp(-2*logC)"
        ),
        "materialized": False,
        "S_bound": _interval_record(
            context,
            _record_or_interval(
                context, source_scale_source, ("S_bound",), "source_scale.S_bound"
            ),
        ),
        "B_bound": _interval_record(context, differentiated["B"]),
        "copied_parent_numerical_enclosure": True,
        "epsilon_power_eight_recomputed": False,
        "absolute_F0_exp_materialized": False,
    }
    family_identity = _lookup(
        capsule, ("family_identity", "family_identity_hashes", "family_identity3hashes")
    )
    preservation = _json_value(
        _lookup(capsule, ("preservation_guards", "guards", "preservation"))
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
        "seed_inputs": {
            "accepted_records": accepted_seed_records,
            "raw": raw_seed,
            "pressure_third_derivative": pressure_seed_records,
            "h_is_shared_tube_Xh_parameter": True,
        },
        "pressure_measure": {
            "B1": _interval_record(context, B1),
            "core_seed_scale": _interval_record(context, core_seed_scale),
            "raw_B1": _json_value(raw_B1),
            "raw_core_seed_scale": _json_value(raw_core_seed),
            "source_function": str(
                _lookup(pressure_measure, ("source_function", "source_definition", "definition"))
                or ""
            ),
            "beta_guard": _json_value(
                _lookup(
                    pressure_measure,
                    ("beta_guard", "beta_range", "beta_bounds", "beta", "beta_support"),
                )
            ),
            "pressure_P0_ZZZ_not_P1_or_m15": True,
        },
        "parent_polynomials": _json_value(parent_polynomials),
        "seed_inputs_derived_CZZ_p3": pressure_seed_records,
        "axial_polynomial": axial_polynomial,
        "axial_rows": axial_rows_packet,
        "average_Q_pressure_gradient": average_packet,
        "physical_velocity_coefficients": physical_ports,
        "source_scale": source_scale_packet,
        "family_identity": _json_value(family_identity),
        "grid": _json_value(grid),
        "preservation_guards": preservation,
        "definition_guards": _json_value(
            _lookup(capsule, ("definition_guards", "definition", "definitions"))
        ),
        "source_recipe_proof_metadata": _json_value(
            _lookup(
                capsule,
                (
                    "source_recipe",
                    "source_recipe_proof_metadata",
                    "source_recipe_literal_identity",
                ),
            )
        ),
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
            "local_ordinary_Z_jet": True,
            "local_Q_source": True,
            "local_velocity_coefficients": True,
            "local_ordinary_Z_jet_enclosed": True,
            "local_Q_source_enclosed": True,
            "local_velocity_coefficients_enclosed": True,
            "global_actual_source_jet_admission": False,
            "global_source_jet_admission": False,
            "physical_materialized": False,
            "physical_field_materialized": False,
            "global_five_moment": False,
            "time_n_recursion": False,
            "full_goal": False,
            "absolute_F0_exp_materialized": False,
            "N_materialized": False,
        },
        "static_limitations": [
            "exact negative anchor-local ordinary Z jet only; no global source jet",
            "Q_source and M_Z are admitted only on the anchor/rho grid",
            "physical outputs are typed coefficient ports with formal F0 products",
            "no absolute F0 or N materialization",
            "no global five-moment closure, time-n recursion, Cartesian NS field, winding, energy, stress, or force claim",
            "physical streamfunction divergence is a conditional structural certificate only",
            "generic exp calls are not counted because interval positive powers may use internal exp",
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
        "local_ordinary_Z_jet": True,
        "local_Q_source": True,
        "local_velocity_coefficients": True,
        "local_ordinary_Z_jet_enclosed": True,
        "local_Q_source_enclosed": True,
        "local_velocity_coefficients_enclosed": True,
        "global_actual_source_jet_admission": False,
        "global_source_jet_admission": False,
        "physical_materialized": False,
        "physical_field_materialized": False,
        "global_five_moment": False,
        "time_n_recursion": False,
        "full_goal": False,
        "full_goal_complete": False,
        "payload": payload,
    }
    _, return_hash, return_bytes = _serialize_result(result_out)
    result_out["return_code_hash"] = return_hash
    result_out["return_code_bytes"] = return_bytes
    result_out["code_hash"] = return_hash
    result_out["code_bytes"] = return_bytes
    result_out["serialization_hash_scope"] = "canonical result before hash metadata"
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





