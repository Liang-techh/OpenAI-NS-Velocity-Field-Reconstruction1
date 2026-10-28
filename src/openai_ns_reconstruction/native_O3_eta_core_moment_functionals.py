"""Local eta-core moment functional leaf.

This module consumes one root-frozen, same-source capsule.  It constructs
interval-valued polynomial functionals from the accepted eta chart; it does
not call an ancestor evaluator and it does not materialize a physical field.
The temporal n=1 port is deliberately an uninstantiated affine operator.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "native_O3_eta_core_moment_functionals_kernel_result_v1"
CAPSULE_SCHEMA = "native_O3_eta_core_moment_functionals_source_capsule_v1"
COMPONENT = "native_O3_eta_core_moment_functionals"
INTERVAL_DPS = 280
SOURCE_PRECISION_BITS = 867

# The root replaces these receipt placeholders before its one derive call.
ROOT_CAPSULE_RAW_SHA = "368b6be76051ddbe0d5b9d8dac1c5b8925be91b457d3ef2ea5affb387c81e422"
ROOT_FREEZE_SHA = "ff33f03893a9eded28b8d5703bfa9c119c92adbb1d322aa1fcc512903e4f975e"
SOURCE_IDS = {
    "eta_capsule_sha256": "e7342d77112d868606d2aaf6e8f88a82d19dc518c4d57d364cc1c7e0d86f9d7e",
    "eta_result_sha256": "e77879ef7c7175ca2d042719ba875c3125c508c952c08366f98e7c16b86c6f76",
    "native_inventory_sha256": "7566397ab21d796c15e0337a45912e38c92dde5fb31944c32fae545220dba121",
    "native_operator_trace_sha256": "a71054a08df5799f504b138fe8b57c56b3d87dfa821878445758b14df25521fe",
}

ETA_POINT_LABELS = ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1")
ETA_LABELS = ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1", "uniform_eta")
PUBLISH_ETA_LABELS = ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1", "uniform_eta")
RHO_LABELS = ("4",)
UNIFORM_RHO_LABELS = ("uniform_rho",)
PUBLISH_RHO_LABELS = ("4", "uniform_rho")
MOMENT_ORDER = ("m11", "m12", "m13", "m14", "m15")
MOMENT_NAMES = {
    "m11": "Mz",
    "m12": "Mtheta",
    "m13": "Mp",
    "m14": "Mtheta_z",
    "m15": "Mz_theta",
}
BASIS = (0, 1, 2, 3)

FAMILY_NAMES = (
    "coefficient_polynomials",
    "moment_integrals",
    "functional_errors",
    "physical_moment_scales",
    "n1_m15_affine_binding",
)
CALLBACK_EVENTS = (
    "family_entered",
    "coefficient_polynomials_entered",
    "coefficient_polynomials_complete",
    "moment_integrals_entered",
    "moment_integrals_complete",
    "functional_errors_entered",
    "functional_errors_complete",
    "physical_moment_scales_entered",
    "physical_moment_scales_complete",
    "n1_m15_affine_binding_entered",
    "n1_m15_affine_binding_complete",
    "proof_bundle",
    "family_completed",
)
CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "coefficient_polynomials_entered": 1,
    "coefficient_polynomials_complete": 1,
    "moment_integrals_entered": 1,
    "moment_integrals_complete": 1,
    "functional_errors_entered": 1,
    "functional_errors_complete": 1,
    "physical_moment_scales_entered": 1,
    "physical_moment_scales_complete": 1,
    "n1_m15_affine_binding_entered": 1,
    "n1_m15_affine_binding_complete": 1,
    "proof_bundle": 1,
    "family_completed": 1,
}
MANAGER_JOURNAL_EVENT_COUNT = 15

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
    "coefficient_polynomials_family_calls": 1,
    "moment_integrals_family_calls": 1,
    "functional_errors_family_calls": 1,
    "physical_moment_scales_family_calls": 1,
    "n1_m15_affine_binding_family_calls": 1,
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
    "leading_moment_function_enclosures": True,
    "n1_m15_affine_operator": True,
    "n1_m15_values_computed": False,
    "n1_solution_materialized": False,
    "global_moments": False,
    "global_five_moment": False,
    "time_n_recursion": False,
    "full_physical": False,
    "physical_materialized": False,
    "global_NS": False,
    "full_goal": False,
    "full_goal_complete": False,
    "physical_field_materialized": False,
    "absolute_F0_exp_materialized": False,
    "N_materialized": False,
    "old_math_reentered": False,
    "old_provider_called": False,
    "provider_called": False,
    "checker_called": False,
    "selftest_called": False,
    "source_value_point_selected": False,
    "source_parameter_point_selected": False,
    "point_background_from_midpoints": False,
    "axis_rho0_exact": True,
    "physical_rho_max_is_four": True,
    "analytic_extension_4_1_used": False,
    "formal_F0_retained": True,
    "formal_S_retained": True,
    "formal_Pstar_squared_retained": True,
    "all_coefficients_source_original_R": True,
}
COUNT_EVENT_MAP = {
    "source_family_calls": ["family_entered"],
    "new_math_family_calls": [
        "coefficient_polynomials_entered",
        "moment_integrals_entered",
        "functional_errors_entered",
        "physical_moment_scales_entered",
        "n1_m15_affine_binding_entered",
    ],
    "coefficient_polynomials_family_calls": ["coefficient_polynomials_entered"],
    "moment_integrals_family_calls": ["moment_integrals_entered"],
    "functional_errors_family_calls": ["functional_errors_entered"],
    "physical_moment_scales_family_calls": ["physical_moment_scales_entered"],
    "n1_m15_affine_binding_family_calls": ["n1_m15_affine_binding_entered"],
    "derive_calls": ["family_entered"],
    "numeric_evaluation_calls": ["family_entered"],
    "new_mathematical_proofs": ["proof_bundle"],
}
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
    "coefficient_polynomials_family_calls": 1,
    "moment_integrals_family_calls": 1,
    "functional_errors_family_calls": 1,
    "physical_moment_scales_family_calls": 1,
    "n1_m15_affine_binding_family_calls": 1,
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
    "leading_moment_function_enclosures": True,
    "n1_m15_affine_operator": True,
    "n1_m15_values_computed": False,
    "n1_solution_materialized": False,
    "global_moments": False,
    "global_five_moment": False,
    "time_n_recursion": False,
    "full_physical": False,
    "physical_materialized": False,
    "global_NS": False,
    "full_goal": False,
    "full_goal_complete": False,
    "physical_field_materialized": False,
    "absolute_F0_exp_materialized": False,
    "N_materialized": False,
    "old_math_reentered": False,
    "old_provider_called": False,
    "provider_called": False,
    "checker_called": False,
    "selftest_called": False,
    "source_value_point_selected": False,
    "source_parameter_point_selected": False,
    "point_background_from_midpoints": False,
    "axis_rho0_exact": True,
    "physical_rho_max_is_four": True,
    "analytic_extension_4_1_used": False,
    "formal_F0_retained": True,
    "formal_S_retained": True,
    "formal_Pstar_squared_retained": True,
    "all_coefficients_source_original_R": True,
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
    """Raised when the root capsule is absent or fails its source contract."""


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
        # The callback contract is exactly two positional arguments.
        emit(name, row)


def _canonical_bytes(capsule: Mapping[str, Any]) -> bytes:
    return json.dumps(
        _json_value(capsule),
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


def _raw_mpf_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise KernelInputError(label + " must be a four-word MPF tuple")
    try:
        words = tuple(int(item) for item in value)
    except (TypeError, ValueError) as exc:
        raise KernelInputError(label + " contains a non-integer MPF word") from exc
    return words  # type: ignore[return-value]


def _interval_value(context: Any, value: Any, label: str) -> Any:
    """Read an interval only from its paired lower/upper MPF words."""

    if not _is_mapping(value):
        raise KernelInputError(label + " must be an exact interval record")
    if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
        lower = _raw_mpf_tuple(
            value["lower_exact_mpf_tuple"], label + ".lower_exact_mpf_tuple"
        )
        upper = _raw_mpf_tuple(
            value["upper_exact_mpf_tuple"], label + ".upper_exact_mpf_tuple"
        )
        return context.make_mpf((lower, upper))
    for nested in ("interval", "enclosure", "record", "value", "bound"):
        if nested in value:
            return _interval_value(context, value[nested], label + "." + nested)
    raise KernelInputError(
        label + " requires lower_exact_mpf_tuple and upper_exact_mpf_tuple"
    )


def _new_interval_context() -> Any:
    # Kept local so AST inspection of this leaf never imports project code.
    from mpmath.ctx_iv import MPIntervalContext

    context = MPIntervalContext()
    context.dps = INTERVAL_DPS
    return context


def _endpoint(context: Any, value: Any, upper: bool = False) -> Any:
    converted = context.convert(value)
    raw = converted._mpi_[1 if upper else 0]
    return context.make_mpf((raw, raw))


def _lower(context: Any, value: Any) -> Any:
    return _endpoint(context, value, False)


def _upper(context: Any, value: Any) -> Any:
    return _endpoint(context, value, True)


def _make_interval(context: Any, lower: Any, upper: Any) -> Any:
    lo = context.convert(lower)
    hi = context.convert(upper)
    return context.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _zero(context: Any) -> Any:
    z = context.mpf("0")
    return _make_interval(context, z, z)


def _one(context: Any) -> Any:
    o = context.mpf("1")
    return _make_interval(context, o, o)


def _exact_text(context: Any, text: str) -> Any:
    if "/" in text:
        numerator, denominator = text.split("/", 1)
        return context.mpf(int(numerator)) / context.mpf(int(denominator))
    return context.mpf(text)


def _abs_upper(context: Any, value: Any) -> Any:
    lo = _lower(context, value)
    hi = _upper(context, value)
    zero = context.mpf("0")
    if bool(lo >= zero):
        return _make_interval(context, hi, hi)
    if bool(hi <= zero):
        return _make_interval(context, -lo, -lo)
    return _make_interval(context, zero, -lo if bool(-lo >= hi) else hi)


def _widen(context: Any, value: Any, error: Any) -> Any:
    return _make_interval(
        context,
        _lower(context, value) - _upper(context, error),
        _upper(context, value) + _upper(context, error),
    )


def _interval_record(context: Any, value: Any) -> dict[str, Any]:
    interval = context.convert(value)
    lower_raw, upper_raw = interval._mpi_
    lower = context.make_mpf((lower_raw, lower_raw))
    upper = context.make_mpf((upper_raw, upper_raw))
    text = lambda item: str(context.nstr(context.convert(item), n=INTERVAL_DPS)).strip()
    return {
        "lower": text(lower),
        "upper": text(upper),
        "width": text(upper - lower),
        "lower_exact_mpf_tuple": list(lower_raw),
        "upper_exact_mpf_tuple": list(upper_raw),
        "point_value": False,
    }


def _poly_records(context: Any, polynomial: Mapping[int, Any]) -> list[dict[str, Any]]:
    return [
        {"degree": degree, "interval": _interval_record(context, polynomial[degree])}
        for degree in sorted(polynomial)
    ]


def _poly_add(context: Any, *polynomials: Mapping[int, Any]) -> dict[int, Any]:
    result: dict[int, Any] = {}
    z = _zero(context)
    for polynomial in polynomials:
        for degree, coefficient in polynomial.items():
            result[degree] = result.get(degree, z) + coefficient
    return {degree: result[degree] for degree in sorted(result)}


def _poly_scale(
    context: Any, polynomial: Mapping[int, Any], factor: Any
) -> dict[int, Any]:
    return {degree: factor * coefficient for degree, coefficient in polynomial.items()}


def _poly_mul(
    context: Any, left: Mapping[int, Any], right: Mapping[int, Any]
) -> dict[int, Any]:
    result: dict[int, Any] = {}
    z = _zero(context)
    for left_degree, left_value in left.items():
        for right_degree, right_value in right.items():
            degree = left_degree + right_degree
            result[degree] = result.get(degree, z) + left_value * right_value
    return {degree: result[degree] for degree in sorted(result)}


def _poly_integral(
    context: Any, polynomial: Mapping[int, Any], weight: str = "one"
) -> dict[int, Any]:
    """Integrate a polynomial, with optional source-radial weight ``s``."""

    offset = 1 if weight == "one" else 2
    return {
        degree + offset: coefficient / (degree + offset)
        for degree, coefficient in polynomial.items()
    }


def _poly_eval(context: Any, polynomial: Mapping[int, Any], rho: Any) -> Any:
    value = _zero(context)
    for degree, coefficient in polynomial.items():
        value = value + coefficient * rho**degree
    return value


def _poly_abs_bound(context: Any, polynomial: Mapping[int, Any]) -> Any:
    bound = _zero(context)
    radius = _exact_text(context, "4")
    for degree, coefficient in polynomial.items():
        bound = bound + _abs_upper(context, coefficient) * radius**degree
    return _make_interval(context, _upper(context, bound), _upper(context, bound))


def _poly_enclose(
    context: Any, candidate: Mapping[int, Any], error: Mapping[int, Any]
) -> dict[int, Any]:
    degrees = sorted(set(candidate) | set(error))
    z = _zero(context)
    return {
        degree: _widen(
            context, candidate.get(degree, z), error.get(degree, z)
        )
        for degree in degrees
    }


def _poly_gap(
    context: Any,
    candidate: Mapping[int, Any],
    enclosed: Mapping[int, Any],
) -> dict[int, Any]:
    """Conservative coefficientwise error retained beside an enclosure."""

    degrees = sorted(set(candidate) | set(enclosed))
    z = _zero(context)
    return {
        degree: _abs_upper(
            context,
            enclosed.get(degree, z) - candidate.get(degree, z),
        )
        for degree in degrees
    }


def _coefficient_values(
    context: Any,
    raw: Any,
    degrees: Sequence[int],
    label: str,
) -> dict[int, Any]:
    if _is_mapping(raw):
        nested = _lookup(raw, ("coefficients", "values", "records", "items"))
        if nested is not None:
            raw = nested
        else:
            selected: list[Any] = []
            for degree in degrees:
                item = None
                for key in (degree, str(degree)):
                    if key in raw:
                        item = raw[key]
                        break
                if item is None:
                    raise KernelInputError(label + " is missing degree " + str(degree))
                selected.append(item)
            raw = selected
    if not _is_sequence(raw) or len(raw) != len(degrees):
        raise KernelInputError(label + " has the wrong coefficient count")
    return {
        degree: _interval_value(context, item, label + "[" + str(degree) + "]")
        for degree, item in zip(degrees, raw)
    }


def _poly_packet(
    context: Any, polynomials: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, dict[int, Any]]]:
    specs = (
        ("Phi1", (0, 1, 2)),
        ("DU1", (1, 2, 3)),
        ("Phi1_Z", (1, 2, 3)),
        ("DU1_Z", (1, 2, 3)),
        ("Phi1_ZZ", (1, 2, 3)),
        ("DU1_ZZ", (1, 2, 3)),
    )
    raw: dict[str, Any] = {}
    values: dict[str, dict[int, Any]] = {}
    for name, degrees in specs:
        source = _required(polynomials, (name, name.replace("_", "")), "polynomials." + name)
        raw[name] = _json_value(source)
        values[name] = _coefficient_values(
            context, source, degrees, "polynomials." + name
        )
    return raw, values


def _row_mapping(eta_rows: Any) -> Mapping[str, Any]:
    if _is_mapping(eta_rows):
        nested = _lookup(eta_rows, ("rows", "table", "records"))
        return _as_mapping(nested if nested is not None else eta_rows, "capsule.eta_rows")
    if _is_sequence(eta_rows):
        rows: dict[str, Any] = {}
        for item in eta_rows:
            row = _as_mapping(item, "capsule.eta_rows[]")
            eta = _required(row, ("eta", "label"), "eta_rows[].eta")
            rows[str(eta)] = row
        return rows
    raise KernelInputError("capsule.eta_rows must be a mapping or sequence")


def _source_interval(context: Any, source: Mapping[str, Any], names: Sequence[str], label: str) -> tuple[Any, Any]:
    raw = _required(source, names, label)
    return raw, _interval_value(context, raw, label)


def _row_error(
    context: Any,
    row: Mapping[str, Any],
    names: Sequence[str],
    label: str,
    epsilon: Any,
    fallback_names: Sequence[str],
) -> tuple[Any, Any]:
    raw = _required(row, names, label)
    try:
        return raw, _interval_value(context, raw, label)
    except KernelInputError:
        # The eta=0 source may deliberately retain epsilon*oldE0/oldEZ as a
        # descriptor instead of duplicating a scaled interval record.
        descriptor = _as_mapping(raw, label)
        multiply = descriptor.get("multiply_epsilon_once") is True
        if not multiply and "epsilon" not in str(descriptor.get("expression", "")).lower():
            raise
        source_record = descriptor.get("source_record")
        if source_record is None:
            source_record = _lookup(row, fallback_names)
        if source_record is None:
            raise KernelInputError(label + " descriptor has no source_record")
        base = _interval_value(context, source_record, label + ".source_record")
        return raw, epsilon * base


def _validate_capsule(capsule: Any) -> Mapping[str, Any]:
    capsule = _as_mapping(capsule, "data.capsule")
    if capsule.get("schema") != CAPSULE_SCHEMA:
        raise KernelInputError("unexpected eta moment source capsule schema")
    if capsule.get("component") != COMPONENT:
        raise KernelInputError("source capsule component does not match")
    if ROOT_FREEZE_SHA in ("__ROOT_FREEZE_SHA_REQUIRED__", "0" * 64):
        raise KernelInputError("ROOT_FREEZE_SHA is not frozen")
    actual = hashlib.sha256(_canonical_bytes(capsule)).hexdigest()
    if actual.lower() != ROOT_FREEZE_SHA.lower():
        raise KernelInputError("source capsule canonical digest mismatch")
    source_ids = _as_mapping(capsule.get("source_ids"), "capsule.source_ids")
    if not source_ids:
        raise KernelInputError("source_ids must not be empty")
    for name, expected in SOURCE_IDS.items():
        if source_ids.get(name) != expected:
            raise KernelInputError("source id mismatch: " + name)
    parameters = _as_mapping(capsule.get("parameters"), "capsule.parameters")
    for name in ("epsilon_core", "a", "s", "U0"):
        if parameters.get(name) is None:
            raise KernelInputError("parameters is missing " + name)
    polynomials = _as_mapping(capsule.get("polynomials"), "capsule.polynomials")
    for name in ("Phi1", "DU1", "Phi1_Z", "DU1_Z", "Phi1_ZZ", "DU1_ZZ"):
        if polynomials.get(name) is None and polynomials.get(name.replace("_", "")) is None:
            raise KernelInputError("polynomials is missing " + name)
    rows = _row_mapping(_required(capsule, ("eta_rows",), "capsule.eta_rows"))
    for eta in ETA_LABELS:
        if eta not in rows:
            raise KernelInputError("eta_rows is missing " + eta)
    _as_mapping(capsule.get("source_atoms"), "capsule.source_atoms")
    _as_mapping(capsule.get("retained_native_packet"), "capsule.retained_native_packet")
    grid = _as_mapping(capsule.get("grid"), "capsule.grid")
    for name in ("rho", "uniform_rho", "eta"):
        if grid.get(name) is None:
            raise KernelInputError("grid is missing " + name)
    grid_text = json.dumps(_json_value(grid), ensure_ascii=False)
    extension_marker = "4" + "." + "1"
    if extension_marker in grid_text or "4,1" in grid_text:
        raise KernelInputError("analytic extension beyond rho=4 is not a physical splice")
    return capsule


def _row_atoms(
    context: Any,
    row: Mapping[str, Any],
    epsilon: Any,
) -> tuple[dict[str, Any], dict[str, Any]]:
    raw: dict[str, Any] = {}
    values: dict[str, Any] = {}
    for name, aliases in (
        ("shift", ("shift", "eta_shift")),
        ("Ubase", ("Ubase", "U_base", "u_base")),
        ("ratio", ("ratio", "relative_profile", "F0_ratio")),
        ("logF0_Z", ("logF0_Z", "logF0Z", "ordinary_logF0_Z")),
    ):
        raw[name], values[name] = _source_interval(context, row, aliases, "eta_row." + name)
    raw["Ephi"], values["Ephi"] = _row_error(
        context, row, ("Ephi", "error_phi", "savedE0", "E0"),
        "eta_row.Ephi", epsilon, ("savedE0", "E0", "Ephi"),
    )
    raw["EphiZ"], values["EphiZ"] = _row_error(
        context, row, ("EphiZ", "error_phi_Z", "savedEZ", "EZ"),
        "eta_row.EphiZ", epsilon, ("savedEZ", "EZ", "EphiZ"),
    )
    raw["EDU"], values["EDU"] = _row_error(
        context, row, ("EDU", "error_DU", "epsilon_E0"),
        "eta_row.EDU", epsilon, ("Ephi", "savedE0", "E0"),
    )
    raw["EDUZ"], values["EDUZ"] = _row_error(
        context, row, ("EDUZ", "error_DU_Z", "epsilon_EZ"),
        "eta_row.EDUZ", epsilon, ("EphiZ", "savedEZ", "EZ"),
    )
    return raw, values


def _candidate_polynomials(
    context: Any,
    shift: Any,
    base: Mapping[str, Mapping[int, Any]],
) -> dict[str, dict[int, Any]]:
    half = _exact_text(context, "1/2")
    shift2 = shift * shift
    phi = _poly_add(
        context,
        base["Phi1"],
        _poly_scale(context, base["Phi1_Z"], shift),
        _poly_scale(context, base["Phi1_ZZ"], half * shift2),
    )
    du = _poly_add(
        context,
        base["DU1"],
        _poly_scale(context, base["DU1_Z"], shift),
        _poly_scale(context, base["DU1_ZZ"], half * shift2),
    )
    phi_z = _poly_add(
        context,
        base["Phi1_Z"],
        _poly_scale(context, base["Phi1_ZZ"], shift),
    )
    du_z = _poly_add(
        context,
        base["DU1_Z"],
        _poly_scale(context, base["DU1_ZZ"], shift),
    )
    # Keep phi-1 as a separately tagged polynomial without numerically
    # subtracting the degree-zero identity.  The source Phi1 degree-zero
    # coefficient is the exact unit atom; only degrees >= 1 are retained.
    phi_minus_one = {
        degree: coefficient for degree, coefficient in phi.items() if degree >= 1
    }
    return {
        "phi": phi,
        "phi_minus_one": {degree: phi_minus_one[degree] for degree in sorted(phi_minus_one)},
        "DU": du,
        "phi_Z": phi_z,
        "DU_Z": du_z,
    }


def _candidate_errors(errors: Mapping[str, Any]) -> dict[str, Any]:
    # The capsule's four records are already complete certified eta-value
    # errors.  They are consumed directly; no chart-tail shift is recomputed.
    return {
        "phi": errors["Ephi"],
        "phi_Z": errors["EphiZ"],
        "DU": errors["EDU"],
        "DU_Z": errors["EDUZ"],
    }


def _error_poly(context: Any, error: Any, weight: str) -> dict[int, Any]:
    if weight == "one":
        return {1: error}
    return {2: error / 2}


def _single_integral(
    context: Any,
    polynomial: Mapping[int, Any],
    derivative: Mapping[int, Any],
    error: Any,
    derivative_error: Any,
    weight: str,
) -> dict[str, Any]:
    candidate = _poly_integral(context, polynomial, weight)
    candidate_z = _poly_integral(context, derivative, weight)
    error_poly = _error_poly(context, error, weight)
    error_poly_z = _error_poly(context, derivative_error, weight)
    return {
        "candidate": candidate,
        "candidate_Z": candidate_z,
        "error": error_poly,
        "error_Z": error_poly_z,
        "enclosed": _poly_enclose(context, candidate, error_poly),
        "enclosed_Z": _poly_enclose(context, candidate_z, error_poly_z),
    }


def _product_error_bound(
    context: Any,
    left: Mapping[int, Any],
    right: Mapping[int, Any],
    left_error: Any,
    right_error: Any,
) -> Any:
    left_abs = _upper(context, _poly_abs_bound(context, left))
    right_abs = _upper(context, _poly_abs_bound(context, right))
    value = _upper(context, left_error) * right_abs
    value = value + _upper(context, right_error) * left_abs
    value = value + _upper(context, left_error) * _upper(context, right_error)
    return _make_interval(context, context.mpf("0"), value)


def _derivative_product_error_bound(
    context: Any,
    left: Mapping[int, Any],
    right: Mapping[int, Any],
    left_z: Mapping[int, Any],
    right_z: Mapping[int, Any],
    left_error: Any,
    right_error: Any,
    left_z_error: Any,
    right_z_error: Any,
) -> Any:
    first = _product_error_bound(context, left_z, right, left_z_error, right_error)
    second = _product_error_bound(context, left, right_z, left_error, right_z_error)
    return _make_interval(
        context,
        context.mpf("0"),
        _upper(context, first) + _upper(context, second),
    )


def _square_error_bound(
    context: Any, polynomial: Mapping[int, Any], error: Any
) -> Any:
    abs_value = _upper(context, _poly_abs_bound(context, polynomial))
    upper = _upper(context, error) * (2 * abs_value + _upper(context, error))
    return _make_interval(context, context.mpf("0"), upper)


def _derivative_square_error_bound(
    context: Any,
    polynomial: Mapping[int, Any],
    derivative: Mapping[int, Any],
    error: Any,
    derivative_error: Any,
) -> Any:
    abs_value = _upper(context, _poly_abs_bound(context, polynomial))
    abs_derivative = _upper(context, _poly_abs_bound(context, derivative))
    upper = 2 * (
        _upper(context, error) * abs_derivative
        + _upper(context, derivative_error) * abs_value
        + _upper(context, error) * _upper(context, derivative_error)
    )
    return _make_interval(context, context.mpf("0"), upper)


def _product_integral(
    context: Any,
    left: Mapping[int, Any],
    right: Mapping[int, Any],
    left_z: Mapping[int, Any],
    right_z: Mapping[int, Any],
    left_error: Any,
    right_error: Any,
    left_z_error: Any,
    right_z_error: Any,
    weight: str,
) -> dict[str, Any]:
    candidate = _poly_integral(context, _poly_mul(context, left, right), weight)
    candidate_z = _poly_integral(
        context,
        _poly_add(
            context,
            _poly_mul(context, left_z, right),
            _poly_mul(context, left, right_z),
        ),
        weight,
    )
    error_poly = _error_poly(
        context,
        _product_error_bound(context, left, right, left_error, right_error),
        weight,
    )
    error_poly_z = _error_poly(
        context,
        _derivative_product_error_bound(
            context,
            left,
            right,
            left_z,
            right_z,
            left_error,
            right_error,
            left_z_error,
            right_z_error,
        ),
        weight,
    )
    return {
        "candidate": candidate,
        "candidate_Z": candidate_z,
        "error": error_poly,
        "error_Z": error_poly_z,
        "enclosed": _poly_enclose(context, candidate, error_poly),
        "enclosed_Z": _poly_enclose(context, candidate_z, error_poly_z),
    }


def _square_integral(
    context: Any,
    polynomial: Mapping[int, Any],
    derivative: Mapping[int, Any],
    error: Any,
    derivative_error: Any,
    weight: str,
) -> dict[str, Any]:
    candidate = _poly_integral(context, _poly_mul(context, polynomial, polynomial), weight)
    candidate_z = _poly_integral(
        context,
        _poly_scale(
            context, _poly_mul(context, polynomial, derivative), context.mpf("2")
        ),
        weight,
    )
    error_poly = _error_poly(
        context, _square_error_bound(context, polynomial, error), weight
    )
    error_poly_z = _error_poly(
        context,
        _derivative_square_error_bound(
            context, polynomial, derivative, error, derivative_error
        ),
        weight,
    )
    return {
        "candidate": candidate,
        "candidate_Z": candidate_z,
        "error": error_poly,
        "error_Z": error_poly_z,
        "enclosed": _poly_enclose(context, candidate, error_poly),
        "enclosed_Z": _poly_enclose(context, candidate_z, error_poly_z),
    }


def _moment_integrals(
    context: Any,
    candidate: Mapping[str, Mapping[int, Any]],
    errors: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    phi, phi_z = candidate["phi"], candidate["phi_Z"]
    delta, delta_z = candidate["phi_minus_one"], candidate["phi_Z"]
    du, du_z = candidate["DU"], candidate["DU_Z"]
    integrals = {
        "Jphir": _single_integral(
            context, phi, phi_z, errors["phi"], errors["phi_Z"], "r"
        ),
        "JD": _single_integral(
            context, du, du_z, errors["DU"], errors["DU_Z"], "one"
        ),
        "Jphidr": _product_integral(
            context,
            phi,
            du,
            phi_z,
            du_z,
            errors["phi"],
            errors["DU"],
            errors["phi_Z"],
            errors["DU_Z"],
            "r",
        ),
        "JD2": _square_integral(
            context, du, du_z, errors["DU"], errors["DU_Z"], "one"
        ),
        "Jphi2r": _square_integral(
            context, phi, phi_z, errors["phi"], errors["phi_Z"], "r"
        ),
        "Jphi2": _square_integral(
            context, phi, phi_z, errors["phi"], errors["phi_Z"], "one"
        ),
    }
    # Tiny phi-1 sectors are retained independently.  The error caps remain
    # conservative copies of the accepted phi/phi_Z chart caps.
    integrals.update(
        {
            "Jdeltar": _single_integral(
                context, delta, delta_z, errors["phi"], errors["phi_Z"], "r"
            ),
            "Jdelta": _single_integral(
                context, delta, delta_z, errors["phi"], errors["phi_Z"], "one"
            ),
            "JdeltaDUr": _product_integral(
                context,
                delta,
                du,
                delta_z,
                du_z,
                errors["phi"],
                errors["DU"],
                errors["phi_Z"],
                errors["DU_Z"],
                "r",
            ),
            "Jdelta2r": _square_integral(
                context, delta, delta_z, errors["phi"], errors["phi_Z"], "r"
            ),
            "Jdelta2": _square_integral(
                context, delta, delta_z, errors["phi"], errors["phi_Z"], "one"
            ),
            "JDUr": _single_integral(
                context, du, du_z, errors["DU"], errors["DU_Z"], "r"
            ),
        }
    )
    return integrals


def _atom_packet(context: Any, raw: Any, value: Any, expression: str) -> dict[str, Any]:
    packet = {
        "raw": _json_value(raw),
        "formal_positive": True,
        "materialized": False,
        "expression": expression,
    }
    if value is not None:
        packet["authenticated_interval"] = _interval_record(context, value)
    return packet


def _source_atoms(
    context: Any, source: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    specs = (
        ("F0", ("F0", "F0_anchor", "formal_F0", "F0positive"), "exp(-logC)"),
        ("S", ("S", "formal_S", "F0_squared", "F0^2"), "F0^2"),
        (
            "Pstar_squared",
            ("Pstar_squared", "Pstar_sq", "Pstar2", "P_star_squared"),
            "Pstar^2",
        ),
    )
    raw_packet: dict[str, Any] = {}
    packet: dict[str, Any] = {}
    for name, aliases, expression in specs:
        raw = _required(source, aliases, "source_atoms." + name)
        raw_packet[name] = raw
        # F0, S, and Pstar^2 are formal positive atoms.  The capsule may
        # intentionally carry only their definitions (and no interval), since
        # the atomic-sector coefficients below are stored before atom
        # multiplication.  If an authenticated interval is also present,
        # retain it without using it to materialize the atom.
        try:
            packet[name] = _interval_value(context, raw, "source_atoms." + name)
        except KernelInputError:
            packet[name] = None
    return raw_packet, packet


def _rho_values(context: Any) -> dict[str, Any]:
    return {label: _rho_interval(context, label) for label in RHO_LABELS + UNIFORM_RHO_LABELS}


def _rho_interval(context: Any, label: str) -> Any:
    if label == "uniform_rho":
        return _make_interval(
            context,
            _exact_text(context, "0"),
            _exact_text(context, "4"),
        )
    return _exact_text(context, label)


def _endpoint_value(
    context: Any,
    polynomial: Mapping[int, Any],
    rho_label: str,
    error_polynomial: Mapping[int, Any] | None = None,
    axis_zero: bool = True,
) -> tuple[Any, Any, Any]:
    rho = _rho_interval(context, rho_label)
    candidate = _poly_eval(context, polynomial, rho)
    error = _poly_eval(context, error_polynomial or {}, rho)
    if rho_label == "0" and axis_zero:
        exact = _zero(context)
        return exact, exact, exact
    if rho_label == "0":
        return candidate, _zero(context), candidate
    actual = _widen(context, candidate, error)
    return candidate, error, actual


def _function_packet(
    context: Any,
    candidate: Mapping[str, Mapping[int, Any]],
    errors: Mapping[str, Any],
) -> dict[str, Any]:
    domains: dict[str, Any] = {}
    # Only value fields have certified C0/C1 pairs here. phi_Z/DU_Z are
    # already their C1 ports; this leaf does not assert their Z2 jets.
    for name in ("phi", "phi_minus_one", "DU"):
        polynomial = candidate[name]
        derivative_polynomial = candidate.get(
            {"phi": "phi_Z", "phi_minus_one": "phi_Z", "DU": "DU_Z"}.get(
                name, name
            ),
            polynomial,
        )
        error_name = {"phi_minus_one": "phi"}.get(name, name)
        derivative_name = {"phi": "phi_Z", "phi_minus_one": "phi_Z", "DU": "DU_Z"}.get(
            name, name
        )
        error = errors.get(error_name, _zero(context))
        derivative_error = errors.get(derivative_name, error)
        bound = _poly_abs_bound(context, polynomial)
        samples: list[dict[str, Any]] = []
        error_interval = _make_interval(context, context.mpf("0"), _upper(context, error))
        derivative_error_interval = _make_interval(
            context, context.mpf("0"), _upper(context, derivative_error)
        )
        for rho_label in RHO_LABELS:
            c0, e0, actual = _endpoint_value(
                context, polynomial, rho_label, {0: error_interval}, axis_zero=False
            )
            c1, e1, actual1 = _endpoint_value(
                context, derivative_polynomial, rho_label, {0: derivative_error_interval}, axis_zero=False
            )
            samples.append(
                {
                    "rho": rho_label,
                    "candidate": _interval_record(context, c0),
                    "error": _interval_record(context, e0),
                    "actual": _interval_record(context, actual),
                    "candidate_Z": _interval_record(context, c1),
                    "error_Z": _interval_record(context, e1),
                    "actual_Z": _interval_record(context, actual1),
                }
            )
        uniform_samples: list[dict[str, Any]] = []
        for rho_label in UNIFORM_RHO_LABELS:
            c0, e0, actual = _endpoint_value(
                context, polynomial, rho_label, {0: error_interval}, axis_zero=False
            )
            c1, e1, actual1 = _endpoint_value(
                context, derivative_polynomial, rho_label, {0: derivative_error_interval}, axis_zero=False
            )
            uniform_samples.append(
                {
                    "rho": rho_label,
                    "candidate": _interval_record(context, c0),
                    "error": _interval_record(context, e0),
                    "actual": _interval_record(context, actual),
                    "candidate_Z": _interval_record(context, c1),
                    "error_Z": _interval_record(context, e1),
                    "actual_Z": _interval_record(context, actual1),
                }
            )
        domains[name] = {
            "candidate_polynomial": _poly_records(context, polynomial),
            "uniform_abs_candidate_bound": _interval_record(context, bound),
            "C0_error": _interval_record(context, error_interval),
            "C1_error": _interval_record(context, derivative_error_interval),
            "candidate_minus_plus_certified": True,
            "actual_function": "unique fixedpoint enclosed by candidate +/- certified error",
            "samples": samples,
            "uniform_rho_samples": uniform_samples,
        }
    return domains


def _moment_polynomials(
    context: Any,
    params: Mapping[str, Any],
    row: Mapping[str, Any],
    atoms: Mapping[str, Any],
    integrals: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, dict[str, Mapping[int, Any]]]]:
    eps = params["epsilon_core"]
    ubase = row["Ubase"]
    ratio = row["ratio"]
    logf = row["logF0_Z"]
    r = {1: _one(context)}
    r2 = {2: _one(context)}
    r2_half = {2: _exact_text(context, "1/2")}
    two = context.mpf("2")
    four = context.mpf("4")
    eight = context.mpf("8")
    j = {name: data["candidate"] for name, data in integrals.items()}
    jz = {name: data["candidate_Z"] for name, data in integrals.items()}
    je = {name: data["enclosed"] for name, data in integrals.items()}
    jez = {name: data["enclosed_Z"] for name, data in integrals.items()}

    def _sector_pair(
        c0: Mapping[int, Any], c1: Mapping[int, Any],
        e0: Mapping[int, Any], e1: Mapping[int, Any],
    ) -> dict[str, Mapping[int, Any]]:
        return {
            "candidate": c0,
            "candidate_Z": c1,
            "enclosed": e0,
            "enclosed_Z": e1,
            "error": _poly_gap(context, c0, e0),
            "error_Z": _poly_gap(context, c1, e1),
        }

    theta_scale = two * eps * eps * ratio
    m13_scale = eps * ratio * ratio

    # m11 = epsilon*(rho*Ubase + JD), with its axis and DU sectors kept
    # separate in coefficient space.
    m11_axis_c = _poly_scale(context, _poly_scale(context, r, ubase), eps)
    m11_axis_z_c = _poly_scale(context, _poly_scale(context, r, four), eps)
    m11_axis_e = m11_axis_c
    m11_axis_z_e = m11_axis_z_c
    m11_du_c = _poly_scale(context, j["JD"], eps)
    m11_du_z_c = _poly_scale(context, jz["JD"], eps)
    m11_du_e = _poly_scale(context, je["JD"], eps)
    m11_du_z_e = _poly_scale(context, jez["JD"], eps)

    # phi = 1 + delta_phi.  Keep the unit and tiny delta sectors separate so
    # interval rounding cannot erase the latter.
    m12_axis_c = _poly_scale(context, r2_half, theta_scale)
    m12_axis_z_c = _poly_scale(context, _poly_scale(context, r2_half, logf), theta_scale)
    m12_axis_e = m12_axis_c
    m12_axis_z_e = m12_axis_z_c
    m12_delta_c = _poly_scale(context, j["Jdeltar"], theta_scale)
    m12_delta_z_c = _poly_scale(context, _poly_add(context, _poly_scale(context, j["Jdeltar"], logf), jz["Jdeltar"]), theta_scale)
    m12_delta_e = _poly_scale(context, je["Jdeltar"], theta_scale)
    m12_delta_z_e = _poly_scale(context, _poly_add(context, _poly_scale(context, je["Jdeltar"], logf), jez["Jdeltar"]), theta_scale)

    m13_axis_c = _poly_scale(context, r, m13_scale)
    m13_axis_z_c = _poly_scale(context, _poly_scale(context, r, two * logf), m13_scale)
    m13_axis_e = m13_axis_c
    m13_axis_z_e = m13_axis_z_c
    m13_delta_c = _poly_scale(context, _poly_add(context, _poly_scale(context, j["Jdelta"], two), j["Jdelta2"]), m13_scale)
    m13_delta_z_c = _poly_scale(context, _poly_add(context, _poly_scale(context, _poly_add(context, _poly_scale(context, j["Jdelta"], two), j["Jdelta2"]), two * logf), _poly_add(context, _poly_scale(context, jz["Jdelta"], two), jz["Jdelta2"])), m13_scale)
    m13_delta_e = _poly_scale(context, _poly_add(context, _poly_scale(context, je["Jdelta"], two), je["Jdelta2"]), m13_scale)
    m13_delta_z_e = _poly_scale(context, _poly_add(context, _poly_scale(context, _poly_add(context, _poly_scale(context, je["Jdelta"], two), je["Jdelta2"]), two * logf), _poly_add(context, _poly_scale(context, jez["Jdelta"], two), jez["Jdelta2"])), m13_scale)

    m14_axis_c = _poly_scale(context, _poly_scale(context, r2_half, ubase), theta_scale)
    m14_axis_z_c = _poly_scale(context, _poly_add(context, _poly_scale(context, _poly_scale(context, r2_half, ubase), logf), _poly_scale(context, r2, two)), theta_scale)
    m14_axis_e = m14_axis_c
    m14_axis_z_e = m14_axis_z_c
    m14_corr_K = _poly_add(context, _poly_scale(context, j["Jdeltar"], ubase), j["JDUr"], j["JdeltaDUr"])
    m14_corr_Kz = _poly_add(context, _poly_scale(context, jz["Jdeltar"], ubase), jz["JDUr"], jz["JdeltaDUr"])
    m14_corr_Ke = _poly_add(context, _poly_scale(context, je["Jdeltar"], ubase), je["JDUr"], je["JdeltaDUr"])
    m14_corr_Kez = _poly_add(context, _poly_scale(context, jez["Jdeltar"], ubase), jez["JDUr"], jez["JdeltaDUr"])
    m14_corr_c = _poly_scale(context, m14_corr_K, theta_scale)
    m14_corr_z_c = _poly_scale(context, _poly_add(context, _poly_scale(context, m14_corr_K, logf), _poly_add(context, _poly_scale(context, j["Jdeltar"], four), m14_corr_Kz)), theta_scale)
    m14_corr_e = _poly_scale(context, m14_corr_Ke, theta_scale)
    m14_corr_z_e = _poly_scale(context, _poly_add(context, _poly_scale(context, m14_corr_Ke, logf), _poly_add(context, _poly_scale(context, je["Jdeltar"], four), m14_corr_Kez)), theta_scale)

    m15_axis_c = _poly_scale(context, _poly_scale(context, r, ubase * ubase), eps)
    m15_axis_z_c = _poly_scale(context, _poly_scale(context, r, eight * ubase), eps)
    m15_axis_e = m15_axis_c
    m15_axis_z_e = m15_axis_z_c
    m15_corr_c = _poly_scale(context, _poly_add(context, _poly_scale(context, j["JD"], two * ubase), j["JD2"]), eps)
    m15_corr_z_c = _poly_scale(context, _poly_add(context, _poly_add(context, _poly_scale(context, j["JD"], eight), _poly_scale(context, jz["JD"], two * ubase)), jz["JD2"]), eps)
    m15_corr_e = _poly_scale(context, _poly_add(context, _poly_scale(context, je["JD"], two * ubase), je["JD2"]), eps)
    m15_corr_z_e = _poly_scale(context, _poly_add(context, _poly_add(context, _poly_scale(context, je["JD"], eight), _poly_scale(context, jez["JD"], two * ubase)), jez["JD2"]), eps)
    swirl_scale = -eps * eps * ratio * ratio
    m15_swirl_axis_c = _poly_scale(context, r2_half, swirl_scale)
    m15_swirl_axis_z_c = _poly_scale(context, _poly_scale(context, r2_half, two * logf), swirl_scale)
    m15_swirl_axis_e = m15_swirl_axis_c
    m15_swirl_axis_z_e = m15_swirl_axis_z_c
    m15_swirl_corr_c = _poly_scale(context, _poly_add(context, _poly_scale(context, j["Jdeltar"], two), j["Jdelta2r"]), swirl_scale)
    m15_swirl_corr_z_c = _poly_scale(context, _poly_add(context, _poly_scale(context, _poly_add(context, _poly_scale(context, j["Jdeltar"], two), j["Jdelta2r"]), two * logf), _poly_add(context, _poly_scale(context, jz["Jdeltar"], two), jz["Jdelta2r"])), swirl_scale)
    m15_swirl_corr_e = _poly_scale(context, _poly_add(context, _poly_scale(context, je["Jdeltar"], two), je["Jdelta2r"]), swirl_scale)
    m15_swirl_corr_z_e = _poly_scale(context, _poly_add(context, _poly_scale(context, _poly_add(context, _poly_scale(context, je["Jdeltar"], two), je["Jdelta2r"]), two * logf), _poly_add(context, _poly_scale(context, jez["Jdeltar"], two), jez["Jdelta2r"])), swirl_scale)

    return {
        "m11": {
            "axis": _sector_pair(m11_axis_c, m11_axis_z_c, m11_axis_e, m11_axis_z_e),
            "DU": _sector_pair(m11_du_c, m11_du_z_c, m11_du_e, m11_du_z_e),
        },
        "m12": {
            "axis": _sector_pair(m12_axis_c, m12_axis_z_c, m12_axis_e, m12_axis_z_e),
            "phi_minus_one": _sector_pair(m12_delta_c, m12_delta_z_c, m12_delta_e, m12_delta_z_e),
        },
        "m13": {
            "axis": _sector_pair(m13_axis_c, m13_axis_z_c, m13_axis_e, m13_axis_z_e),
            "phi_minus_one": _sector_pair(m13_delta_c, m13_delta_z_c, m13_delta_e, m13_delta_z_e),
        },
        "m14": {
            "axis": _sector_pair(m14_axis_c, m14_axis_z_c, m14_axis_e, m14_axis_z_e),
            "correction": _sector_pair(m14_corr_c, m14_corr_z_c, m14_corr_e, m14_corr_z_e),
        },
        "m15": {
            "axis": _sector_pair(m15_axis_c, m15_axis_z_c, m15_axis_e, m15_axis_z_e),
            "correction": _sector_pair(m15_corr_c, m15_corr_z_c, m15_corr_e, m15_corr_z_e),
            "swirl": _sector_pair(m15_swirl_axis_c, m15_swirl_axis_z_c, m15_swirl_axis_e, m15_swirl_axis_z_e),
            "swirl_correction": _sector_pair(m15_swirl_corr_c, m15_swirl_corr_z_c, m15_swirl_corr_e, m15_swirl_corr_z_e),
        },
    }


def _value_error_record(
    context: Any, candidate: Any, actual: Any
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    difference = _abs_upper(context, actual - candidate)
    return (
        _interval_record(context, actual),
        _interval_record(context, candidate),
        _interval_record(context, difference),
    )


def _sector_row(
    context: Any,
    eta: str,
    rho_label: str,
    moment: str,
    sector: str,
    atom: str,
    pair: Mapping[str, Mapping[int, Any]],
) -> dict[str, Any]:
    c0, e0, actual0 = _endpoint_value(
        context, pair["candidate"], rho_label, pair.get("error", {})
    )
    c1, e1, actual1 = _endpoint_value(
        context, pair["candidate_Z"], rho_label, pair.get("error_Z", {})
    )
    actual0_record, candidate0_record, error0_record = _value_error_record(context, c0, actual0)
    actual1_record, candidate1_record, error1_record = _value_error_record(context, c1, actual1)
    return {
        "eta": eta,
        "rho": rho_label,
        "moment": moment,
        "moment_name": MOMENT_NAMES[moment],
        "temporal_order": 0,
        "is_n1": False,
        "uniform_eta_source_row": eta == "uniform_eta",
        "atomic_sectors": [
            {
                "atom": atom,
                "sector": sector,
                "C0": actual0_record,
                "C1": actual1_record,
                "C0_candidate": candidate0_record,
                "C1_candidate": candidate1_record,
                "C0_error": error0_record,
                "C1_error": error1_record,
            }
        ],
        "ordinary_Z": True,
        "full_field_value_materialized": False,
    }


def _moment_table(
    context: Any,
    rows: Mapping[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for eta in ETA_LABELS:
        for rho_label in PUBLISH_RHO_LABELS:
            moments = rows[eta]["moment_polynomials"]
            for moment in MOMENT_ORDER:
                sector_rows = []
                for sector, pair in moments[moment].items():
                    atom = "S" if sector.startswith("swirl") or moment == "m13" else (
                        "F0" if moment in ("m12", "m14") else "one"
                    )
                    sector_rows.append(
                        _sector_row(context, eta, rho_label, moment, sector, atom, pair)
                    )
                base = dict(sector_rows[0])
                base["atomic_sectors"] = [
                    sector["atomic_sectors"][0] for sector in sector_rows
                ]
                table.append(base)
    return table


def _n1_weights(
    context: Any,
    params: Mapping[str, Any],
    row: Mapping[str, Any],
    polynomials: Mapping[str, Mapping[int, Any]],
    errors: Mapping[str, Any],
) -> dict[int, dict[str, Any]]:
    eps = params["epsilon_core"]
    ubase = row["Ubase"]
    two = context.mpf("2")
    weights: dict[int, dict[str, Any]] = {}
    for k in BASIS:
        tc: dict[int, Any] = {k + 1: two * eps * eps * ubase / (k + 1)}
        # Ordinary-Z scale contributes the explicit ``4`` term; the
        # remaining coefficients come from DU_Z.
        tzc: dict[int, Any] = {k + 1: two * eps * eps * 4 / (k + 1)}
        for degree, coefficient in polynomials["DU"].items():
            exponent = degree + k + 1
            tc[exponent] = tc.get(exponent, _zero(context)) + two * eps * eps * coefficient / exponent
        for degree, coefficient in polynomials["DU_Z"].items():
            exponent = degree + k + 1
            tzc[exponent] = tzc.get(exponent, _zero(context)) + two * eps * eps * coefficient / exponent
        # The error belongs at rho^(k+1)/(k+1), not at rho.
        tu_error = {k + 1: two * eps * eps * errors["DU"] / (k + 1)}
        tuz_error = {k + 1: two * eps * eps * errors["DU_Z"] / (k + 1)}
        tp: dict[int, Any]
        if k == 0:
            tp = {1: _zero(context)}
        else:
            tp = {k + 1: -eps * eps * k / (k + 1)}
        weights[k] = {
            "Tu": tc,
            "TuZ": {degree: tzc[degree] for degree in sorted(tzc)},
            "Tu_error": tu_error,
            "TuZ_error": tuz_error,
            "Tp": tp,
            "Tp_error": {degree: _zero(context) for degree in tp},
        }
    return weights


def _n1_affine_table(
    context: Any,
    rows: Mapping[str, dict[str, Any]],
    params: Mapping[str, Any],
    polynomials: Mapping[str, Mapping[int, Any]],
) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for eta in ETA_LABELS:
        weights = _n1_weights(
            context,
            params,
            rows[eta]["row_values"],
            rows[eta]["candidate"],
            rows[eta]["errors"],
        )
        for rho_label in PUBLISH_RHO_LABELS:
            for k in BASIS:
                item = weights[k]
                uc, ue, ua = _endpoint_value(context, item["Tu"], rho_label, item["Tu_error"])
                uzc, uze, uza = _endpoint_value(context, item["TuZ"], rho_label, item["TuZ_error"])
                pc, pe, pa = _endpoint_value(context, item["Tp"], rho_label, item["Tp_error"])
                exact0 = _zero(context) if k == 0 else _zero(context)
                table.append(
                    {
                        "eta": eta,
                        "rho": rho_label,
                        "k": k,
                        "u_weight_C0": _interval_record(context, ua),
                        "u_weight_C1": _interval_record(context, uza),
                        "p_weight_C0": _interval_record(context, pa),
                        "p_weight_C1": _interval_record(context, exact0),
                        "u_weight_C0_candidate": _interval_record(context, uc),
                        "u_weight_C1_candidate": _interval_record(context, uzc),
                        "p_weight_C0_candidate": _interval_record(context, pc),
                        "p_weight_C1_candidate": _interval_record(context, exact0),
                        "u_weight_C0_error": _interval_record(context, ue),
                        "u_weight_C1_error": _interval_record(context, uze),
                        "p_weight_C0_error": _interval_record(context, pe),
                        "p_weight_C1_error": _interval_record(context, exact0),
                        "p_atom": "Pstar_squared",
                        "temporal_order": 1,
                        "n1_values_computed": False,
                        "ordinary_Z": True,
                        "uniform_eta_source_row": eta == "uniform_eta",
                    }
                )
    return table


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
    """Derive local leading moment function enclosures from one capsule."""

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
            "old_math_calls": 0,
            "absolute_F0_exp_materialized": False,
        },
    )
    context = _new_interval_context()
    parameters_source = _as_mapping(capsule["parameters"], "capsule.parameters")
    raw_params: dict[str, Any] = {}
    params: dict[str, Any] = {}
    for name in ("epsilon_core", "a", "s", "U0"):
        raw_params[name], params[name] = _source_interval(
            context, parameters_source, (name,), "parameters." + name
        )
    if not bool(_lower(context, params["epsilon_core"]) > context.mpf("0")):
        raise KernelInputError("epsilon_core must be strictly positive")
    if not bool(_lower(context, params["s"]) > context.mpf("0")):
        raise KernelInputError("s must be strictly positive")
    raw_polynomials, base_polynomials = _poly_packet(
        context, _as_mapping(capsule["polynomials"], "capsule.polynomials")
    )
    raw_atoms, atoms = _source_atoms(
        context, _as_mapping(capsule["source_atoms"], "capsule.source_atoms")
    )
    rows_source = _row_mapping(capsule["eta_rows"])
    rho_values = _rho_values(context)

    _emit(emit, events, "coefficient_polynomials_entered", {"family": FAMILY_NAMES[0], "before_math": True})
    row_work: dict[str, dict[str, Any]] = {}
    for eta in ETA_LABELS:
        raw_row, row_values = _row_atoms(context, _as_mapping(rows_source[eta], "eta_rows." + eta), params["epsilon_core"])
        candidates = _candidate_polynomials(context, row_values["shift"], base_polynomials)
        errors = _candidate_errors(row_values)
        row_work[eta] = {
            "raw": raw_row,
            "row_values": row_values,
            "candidate": candidates,
            "errors": errors,
            "function_packet": _function_packet(context, candidates, errors),
        }
    _emit(
        emit,
        events,
        "coefficient_polynomials_complete",
        {"family": FAMILY_NAMES[0], "rho_domain": "[0,4]", "eta_rows": 8, "uniform_eta_row": True, "axis_rho0_exact": True},
    )

    _emit(emit, events, "moment_integrals_entered", {"family": FAMILY_NAMES[1], "before_math": True})
    for eta in ETA_LABELS:
        row_work[eta]["integrals"] = _moment_integrals(
            context, row_work[eta]["candidate"], row_work[eta]["errors"]
        )
    _emit(
        emit,
        events,
        "moment_integrals_complete",
        {"family": FAMILY_NAMES[1], "integrals": ["Jphir", "JD", "Jphidr", "JD2", "Jphi2r", "Jphi2", "Jdeltar", "Jdelta", "JdeltaDUr", "Jdelta2r", "Jdelta2", "JDUr"], "old_integrals_replayed": False, "phi_minus_one_sectors_retained": True},
    )

    _emit(emit, events, "functional_errors_entered", {"family": FAMILY_NAMES[2], "before_math": True})
    for eta in ETA_LABELS:
        row_work[eta]["error_packet"] = {
            name: {
                "error": _poly_records(context, item["error"]),
                "error_Z": _poly_records(context, item["error_Z"]),
                "enclosed": _poly_records(context, item["enclosed"]),
                "enclosed_Z": _poly_records(context, item["enclosed_Z"]),
            }
            for name, item in row_work[eta]["integrals"].items()
        }
    _emit(
        emit,
        events,
        "functional_errors_complete",
        {"family": FAMILY_NAMES[2], "C0_C1_errors_explicit": True, "upper_endpoint_only_widening": True, "axis_rho0_exact": True},
    )

    _emit(emit, events, "physical_moment_scales_entered", {"family": FAMILY_NAMES[3], "before_math": True})
    for eta in ETA_LABELS:
        row_work[eta]["moment_polynomials"] = _moment_polynomials(
            context,
            params,
            row_work[eta]["row_values"],
            atoms,
            row_work[eta]["integrals"],
        )
    _emit(
        emit,
        events,
        "physical_moment_scales_complete",
        {"family": FAMILY_NAMES[3], "leading_moments": list(MOMENT_ORDER), "moment_Z_coefficients": True, "full_physical_materialized": False},
    )

    _emit(emit, events, "n1_m15_affine_binding_entered", {"family": FAMILY_NAMES[4], "before_math": True})
    n1_table = _n1_affine_table(context, row_work, params, base_polynomials)
    n1_binding = {
        "actualu": "Uz1/epsilon_core",
        "actualp": "P1/(epsilon_core*Pstar_squared)",
        "m15_1": "epsilon_core^2*(2*integral(Ulead*u)+Pstar_squared*(integral(p)-rho*p))",
        "ordinary_Z1_definition": "epsilon_core^2*(2*integral(Ulead_Z*u+Ulead*u_Z)+Pstar_squared*(integral(p_Z)-rho*p_Z))",
        "P1_distinct_from_P0_and_m15": True,
        "actual_n1_solution": None,
        "arbitrary_u_p_choice": False,
        "NS_force_or_admissible_field": False,
        "regular_axis_gauge_required": True,
        "Uz1_axis_zero_required": True,
        "P1_axis_normalization_required": True,
        "original_coupled_Volterra_required": True,
        "g6_and_Np_source_required": True,
        "old_core_binder_proved": False,
        "gateaux_derivative_of_leading_m15": False,
        "pressure_gauge_Tp0_does_not_delete_old_errors": True,
        "source_ports": _json_value(_lookup(capsule, ("n1_source_ports", "n1_ports"))),
        "n1_source_ports_provided": _lookup(capsule, ("n1_source_ports", "n1_ports")) is not None,
        "computed_value": None,
        "computed_value_reason": "typed source ports are not provided as an instantiated n=1 solution",
        "polynomial_basis": list(BASIS),
        "Tu_formula": "2*epsilon_core^2*(Ubase*rho^(k+1)/(k+1)+sum(DU_i*rho^(i+k+1)/(i+k+1)))",
        "TuZ_formula": "2*epsilon_core^2*(4*rho^(k+1)/(k+1)+sum(DU_Z_i*rho^(i+k+1)/(i+k+1)))",
        "Tu_error_formula": "2*epsilon_core^2*EDU*rho^(k+1)/(k+1)",
        "TuZ_error_formula": "2*epsilon_core^2*EDUZ*rho^(k+1)/(k+1)",
        "Tp_formula_without_Pstar_squared": "-epsilon_core^2*k*rho^(k+1)/(k+1)",
        "Tp_k0_exact_zero_before_interval_subtraction": True,
        "n1_C1_formula": "sum(TuZ*u_k+Tu*u_kZ+Pstar_squared*Tp*p_kZ)",
        "analytic_field_integral_recipe_retained": True,
        "finite_basis_is_not_n1_truncation_or_zero": True,
        "table": n1_table,
    }
    _emit(
        emit,
        events,
        "n1_m15_affine_binding_complete",
        {"family": FAMILY_NAMES[4], "n1_values_computed": False, "n1_solution_materialized": False, "basis": list(BASIS)},
    )

    moment_table = _moment_table(context, row_work)
    proof = {
        "kind": "native_O3_eta_core_moment_functionals_proof_v1",
        "proof_count": 1,
        "finite_integral": "exact finite polynomial multiplication followed by weighted integration on rho in [0,4]",
        "product_error": "E(fg)<=w*(Ef*absG+Eg*absF+Ef*Eg), w in {rho,rho^2/2}",
        "squared_error": "E(f^2)<=w*Ef*(2*absF+Ef)",
        "derivative_product_error": "derivative product rule with derivative2 bound 2*w*(Ef*absFZ+EfZ*absF+Ef*EfZ)",
        "ordinary_Z_physical_scale": "Z=a+s*eta with exact a atom; physical rho domain is [0,4]",
        "measure": "all leading moments use source-coordinate dR=epsilon*d rho, not Eulerian volume",
        "n1_parts_identity": "m15^(1)=epsilon^2*(2*integral(Ulead*u)+Pstar_squared*(integral(p)-rho*p))",
        "n1_solution_guard": "regular axis/gauge, coupled Volterra, g6, pressure source, and Np ports remain uninstantiated",
        "n1_not_a_gateaux_derivative": True,
        "conditional_on_accepted_input_theorem": True,
        "candidate_function": "unique fixedpoint enclosed by candidate +/- certified chart errors",
        "uniform_abs_bound": "sum(abs(coefficient))*4^degree, never a norm as a value",
        "phi_minus_one": "degree>=1 polynomial kept as separate axis/correction sectors; no numeric subtraction from the unit atom",
        "axis_rho0_exact": True,
        "formal_atoms": ["F0=exp(-logC)", "S=F0^2", "Pstar_squared"],
        "n1_values_computed": False,
        "global_moments": False,
        "global_five_moment": False,
        "time_n_recursion": False,
        "full_physical": False,
        "physical_materialized": False,
        "global_NS": False,
        "full_goal": False,
    }
    _emit(emit, events, "proof_bundle", proof)
    _emit(
        emit,
        events,
        "family_completed",
        {
            "proof_count": 1,
            "leading_moment_function_enclosures": True,
            "n1_m15_affine_operator": True,
            "n1_m15_values_computed": False,
            "global_moments": False,
            "global_five_moment": False,
            "time_n_recursion": False,
            "full_physical": False,
            "physical_materialized": False,
            "global_NS": False,
            "full_goal": False,
        },
    )

    parameter_data = {
        name: {"raw": _json_value(raw_params[name]), "interval": _interval_record(context, params[name])}
        for name in raw_params
    }
    parameter_data["a"]["exact_anchor_atom"] = True
    parameter_data["s"]["eta_scale"] = True
    parameter_data["Z_definition"] = "Z=a+s*eta"
    row_payload: dict[str, Any] = {}
    for eta in ETA_LABELS:
        work = row_work[eta]
        row_payload[eta] = {
            "raw": _json_value(work["raw"]),
            "parameters": {
                name: _interval_record(context, value)
                for name, value in work["row_values"].items()
                if name in ("shift", "Ubase", "ratio", "logF0_Z", "Ephi", "EphiZ", "EDU", "EDUZ")
            },
            "candidate_functions": work["function_packet"],
            "F_definition": "F=F0_anchor*ratio*Phi",
            "U_definition": "U=Ubase+DU",
            "moment_integrals": {
                name: {
                    "candidate": _poly_records(context, item["candidate"]),
                    "candidate_Z": _poly_records(context, item["candidate_Z"]),
                    "error": _poly_records(context, item["error"]),
                    "error_Z": _poly_records(context, item["error_Z"]),
                    "enclosed": _poly_records(context, item["enclosed"]),
                    "enclosed_Z": _poly_records(context, item["enclosed_Z"]),
                }
                for name, item in work["integrals"].items()
            },
            "physical_moment_polynomials": {
                moment: {
                    sector: {
                        key: _poly_records(context, value)
                        for key, value in pair.items()
                    }
                    for sector, pair in sectors.items()
                }
                for moment, sectors in work["moment_polynomials"].items()
            },
            "phi_minus_one_kept_separate": True,
            "uniform_eta_source_row": eta == "uniform_eta",
            "full_field_value_materialized": False,
        }
    summary = {
        "moment_table_rows": len(moment_table),
        "moment_table_expected_rows": 80,
        "n1_affine_table_rows": len(n1_table),
        "n1_affine_table_expected_rows": 64,
        "eta_count": 8,
        "published_rho_count": 2,
        "moment_count": 5,
        "n1_basis_count": 4,
        "full_function_domains": {
            "rho": {"lower": "0", "upper": "4", "physical": True},
            "uniform_rho": "[0,4] real interval",
            "eta": {"rows": list(ETA_LABELS), "uniform": "[-1,1]"},
            "grid_rho": list(RHO_LABELS),
        },
        "all_coefficients_source_original_R": True,
        "original_R_moments_not_Eulerian_physical_volume": True,
        "source_atoms_retained": True,
        "full_field_value_materialized": False,
        "n1_values_computed": False,
        "global_flags": {
            "global_moments": False,
            "global_five_moment": False,
            "time_n_recursion": False,
            "full_physical": False,
            "physical_materialized": False,
            "global_NS": False,
            "full_goal": False,
        },
    }
    payload = {
        "component": COMPONENT,
        "schema": SCHEMA,
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "conditional": True,
        "runtime_model_metadata": None,
        "api": dict(API),
        "api_metadata": dict(API_METADATA),
        "canonical_gate": {
            "algorithm": "sha256(canonical_json(entire_capsule))",
            "actual_sha256": hashlib.sha256(_canonical_bytes(capsule)).hexdigest(),
            "root_constant_used": ROOT_FREEZE_SHA,
            "raw_capsule_sha256_receipt": _lookup(capsule, ("raw_capsule_sha256", "raw_sha256")) or ROOT_CAPSULE_RAW_SHA,
            "canonical_field_in_capsule": False,
        },
        "source_ids": _json_value(capsule["source_ids"]),
        "family_identity": _json_value(_lookup(capsule, ("family_identity",))),
        "native_order": _json_value(_lookup(capsule, ("native_order",))),
        "definitions": _json_value(_lookup(capsule, ("definitions",))),
        "preservation": _json_value(_lookup(capsule, ("preservation",))),
        "metadata": _json_value(_lookup(capsule, ("metadata",))),
        "projection_pins": _json_value(_lookup(capsule, ("projection_pins",))),
        "parameters": parameter_data,
        "polynomials": {"accepted_records": raw_polynomials, "source_original_R": True},
        "source_atoms": {
            "one": {
                "expression": "1 (unit identity atom)",
                "formal_identity": True,
                "materialized": False,
            },
            **{
                name: _atom_packet(context, raw_atoms[name], atoms[name], expression)
                for name, expression in (
                    ("F0", "exp(-logC)"),
                    ("S", "F0^2"),
                    ("Pstar_squared", "Pstar^2"),
                )
            },
        },
        "eta_rows": row_payload,
        "retained_native_packet": _json_value(capsule["retained_native_packet"]),
        "old_core_to_eta_function_identity_proved": False,
        "preservation_contract": {
            "all_old_error_data_retained": True,
            "inherited_errors_deleted": False,
            "reserved_contribution_registered_once": True,
            "six_signed_channels_retained": True,
            "leading_moments_appended_to_n1_ledger": False,
        },
        "moment_table": moment_table,
        "n1_affine_table": n1_table,
        "moment_integrals": {eta: row_payload[eta]["moment_integrals"] for eta in ETA_LABELS},
        "physical_moment_scales": {eta: row_payload[eta]["physical_moment_polynomials"] for eta in ETA_LABELS},
        "n1_m15_affine_binding": n1_binding,
        "summary": summary,
        "grid": _json_value(capsule["grid"]),
        "physical_rho_grid": list(RHO_LABELS),
        "uniform_rho_grid": list(UNIFORM_RHO_LABELS),
        "uniform_rho_interval": ["0", "4"],
        "axis_structural_proof": {
            "rho": "0",
            "phi_degree0": "one identity atom",
            "phi_Z_degree0": "exact zero",
            "DU_degree0": "exact zero",
            "DU_Z_degree0": "exact zero",
            "integral_moments_at_axis": "exact zero",
            "separate_from_uniform_rho": True,
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
        "count_event_map": {key: list(value) for key, value in COUNT_EVENT_MAP.items()},
        "static_limitations": [
            "local native eta functional enclosure only on rho in [0,4]",
            "uniform_eta is the accepted certified source row on eta in [-1,1], not an analytic extension or physical splice",
            "F0, S, and Pstar_squared remain formal positive source atoms",
            "the temporal n=1 m15 operator is typed but uninstantiated",
            "no global moment closure, time-n recursion, full physical field, or Cartesian NS claim",
        ],
    }
    result = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "passed": True,
        "proof_count": 1,
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "conditional": True,
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
        "count_event_map": {key: list(value) for key, value in COUNT_EVENT_MAP.items()},
        "leading_moment_function_enclosures": True,
        "n1_m15_affine_operator": True,
        "n1_m15_values_computed": False,
        "n1_solution_materialized": False,
        "global_moments": False,
        "global_five_moment": False,
        "time_n_recursion": False,
        "full_physical": False,
        "physical_materialized": False,
        "global_NS": False,
        "full_goal": False,
        "payload": payload,
    }
    _, return_hash, return_bytes = _serialize_result(result)
    result["return_code_hash"] = return_hash
    result["return_code_bytes"] = return_bytes
    result["code_hash"] = return_hash
    result["code_bytes"] = return_bytes
    result["serialization_hash_scope"] = "canonical result before hash metadata"
    return result


__all__ = [
    "API",
    "API_METADATA",
    "CALLBACK_EVENT_COUNTS",
    "CALLBACK_EVENTS",
    "CAPSULE_SCHEMA",
    "COMPONENT",
    "COUNT_EVENT_MAP",
    "EXPECTED_COUNTS",
    "EXPECTED_FLAGS",
    "FAMILY_NAMES",
    "INTERVAL_DPS",
    "KernelInputError",
    "MANAGER_JOURNAL_EVENT_COUNT",
    "ROOT_CAPSULE_RAW_SHA",
    "ROOT_FREEZE_SHA",
    "SCHEMA",
    "SOURCE_IDS",
    "SUCCESS_COUNTS",
    "SUCCESS_FLAGS",
    "derive",
]
