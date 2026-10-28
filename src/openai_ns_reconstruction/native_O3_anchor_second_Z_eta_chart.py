"""Standalone second ordinary-Z jet and analytic eta chart leaf.

The root supplies one authenticated, root-frozen capsule.  All source
numbers are directed MPF records; source decimal display strings are never
parsed.  The leaf consumes accepted first-Picard/axial records and performs
only the new pressure-fourth-derivative, second-Z, eta-chart, compatible-Q,
and physical-coefficient constructions described by the source traces.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "native_O3_anchor_second_Z_eta_chart_kernel_result_v1"
CAPSULE_SCHEMA = "native_O3_anchor_second_Z_eta_chart_source_capsule_v1"
COMPONENT = "native_O3_anchor_second_Z_eta_chart"
INTERVAL_DPS = 280
SOURCE_PRECISION_BITS = 867

# The raw capsule receipt and the canonical JSON digest are distinct.
ROOT_CAPSULE_RAW_SHA = (
    "e7342d77112d868606d2aaf6e8f88a82d19dc518c4d57d364cc1c7e0d86f9d7e"
)
ROOT_FREEZE_SHA = (
    "27db749374f496abcb3d50c68e87bf9076e8def0180c4b25d068d01837833eea"
)

SOURCE_IDS = {
    "axial_Q_result_sha256": (
        "f543d5d7064bf9c95fdc10a3f54b3de2c584907045ccfd41915190d24007926e"
    ),
    "anchor_axis_result_sha256": (
        "c5a6f8273105a8f8710a44455de072ca93b5cd02b60f382249b8a0e5a3b016e8"
    ),
    "common_capsule_sha256": (
        "012511ca881a0207e04b7f494ce8ca275d7ea8a0a04d2914a824b5f846d64e0e"
    ),
}

FAMILY_NAMES = (
    "pressure4_C3",
    "second_Z_operator",
    "analytic_eta_chart",
    "compatible_pressure_Q",
    "physical_eta_velocity",
)
CALLBACK_EVENTS = (
    "family_entered",
    "pressure4_C3_entered",
    "pressure4_C3_complete",
    "second_Z_operator_entered",
    "second_Z_operator_complete",
    "analytic_eta_chart_entered",
    "analytic_eta_chart_complete",
    "compatible_pressure_Q_entered",
    "compatible_pressure_Q_complete",
    "physical_eta_velocity_entered",
    "physical_eta_velocity_complete",
    "proof_bundle",
    "family_completed",
)
CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "pressure4_C3_entered": 1,
    "pressure4_C3_complete": 1,
    "second_Z_operator_entered": 1,
    "second_Z_operator_complete": 1,
    "analytic_eta_chart_entered": 1,
    "analytic_eta_chart_complete": 1,
    "compatible_pressure_Q_entered": 1,
    "compatible_pressure_Q_complete": 1,
    "physical_eta_velocity_entered": 1,
    "physical_eta_velocity_complete": 1,
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
    "pressure4_C3_family_calls": 1,
    "second_Z_operator_family_calls": 1,
    "analytic_eta_chart_family_calls": 1,
    "compatible_pressure_Q_family_calls": 1,
    "physical_eta_velocity_family_calls": 1,
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
    "local_eta_functional_enclosures": True,
    "local_second_Z_jet": True,
    "local_Q": True,
    "local_coefficients": True,
    "local_Q_source": True,
    "local_velocity_coefficients": True,
    "local_eta_functional_enclosures_enclosed": True,
    "local_second_Z_jet_enclosed": True,
    "local_Q_enclosed": True,
    "local_coefficients_enclosed": True,
    "global_jet_admission": False,
    "global_Q": False,
    "global_source_jet_admission": False,
    "full_physical": False,
    "global_five_moment": False,
    "time_n_recursion": False,
    "full_goal": False,
    "full_goal_complete": False,
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
    "interval_power_internal_exp_allowed": True,
}
COUNT_EVENT_MAP = {
    "source_family_calls": ["family_entered"],
    "new_math_family_calls": [
        "pressure4_C3_entered",
        "second_Z_operator_entered",
        "analytic_eta_chart_entered",
        "compatible_pressure_Q_entered",
        "physical_eta_velocity_entered",
    ],
    "pressure4_C3_family_calls": ["pressure4_C3_entered"],
    "second_Z_operator_family_calls": ["second_Z_operator_entered"],
    "analytic_eta_chart_family_calls": ["analytic_eta_chart_entered"],
    "compatible_pressure_Q_family_calls": ["compatible_pressure_Q_entered"],
    "physical_eta_velocity_family_calls": ["physical_eta_velocity_entered"],
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
    "pressure4_C3_family_calls": 1,
    "second_Z_operator_family_calls": 1,
    "analytic_eta_chart_family_calls": 1,
    "compatible_pressure_Q_family_calls": 1,
    "physical_eta_velocity_family_calls": 1,
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
    "local_eta_functional_enclosures": True,
    "local_second_Z_jet": True,
    "local_Q": True,
    "local_coefficients": True,
    "local_Q_source": True,
    "local_velocity_coefficients": True,
    "local_eta_functional_enclosures_enclosed": True,
    "local_second_Z_jet_enclosed": True,
    "local_Q_enclosed": True,
    "local_coefficients_enclosed": True,
    "global_jet_admission": False,
    "global_Q": False,
    "global_source_jet_admission": False,
    "full_physical": False,
    "global_five_moment": False,
    "time_n_recursion": False,
    "full_goal": False,
    "full_goal_complete": False,
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
    return json.dumps(
        _json_value(capsule),
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
        return any(
            _contains_text(str(key), needles) or _contains_text(child, needles)
            for key, child in value.items()
        )
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
    """Construct an interval from its authenticated paired MPF words only."""

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
    return {
        "lower": _display_endpoint(context, lower),
        "upper": _display_endpoint(context, upper),
        "width": _display_endpoint(context, upper - lower),
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
        return _make_interval(context, hi, hi)
    if bool(hi <= zero):
        return _make_interval(context, -lo, -lo)
    return _make_interval(context, zero, -lo if bool(-lo >= hi) else hi)


def _widen(context: Any, value: Any, error: Any) -> Any:
    err = _upper(context, error)
    return _make_interval(
        context,
        _lower(context, value) - err,
        _upper(context, value) + err,
    )


def _parameter_record(
    context: Any, parameters: Mapping[str, Any], name: str
) -> tuple[Any, Any]:
    raw = parameters.get(name)
    if raw is None:
        raise KernelInputError("missing source parameter " + name)
    return raw, _interval_value(context, raw, "parameters." + name)


def _capsule_record(
    context: Any, mapping: Mapping[str, Any], name: str, label: str
) -> tuple[Any, Any]:
    raw = mapping.get(name)
    if raw is None:
        raise KernelInputError("missing source record " + label)
    return raw, _interval_value(context, raw, label)


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


def _new_interval_context() -> Any:
    # The root invokes this leaf once; the interval context is private.
    from mpmath.ctx_iv import MPIntervalContext

    context = MPIntervalContext()
    context.dps = INTERVAL_DPS
    return context


def _exact_text(context: Any, text: str) -> Any:
    if "/" in text:
        numerator, denominator = text.split("/", 1)
        return context.mpf(int(numerator)) / context.mpf(int(denominator))
    return context.mpf(text)


def _poly_eval(context: Any, coefficients: Mapping[int, Any], rho: Any) -> Any:
    value = context.mpf("0")
    for degree, coefficient in coefficients.items():
        value = value + coefficient * rho**degree
    return value


def _poly_add(
    context: Any, *polynomials: Mapping[int, Any]
) -> dict[int, Any]:
    result: dict[int, Any] = {}
    for polynomial in polynomials:
        for degree, coefficient in polynomial.items():
            result[degree] = result.get(degree, context.mpf("0")) + coefficient
    return {degree: result[degree] for degree in sorted(result)}


def _poly_scale(
    context: Any, polynomial: Mapping[int, Any], factor: Any
) -> dict[int, Any]:
    return {degree: factor * coefficient for degree, coefficient in polynomial.items()}


def _poly_mul(
    context: Any, left: Mapping[int, Any], right: Mapping[int, Any]
) -> dict[int, Any]:
    result: dict[int, Any] = {}
    for left_degree, left_value in left.items():
        for right_degree, right_value in right.items():
            degree = left_degree + right_degree
            result[degree] = result.get(degree, context.mpf("0")) + left_value * right_value
    return {degree: result[degree] for degree in sorted(result)}


def _poly_integral(context: Any, polynomial: Mapping[int, Any]) -> dict[int, Any]:
    """V f = integral_0^rho f(s) ds, with no constant of integration."""

    return {
        degree + 1: coefficient / (degree + 1)
        for degree, coefficient in polynomial.items()
    }


def _poly_average(context: Any, polynomial: Mapping[int, Any]) -> dict[int, Any]:
    """rho^-1 V f, represented continuously at rho = 0."""

    return {
        degree: coefficient / (degree + 1)
        for degree, coefficient in polynomial.items()
    }


def _jnu(
    context: Any, polynomial: Mapping[int, Any], nu: int
) -> dict[int, Any]:
    """J_nu(rho^n) = rho^(n+1) / ((n+1)(n+nu))."""

    return {
        degree + 1: coefficient / ((degree + 1) * (degree + nu))
        for degree, coefficient in polynomial.items()
    }


def _poly_value(context: Any, polynomial: Mapping[int, Any], rho: Any) -> Any:
    value = context.mpf("0")
    for degree, coefficient in polynomial.items():
        value = value + coefficient * rho**degree
    return value


def _poly_abs_bound(
    context: Any, polynomial: Mapping[int, Any], radius: Any
) -> Any:
    bound = context.mpf("0")
    for degree, coefficient in polynomial.items():
        bound = bound + _abs_upper(context, coefficient) * radius**degree
    return _make_interval(context, bound, bound)


def _poly_records(
    context: Any, polynomial: Mapping[int, Any]
) -> list[dict[str, Any]]:
    return [
        {"degree": degree, "interval": _interval_record(context, polynomial[degree])}
        for degree in sorted(polynomial)
    ]


def _coefficient_values(
    context: Any,
    mapping: Mapping[str, Any],
    names: Sequence[str],
    label: str,
    count: int,
    first_degree: int,
) -> dict[int, Any]:
    raw = _lookup(mapping, names)
    if raw is None:
        raise KernelInputError("missing coefficient packet " + label)
    if _is_mapping(raw):
        nested = _lookup(raw, ("coefficients", "values", "records", "items"))
        if nested is not None:
            raw = nested
        else:
            selected: list[Any] = []
            for degree in range(first_degree, first_degree + count):
                item = None
                for key in (str(degree), degree, str(degree - first_degree)):
                    if key in raw:
                        item = raw[key]
                        break
                if item is None:
                    raise KernelInputError(
                        label + " is missing coefficient degree " + str(degree)
                    )
                selected.append(item)
            raw = selected
    if not _is_sequence(raw) or len(raw) != count:
        raise KernelInputError(
            label + " must contain exactly " + str(count) + " records"
        )
    return {
        first_degree + index: _interval_value(
            context, value, label + "[" + str(first_degree + index) + "]"
        )
        for index, value in enumerate(raw)
    }


def _required_mapping(
    capsule: Mapping[str, Any], names: Sequence[str], label: str
) -> Mapping[str, Any]:
    value = _lookup(capsule, names)
    if value is None:
        raise KernelInputError("missing source mapping " + label)
    return _as_mapping(value, label)


def _rows_mapping(mapping: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    rows = _lookup(mapping, ("rows", "table", "records"))
    if rows is None:
        rows = mapping
    return _as_mapping(rows, label)


def _physical_mapping(mapping: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    slices = _lookup(mapping, ("time_slices", "physical_time_slices", "slices"))
    if slices is None:
        slices = mapping
    return _as_mapping(slices, label)


def _validate_capsule(capsule: Any) -> Mapping[str, Any]:
    capsule = _as_mapping(capsule, "data.capsule")
    if capsule.get("schema") != CAPSULE_SCHEMA:
        raise KernelInputError("unexpected second-Z eta source capsule schema")
    if capsule.get("component") != COMPONENT:
        raise KernelInputError("source capsule component does not match")
    if (
        not _valid_sha(ROOT_FREEZE_SHA)
        or ROOT_FREEZE_SHA == "0" * 64
        or ROOT_FREEZE_SHA == "__ROOT_FREEZE_SHA_REQUIRED__"
    ):
        raise KernelInputError("ROOT_FREEZE_SHA is not frozen")
    actual = hashlib.sha256(_canonical_bytes(capsule)).hexdigest()
    if actual.lower() != ROOT_FREEZE_SHA.lower():
        raise KernelInputError("source capsule canonical digest mismatch")

    source_ids = _as_mapping(capsule.get("source_ids"), "capsule.source_ids")
    for name, expected in SOURCE_IDS.items():
        if source_ids.get(name) != expected:
            raise KernelInputError("source id mismatch: " + name)

    parameters = _required_mapping(capsule, ("parameters",), "capsule.parameters")
    for name in ("delta", "epsilon_core", "j", "sigma", "h", "s"):
        if parameters.get(name) is None:
            raise KernelInputError("parameters is missing " + name)

    seed = _required_mapping(capsule, ("seed", "seed_inputs"), "capsule.seed")
    seed_names = (
        "a",
        "d",
        "L",
        "U0",
        "U0_Z",
        "H1",
        "H2",
        "A_ZZ",
        "C",
        "C1",
        "C2",
        "beta",
        "K",
        "K1",
        "L1",
        "chi2",
        "s2",
        "Qtheta",
        "Qhat",
        "D",
        "Dstar",
        "Q0",
        "p0",
        "p1",
        "p2",
        "p3",
        "I2",
        "I3",
        "qP",
        "B1",
        "core_seed_scale",
    )
    for name in seed_names:
        if seed.get(name) is None:
            raise KernelInputError("seed is missing " + name)

    polynomials = _required_mapping(
        capsule, ("polynomials", "parent_polynomials"), "capsule.polynomials"
    )
    for name in ("Phi1", "DeltaU1", "Phi1_Z", "DeltaU1_Z"):
        if polynomials.get(name) is None:
            raise KernelInputError("polynomials is missing " + name)

    parent = _required_mapping(capsule, ("parent", "parent_data"), "capsule.parent")
    axial = _rows_mapping(
        _as_mapping(
            _lookup(parent, ("axial_rows", "axial", "rows_axial")),
            "parent.axial_rows",
        ),
        "parent.axial_rows",
    )
    average = _rows_mapping(
        _as_mapping(
            _lookup(parent, ("average_rows", "average_Q", "average")),
            "parent.average_rows",
        ),
        "parent.average_rows",
    )
    physical = _physical_mapping(
        _as_mapping(
            _lookup(
                parent,
                ("physical_time_slices", "physical_velocity", "physical"),
            ),
            "parent.physical_time_slices",
        ),
        "parent.physical_time_slices",
    )
    rho_names = ("0", "1/4", "1", "2", "4.1", "uniform_rho")
    for label in rho_names:
        if label not in axial:
            raise KernelInputError("parent.axial_rows is missing " + label)
        if label not in average:
            raise KernelInputError("parent.average_rows is missing " + label)
    for label in ("1", "1/4", "1/16", "1/64"):
        if label not in physical:
            raise KernelInputError(
                "parent.physical_time_slices is missing " + label
            )

    remainder = _required_mapping(capsule, ("remainder",), "capsule.remainder")
    for name in ("E_norm", "E0", "EZ", "M_Phi", "M_Psi", "uniform_abs_Phi1"):
        if remainder.get(name) is None:
            raise KernelInputError("remainder is missing " + name)

    source_scale = _required_mapping(
        capsule, ("source_scale",), "capsule.source_scale"
    )
    for name in ("S_bound", "B_bound", "logC_record"):
        if source_scale.get(name) is None:
            raise KernelInputError("source_scale is missing " + name)
    if source_scale.get("formal_positive") is False or source_scale.get(
        "formal_positiveTrue"
    ) is False:
        raise KernelInputError("source_scale must retain formal positive atoms")
    if source_scale.get("materialized") is True or source_scale.get(
        "absolute_exponential_materialized"
    ) is True:
        raise KernelInputError("source scale materialization is forbidden")

    f0_chart = _required_mapping(capsule, ("F0_chart", "f0_chart"), "capsule.F0_chart")
    f0_rows = _rows_mapping(f0_chart, "capsule.F0_chart.rows")
    for label in ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1"):
        if label not in f0_rows:
            raise KernelInputError("F0_chart.rows is missing " + label)
        row = _as_mapping(f0_rows[label], "F0_chart.rows." + label)
        if _lookup(row, ("ratio", "relative_profile", "relative_F0")) is None:
            raise KernelInputError("F0_chart ratio is missing " + label)
        if _lookup(row, ("logF0_Z", "ordinary_logF0_Z", "logF0Z")) is None:
            raise KernelInputError("F0_chart logF0_Z is missing " + label)
    if _lookup(f0_chart, ("B_all", "B_all_record", "positive_B")) is None:
        raise KernelInputError("F0_chart B_all is required")
    if not str(
        _lookup(f0_chart, ("exact_G_definition", "G_definition", "definition"))
        or ""
    ).strip():
        raise KernelInputError("F0_chart exact G definition is required")

    family_identity = _lookup(
        capsule, ("family_identity", "family_identity_hashes")
    )
    if not _is_mapping(family_identity) or not family_identity:
        raise KernelInputError("family identity hashes are required")
    for name, value in family_identity.items():
        if not _valid_sha(value):
            raise KernelInputError("family identity hash is invalid: " + str(name))

    grid = _required_mapping(capsule, ("grid",), "capsule.grid")
    expected_grid = {
        "rho": ["0", "1/4", "1", "2", "4.1"],
        "eta": ["-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1"],
        "uniform_rho": ["0", "4.1"],
        "uniform_eta": ["-1", "1"],
        "physical_tau": ["1", "1/4", "1/16", "1/64"],
        "physical_azimuth_labels": ["0", "pi/2", "pi", "3pi/2"],
    }
    for name, expected in expected_grid.items():
        if list(grid.get(name, ())) != expected:
            raise KernelInputError("capsule grid does not match " + name)

    guards = _lookup(capsule, ("preservation", "preservation_guards", "guards"))
    if not _is_mapping(guards):
        raise KernelInputError("preservation guards are required")
    for name in (
        "global_jet_admission",
        "global_actual_source_jet_admission",
        "global_actual_source_Z_jets_admitted",
        "global_Q",
        "global_source_jet_admission",
        "full_physical",
        "physical_materialized",
        "physical_field_materialized",
        "full_physical_field_materialized",
        "global_five_moment",
        "time_n_recursion",
        "full_goal",
        "full_goal_complete",
        "root_midpoint_selected",
        "old_fixedpoint_norm_error_recomputed",
        "old_source_fields_recomputed",
    ):
        if name in guards and guards[name] is not False:
            raise KernelInputError("preservation guard is not false: " + name)
    for name in (
        "N_materialization_calls",
        "absolute_F0_exp_materialization_calls",
        "F0_absolute_exp_calls",
        "ancestor_calls",
        "old_math_calls",
        "old_provider_calls",
        "provider_calls",
        "checker_calls",
        "selftest_calls",
        "source_value_reducer_calls",
    ):
        if name in guards and guards[name] != 0:
            raise KernelInputError("preservation count is not zero: " + name)

    definitions = _lookup(
        capsule, ("definition_guards", "definitions", "definition")
    )
    if not _is_mapping(definitions):
        raise KernelInputError("definition guards are required")
    if not _contains_text(
        definitions, ("ordinary Z", "ordinary_Z", "fixed rho")
    ):
        raise KernelInputError("ordinary-Z definition guard is required")
    if _contains_text(definitions, ("midpoint selected", "midpoint_selected=true")):
        raise KernelInputError("midpoint selection is forbidden")
    return capsule


def _parameter_packet(
    context: Any, parameters: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    raw_packet: dict[str, Any] = {}
    packet: dict[str, Any] = {}
    for name in ("delta", "epsilon_core", "j", "sigma", "h", "s"):
        raw, value = _parameter_record(context, parameters, name)
        raw_packet[name] = _json_value(raw)
        packet[name] = value
    for name in ("epsilon_core", "sigma", "h", "s"):
        if not _strict_positive(context, packet[name]):
            raise KernelInputError(name + " must be strictly positive")
    if not bool(_lower(context, packet["delta"]) > context.mpf("0")):
        raise KernelInputError("delta must be positive")
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
        "H2",
        "A_ZZ",
        "C",
        "C1",
        "C2",
        "beta",
        "K",
        "K1",
        "L1",
        "chi2",
        "s2",
        "Qtheta",
        "Qhat",
        "D",
        "Dstar",
        "Q0",
        "p0",
        "p1",
        "p2",
        "p3",
        "I2",
        "I3",
        "qP",
        "B1",
        "core_seed_scale",
    )
    raw_packet: dict[str, Any] = {}
    packet: dict[str, Any] = {}
    for name in names:
        raw, value = _capsule_record(context, seed, name, "seed." + name)
        raw_packet[name] = _json_value(raw)
        packet[name] = value
    if not bool(_upper(context, packet["a"]) < context.mpf("0")):
        raise KernelInputError("exact anchor must be strictly negative")
    for name in ("d", "L", "qP", "B1", "core_seed_scale"):
        if not _strict_positive(context, packet[name]):
            raise KernelInputError("seed " + name + " must be strictly positive")
    return raw_packet, packet


def _polynomial_packet(
    context: Any, polynomials: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, dict[int, Any]]]:
    specs = (
        ("Phi1", ("Phi1", "Phi1_coefficients"), 3, 0),
        ("DeltaU1", ("DeltaU1", "DeltaU1_coefficients"), 3, 1),
        ("Phi1_Z", ("Phi1_Z", "Phi1_Z_coefficients"), 3, 1),
        ("DeltaU1_Z", ("DeltaU1_Z", "DeltaU1_Z_coefficients"), 3, 1),
    )
    raw: dict[str, Any] = {}
    values: dict[str, dict[int, Any]] = {}
    for name, aliases, count, first_degree in specs:
        raw[name] = _json_value(_lookup(polynomials, aliases))
        values[name] = _coefficient_values(
            context, polynomials, aliases, "polynomials." + name, count, first_degree
        )
    return raw, values


def _pressure4_C3(
    context: Any,
    parameters: Mapping[str, Any],
    seed: Mapping[str, Any],
) -> dict[str, Any]:
    delta = parameters["delta"]
    epsilon = parameters["epsilon_core"]
    a = seed["a"]
    d = seed["d"]
    L = seed["L"]
    B1 = seed["B1"]
    qP = seed["qP"]
    core_seed_scale = seed["core_seed_scale"]
    if not bool(_lower(context, qP) >= context.mpf("1")):
        raise KernelInputError("pressure qP must satisfy qP >= 1")
    I2 = seed["I2"]
    I3 = seed["I3"]
    I4 = _make_interval(
        context,
        6 * _lower(context, B1) / _upper(context, qP) ** 6,
        60 * _upper(context, B1) / _lower(context, qP) ** 4,
    )
    Pnorm4 = -12 * I2 + 48 * a**2 * I3 - 16 * a**4 * I4
    p4 = core_seed_scale * Pnorm4
    H3 = context.mpf("-24")
    N3 = (
        epsilon * (12 * seed["A_ZZ"] - 4 * H3)
        + (12 + 6 * delta) * seed["p2"]
        + (8 + 2 * delta) * a * seed["p3"]
        - d * p4
    )
    C3 = N3 / (2 * L) + 6 * delta * a * seed["C2"] / L + 6 * delta * seed["C1"] / L
    zmax_upper = _upper(context, _abs_upper(context, a)) + _upper(
        context, _abs_upper(context, parameters["s"])
    )
    uniform_p4_sup = _make_interval(
        context,
        context.mpf("0"),
        _upper(context, core_seed_scale)
        * _upper(context, B1)
        * (
            context.mpf("36")
            + context.mpf("576") * zmax_upper**2
            + context.mpf("960") * zmax_upper**4
        ),
    )
    return {
        "H3": H3,
        "I4": I4,
        "Pnorm4": Pnorm4,
        "p4": p4,
        "N3": N3,
        "C3": C3,
        "uniformP4sup": uniform_p4_sup,
        "Zmax_upper": zmax_upper,
        "definitions": {
            "I4": "integral beta*(beta+1)*(beta+2)*(beta+3)*q^(-beta-4) dm(beta)",
            "I4_enclosure": "6*B1.lower/qP.upper^6 <= I4 <= 60*B1.upper/qP.lower^4",
            "Pnorm4": "-12*I2+48*a^2*I3-16*a^4*I4",
            "p4": "core_seed_scale*Pnorm4",
            "N3": "epsilon*(12*A_ZZ-4*H3)+(12+6*delta)*p2+(8+2*delta)*a*p3-d*p4",
            "C3": "N3/(2*L)+6*delta*a*C2/L+6*delta*C1/L",
            "uniformP4sup": "core_seed_scale.upper*B1.upper*(36+576*Zmax^2+960*Zmax^4)",
        },
    }


def _second_Z_operator(
    context: Any,
    parameters: Mapping[str, Any],
    seed: Mapping[str, Any],
    pressure4: Mapping[str, Any],
    source_scale: Mapping[str, Any],
    polynomials: Mapping[str, dict[int, Any]],
) -> dict[str, Any]:
    delta = parameters["delta"]
    epsilon = parameters["epsilon_core"]
    sigma = parameters["sigma"]
    a = seed["a"]
    d = seed["d"]
    L = seed["L"]
    U0 = seed["U0"]
    H1 = seed["H1"]
    H2 = seed["H2"]
    C = seed["C"]
    C1 = seed["C1"]
    C2 = seed["C2"]
    C3 = pressure4["C3"]
    B = _record_or_interval(
        context, source_scale, ("B_bound", "B"), "source_scale.B_bound"
    )
    L1 = seed["L1"]
    L2 = -2 * delta
    d1 = -2 * a
    d2 = context.mpf("-2")
    chi2 = seed["chi2"]
    chi3 = 6 * H1 * H2 / sigma**2
    s3 = -2 * (2 * L1 * H1 + L * H2) / (epsilon * sigma**2)

    J = (1 - delta) * a * C + d * C1
    J1 = (1 - delta) * C - (1 + delta) * a * C1 + d * C2
    J2_comb = -2 * delta * C1 - (3 + delta) * a * C2 + d * C3
    Hs0 = d * C
    Hs1 = d1 * C + d * C1
    Hs2 = d2 * C + 2 * d1 * C1 + d * C2
    Qtheta2 = 2 * delta * C1 - (3 - delta) * a * C2 + d * C3
    Wb = 1 - (1 - delta) * a * U0 - 4 * d
    Athetab = delta / 2 * (1 - 2 * a * U0)
    g1 = d * H1 / sigma**2
    g2 = (d * H2 + 2 * d1 * H1) / sigma**2

    Ttheta0 = {
        0: -L * seed["beta"],
        1: -seed["Qtheta"] / 2,
    }
    Ttheta1 = {
        0: -parameters["j"],
        1: -seed["Qhat"] / 2,
        2: -d * C * chi2 / 4,
    }
    Ttheta2 = {
        1: -Qtheta2 / 2 - chi2 * (Wb / 2 + Athetab / 4 + H1 / 2),
        2: -chi2 * (-J / 4 - delta * a * C / 4 + Hs1 / 2)
        - chi3 * d * C / 4,
    }
    quotient_factor = 2 * L1**2 / L**3 - L2 / L**2
    Ebar_theta2 = _poly_add(
        context,
        _poly_scale(context, Ttheta2, epsilon / L),
        _poly_scale(context, Ttheta1, -2 * epsilon * L1 / L**2),
        _poly_scale(context, Ttheta0, epsilon * quotient_factor),
        {1: -(g2 * C + 2 * g1 * C1)},
    )
    J2_theta = _jnu(context, Ebar_theta2, 2)
    phi_minus_one = {
        degree: coefficient
        for degree, coefficient in polynomials["Phi1"].items()
        if degree >= 1
    }
    J2_phi_minus_one = _jnu(context, phi_minus_one, 2)
    delta_phi_term_theta = _poly_scale(context, J2_theta, context.mpf("0.5"))
    delta_phi_term_resolvent = _poly_scale(
        context, J2_phi_minus_one, -chi2 / 2
    )
    delta_phi_zz = _poly_add(
        context, delta_phi_term_theta, delta_phi_term_resolvent
    )
    phi_zz = _poly_add(
        context, {1: -chi2 / 4}, delta_phi_zz
    )

    Tbar0 = {
        1: epsilon * seed["K"] * C - 2 * (2 + delta) * a * B,
        2: epsilon * seed["D"] / 2,
    }
    Tbar1 = {
        1: epsilon * (seed["K1"] * C + (seed["K"] + H1) * C1)
        + d * B * seed["s2"]
        - 2 * (2 + delta) * B,
        2: epsilon * seed["Dstar"] / 2 - d * B * chi2 / 4,
    }
    Tbar2 = {
        1: epsilon
        * (
            -8 * (3 + delta) * C
            + 2 * seed["K1"] * C1
            + seed["K"] * C2
            + H2 * C1
            + 2 * H1 * C2
        )
        + B * (d * s3 - 2 * (4 + delta) * a * seed["s2"]),
        2: -epsilon
        / 2
        * (J2_comb * C + 2 * J1 * C1 + J * C2)
        - epsilon
        * (1 + delta)
        * (4 * C * C1 + 2 * a * (C1**2 + C * C2))
        + epsilon * (Hs2 * C1 + 2 * Hs1 * C2 + Hs0 * C3)
        + B * ((5 + delta) * a * chi2 / 2 - d * chi3 / 4),
    }
    Ebar_z2 = _poly_add(
        context,
        _poly_scale(context, Tbar2, 1 / L),
        _poly_scale(context, Tbar1, -2 * L1 / L**2),
        _poly_scale(context, Tbar0, quotient_factor),
    )
    delta_u_zz = _poly_add(
        context,
        {1: C2},
        _poly_scale(context, _jnu(context, Ebar_z2, 1), context.mpf("0.5")),
    )

    return {
        "combinations": {
            "J": J,
            "J1": J1,
            "J2": J2_comb,
            "Hs0": Hs0,
            "Hs1": Hs1,
            "Hs2": Hs2,
            "Qtheta2": Qtheta2,
            "Wb": Wb,
            "Athetab": Athetab,
            "g1": g1,
            "g2": g2,
            "chi3": chi3,
            "s3": s3,
            "K2": -8 * (3 + delta),
            "L2": L2,
            "d1": d1,
            "d2": d2,
        },
        "Ttheta0": Ttheta0,
        "Ttheta1": Ttheta1,
        "Ttheta2": Ttheta2,
        "Ebar_theta2": Ebar_theta2,
        "J2_Ebar_theta2": J2_theta,
        "J2_Phi1_minus_one": J2_phi_minus_one,
        "DeltaPhi1ZZ_theta_term": delta_phi_term_theta,
        "DeltaPhi1ZZ_resolvent_term": delta_phi_term_resolvent,
        "DeltaPhi1ZZ": delta_phi_zz,
        "Phi1ZZ": phi_zz,
        "Tbar0": Tbar0,
        "Tbar1": Tbar1,
        "Tbar2": Tbar2,
        "Ebar_z2": Ebar_z2,
        "DeltaU1ZZ": delta_u_zz,
        "definitions": {
            "ordinary_Z": "differentiate the full global map at fixed rho before H0(a)=0",
            "Jnu": "J_nu(rho^n)=rho^(n+1)/((n+1)*(n+nu))",
            "Phi0_ZZ": "-chi2*rho/4",
            "DeltaPhi1ZZ": "(1/2)*J2(Ebar_theta2)-(chi2/2)*J2(Phi1-1 accepted DATA)",
            "Phi1ZZ": "-chi2*rho/4+DeltaPhi1ZZ",
            "DeltaU1ZZ": "C2*rho+(1/2)*J1(Ebar_z2)",
            "Pcal_ZZ": "S*s2*rho-S*chi2*rho^2/4",
            "Pcal_ZZZ": "S*s3*rho-S*chi3*rho^2/4",
        },
        "source_units": "B=epsilon_core^2*S is formal positive source DATA",
    }


def _anchor_second_rows(
    context: Any,
    parent_axial: Mapping[str, Any],
    second: Mapping[str, Any],
    ezz_bound: Any,
    epsilon: Any,
) -> tuple[dict[str, Any], dict[str, Any]]:
    rho_names = ("0", "1/4", "1", "2", "4.1", "uniform_rho")
    rows: dict[str, Any] = {}
    states: dict[str, Any] = {}
    delta_m_zz_poly = _poly_average(context, second["DeltaU1ZZ"])
    for name in rho_names:
        rho = (
            _make_interval(context, context.mpf("0"), _exact_text(context, "4.1"))
            if name == "uniform_rho"
            else _exact_text(context, name)
        )
        phi_candidate = _poly_value(context, second["Phi1ZZ"], rho)
        du_candidate = _poly_value(context, second["DeltaU1ZZ"], rho)
        dmzz_candidate = _poly_value(context, delta_m_zz_poly, rho)
        axis = _exact_zero(context, rho)
        phi_actual = context.mpf("0") if axis else _widen(
            context, phi_candidate, ezz_bound
        )
        du_actual = context.mpf("0") if axis else _widen(
            context, du_candidate, epsilon * ezz_bound
        )
        dmzz_actual = context.mpf("0") if axis else _widen(
            context, dmzz_candidate, epsilon * ezz_bound
        )
        row = dict(_json_value(parent_axial[name]))
        row.update(
            {
                "rho": _interval_record(context, rho),
                "Phi1_ZZ_candidate": _interval_record(context, phi_candidate),
                "Phi1_ZZ_actual": _interval_record(context, phi_actual),
                "DeltaU1_ZZ_candidate": _interval_record(context, du_candidate),
                "DeltaU1_ZZ_actual": _interval_record(context, du_actual),
                "DeltaM_average_ZZ_candidate": _interval_record(
                    context, dmzz_candidate
                ),
                "DeltaM_average_ZZ_actual": _interval_record(
                    context, dmzz_actual
                ),
                "Phi_ZZ_actual": _interval_record(context, phi_actual),
                "DeltaU_ZZ_actual": _interval_record(context, du_actual),
                "Phi1_ZZ_error_bound": _interval_record(
                    context, context.mpf("0") if axis else ezz_bound
                ),
                "DeltaU1_ZZ_error_bound": _interval_record(
                    context, context.mpf("0") if axis else epsilon * ezz_bound
                ),
                "axis_exact": axis,
                "ordinary_Z_frame": "fixed rho",
                "point_value": False,
            }
        )
        rows[name] = row
        states[name] = {
            "rho": rho,
            "phi_candidate": phi_candidate,
            "du_candidate": du_candidate,
            "phi_actual": phi_actual,
            "du_actual": du_actual,
            "dmzz_candidate": dmzz_candidate,
            "dmzz_actual": dmzz_actual,
            "axis": axis,
            "range_row": name == "uniform_rho",
        }
    return rows, states


def _f0_records(
    context: Any, f0_chart: Mapping[str, Any], s: Any
) -> tuple[dict[str, Any], Any]:
    rows_source = _rows_mapping(f0_chart, "F0_chart.rows")
    labels = ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1")
    rows: dict[str, Any] = {}
    for label in labels:
        source = _as_mapping(rows_source[label], "F0_chart.rows." + label)
        ratio = _record_or_interval(
            context,
            source,
            ("ratio", "relative_profile", "relative_F0"),
            "F0_chart.ratio." + label,
        )
        logf = _record_or_interval(
            context,
            source,
            ("logF0_Z", "ordinary_logF0_Z", "logF0Z"),
            "F0_chart.logF0_Z." + label,
        )
        if not _strict_positive(context, ratio):
            raise KernelInputError("F0 ratio must be positive: " + label)
        rows[label] = {
            "source": _json_value(source),
            "ratio": ratio,
            "logF0_Z": logf,
        }
    endpoint_lower = min(
        _lower(context, rows["-1"]["ratio"]),
        _lower(context, rows["1"]["ratio"]),
    )
    uniform_ratio = _make_interval(context, endpoint_lower, context.mpf("1"))
    B_all = _record_or_interval(
        context, f0_chart, ("B_all", "B_all_record", "positive_B"), "F0_chart.B_all"
    )
    uniform_logf = _make_interval(
        context,
        -_upper(context, B_all) / _lower(context, s),
        _upper(context, B_all) / _lower(context, s),
    )
    rows["uniform_eta"] = {
        "source": "uniform signed eta segment enclosure",
        "ratio": uniform_ratio,
        "logF0_Z": uniform_logf,
    }
    return rows, B_all


def _chart_error_records(
    context: Any,
    shift: Any,
    h: Any,
    t: Any,
    E0: Any,
    EZ: Any,
    EZZ: Any,
    Mphi: Any,
    Mpsi: Any,
    qeta: Any,
) -> dict[str, Any]:
    one = context.mpf("1")
    denom = one - qeta
    tail_phi = _make_interval(
        context,
        context.mpf("0"),
        _upper(context, Mphi) / (one - t) * qeta**3 / (16 * denom),
    )
    tail_psi = _make_interval(
        context,
        context.mpf("0"),
        _upper(context, Mpsi) / (one - t) * qeta**3 / (16 * denom),
    )
    tail_phi_z = _make_interval(
        context,
        context.mpf("0"),
        3 * _upper(context, Mphi) / (16 * _lower(context, h) * (one - t) ** 2)
        * qeta**2
        / denom,
    )
    tail_psi_z = _make_interval(
        context,
        context.mpf("0"),
        3 * _upper(context, Mpsi) / (16 * _lower(context, h) * (one - t) ** 2)
        * qeta**2
        / denom,
    )
    abs_shift = _abs_upper(context, shift)
    common = E0 + abs_shift * EZ + shift**2 * EZZ / 2
    return {
        "tailPhi": tail_phi,
        "tailPsi": tail_psi,
        "tailPhiZ": tail_phi_z,
        "tailPsiZ": tail_psi_z,
        "errorPhi": common + tail_phi,
        "errorDeltaU": common + tail_psi,
        "errorPhiZ": EZ + abs_shift * EZZ + tail_phi_z,
        "errorDeltaUZ": EZ + abs_shift * EZZ + tail_psi_z,
        "qeta": qeta,
        "common_anchor_error": common,
    }


def _analytic_eta_chart(
    context: Any,
    parameters: Mapping[str, Any],
    seed: Mapping[str, Any],
    polynomials: Mapping[str, dict[int, Any]],
    parent: Mapping[str, Any],
    remainder: Mapping[str, Any],
    f0_chart: Mapping[str, Any],
    anchor_states: Mapping[str, Any],
) -> dict[str, Any]:
    epsilon = parameters["epsilon_core"]
    delta = parameters["delta"]
    a = seed["a"]
    h = parameters["h"]
    s = parameters["s"]
    t = _exact_text(context, "4.1") / 20
    E_norm = _record_or_interval(context, remainder, ("E_norm",), "remainder.E_norm")
    E0 = _record_or_interval(context, remainder, ("E0",), "remainder.E0")
    EZ = _record_or_interval(context, remainder, ("EZ",), "remainder.EZ")
    Mphi = _record_or_interval(context, remainder, ("M_Phi",), "remainder.M_Phi")
    Mpsi = _record_or_interval(context, remainder, ("M_Psi",), "remainder.M_Psi")
    uniform_abs_phi = _record_or_interval(
        context, remainder, ("uniform_abs_Phi1",), "remainder.uniform_abs_Phi1"
    )
    ezz_bound = _make_interval(
        context,
        context.mpf("0"),
        2 * _upper(context, E_norm)
        / (9 * _lower(context, h) ** 2 * (1 - t)),
    )
    f0_rows, B_all = _f0_records(context, f0_chart, s)
    eta_labels = ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1")
    rho_labels = ("0", "1/4", "1", "2", "4.1", "uniform_rho")
    axial_parent = _rows_mapping(
        _as_mapping(_lookup(parent, ("axial_rows", "axial")), "parent.axial_rows"),
        "parent.axial_rows",
    )
    average_parent = _rows_mapping(
        _as_mapping(
            _lookup(parent, ("average_rows", "average_Q", "average")),
            "parent.average_rows",
        ),
        "parent.average_rows",
    )
    chart_rows: dict[str, Any] = {}
    chart_states: dict[str, Any] = {}
    old_phi = polynomials["Phi1"]
    old_phi_small = {
        degree: coefficient for degree, coefficient in old_phi.items() if degree >= 1
    }
    old_du = polynomials["DeltaU1"]
    old_phi_z = polynomials["Phi1_Z"]
    old_du_z = polynomials["DeltaU1_Z"]
    radius = _exact_text(context, "4.1")

    for eta_label in eta_labels + ("uniform_eta",):
        zero_eta = eta_label == "0"
        if zero_eta:
            chart_rows[eta_label] = {
                "eta": eta_label,
                "rows": _json_value(axial_parent),
                "average_rows": _json_value(average_parent),
                "parent_eta0_copied_without_recomputation": True,
                "new_second_Z_rows_ref": "anchor_second_Z_rows",
                "point_value": False,
            }
            chart_states[eta_label] = {
                "eta": eta_label,
                "zero": True,
                "rows": None,
                "average_rows": None,
                "coordinate_Z": a,
                "shift": context.mpf("0"),
                "ratio": f0_rows[eta_label]["ratio"],
                "logF0_Z": f0_rows[eta_label]["logF0_Z"],
                "errorPhi": context.mpf("0"),
                "errorDeltaU": context.mpf("0"),
                "errorPhiZ": context.mpf("0"),
                "errorDeltaUZ": context.mpf("0"),
                "uniform_abs_phi_candidate": uniform_abs_phi,
                "uniform_abs_phi_z_candidate": context.mpf("0"),
            }
            continue

        eta = (
            _make_interval(context, context.mpf("-1"), context.mpf("1"))
            if eta_label == "uniform_eta"
            else _exact_text(context, eta_label)
        )
        shift = s * eta
        coordinate_Z = a + shift
        dZ = seed["d"] - 2 * a * shift - shift**2
        LZ = seed["L"] - 2 * delta * a * shift - delta * shift**2
        U0Z = seed["U0"] + 4 * shift
        if not _strict_positive(context, dZ):
            raise KernelInputError("eta chart d_chart must be positive")
        if not _strict_positive(context, LZ):
            raise KernelInputError("eta chart L_chart must be positive")
        abs_shift_upper = _upper(context, _abs_upper(context, shift))
        qeta_upper = abs_shift_upper / (
            _lower(context, h) * (1 - t)
        )
        if not bool(qeta_upper < context.mpf("1")):
            raise KernelInputError("eta Taylor q_eta condition is not certified")
        qeta = _make_interval(context, context.mpf("0"), qeta_upper)
        error_data = _chart_error_records(
            context,
            shift,
            h,
            t,
            E0,
            EZ,
            ezz_bound,
            Mphi,
            Mpsi,
            qeta,
        )
        error_data["errorDeltaU"] = epsilon * error_data["errorDeltaU"]
        error_data["errorDeltaUZ"] = epsilon * error_data["errorDeltaUZ"]
        phi_poly = _poly_add(
            context,
            old_phi,
            _poly_scale(context, old_phi_z, shift),
            _poly_scale(context, anchor_states["__poly_phiZZ"], shift**2 / 2),
        )
        delta_phi_poly = _poly_add(
            context,
            old_phi_small,
            _poly_scale(context, old_phi_z, shift),
            _poly_scale(context, anchor_states["__poly_phiZZ"], shift**2 / 2),
        )
        du_poly = _poly_add(
            context,
            old_du,
            _poly_scale(context, old_du_z, shift),
            _poly_scale(context, anchor_states["__poly_duZZ"], shift**2 / 2),
        )
        phi_z_poly = _poly_add(
            context,
            old_phi_z,
            _poly_scale(context, anchor_states["__poly_phiZZ"], shift),
        )
        du_z_poly = _poly_add(
            context,
            old_du_z,
            _poly_scale(context, anchor_states["__poly_duZZ"], shift),
        )
        abs_phi_candidate = _poly_abs_bound(context, phi_poly, radius)
        abs_phi_z_candidate = _poly_abs_bound(context, phi_z_poly, radius)
        eta_row_states: dict[str, Any] = {}
        eta_rows: dict[str, Any] = {}
        for rho_label in rho_labels:
            rho = (
                _make_interval(context, context.mpf("0"), radius)
                if rho_label == "uniform_rho"
                else _exact_text(context, rho_label)
            )
            axis = _exact_zero(context, rho)
            delta_phi_candidate = _poly_value(context, delta_phi_poly, rho)
            phi_candidate = _poly_value(context, phi_poly, rho)
            du_candidate = _poly_value(context, du_poly, rho)
            phi_z_candidate = _poly_value(context, phi_z_poly, rho)
            du_z_candidate = _poly_value(context, du_z_poly, rho)
            if axis:
                delta_phi_actual = context.mpf("0")
                phi_actual = context.mpf("1")
                du_actual = context.mpf("0")
                phi_z_actual = context.mpf("0")
                du_z_actual = context.mpf("0")
                u_actual = U0Z
                phi_error = context.mpf("0")
                du_error = context.mpf("0")
                phi_z_error = context.mpf("0")
                du_z_error = context.mpf("0")
            else:
                delta_phi_actual = _widen(
                    context, delta_phi_candidate, error_data["errorPhi"]
                )
                phi_actual = 1 + delta_phi_actual
                du_actual = _widen(
                    context, du_candidate, error_data["errorDeltaU"]
                )
                phi_z_actual = _widen(
                    context, phi_z_candidate, error_data["errorPhiZ"]
                )
                du_z_actual = _widen(
                    context, du_z_candidate, error_data["errorDeltaUZ"]
                )
                u_actual = U0Z + du_actual
                phi_error = error_data["errorPhi"]
                du_error = error_data["errorDeltaU"]
                phi_z_error = error_data["errorPhiZ"]
                du_z_error = error_data["errorDeltaUZ"]
            eta_rows[rho_label] = {
                "eta": eta_label,
                "rho": _interval_record(context, rho),
                "coordinate_Z": _interval_record(context, coordinate_Z),
                "Z_offset": _interval_record(context, shift),
                "d_chart": _interval_record(context, dZ),
                "L_chart": _interval_record(context, LZ),
                "U0_chart": _interval_record(context, U0Z),
                "Phi1_candidate": _interval_record(context, phi_candidate),
                "DeltaPhi1_candidate": _interval_record(
                    context, delta_phi_candidate
                ),
                "Phi1_actual": _interval_record(context, phi_actual),
                "DeltaPhi1_actual": _interval_record(
                    context, delta_phi_actual
                ),
                "DeltaU1_candidate": _interval_record(context, du_candidate),
                "DeltaU1_actual": _interval_record(context, du_actual),
                "U_actual": _interval_record(context, u_actual),
                "U1_actual": _interval_record(context, u_actual),
                "Phi1_Z_candidate": _interval_record(
                    context, phi_z_candidate
                ),
                "Phi1_Z_actual": _interval_record(context, phi_z_actual),
                "DeltaU1_Z_candidate": _interval_record(
                    context, du_z_candidate
                ),
                "DeltaU1_Z_actual": _interval_record(context, du_z_actual),
                "Phi_Z_actual": _interval_record(context, phi_z_actual),
                "DeltaU_Z_actual": _interval_record(context, du_z_actual),
                "U_Z_actual": _interval_record(
                    context, 4 if axis else 4 + du_z_actual
                ),
                "Phi1_error_bound": _interval_record(context, phi_error),
                "DeltaU1_error_bound": _interval_record(context, du_error),
                "Phi1_Z_error_bound": _interval_record(context, phi_z_error),
                "DeltaU1_Z_error_bound": _interval_record(context, du_z_error),
                "axis_exact": axis,
                "point_value": False,
            }
            eta_row_states[rho_label] = {
                "rho": rho,
                "coordinate_Z": coordinate_Z,
                "shift": shift,
                "dZ": dZ,
                "LZ": LZ,
                "U0Z": U0Z,
                "phi_poly": phi_poly,
                "delta_phi_poly": delta_phi_poly,
                "du_poly": du_poly,
                "phi_z_poly": phi_z_poly,
                "du_z_poly": du_z_poly,
                "phi_candidate": phi_candidate,
                "delta_phi_candidate": delta_phi_candidate,
                "du_candidate": du_candidate,
                "phi_z_candidate": phi_z_candidate,
                "du_z_candidate": du_z_candidate,
                "phi_actual": phi_actual,
                "delta_phi_actual": delta_phi_actual,
                "du_actual": du_actual,
                "phi_z_actual": phi_z_actual,
                "du_z_actual": du_z_actual,
                "u_actual": u_actual,
                "phi_error": phi_error,
                "du_error": du_error,
                "phi_z_error": phi_z_error,
                "du_z_error": du_z_error,
                "axis": axis,
            }
        chart_rows[eta_label] = {
            "eta": eta_label,
            "eta_interval": _interval_record(context, eta),
            "shift": _interval_record(context, shift),
            "coordinate_Z": _interval_record(context, coordinate_Z),
            "d_chart": _interval_record(context, dZ),
            "L_chart": _interval_record(context, LZ),
            "U0_chart": _interval_record(context, U0Z),
            "ratio": _interval_record(context, f0_rows[eta_label]["ratio"]),
            "logF0_Z": _interval_record(context, f0_rows[eta_label]["logF0_Z"]),
            "q_eta": _interval_record(context, qeta),
            "tailPhi": _interval_record(context, error_data["tailPhi"]),
            "tailPsi": _interval_record(context, error_data["tailPsi"]),
            "tailPhiZ": _interval_record(context, error_data["tailPhiZ"]),
            "tailPsiZ": _interval_record(context, error_data["tailPsiZ"]),
            "errorPhi": _interval_record(context, error_data["errorPhi"]),
            "errorDeltaU": _interval_record(context, error_data["errorDeltaU"]),
            "errorPhiZ": _interval_record(context, error_data["errorPhiZ"]),
            "errorDeltaUZ": _interval_record(context, error_data["errorDeltaUZ"]),
            "uniform_abs_Phi1_candidate": _interval_record(
                context, abs_phi_candidate
            ),
            "uniform_abs_Phi1_Z_candidate": _interval_record(
                context, abs_phi_z_candidate
            ),
            "rows": eta_rows,
            "point_value": False,
        }
        chart_states[eta_label] = {
            "eta": eta_label,
            "zero": False,
            "eta_interval": eta,
            "ratio": f0_rows[eta_label]["ratio"],
            "logF0_Z": f0_rows[eta_label]["logF0_Z"],
            "shift": shift,
            "coordinate_Z": coordinate_Z,
            "dZ": dZ,
            "LZ": LZ,
            "U0Z": U0Z,
            "qeta": qeta,
            "rows": eta_row_states,
            "phi_poly": phi_poly,
            "delta_phi_poly": delta_phi_poly,
            "du_poly": du_poly,
            "phi_z_poly": phi_z_poly,
            "du_z_poly": du_z_poly,
            "uniform_abs_phi_candidate": abs_phi_candidate,
            "uniform_abs_phi_z_candidate": abs_phi_z_candidate,
            "errorPhi": error_data["errorPhi"],
            "errorDeltaU": error_data["errorDeltaU"],
            "errorPhiZ": error_data["errorPhiZ"],
            "errorDeltaUZ": error_data["errorDeltaUZ"],
        }

    return {
        "payload": {
            "certificate": {
                "coordinate": "Z=a+s*eta",
                "ordinary_Z_is_not_eta": True,
                "a_is_exact_saved_anchor": True,
                "rho_domain": "[0,4.1]",
                "eta_domain": "[-1,1]",
                "all_rho_eta_functional_enclosures": True,
                "M_Phi_and_M_Psi_are_norm_majorants_not_values": True,
                "anchor_second_Z_error": "2*E_norm.upper/(9*h.lower^2*(1-4.1/20))",
                "tailPhi": "M_Phi/(1-t)*q_eta^3/(16*(1-q_eta))",
                "tailPhiZ": "3*M_Phi/(16*h*(1-t)^2)*q_eta^2/(1-q_eta)",
                "eta0_parent_rows_copied": True,
            },
            "eta_rows": chart_rows,
            "table": chart_rows,
            "function_definition": {
                "Z": "a+s*eta",
                "dZ": "d-2*a*shift-shift^2",
                "LZ": "L-2*delta*a*shift-delta*shift^2",
                "U0Z": "U0+4*shift",
                "Phi1": "Phi1_old+shift*Phi1_Z+shift^2*Phi1_ZZ/2",
                "DeltaPhi1": "(Phi1-1)_old+shift*Phi1_Z+shift^2*Phi1_ZZ/2",
                "DeltaU1": "DeltaU1_old+shift*DeltaU1_Z+shift^2*DeltaU1_ZZ/2",
                "Phi1_Z": "Phi1_Z_old+shift*Phi1_ZZ",
                "DeltaU1_Z": "DeltaU1_Z_old+shift*DeltaU1_ZZ",
                "F0": "F0_anchor*relative_profile_ratio_eta",
                "F0_anchor": "exp(-logC), formal positive atom",
            },
            "f0_rows": {
                name: {
                    "ratio": _interval_record(context, row["ratio"]),
                    "logF0_Z": _interval_record(context, row["logF0_Z"]),
                    "source": row["source"],
                }
                for name, row in f0_rows.items()
            },
            "B_all": _interval_record(context, B_all),
            "exact_G_definition": str(
                _lookup(
                    f0_chart,
                    ("exact_G_definition", "G_definition", "definition"),
                )
                or ""
            ),
            "point_value": False,
        },
        "states": chart_states,
        "ezz_bound": ezz_bound,
        "f0_rows": f0_rows,
        "B_all": B_all,
        "uniform_abs_phi": uniform_abs_phi,
        "parent_average": average_parent,
    }


def _compatible_pressure_Q(
    context: Any,
    parameters: Mapping[str, Any],
    seed: Mapping[str, Any],
    source_scale: Mapping[str, Any],
    pressure4: Mapping[str, Any],
    chart: Mapping[str, Any],
    parent: Mapping[str, Any],
) -> dict[str, Any]:
    epsilon = parameters["epsilon_core"]
    delta = parameters["delta"]
    B_anchor = _record_or_interval(
        context, source_scale, ("B_bound", "B"), "source_scale.B_bound"
    )
    S_anchor = _record_or_interval(
        context, source_scale, ("S_bound", "S"), "source_scale.S_bound"
    )
    parent_average = _rows_mapping(
        _as_mapping(
            _lookup(parent, ("average_rows", "average_Q", "average")),
            "parent.average_rows",
        ),
        "parent.average_rows",
    )
    result_rows: dict[str, Any] = {}
    q_states: dict[str, Any] = {}
    for eta_label, state in chart["states"].items():
        if state["zero"]:
            result_rows[eta_label] = {
                "eta": eta_label,
                "rows": _json_value(parent_average),
                "parent_eta0_copied_without_recomputation": True,
                "point_value": False,
            }
            q_states[eta_label] = {
                "zero": True,
                "rows": None,
                "ratio": state["ratio"],
                "logF0_Z": state["logF0_Z"],
                "S_eta": S_anchor,
                "B_eta": B_anchor,
            }
            continue

        ratio = state["ratio"]
        ratio2 = ratio**2
        S_eta = S_anchor * ratio2
        B_eta = B_anchor * ratio2
        phi_poly = state["phi_poly"]
        phi_z_poly = state["phi_z_poly"]
        abs_phi = state["uniform_abs_phi_candidate"]
        abs_phi_z = state["uniform_abs_phi_z_candidate"]
        chart_q_rows: dict[str, Any] = {}
        q_row_states: dict[str, Any] = {}
        shift = state["shift"]
        logf = state["logF0_Z"]
        p0_candidate = (
            seed["p0"]
            + seed["p1"] * shift
            + seed["p2"] * shift**2 / 2
            + seed["p3"] * shift**3 / 6
        )
        p0_error = pressure4["uniformP4sup"] * _abs_upper(
            context, shift
        ) ** 4 / 24
        p0_actual = _widen(context, p0_candidate, p0_error)
        p1_candidate = seed["p1"] + seed["p2"] * shift + seed["p3"] * shift**2 / 2
        p1_error = pressure4["uniformP4sup"] * _abs_upper(
            context, shift
        ) ** 3 / 6
        p1_actual = _widen(context, p1_candidate, p1_error)
        delta_N = (
            ((1 + delta) * seed["U0"] + (12 + 4 * delta) * seed["a"]) * shift
            + (8 + 4 * delta) * shift**2
        )
        delta_L = -2 * delta * seed["a"] * shift - delta * shift**2
        delta_Q0 = (delta_N - seed["Q0"] * delta_L) / state["LZ"]
        Q0_chart = seed["Q0"] + delta_Q0
        for rho_label, row_state in state["rows"].items():
            rho = row_state["rho"]
            axis = row_state["axis"]
            phi = row_state["phi_candidate"]
            phi_z = row_state["phi_z_candidate"]
            du = row_state["du_candidate"]
            du_z = row_state["du_z_candidate"]
            delta_m_poly = _poly_average(context, state["du_poly"])
            delta_mz_poly = _poly_average(context, state["du_z_poly"])
            dm_candidate = _poly_value(context, delta_m_poly, rho)
            dmz_candidate = _poly_value(context, delta_mz_poly, rho)
            if axis:
                dm_actual = context.mpf("0")
                dmz_actual = context.mpf("0")
                m_actual = row_state["U0Z"]
                mz_actual = context.mpf("4")
                delta_q = context.mpf("0")
                q_nonlinear = context.mpf("0")
                q_source = Q0_chart
            else:
                dm_actual = _widen(context, dm_candidate, row_state["du_error"])
                dmz_actual = _widen(
                    context, dmz_candidate, row_state["du_z_error"]
                )
                m_actual = row_state["U0Z"] + dm_actual
                mz_actual = 4 + dmz_actual
                delta_q = (
                    2 * state["coordinate_Z"] * row_state["du_actual"]
                    - (1 - delta) * state["coordinate_Z"] * dm_actual
                    - state["dZ"] * dmz_actual
                ) / state["LZ"]
                q_nonlinear = delta_q
                q_source = Q0_chart + q_nonlinear
            q0_actual = Q0_chart
            q_source_candidate = (
                Q0_chart
                + (
                    2 * state["coordinate_Z"] * du
                    - (1 - delta) * state["coordinate_Z"] * dm_candidate
                    - state["dZ"] * dmz_candidate
                )
                / state["LZ"]
                if not axis
                else Q0_chart
            )
            q_source_actual = q_source
            phi_sq = _poly_mul(context, phi_poly, phi_poly)
            phi_phiz = _poly_mul(context, phi_poly, phi_z_poly)
            i_poly = _poly_integral(context, phi_sq)
            iz_poly = _poly_scale(
                context, _poly_integral(context, phi_phiz), context.mpf("2")
            )
            i_candidate = _poly_value(context, i_poly, rho)
            iz_candidate = _poly_value(context, iz_poly, rho)
            rho_upper = _upper(context, rho)
            ephi = row_state["phi_error"]
            ephiz = row_state["phi_z_error"]
            i_error = rho_upper * ephi * (2 * _upper(context, abs_phi) + ephi)
            iz_error = 2 * rho_upper * (
                ephi * _upper(context, abs_phi_z)
                + ephiz * _upper(context, abs_phi)
                + ephi * ephiz
            )
            i_actual = (
                context.mpf("0")
                if axis
                else _widen(
                    context,
                    i_candidate,
                    _make_interval(context, context.mpf("0"), i_error),
                )
            )
            iz_actual = (
                context.mpf("0")
                if axis
                else _widen(
                    context,
                    iz_candidate,
                    _make_interval(context, context.mpf("0"), iz_error),
                )
            )
            epsP_candidate = p0_candidate + B_eta * i_candidate
            epsP_actual = p0_actual + B_eta * i_actual
            epsPz_candidate = p1_candidate + B_eta * (
                2 * logf * i_candidate + iz_candidate
            )
            epsPz_actual = p1_actual + B_eta * (
                2 * logf * i_actual + iz_actual
            )
            if axis:
                epsP_actual = p0_actual
                epsPz_actual = p1_actual
                epsP_candidate = p0_candidate
                epsPz_candidate = p1_candidate
            chart_q_rows[rho_label] = {
                "eta": eta_label,
                "rho": _interval_record(context, rho),
                "coordinate_Z": _interval_record(context, state["coordinate_Z"]),
                "d_chart": _interval_record(context, state["dZ"]),
                "L_chart": _interval_record(context, state["LZ"]),
                "ratio_eta": _interval_record(context, ratio),
                "S_eta": _interval_record(context, S_eta),
                "B_eta": _interval_record(context, B_eta),
                "DeltaM_average_candidate": _interval_record(context, dm_candidate),
                "DeltaM_average_actual": _interval_record(context, dm_actual),
                "M_average": _interval_record(context, m_actual),
                "DeltaM_average_Z_candidate": _interval_record(
                    context, dmz_candidate
                ),
                "DeltaM_average_Z_actual": _interval_record(context, dmz_actual),
                "M_average_Z": _interval_record(context, mz_actual),
                "Q0_chart": _interval_record(context, q0_actual),
                "DeltaQ0": _interval_record(context, delta_Q0),
                "Q_nonlinear_candidate": _interval_record(
                    context, q_source_candidate - Q0_chart
                ),
                "Q_nonlinear_actual": _interval_record(context, q_nonlinear),
                "DeltaQ": _interval_record(context, q_source_actual - Q0_chart),
                "Q_source_candidate": _interval_record(
                    context, q_source_candidate
                ),
                "Q_source": _interval_record(context, q_source_actual),
                "Q_chart": _interval_record(context, q_source_actual),
                "I_candidate": _interval_record(context, i_candidate),
                "I_actual": _interval_record(context, i_actual),
                "I_error_bound": _interval_record(
                    context, _make_interval(context, context.mpf("0"), i_error)
                ),
                "I_Z_candidate": _interval_record(context, iz_candidate),
                "I_Z_actual": _interval_record(context, iz_actual),
                "I_Z_error_bound": _interval_record(
                    context, _make_interval(context, context.mpf("0"), iz_error)
                ),
                "p0_chart_candidate": _interval_record(context, p0_candidate),
                "p0_chart": _interval_record(context, p0_actual),
                "p0_chart_error_bound": _interval_record(context, p0_error),
                "p1_chart_candidate": _interval_record(context, p1_candidate),
                "p1_chart": _interval_record(context, p1_actual),
                "p1_chart_error_bound": _interval_record(context, p1_error),
                "epsilonP_candidate": _interval_record(context, epsP_candidate),
                "epsilonP": _interval_record(context, epsP_actual),
                "epsilonP_Z_candidate": _interval_record(
                    context, epsPz_candidate
                ),
                "epsilonP_Z": _interval_record(context, epsPz_actual),
                "rawP": _interval_record(context, epsP_actual / epsilon),
                "rawP_Z": _interval_record(context, epsPz_actual / epsilon),
                "formal_S_identity": "S_eta=S_anchor*ratio_eta^2; positive formal atom",
                "S_Z_over_S": _interval_record(context, 2 * logf),
                "Q_source_definition": "[2*Z*DeltaU-(1-delta)*Z*DeltaM-dZ*DeltaM_Z]/LZ",
                "M_definition": "U0_chart+(1/rho)*integral_0^rho DeltaU",
                "M_Z_definition": "4+(1/rho)*integral_0^rho DeltaU_Z",
                "axis_exact": axis,
                "parent_eta0_copied_without_recomputation": False,
                "point_value": False,
            }
            q_row_states[rho_label] = {
                "rho": rho,
                "Q_source": q_source_actual,
                "Q_candidate": q_source_candidate,
                "M_average": m_actual,
                "M_average_Z": mz_actual,
                "DeltaM": dm_actual,
                "DeltaM_Z": dmz_actual,
                "axis": axis,
            }
        result_rows[eta_label] = {
            "eta": eta_label,
            "rows": chart_q_rows,
            "coordinate_Z": _interval_record(context, state["coordinate_Z"]),
            "ratio": _interval_record(context, ratio),
            "logF0_Z": _interval_record(context, logf),
            "S_eta": _interval_record(context, S_eta),
            "B_eta": _interval_record(context, B_eta),
            "Q0_chart": _interval_record(context, Q0_chart),
            "DeltaQ0": _interval_record(context, delta_Q0),
            "pressure_units": "epsilon_P=epsilon_core*P in retained source units",
            "point_value": False,
        }
        q_states[eta_label] = {
            "zero": False,
            "rows": q_row_states,
            "ratio": ratio,
            "logF0_Z": logf,
            "S_eta": S_eta,
            "B_eta": B_eta,
        }
    return {
        "payload": {
            "rows": result_rows,
            "table": result_rows,
            "function_definition": {
                "Q0_chart": "Q0a+(DeltaN-Q0a*DeltaL)/LZ",
                "DeltaN": "[(1+delta)*U0+(12+4*delta)*a]*shift+(8+4*delta)*shift^2",
                "DeltaL": "-2*delta*a*shift-delta*shift^2",
                "Q_nonlinear": "[2*Z*DeltaU-(1-delta)*Z*DeltaM-dZ*DeltaM_Z]/LZ",
                "I": "V(Phi_eta^2)",
                "I_Z": "2*V(Phi_eta*Phi_Z_eta)",
                "epsilonP": "p0_chart+B_eta*I",
                "epsilonP_Z": "p1_chart+B_eta*(2*logF0Z*I+I_Z)",
                "F0_ratio": "saved seven-row DATA; uniform endpoint enclosure",
            },
            "eta0_parent_average_copied_without_recomputation": True,
            "source_scale_formal_positive": True,
            "point_value": False,
        },
        "payload": result_rows,
        "states": q_states,
    }


def _affine_reference(
    expression: str, base_field: str, base_sign: int, f0_field: str, f0_sign: int
) -> dict[str, Any]:
    return {
        "formal_expression": expression,
        "constant_term_ref": {"field": base_field, "sign": base_sign},
        "F0_coefficient_ref": {"field": f0_field, "sign": f0_sign},
        "canonical_numeric_coefficients": True,
        "absolute_F0_materialized": False,
        "point_value": False,
    }


def _physical_eta_velocity(
    context: Any,
    parameters: Mapping[str, Any],
    seed: Mapping[str, Any],
    chart: Mapping[str, Any],
    q_data: Mapping[str, Any],
    parent: Mapping[str, Any],
) -> dict[str, Any]:
    delta = parameters["delta"]
    epsilon = parameters["epsilon_core"]
    tau_labels = ("1", "1/4", "1/16", "1/64")
    rho_labels = ("0", "1/4", "1", "2", "4.1", "uniform_rho")
    angle_labels = ("0", "pi/2", "pi", "3pi/2")
    rotations = {
        "0": ("1", "0"),
        "pi/2": ("0", "1"),
        "pi": ("-1", "0"),
        "3pi/2": ("0", "-1"),
    }
    parent_physical = _physical_mapping(
        _as_mapping(
            _lookup(parent, ("physical_time_slices", "physical_velocity", "physical")),
            "parent.physical_time_slices",
        ),
        "parent.physical_time_slices",
    )
    eta_slices: dict[str, Any] = {}
    eta_states: dict[str, Any] = {}
    for eta_label, state in chart["states"].items():
        if state["zero"]:
            eta_slices[eta_label] = {
                "eta": eta_label,
                "time_slices": _json_value(parent_physical),
                "parent_eta0_physical_copied_without_recomputation": True,
                "point_value": False,
            }
            eta_states[eta_label] = {"zero": True}
            continue
        q_state = q_data["states"][eta_label]
        time_slices: dict[str, Any] = {}
        for tau_label in tau_labels:
            tau = _exact_text(context, tau_label)
            tphys = 1 - tau
            dZ = state["dZ"]
            if not _strict_positive(context, dZ):
                raise KernelInputError("physical chart d_chart must be positive")
            lam = context.sqrt(tau / dZ)
            kappa = lam ** (1 - delta)
            rho_slices: dict[str, Any] = {}
            for rho_label in rho_labels:
                row_state = state["rows"][rho_label]
                q_row = q_state["rows"][rho_label]
                rho = row_state["rho"]
                R = epsilon * rho
                r_phys = lam * context.sqrt(2 * R)
                z_phys = kappa * state["coordinate_Z"]
                Q = q_row["Q_source"]
                U = row_state["u_actual"]
                phi = row_state["phi_actual"]
                if row_state["axis"]:
                    ur = context.mpf("0")
                    swirl = context.mpf("0")
                    uz = lam ** (-1 - delta) * state["U0Z"]
                else:
                    ur = lam ** (-1) * context.sqrt(R / 2) * Q
                    swirl = (
                        lam ** (-1 - delta)
                        * context.sqrt(2 * R)
                        * q_state["ratio"]
                        * phi
                    )
                    uz = lam ** (-1 - delta) * U
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
                        "u": _affine_reference(
                            "ur*cos(theta)-F0_anchor*swirlcoef*sin(theta)",
                            "ur",
                            base_u_sign,
                            "swirlcoef",
                            f0_u_sign,
                        ),
                        "v": _affine_reference(
                            "ur*sin(theta)+F0_anchor*swirlcoef*cos(theta)",
                            "ur",
                            base_v_sign,
                            "swirlcoef",
                            f0_v_sign,
                        ),
                        "w": {
                            "canonical_field_ref": "uz",
                            "point_value": False,
                        },
                        "cartesian_formula": {
                            "u": "ur*cos(theta)-F0_anchor*swirlcoef*sin(theta)",
                            "v": "ur*sin(theta)+F0_anchor*swirlcoef*cos(theta)",
                            "w": "uz",
                        },
                        "point_value": False,
                    }
                rho_slices[rho_label] = {
                    "rho": _interval_record(context, rho),
                    "R": _interval_record(context, R),
                    "coordinate_Z": _interval_record(context, state["coordinate_Z"]),
                    "Z_offset": _interval_record(context, state["shift"]),
                    "d_chart": _interval_record(context, dZ),
                    "lambda": _interval_record(context, lam),
                    "kappa": _interval_record(context, kappa),
                    "r_phys": _interval_record(context, r_phys),
                    "z_phys": _interval_record(context, z_phys),
                    "t_phys": _interval_record(context, tphys),
                    "Q_source": _interval_record(context, Q),
                    "U": _interval_record(context, U),
                    "ur": _interval_record(context, ur),
                    "uz": _interval_record(context, uz),
                    "swirlcoef": _interval_record(context, swirl),
                    "actualutheta": {
                        "formal_expression": "F0_anchor*swirlcoef",
                        "F0_anchor_formal": "exp(-logC)",
                        "F0_coefficient": _interval_record(context, swirl),
                        "absolute_F0_materialized": False,
                        "positive_formal_factor": True,
                        "point_value": False,
                    },
                    "angles": angle_ports,
                    "streamfunction": {
                        "psi": _interval_record(
                            context, kappa * R * q_row["M_average"]
                        ),
                        "formal_expression": "psi=kappa*R*M_average",
                        "ur_pullback": "ur=-psi_z/r_phys",
                        "uz_pullback": "uz=psi_r/r_phys",
                        "mixed_partials_divergence_zero": True,
                        "numeric_cartesian_NS_check": False,
                        "point_value": False,
                    },
                    "axis_regularization": {
                        "r_zero": row_state["axis"],
                        "ur_zero": row_state["axis"],
                        "swirlcoef_zero": row_state["axis"],
                        "Phi1_axis_exact": row_state["axis"],
                        "DeltaPhi_axis_exact": row_state["axis"],
                        "DeltaU_axis_exact": row_state["axis"],
                        "Phi1_Z_axis_exact": row_state["axis"],
                        "DeltaU_Z_axis_exact": row_state["axis"],
                    },
                    "point_value": False,
                }
            time_slices[tau_label] = {
                "tau": _interval_record(context, tau),
                "t_phys": _interval_record(context, tphys),
                "lambda": _interval_record(context, lam),
                "kappa": _interval_record(context, kappa),
                "rho_slices": rho_slices,
                "positive_tau": True,
                "positive_power_operations_allowed": True,
                "point_value": False,
            }
        eta_slices[eta_label] = {
            "eta": eta_label,
            "coordinate_Z": _interval_record(context, state["coordinate_Z"]),
            "Z_offset": _interval_record(context, state["shift"]),
            "time_slices": time_slices,
            "parent_eta0_physical_copied_without_recomputation": False,
            "point_value": False,
        }
        eta_states[eta_label] = {"zero": False, "time_slices": time_slices}
    return {
        "payload": {
            "eta_rows": eta_slices,
            "physical_time_slices": eta_slices,
            "domain": {
                "Z": "exact atom a+s*eta; floating overlap is not identity",
                "rho": list(rho_labels[:5]),
                "uniform_rho": ["0", "4.1"],
                "eta": ["-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1"],
                "physical_tau": list(tau_labels),
                "physical_azimuth_labels": list(angle_labels),
                "source_phase_is_not_physical_azimuth": True,
                "source_profile_scale": "Lambda_source=1/epsilon_core",
                "physical_inverse_scale": "lambda=sqrt(tau/d_chart)",
                "coordinate_Z_offset_retained": True,
            },
            "conditional_streamfunction_divergence": {
                "frame": "physical cylindrical",
                "streamfunction": "psi=kappa*R*M_average",
                "certificate": "mixed partial cancellation; structural conditional only",
                "actual_field_source": "same fixed analytic source; no fitted independent components",
                "numeric_cartesian_NS_claim": False,
                "full_physical_field_materialized": False,
            },
            "viscosity": {"paper_unit": 1, "residual_defined_force": False},
            "point_value": False,
        },
        "states": eta_states,
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
    """Evaluate one bounded second-Z/eta source leaf."""

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
    seed_source = _required_mapping(capsule, ("seed", "seed_inputs"), "capsule.seed")
    polynomials_source = _required_mapping(
        capsule, ("polynomials", "parent_polynomials"), "capsule.polynomials"
    )
    parent = _as_mapping(capsule["parent"], "capsule.parent")
    remainder = _as_mapping(capsule["remainder"], "capsule.remainder")
    source_scale = _as_mapping(capsule["source_scale"], "capsule.source_scale")
    f0_chart = _as_mapping(capsule["F0_chart"], "capsule.F0_chart")

    _emit(emit, events, "pressure4_C3_entered", {"family": "pressure4_C3", "before_math": True})
    raw_parameters, parameters = _parameter_packet(context, parameters_source)
    raw_seed, seed = _seed_packet(context, seed_source)
    raw_polynomials, polynomials = _polynomial_packet(context, polynomials_source)
    pressure4 = _pressure4_C3(context, parameters, seed)
    _emit(
        emit,
        events,
        "pressure4_C3_complete",
        {
            "family": "pressure4_C3",
            "p4_new_pressure_datum": True,
            "positive_kernel_I4_enclosure": True,
            "uniformP4sup_is_bound_not_function_value": True,
        },
    )

    _emit(emit, events, "second_Z_operator_entered", {"family": "second_Z_operator", "before_math": True})
    second = _second_Z_operator(
        context, parameters, seed, pressure4, source_scale, polynomials
    )
    _emit(
        emit,
        events,
        "second_Z_operator_complete",
        {
            "family": "second_Z_operator",
            "ordinary_Z_before_anchor_restriction": True,
            "chi3_and_s3_retained": True,
            "resolvent_second_term_separate": True,
            "new_polynomial_degree": 3,
        },
    )

    _emit(emit, events, "analytic_eta_chart_entered", {"family": "analytic_eta_chart", "before_math": True})
    e_norm = _record_or_interval(context, remainder, ("E_norm",), "remainder.E_norm")
    t = _exact_text(context, "4.1") / 20
    ezz_bound = _make_interval(
        context,
        context.mpf("0"),
        2 * _upper(context, e_norm)
        / (9 * _lower(context, parameters["h"]) ** 2 * (1 - t)),
    )
    parent_axial = _rows_mapping(
        _as_mapping(_lookup(parent, ("axial_rows", "axial")), "parent.axial_rows"),
        "parent.axial_rows",
    )
    # Internal polynomial ports are kept typed; they are not duplicate output
    # records in every eta row.
    anchor_states: dict[str, Any] = {
        "__poly_phiZZ": second["Phi1ZZ"],
        "__poly_duZZ": second["DeltaU1ZZ"],
    }
    anchor_rows, anchor_row_states = _anchor_second_rows(
        context,
        parent_axial,
        second,
        ezz_bound,
        parameters["epsilon_core"],
    )
    anchor_states.update(anchor_row_states)
    chart = _analytic_eta_chart(
        context,
        parameters,
        seed,
        polynomials,
        parent,
        remainder,
        f0_chart,
        anchor_states,
    )
    _emit(
        emit,
        events,
        "analytic_eta_chart_complete",
        {
            "family": "analytic_eta_chart",
            "eta_rows": 7,
            "uniform_eta_row": True,
            "anchor_second_Z_rows": 6,
            "eta0_parent_copied": True,
            "functional_chart_certificate": True,
        },
    )

    _emit(emit, events, "compatible_pressure_Q_entered", {"family": "compatible_pressure_Q", "before_math": True})
    q_data = _compatible_pressure_Q(
        context, parameters, seed, source_scale, pressure4, chart, parent
    )
    _emit(
        emit,
        events,
        "compatible_pressure_Q_complete",
        {
            "family": "compatible_pressure_Q",
            "Q_background_kept_separate": True,
            "M_average_not_raw_Mz": True,
            "old_eta0_pressure_Q_copied": True,
            "S_Z_over_S_retained": True,
        },
    )

    _emit(emit, events, "physical_eta_velocity_entered", {"family": "physical_eta_velocity", "before_math": True})
    physical = _physical_eta_velocity(
        context, parameters, seed, chart, q_data, parent
    )
    _emit(
        emit,
        events,
        "physical_eta_velocity_complete",
        {
            "family": "physical_eta_velocity",
            "eta_rows": 7,
            "tau_slices": 4,
            "rho_slices": 6,
            "azimuth_labels": 4,
            "eta0_parent_physical_copied": True,
            "formal_F0_products_retained": True,
            "physical_field_materialized": False,
        },
    )

    proof = {
        "kind": "native_O3_second_Z_eta_chart_proof_v1",
        "proof_count": 1,
        "pressure4_C3": "positive-kernel I4 enclosure and LC=N/2 C3 seed",
        "second_Z_operator": "full ordinary-Z Picard map differentiated before H0(a)=0",
        "Jnu": "J_nu(rho^n)=rho^(n+1)/((n+1)*(n+nu))",
        "resolvent": "DeltaPhi1ZZ keeps theta and -(chi2/2)J2(Phi1-1) terms separate",
        "eta_chart": "Z=a+s*eta with common analytic h Taylor tails on rho<=4.1",
        "q_eta_condition": "abs(shift)/(h*(1-rho/20))<1",
        "pressure_Q": "Q0 background and nonlinear correction remain separate",
        "physical_velocity": "lambda=sqrt(tau/d_chart), kappa=lambda^(1-delta)",
        "axis": "rho=0 Phi1/DeltaPhi/DeltaU/PhiZ/DeltaUZ exact atoms",
        "eta0": "accepted parent rows, Q, pressure, and physical time slices copied",
        "source_scale": "B=epsilon^2*S is formal positive source DATA",
        "point_value_selection": False,
        "global_jet_admission": False,
        "global_Q": False,
        "full_physical": False,
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
            "local_eta_functional_enclosures": True,
            "local_second_Z_jet": True,
            "local_Q": True,
            "local_coefficients": True,
            "global_jet_admission": False,
            "global_Q": False,
            "full_physical": False,
            "global_five_moment": False,
            "time_n_recursion": False,
            "full_goal": False,
        },
    )

    parameter_data = {
        name: {
            "raw": raw_parameters[name],
            "interval": _interval_record(context, parameters[name]),
        }
        for name in ("delta", "epsilon_core", "j", "sigma", "h", "s")
    }
    seed_data = {
        name: _interval_record(context, seed[name]) for name in raw_seed
    }
    pressure_payload = {
        name: _interval_record(context, pressure4[name])
        for name in ("H3", "I4", "Pnorm4", "p4", "N3", "C3", "uniformP4sup")
    }
    second_payload = {
        "combinations": {
            name: _interval_record(context, value)
            for name, value in second["combinations"].items()
        },
        "Ttheta0": _poly_records(context, second["Ttheta0"]),
        "Ttheta1": _poly_records(context, second["Ttheta1"]),
        "Ttheta2": _poly_records(context, second["Ttheta2"]),
        "Ebar_theta2": _poly_records(context, second["Ebar_theta2"]),
        "J2_Ebar_theta2": _poly_records(context, second["J2_Ebar_theta2"]),
        "J2_Phi1_minus_one": _poly_records(
            context, second["J2_Phi1_minus_one"]
        ),
        "DeltaPhi1ZZ_theta_term": _poly_records(
            context, second["DeltaPhi1ZZ_theta_term"]
        ),
        "DeltaPhi1ZZ_resolvent_term": _poly_records(
            context, second["DeltaPhi1ZZ_resolvent_term"]
        ),
        "RchiZZterm": _poly_records(
            context, second["DeltaPhi1ZZ_resolvent_term"]
        ),
        "DeltaPhi1ZZ": _poly_records(context, second["DeltaPhi1ZZ"]),
        "Phi1ZZ": _poly_records(context, second["Phi1ZZ"]),
        "Tbar0": _poly_records(context, second["Tbar0"]),
        "Tbar1": _poly_records(context, second["Tbar1"]),
        "Tbar2": _poly_records(context, second["Tbar2"]),
        "Tbar2_linear": [
            row for row in _poly_records(context, second["Tbar2"]) if row["degree"] == 1
        ],
        "Tbar2_quadratic": [
            row for row in _poly_records(context, second["Tbar2"]) if row["degree"] == 2
        ],
        "Ebar_z2": _poly_records(context, second["Ebar_z2"]),
        "DeltaU1ZZ": _poly_records(context, second["DeltaU1ZZ"]),
        "definitions": _json_value(second["definitions"]),
        "source_units": second["source_units"],
        "point_value": False,
    }
    capsule_family = _lookup(capsule, ("family_identity", "family_identity_hashes"))
    pressure4_packet = {
        "records": pressure_payload,
        "definitions": _json_value(pressure4["definitions"]),
        "point_value": False,
    }
    anchor_second_packet = {
        "rows": anchor_rows,
        "error_bound": _interval_record(context, ezz_bound),
        "MavgZZ_definition": "averageDeltaUZZpoly; coefficient_n divided by n+1",
        "MavgZZ_error_bound": "epsilon_core*EZZ",
        "rows_are_independent_new_records": True,
        "axis_rho0_exact": True,
        "point_value": False,
    }
    payload = {
        "component": COMPONENT,
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "conditional": True,
        "full_physical": False,
        "global_jet_admission": False,
        "global_Q": False,
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
            "raw_capsule_sha256_receipt": ROOT_CAPSULE_RAW_SHA,
            "canonical_field_in_capsule": False,
        },
        "source_ids": _json_value(capsule["source_ids"]),
        "parameters": parameter_data,
        "seed": {
            "accepted_records": seed_data,
            "raw": raw_seed,
            "parent_seed_consumed_without_recomputation": True,
        },
        "polynomials": {
            "accepted_records": _json_value(raw_polynomials),
            "Phi1": _poly_records(context, polynomials["Phi1"]),
            "DeltaU1": _poly_records(context, polynomials["DeltaU1"]),
            "Phi1_Z": _poly_records(context, polynomials["Phi1_Z"]),
            "DeltaU1_Z": _poly_records(context, polynomials["DeltaU1_Z"]),
            "parent_polynomials_consumed_without_recomputation": True,
        },
        "pressure4_C3": pressure4_packet,
        "pressure4_seed": pressure4_packet,
        "second_Z_operator": second_payload,
        "second_Z_polynomial": second_payload,
        "anchor_second_Z_rows": anchor_second_packet,
        "second_Z_anchor_rows": anchor_second_packet,
        "analytic_eta_chart": chart["payload"],
        "compatible_pressure_Q": q_data["payload"],
        "physical_eta_velocity": physical["payload"],
        "source_scale": {
            "S_bound": _interval_record(
                context,
                _record_or_interval(context, source_scale, ("S_bound", "S"), "source_scale.S_bound"),
            ),
            "B_bound": _interval_record(
                context,
                _record_or_interval(context, source_scale, ("B_bound", "B"), "source_scale.B_bound"),
            ),
            "logC_record": _json_value(source_scale["logC_record"]),
            "formal_positive": True,
            "formal_positiveTrue": True,
            "materialized": False,
            "formalS": "S=F0_anchor^2=exp(-2*logC), retained formal positive atom",
            "absolute_F0_exp_materialized": False,
            "S_Z_over_S": "2*logF0_Z away from eta=0; saved eta0 DATA copied",
        },
        "family_identity": _json_value(capsule_family),
        "grid": _json_value(capsule["grid"]),
        "preservation": _json_value(
            _lookup(capsule, ("preservation", "preservation_guards", "guards"))
        ),
        "definition_guards": _json_value(
            _lookup(capsule, ("definition_guards", "definitions", "definition"))
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
        "static_limitations": [
            "local analytic eta chart only; no global source jet admission",
            "M_Phi and M_Psi are norm majorants, never function values",
            "eta=0 parent Q/pressure/physical data are copied without recomputation",
            "F0_anchor and S remain positive formal atoms; no absolute exponential materialization",
            "no global Q closure, five-moment closure, time-n recursion, residual force, or Cartesian NS claim",
            "positive physical powers are local chart operations and do not assert generic exp count zero",
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
        "local_eta_functional_enclosures": True,
        "local_second_Z_jet": True,
        "local_Q": True,
        "local_coefficients": True,
        "local_Q_source": True,
        "local_velocity_coefficients": True,
        "global_jet_admission": False,
        "global_Q": False,
        "global_source_jet_admission": False,
        "physical_materialized": False,
        "physical_field_materialized": False,
        "global_five_moment": False,
        "time_n_recursion": False,
        "full_goal": False,
        "full_goal_complete": False,
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
    "SOURCE_PRECISION_BITS",
    "SUCCESS_COUNTS",
    "SUCCESS_FLAGS",
    "derive",
]
